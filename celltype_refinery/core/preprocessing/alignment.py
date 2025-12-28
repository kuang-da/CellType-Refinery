"""Cross-sample alignment (Stage D).

Provides percentile-based alignment to harmonize marker intensities
across samples within a study.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from .config import AlignmentConfig


@dataclass
class AlignmentParams:
    """Alignment parameters for a single marker.

    Attributes
    ----------
    marker : str
        Marker name
    slope : float
        Alignment slope
    intercept : float
        Alignment intercept
    local_lower : float
        Local lower percentile value
    local_upper : float
        Local upper percentile value
    target_lower : float
        Global target lower value
    target_upper : float
        Global target upper value
    """

    marker: str
    slope: float
    intercept: float
    local_lower: float
    local_upper: float
    target_lower: float
    target_upper: float


@dataclass
class AlignmentResult:
    """Result from aligning a single sample.

    Attributes
    ----------
    sample_id : str
        Sample identifier
    aligned_matrix : pd.DataFrame
        Aligned intensity matrix
    params : List[AlignmentParams]
        Per-marker alignment parameters
    """

    sample_id: str
    aligned_matrix: Optional[pd.DataFrame] = None
    params: List[AlignmentParams] = field(default_factory=list)

    def params_to_dataframe(self) -> pd.DataFrame:
        """Convert params to DataFrame."""
        records = []
        for p in self.params:
            records.append({
                "sample_id": self.sample_id,
                "marker": p.marker,
                "slope": p.slope,
                "intercept": p.intercept,
                "local_lower": p.local_lower,
                "local_upper": p.local_upper,
                "target_lower": p.target_lower,
                "target_upper": p.target_upper,
            })
        return pd.DataFrame(records)


class CrossSampleAligner:
    """Cross-sample intensity aligner.

    Aligns marker intensities across samples using percentile-based
    linear transformation to global target values.

    Parameters
    ----------
    config : AlignmentConfig
        Alignment configuration

    Example
    -------
    >>> from celltype_refinery.core.preprocessing import CrossSampleAligner, AlignmentConfig
    >>> config = AlignmentConfig(lower_percentile=5.0, upper_percentile=95.0)
    >>> aligner = CrossSampleAligner(config)
    >>> # First, compute global targets from all samples
    >>> targets = aligner.compute_global_targets(sample_matrices)
    >>> # Then align each sample
    >>> result = aligner.align_sample(matrix, "sample_01", targets)
    """

    EPSILON = 1e-6

    def __init__(self, config: Optional[AlignmentConfig] = None):
        self.config = config or AlignmentConfig()

    def compute_local_percentiles(
        self,
        intensity_df: pd.DataFrame,
    ) -> pd.DataFrame:
        """Compute local percentiles for a sample.

        Parameters
        ----------
        intensity_df : pd.DataFrame
            Marker intensity matrix

        Returns
        -------
        pd.DataFrame
            Percentiles with columns 'lower' and 'upper', indexed by marker
        """
        lower_p = self.config.lower_percentile
        upper_p = self.config.upper_percentile

        records = []
        for marker in intensity_df.columns:
            values = pd.to_numeric(intensity_df[marker], errors="coerce").dropna()
            if values.empty:
                q_lower = float("nan")
                q_upper = float("nan")
            else:
                q_lower = float(np.percentile(values, lower_p))
                q_upper = float(np.percentile(values, upper_p))
            records.append({
                "marker": marker,
                "lower": q_lower,
                "upper": q_upper,
            })

        return pd.DataFrame(records).set_index("marker")

    def compute_global_targets(
        self,
        sample_matrices: Dict[str, pd.DataFrame],
        cell_id_col: str = "cell_mask_id",
    ) -> Dict[str, Dict[str, float]]:
        """Compute global target percentiles from all samples.

        Parameters
        ----------
        sample_matrices : Dict[str, pd.DataFrame]
            Map of sample_id to normalized intensity matrix
        cell_id_col : str
            Cell ID column name

        Returns
        -------
        Dict[str, Dict[str, float]]
            Map of marker to {'target_lower': float, 'target_upper': float}
        """
        lower_p = self.config.lower_percentile
        upper_p = self.config.upper_percentile

        # Collect all values per marker
        marker_values: Dict[str, List[np.ndarray]] = {}

        for sample_id, matrix in sample_matrices.items():
            marker_cols = [c for c in matrix.columns if c != cell_id_col]
            for marker in marker_cols:
                if marker not in marker_values:
                    marker_values[marker] = []
                values = pd.to_numeric(matrix[marker], errors="coerce").dropna()
                if len(values) > 0:
                    marker_values[marker].append(values.to_numpy(dtype=float))

        # Compute global percentiles
        targets: Dict[str, Dict[str, float]] = {}
        for marker, value_arrays in marker_values.items():
            if not value_arrays:
                continue
            combined = np.concatenate(value_arrays)
            targets[marker] = {
                "target_lower": float(np.percentile(combined, lower_p)),
                "target_upper": float(np.percentile(combined, upper_p)),
            }

        return targets

    def align_sample(
        self,
        cell_matrix: pd.DataFrame,
        sample_id: str,
        global_targets: Dict[str, Dict[str, float]],
        cell_id_col: str = "cell_mask_id",
    ) -> AlignmentResult:
        """Align a single sample to global targets.

        Parameters
        ----------
        cell_matrix : pd.DataFrame
            Normalized intensity matrix
        sample_id : str
            Sample identifier
        global_targets : Dict[str, Dict[str, float]]
            Global target percentiles per marker
        cell_id_col : str
            Cell ID column name

        Returns
        -------
        AlignmentResult
            Alignment result with transformed matrix
        """
        result = AlignmentResult(sample_id=sample_id)

        marker_cols = [c for c in cell_matrix.columns if c != cell_id_col]
        intensity_df = cell_matrix[marker_cols].apply(pd.to_numeric, errors="coerce")

        # Compute local percentiles
        local_percentiles = self.compute_local_percentiles(intensity_df)

        # Build aligned matrix
        aligned = pd.DataFrame(index=cell_matrix.index)
        aligned[cell_id_col] = cell_matrix[cell_id_col]

        for marker in marker_cols:
            values = intensity_df[marker].to_numpy(dtype=float)

            # Get local and global values
            local_lower = float(local_percentiles.loc[marker, "lower"])
            local_upper = float(local_percentiles.loc[marker, "upper"])

            target = global_targets.get(marker)
            if target is None:
                aligned[marker] = values
                continue

            target_lower = target["target_lower"]
            target_upper = target["target_upper"]

            # Compute linear transformation
            denom = max(local_upper - local_lower, self.EPSILON)
            slope = (target_upper - target_lower) / denom
            intercept = target_lower - slope * local_lower

            # Apply transformation with optional clipping
            aligned_values = (values - local_lower) * slope + target_lower

            if self.config.clip_to_targets:
                aligned_values = np.clip(aligned_values, target_lower, target_upper)

            aligned[marker] = aligned_values

            # Record params
            result.params.append(AlignmentParams(
                marker=marker,
                slope=slope,
                intercept=intercept,
                local_lower=local_lower,
                local_upper=local_upper,
                target_lower=target_lower,
                target_upper=target_upper,
            ))

        result.aligned_matrix = aligned
        return result

    @staticmethod
    def ks_distance(a: np.ndarray, b: np.ndarray) -> float:
        """Compute Kolmogorov-Smirnov distance between distributions.

        Parameters
        ----------
        a : np.ndarray
            First distribution
        b : np.ndarray
            Second distribution

        Returns
        -------
        float
            KS distance
        """
        from scipy import stats

        a_clean = a[np.isfinite(a)]
        b_clean = b[np.isfinite(b)]

        if a_clean.size < 2 or b_clean.size < 2:
            return float("nan")

        stat, _ = stats.ks_2samp(a_clean, b_clean)
        return float(stat)

    @staticmethod
    def earth_movers_distance(a: np.ndarray, b: np.ndarray) -> float:
        """Compute Earth Mover's Distance between distributions.

        Parameters
        ----------
        a : np.ndarray
            First distribution
        b : np.ndarray
            Second distribution

        Returns
        -------
        float
            EMD
        """
        from scipy import stats

        a_clean = a[np.isfinite(a)]
        b_clean = b[np.isfinite(b)]

        if a_clean.size < 2 or b_clean.size < 2:
            return float("nan")

        return float(stats.wasserstein_distance(a_clean, b_clean))
