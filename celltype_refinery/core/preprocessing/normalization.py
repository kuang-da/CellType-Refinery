"""Within-sample normalization (Stage C).

Provides background correction and variance stabilization transforms
for marker intensity normalization.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

from .config import NormalizationConfig


@dataclass
class TransformSpec:
    """Specification for a variance stabilization transform.

    Attributes
    ----------
    name : str
        Transform name (raw, log1p, asinh_c5, asinh_c10)
    label : str
        Human-readable label
    """

    name: str
    label: str


# Standard transforms
TRANSFORMS = {
    "raw": TransformSpec("raw", "raw"),
    "log1p": TransformSpec("log1p", "log1p"),
    "asinh_c5": TransformSpec("asinh_c5", "asinh(x/5)"),
    "asinh_c10": TransformSpec("asinh_c10", "asinh(x/10)"),
}


@dataclass
class NormalizationResult:
    """Result from normalizing a single sample.

    Attributes
    ----------
    sample_id : str
        Sample identifier
    normalized_matrix : pd.DataFrame
        Normalized intensity matrix
    quantiles : pd.DataFrame
        Quantile statistics per marker
    bg_params : Dict[str, float]
        Background parameters used
    """

    sample_id: str
    normalized_matrix: Optional[pd.DataFrame] = None
    quantiles: Optional[pd.DataFrame] = None
    bg_params: Dict[str, float] = field(default_factory=dict)


class Normalizer:
    """Within-sample intensity normalizer.

    Parameters
    ----------
    config : NormalizationConfig
        Normalization configuration

    Example
    -------
    >>> from celltype_refinery.core.preprocessing import Normalizer, NormalizationConfig
    >>> config = NormalizationConfig(default_bg_mode="clip", default_p=3)
    >>> normalizer = Normalizer(config)
    >>> result = normalizer.normalize_sample(cell_matrix, "sample_01")
    """

    def __init__(self, config: Optional[NormalizationConfig] = None):
        self.config = config or NormalizationConfig()

    def compute_quantiles(
        self,
        intensity_df: pd.DataFrame,
        percentiles: Optional[List[int]] = None,
    ) -> pd.DataFrame:
        """Compute quantile statistics per marker.

        Parameters
        ----------
        intensity_df : pd.DataFrame
            Marker intensity matrix (cells x markers)
        percentiles : List[int], optional
            Percentiles to compute (default: [1, 3, 5, 7, 10, 50, 90, 95, 99])

        Returns
        -------
        pd.DataFrame
            Quantiles indexed by marker
        """
        if percentiles is None:
            percentiles = [1, 3, 5, 7, 10, 50, 90, 95, 99]

        records = []
        for marker in intensity_df.columns:
            values = pd.to_numeric(intensity_df[marker], errors="coerce").dropna()
            if values.empty:
                record = {"marker": marker}
                for p in percentiles:
                    record[f"Q{p}"] = float("nan")
            else:
                record = {"marker": marker}
                for p in percentiles:
                    record[f"Q{p}"] = float(np.percentile(values, p))
            records.append(record)

        return pd.DataFrame(records).set_index("marker")

    def apply_background_correction(
        self,
        values: np.ndarray,
        mode: str,
        threshold: float,
    ) -> np.ndarray:
        """Apply background correction to values.

        Parameters
        ----------
        values : np.ndarray
            Raw intensity values
        mode : str
            Correction mode: 'clip' (set below threshold to 0) or 'center' (subtract threshold)
        threshold : float
            Background threshold

        Returns
        -------
        np.ndarray
            Background-corrected values
        """
        values = np.asarray(values, dtype=float)

        if mode == "clip":
            # Clip values below threshold to 0
            return np.where(values < threshold, 0.0, values - threshold)
        elif mode == "center":
            # Center by subtracting threshold (allows negative values)
            return values - threshold
        else:
            raise ValueError(f"Unknown background mode: {mode}")

    def apply_transform(
        self,
        values: np.ndarray,
        transform: str,
    ) -> np.ndarray:
        """Apply variance stabilization transform.

        Parameters
        ----------
        values : np.ndarray
            Background-corrected values
        transform : str
            Transform name: 'raw', 'log1p', 'asinh_c5', 'asinh_c10'

        Returns
        -------
        np.ndarray
            Transformed values
        """
        values = np.asarray(values, dtype=float)

        if transform == "raw":
            return values
        elif transform == "log1p":
            return np.log1p(np.maximum(values, 0))
        elif transform == "asinh_c5":
            return np.arcsinh(values / 5.0)
        elif transform == "asinh_c10":
            return np.arcsinh(values / 10.0)
        else:
            raise ValueError(f"Unknown transform: {transform}")

    def normalize_sample(
        self,
        cell_matrix: pd.DataFrame,
        sample_id: str,
        cell_id_col: str = "cell_mask_id",
        bg_mode: Optional[str] = None,
        bg_percentile: Optional[int] = None,
        transform: Optional[str] = None,
    ) -> NormalizationResult:
        """Normalize a single sample.

        Parameters
        ----------
        cell_matrix : pd.DataFrame
            Cell-by-marker intensity matrix
        sample_id : str
            Sample identifier
        cell_id_col : str
            Cell ID column name
        bg_mode : str, optional
            Background mode (default from config)
        bg_percentile : int, optional
            Background percentile (default from config)
        transform : str, optional
            Transform to apply (default from config)

        Returns
        -------
        NormalizationResult
            Normalization result with transformed matrix
        """
        result = NormalizationResult(sample_id=sample_id)

        # Use config defaults if not specified
        bg_mode = bg_mode or self.config.default_bg_mode
        bg_percentile = bg_percentile or self.config.default_p
        transform = transform or self.config.default_transform

        # Get marker columns
        marker_cols = [c for c in cell_matrix.columns if c != cell_id_col]
        intensity_df = cell_matrix[marker_cols].apply(pd.to_numeric, errors="coerce")

        # Compute quantiles
        result.quantiles = self.compute_quantiles(intensity_df)

        # Build normalized matrix
        normalized = pd.DataFrame(index=cell_matrix.index)
        normalized[cell_id_col] = cell_matrix[cell_id_col]

        for marker in marker_cols:
            values = intensity_df[marker].to_numpy(dtype=float)

            # Get background threshold
            q_label = f"Q{bg_percentile}"
            threshold = float(result.quantiles.loc[marker, q_label])

            if np.isfinite(threshold):
                # Apply background correction
                corrected = self.apply_background_correction(values, bg_mode, threshold)
                # Apply transform
                transformed = self.apply_transform(corrected, transform)
            else:
                transformed = values

            normalized[marker] = transformed
            result.bg_params[marker] = threshold

        result.normalized_matrix = normalized
        return result

    def normalize_matrix(
        self,
        intensity_df: pd.DataFrame,
        quantiles: pd.DataFrame,
        bg_mode: str,
        bg_percentile: int,
        transform: str,
    ) -> pd.DataFrame:
        """Normalize intensity matrix using pre-computed quantiles.

        Parameters
        ----------
        intensity_df : pd.DataFrame
            Marker intensity matrix
        quantiles : pd.DataFrame
            Pre-computed quantiles
        bg_mode : str
            Background correction mode
        bg_percentile : int
            Background percentile
        transform : str
            Transform to apply

        Returns
        -------
        pd.DataFrame
            Normalized matrix
        """
        normalized = pd.DataFrame(index=intensity_df.index)
        q_label = f"Q{bg_percentile}"

        for marker in intensity_df.columns:
            values = intensity_df[marker].to_numpy(dtype=float)

            if marker in quantiles.index:
                threshold = float(quantiles.loc[marker, q_label])
            else:
                threshold = float("nan")

            if np.isfinite(threshold):
                corrected = self.apply_background_correction(values, bg_mode, threshold)
                transformed = self.apply_transform(corrected, transform)
            else:
                transformed = values

            normalized[marker] = transformed

        return normalized

    @staticmethod
    def zero_fraction(values: np.ndarray) -> float:
        """Compute fraction of zero values.

        Parameters
        ----------
        values : np.ndarray
            Values to check

        Returns
        -------
        float
            Fraction of zeros
        """
        values = np.asarray(values, dtype=float)
        clean = values[np.isfinite(values)]
        if clean.size == 0:
            return float("nan")
        return float(np.sum(clean == 0) / clean.size)

    @staticmethod
    def negative_fraction(values: np.ndarray) -> float:
        """Compute fraction of negative values.

        Parameters
        ----------
        values : np.ndarray
            Values to check

        Returns
        -------
        float
            Fraction of negatives
        """
        values = np.asarray(values, dtype=float)
        clean = values[np.isfinite(values)]
        if clean.size == 0:
            return float("nan")
        return float(np.sum(clean < 0) / clean.size)

    def compute_cv(self, values: np.ndarray) -> float:
        """Compute coefficient of variation.

        Parameters
        ----------
        values : np.ndarray
            Values to analyze

        Returns
        -------
        float
            CV (std / mean)
        """
        values = np.asarray(values, dtype=float)
        clean = values[np.isfinite(values)]
        if clean.size == 0:
            return float("nan")
        mean = np.mean(clean)
        if np.isclose(mean, 0):
            return float("nan")
        return float(np.std(clean) / abs(mean))
