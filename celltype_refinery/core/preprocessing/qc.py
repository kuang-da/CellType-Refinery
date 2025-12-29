"""Cell-level quality control (Stage B).

Provides percentile-based cell filtering for removing
outlier cells based on morphology and intensity metrics.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from .config import QCConfig


# Reason columns for tracking removal causes
REASON_COLUMNS = [
    "low_area",
    "high_area",
    "low_nucleus_ratio",
    "high_nucleus_ratio",
    "low_intensity",
    "high_autofluorescence",
]


@dataclass
class QCResult:
    """Result from QC filtering a single sample.

    Attributes
    ----------
    sample_id : str
        Sample identifier
    cells_total : int
        Total cells before filtering
    cells_removed : int
        Number of cells removed
    removal_fraction : float
        Fraction of cells removed
    capped_by_max : bool
        Whether removal was capped by max_removal_fraction
    reason_counts : Dict[str, int]
        Counts per removal reason
    filtered_matrix : pd.DataFrame
        Filtered cell matrix
    filtered_metadata : pd.DataFrame
        Filtered cell metadata
    removal_records : List[Dict]
        Details of removed cells
    """

    sample_id: str
    cells_total: int = 0
    cells_removed: int = 0
    removal_fraction: float = 0.0
    capped_by_max: bool = False
    reason_counts: Dict[str, int] = field(default_factory=dict)
    filtered_matrix: Optional[pd.DataFrame] = None
    filtered_metadata: Optional[pd.DataFrame] = None
    removal_records: List[Dict[str, Any]] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for reporting."""
        result = {
            "sample_id": self.sample_id,
            "cells_total": self.cells_total,
            "cells_removed": self.cells_removed,
            "removal_fraction": round(self.removal_fraction, 4),
            "capped_by_max": self.capped_by_max,
        }
        for reason in REASON_COLUMNS:
            result[f"removed_{reason}"] = self.reason_counts.get(reason, 0)
        return result


class CellQC:
    """Cell-level quality control filter.

    Parameters
    ----------
    config : QCConfig
        QC configuration

    Example
    -------
    >>> from celltype_refinery.core.preprocessing import CellQC, QCConfig
    >>> config = QCConfig(max_removal_fraction=0.15)
    >>> qc = CellQC(config)
    >>> result = qc.filter_sample(cell_matrix, cell_metadata, "sample_01")
    """

    def __init__(self, config: Optional[QCConfig] = None):
        self.config = config or QCConfig()

    def _compute_percentile_bounds(
        self, series: pd.Series, lower: float, upper: float
    ) -> Tuple[float, float]:
        """Compute percentile bounds for a series.

        Parameters
        ----------
        series : pd.Series
            Numeric series
        lower : float
            Lower percentile (0-100)
        upper : float
            Upper percentile (0-100)

        Returns
        -------
        Tuple[float, float]
            (lower_bound, upper_bound)
        """
        clean = pd.to_numeric(series, errors="coerce").dropna()
        if clean.empty:
            return float("nan"), float("nan")
        return (
            float(np.percentile(clean, lower)),
            float(np.percentile(clean, upper)),
        )

    def _safe_ratio(
        self, numerator: pd.Series, denominator: pd.Series
    ) -> pd.Series:
        """Compute safe ratio handling division by zero.

        Parameters
        ----------
        numerator : pd.Series
            Numerator values
        denominator : pd.Series
            Denominator values

        Returns
        -------
        pd.Series
            Ratio with NaN for invalid entries
        """
        num = pd.to_numeric(numerator, errors="coerce")
        den = pd.to_numeric(denominator, errors="coerce")
        ratio = pd.Series(np.nan, index=num.index, dtype=float)
        mask = den > 0
        ratio.loc[mask] = num.loc[mask] / den.loc[mask]
        return ratio

    def filter_sample(
        self,
        cell_matrix: pd.DataFrame,
        cell_metadata: pd.DataFrame,
        sample_id: str,
        cell_id_col: str = "cell_mask_id",
    ) -> QCResult:
        """Filter cells based on QC criteria.

        Parameters
        ----------
        cell_matrix : pd.DataFrame
            Cell-by-marker intensity matrix
        cell_metadata : pd.DataFrame
            Cell metadata with morphology columns
        sample_id : str
            Sample identifier
        cell_id_col : str
            Cell ID column name

        Returns
        -------
        QCResult
            Filtering result with cleaned data
        """
        result = QCResult(sample_id=sample_id)

        # Align on common cell IDs
        common_ids = sorted(
            set(cell_matrix[cell_id_col]) & set(cell_metadata[cell_id_col])
        )
        cell_matrix = cell_matrix[cell_matrix[cell_id_col].isin(common_ids)].copy()
        cell_metadata = cell_metadata[cell_metadata[cell_id_col].isin(common_ids)].copy()

        cell_matrix.set_index(cell_id_col, inplace=True)
        cell_metadata.set_index(cell_id_col, inplace=True)

        # Get marker columns and compute intensities (exclude non-marker columns)
        exclude_cols = set(self.config.exclude_from_markers)
        marker_cols = [c for c in cell_matrix.columns if c not in exclude_cols]
        intensity_matrix = cell_matrix[marker_cols].astype(float)
        total_intensity = intensity_matrix.sum(axis=1)
        median_intensity = intensity_matrix.median(axis=1)

        # Get area columns
        area_col = self.config.area_col
        nucleus_area_col = self.config.nucleus_area_col

        # Compute area bounds
        cell_area_series = pd.to_numeric(
            cell_metadata.get(area_col), errors="coerce"
        )
        area_low, area_high = self._compute_percentile_bounds(
            cell_area_series,
            self.config.area_percentile_low,
            self.config.area_percentile_high,
        )

        # Compute nucleus ratio
        ratio_series = self._safe_ratio(
            cell_metadata.get(nucleus_area_col),
            cell_metadata.get(area_col),
        )

        # Build reason flags
        reasons = pd.DataFrame(index=intensity_matrix.index)

        if not np.isnan(area_low):
            reasons["low_area"] = cell_area_series < area_low
        else:
            reasons["low_area"] = False

        if not np.isnan(area_high):
            reasons["high_area"] = cell_area_series > area_high
        else:
            reasons["high_area"] = False

        reasons["low_nucleus_ratio"] = ratio_series < self.config.nucleus_ratio_min
        reasons["high_nucleus_ratio"] = ratio_series > self.config.nucleus_ratio_max

        # Intensity-based filtering
        if not total_intensity.empty:
            intensity_low = float(
                np.percentile(total_intensity, self.config.intensity_percentile_low)
            )
            reasons["low_intensity"] = total_intensity <= intensity_low
        else:
            reasons["low_intensity"] = False

        if not median_intensity.empty:
            autofluorescence_high = float(
                np.percentile(
                    median_intensity, self.config.autofluorescence_percentile_high
                )
            )
            reasons["high_autofluorescence"] = median_intensity >= autofluorescence_high
        else:
            reasons["high_autofluorescence"] = False

        reasons = reasons.fillna(False)

        # Flag cells with any issue
        flagged = reasons.any(axis=1)
        severity = reasons.sum(axis=1)
        reasons["severity"] = severity
        reasons["total_intensity"] = total_intensity

        # Compute removal counts
        result.cells_total = len(reasons)
        flagged_count = int(flagged.sum())
        max_remove = int(result.cells_total * self.config.max_removal_fraction)
        remove_count = min(flagged_count, max_remove)
        result.capped_by_max = flagged_count > remove_count

        # Select cells to remove (worst first)
        removal_ids: List[str] = []
        if remove_count > 0:
            removal_ids = (
                reasons.loc[flagged]
                .sort_values(by=["severity", "total_intensity"], ascending=[False, True])
                .head(remove_count)
                .index.tolist()
            )

        # Record removal details
        if removal_ids:
            reasons_for_removed = reasons.loc[removal_ids]
            for cell_id, row_values in reasons_for_removed.iterrows():
                cell_reasons = [
                    name for name in REASON_COLUMNS if bool(row_values.get(name))
                ]
                if not cell_reasons:
                    cell_reasons = ["other"]
                result.removal_records.append(
                    {
                        "sample_id": sample_id,
                        "cell_id": cell_id,
                        "reasons": ";".join(sorted(cell_reasons)),
                    }
                )
                for reason in cell_reasons:
                    result.reason_counts[reason] = (
                        result.reason_counts.get(reason, 0) + 1
                    )

        # Create filtered dataframes
        result.filtered_matrix = cell_matrix.drop(index=removal_ids).reset_index()
        result.filtered_metadata = cell_metadata.drop(index=removal_ids).reset_index()

        result.cells_removed = len(removal_ids)
        result.removal_fraction = (
            result.cells_removed / result.cells_total
            if result.cells_total > 0
            else 0.0
        )

        return result

    def get_area_distribution(
        self, cell_metadata: pd.DataFrame
    ) -> Optional[pd.Series]:
        """Get cell area distribution for diagnostics.

        Parameters
        ----------
        cell_metadata : pd.DataFrame
            Cell metadata

        Returns
        -------
        pd.Series or None
            Area values
        """
        area_col = self.config.area_col
        if area_col not in cell_metadata.columns:
            return None
        return pd.to_numeric(cell_metadata[area_col], errors="coerce")
