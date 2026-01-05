"""Abstract biology metrics for tissue-specific analysis.

This module provides a base class for tissue-specific biology metrics.
Each tissue type can implement its own metrics by subclassing
TissueBiologyMetrics.

Example tissue-specific implementations:
- Fallopian tube: Epithelial:Stromal ratio, Ciliated:Secretory ratio
- Colon: Goblet cell percentage, Crypt:Villus ratio
- Lung: Alveolar type I:II ratio, Immune infiltration
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd

from .config import PatternConfig


@dataclass
class BiologyMetric:
    """A computed biology metric.

    Attributes
    ----------
    name : str
        Metric name
    value : float
        Computed value
    category : str
        Metric category (e.g., "ratio", "percentage", "index")
    description : str
        Human-readable description
    is_valid : bool
        Whether the metric could be computed (enough cells)
    details : Dict[str, Any]
        Additional details about the computation
    """

    name: str
    value: float
    category: str
    description: str
    is_valid: bool = True
    details: Dict[str, Any] = field(default_factory=dict)


@dataclass
class BiologyResult:
    """Result of biology metrics computation.

    Attributes
    ----------
    metrics : List[BiologyMetric]
        List of computed metrics
    group_id : str
        Group identifier (sample_id, region, etc.)
    n_cells : int
        Total cells analyzed
    """

    metrics: List[BiologyMetric]
    group_id: str
    n_cells: int

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        result = {
            "group_id": self.group_id,
            "n_cells": self.n_cells,
        }
        for m in self.metrics:
            result[m.name] = m.value
            result[f"{m.name}_valid"] = m.is_valid
        return result


def classify_cell_type(
    label: str,
    patterns: PatternConfig,
) -> Optional[str]:
    """Classify a cell type label into a category.

    Parameters
    ----------
    label : str
        Cell type label
    patterns : PatternConfig
        Pattern configuration

    Returns
    -------
    str or None
        Category name or None if no match
    """
    label_lower = label.lower()

    for category, category_patterns in patterns.get_all_categories().items():
        for pattern in category_patterns:
            if pattern.lower() in label_lower:
                return category

    return None


def count_by_pattern(
    cell_types: pd.Series,
    patterns: List[str],
) -> int:
    """Count cells matching any of the patterns.

    Parameters
    ----------
    cell_types : pd.Series
        Series of cell type labels
    patterns : List[str]
        List of patterns to match

    Returns
    -------
    int
        Number of cells matching any pattern
    """
    if not patterns:
        return 0

    mask = pd.Series(False, index=cell_types.index)
    for pattern in patterns:
        mask |= cell_types.str.contains(pattern, case=False, na=False)

    return int(mask.sum())


def compute_ratio(
    numerator: int,
    denominator: int,
    default: float = np.nan,
) -> float:
    """Compute ratio with division-by-zero handling.

    Parameters
    ----------
    numerator : int
        Numerator count
    denominator : int
        Denominator count
    default : float
        Value to return if denominator is zero

    Returns
    -------
    float
        Ratio or default value
    """
    if denominator == 0:
        return default
    return numerator / denominator


def compute_percentage(
    count: int,
    total: int,
    default: float = np.nan,
) -> float:
    """Compute percentage with handling for zero total.

    Parameters
    ----------
    count : int
        Count of interest
    total : int
        Total count
    default : float
        Value to return if total is zero

    Returns
    -------
    float
        Percentage or default value
    """
    if total == 0:
        return default
    return 100.0 * count / total


class TissueBiologyMetrics(ABC):
    """Abstract base class for tissue-specific biology metrics.

    Subclass this to implement tissue-specific metrics.

    Example
    -------
    >>> class FallopianTubeMetrics(TissueBiologyMetrics):
    ...     def compute_metrics(self, df, cell_type_col, patterns):
    ...         # Compute FT-specific metrics
    ...         pass
    """

    @property
    @abstractmethod
    def tissue_name(self) -> str:
        """Return tissue name."""
        pass

    @property
    @abstractmethod
    def metric_names(self) -> List[str]:
        """Return list of metric names computed by this class."""
        pass

    @abstractmethod
    def compute_metrics(
        self,
        df: pd.DataFrame,
        cell_type_col: str,
        patterns: PatternConfig,
    ) -> List[BiologyMetric]:
        """Compute tissue-specific biology metrics.

        Parameters
        ----------
        df : pd.DataFrame
            DataFrame with cell data
        cell_type_col : str
            Column with cell type labels
        patterns : PatternConfig
            Pattern configuration

        Returns
        -------
        List[BiologyMetric]
            List of computed metrics
        """
        pass

    def compute_for_group(
        self,
        df: pd.DataFrame,
        group_id: str,
        cell_type_col: str,
        patterns: PatternConfig,
    ) -> BiologyResult:
        """Compute metrics for a single group.

        Parameters
        ----------
        df : pd.DataFrame
            DataFrame with cell data for one group
        group_id : str
            Group identifier
        cell_type_col : str
            Column with cell type labels
        patterns : PatternConfig
            Pattern configuration

        Returns
        -------
        BiologyResult
            Computed metrics for this group
        """
        metrics = self.compute_metrics(df, cell_type_col, patterns)
        return BiologyResult(
            metrics=metrics,
            group_id=group_id,
            n_cells=len(df),
        )

    def compute_by_group(
        self,
        df: pd.DataFrame,
        group_col: str,
        cell_type_col: str,
        patterns: PatternConfig,
    ) -> pd.DataFrame:
        """Compute metrics for each group.

        Parameters
        ----------
        df : pd.DataFrame
            DataFrame with cell data
        group_col : str
            Column to group by
        cell_type_col : str
            Column with cell type labels
        patterns : PatternConfig
            Pattern configuration

        Returns
        -------
        pd.DataFrame
            DataFrame with metrics for each group
        """
        results = []
        for group_id, group_df in df.groupby(group_col, observed=True):
            result = self.compute_for_group(
                group_df,
                str(group_id),
                cell_type_col,
                patterns,
            )
            results.append(result.to_dict())

        return pd.DataFrame(results)


class GenericBiologyMetrics(TissueBiologyMetrics):
    """Generic biology metrics applicable to most tissues.

    Computes common metrics like:
    - Epithelial percentage
    - Stromal percentage
    - Immune infiltration
    - Endothelial percentage
    - Epithelial:Stromal ratio
    """

    @property
    def tissue_name(self) -> str:
        return "generic"

    @property
    def metric_names(self) -> List[str]:
        return [
            "epithelial_pct",
            "stromal_pct",
            "immune_pct",
            "endothelial_pct",
            "epithelial_stromal_ratio",
        ]

    def compute_metrics(
        self,
        df: pd.DataFrame,
        cell_type_col: str,
        patterns: PatternConfig,
    ) -> List[BiologyMetric]:
        """Compute generic biology metrics."""
        if len(df) == 0:
            return [
                BiologyMetric(name=m, value=np.nan, category="percentage",
                             description="", is_valid=False)
                for m in self.metric_names
            ]

        cell_types = df[cell_type_col]
        total = len(df)

        # Count by category
        epithelial_count = count_by_pattern(cell_types, patterns.epithelial)
        stromal_count = count_by_pattern(cell_types, patterns.stromal)
        immune_count = count_by_pattern(cell_types, patterns.immune)
        endothelial_count = count_by_pattern(cell_types, patterns.endothelial)

        # Compute metrics
        metrics = [
            BiologyMetric(
                name="epithelial_pct",
                value=compute_percentage(epithelial_count, total),
                category="percentage",
                description="Percentage of epithelial cells",
                details={"count": epithelial_count, "total": total},
            ),
            BiologyMetric(
                name="stromal_pct",
                value=compute_percentage(stromal_count, total),
                category="percentage",
                description="Percentage of stromal cells",
                details={"count": stromal_count, "total": total},
            ),
            BiologyMetric(
                name="immune_pct",
                value=compute_percentage(immune_count, total),
                category="percentage",
                description="Percentage of immune cells (immune infiltration)",
                details={"count": immune_count, "total": total},
            ),
            BiologyMetric(
                name="endothelial_pct",
                value=compute_percentage(endothelial_count, total),
                category="percentage",
                description="Percentage of endothelial cells",
                details={"count": endothelial_count, "total": total},
            ),
            BiologyMetric(
                name="epithelial_stromal_ratio",
                value=compute_ratio(epithelial_count, stromal_count),
                category="ratio",
                description="Epithelial to stromal cell ratio",
                details={"epithelial": epithelial_count, "stromal": stromal_count},
            ),
        ]

        return metrics


def get_default_biology_metrics() -> TissueBiologyMetrics:
    """Get default biology metrics implementation.

    Returns
    -------
    TissueBiologyMetrics
        GenericBiologyMetrics instance
    """
    return GenericBiologyMetrics()


# =============================================================================
# Organ Registry
# =============================================================================
# Maps organ names to their TissueBiologyMetrics implementations.
# Use get_organ_metrics() to get the appropriate class for an organ.

ORGAN_METRICS_REGISTRY: Dict[str, type] = {
    # Default/generic
    "generic": GenericBiologyMetrics,
    # Add organ-specific implementations here as they are created
    # "fallopian_tube": FallopianTubeMetrics,  # Added via register_organ_metrics()
}

# Aliases for common organ names
ORGAN_ALIASES: Dict[str, str] = {
    # Fallopian tube aliases
    "ft": "fallopian_tube",
    "fallopian": "fallopian_tube",
    "oviduct": "fallopian_tube",
    # Add more aliases as needed
}


def register_organ_metrics(
    organ: str,
    metrics_class: type,
    aliases: Optional[List[str]] = None,
) -> None:
    """Register an organ-specific metrics class.

    Parameters
    ----------
    organ : str
        Canonical organ name (lowercase, underscores)
    metrics_class : type
        TissueBiologyMetrics subclass
    aliases : List[str], optional
        Alternative names for this organ

    Example
    -------
    >>> from celltype_refinery.core.composition.biology import register_organ_metrics
    >>> from celltype_refinery.core.composition.biology_ft import FallopianTubeMetrics
    >>> register_organ_metrics("fallopian_tube", FallopianTubeMetrics, aliases=["ft"])
    """
    organ_lower = organ.lower().replace(" ", "_").replace("-", "_")
    ORGAN_METRICS_REGISTRY[organ_lower] = metrics_class

    if aliases:
        for alias in aliases:
            ORGAN_ALIASES[alias.lower()] = organ_lower


def get_organ_metrics(organ: Optional[str] = None) -> TissueBiologyMetrics:
    """Get biology metrics implementation for an organ.

    Parameters
    ----------
    organ : str, optional
        Organ name (case-insensitive). If None or "generic", returns GenericBiologyMetrics.

    Returns
    -------
    TissueBiologyMetrics
        Organ-specific or generic metrics instance

    Raises
    ------
    ValueError
        If organ is not registered

    Example
    -------
    >>> metrics = get_organ_metrics("fallopian_tube")
    >>> print(metrics.tissue_name)
    'fallopian_tube'
    """
    if organ is None:
        return GenericBiologyMetrics()

    # Normalize organ name
    organ_lower = organ.lower().replace(" ", "_").replace("-", "_")

    # Check aliases
    if organ_lower in ORGAN_ALIASES:
        organ_lower = ORGAN_ALIASES[organ_lower]

    # Look up in registry
    if organ_lower not in ORGAN_METRICS_REGISTRY:
        available = list(ORGAN_METRICS_REGISTRY.keys())
        available_aliases = [f"{k} -> {v}" for k, v in ORGAN_ALIASES.items()]
        raise ValueError(
            f"Unknown organ: '{organ}'. "
            f"Available organs: {available}. "
            f"Aliases: {available_aliases}"
        )

    metrics_class = ORGAN_METRICS_REGISTRY[organ_lower]
    return metrics_class()


def list_available_organs() -> List[str]:
    """List all available organ names.

    Returns
    -------
    List[str]
        List of registered organ names and their aliases
    """
    organs = list(ORGAN_METRICS_REGISTRY.keys())
    return organs


def compute_biology_metrics(
    df: pd.DataFrame,
    group_col: str,
    cell_type_col: str,
    patterns: PatternConfig,
    metrics_impl: Optional[TissueBiologyMetrics] = None,
) -> pd.DataFrame:
    """Convenience function to compute biology metrics.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with cell data
    group_col : str
        Column to group by
    cell_type_col : str
        Column with cell type labels
    patterns : PatternConfig
        Pattern configuration
    metrics_impl : TissueBiologyMetrics, optional
        Custom metrics implementation. Uses GenericBiologyMetrics if None.

    Returns
    -------
    pd.DataFrame
        DataFrame with metrics for each group
    """
    if metrics_impl is None:
        metrics_impl = get_default_biology_metrics()

    return metrics_impl.compute_by_group(df, group_col, cell_type_col, patterns)
