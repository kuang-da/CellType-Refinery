"""Fallopian Tube-specific biology metrics.

This module provides FT-specific metrics:
- Epithelial:Stromal ratio
- Ciliated:Secretory ratio (FT-specific)
- Immune cell infiltration percentage
- Smooth muscle content
- Endothelial content

Usage:
    from celltype_refinery.core.composition.biology_ft import FallopianTubeMetrics

    engine = CompositionEngine(biology_metrics=FallopianTubeMetrics())
    result = engine.execute(adata, cell_type_col="cell_type_phenocycler")
"""

from typing import List

import numpy as np
import pandas as pd

from .biology import (
    BiologyMetric,
    TissueBiologyMetrics,
    count_by_pattern,
    compute_ratio,
    compute_percentage,
    register_organ_metrics,
)
from .config import PatternConfig


class FallopianTubeMetrics(TissueBiologyMetrics):
    """Fallopian Tube-specific biology metrics.

    Computes:
    - epithelial_pct: Percentage of epithelial cells
    - stromal_pct: Percentage of stromal cells
    - smooth_muscle_pct: Percentage of smooth muscle cells
    - immune_pct: Immune cell infiltration percentage
    - endothelial_pct: Percentage of endothelial cells
    - ciliated_pct: Percentage of ciliated epithelial cells
    - secretory_pct: Percentage of secretory epithelial cells
    - epithelial_stromal_ratio: Epithelial to stromal ratio
    - ciliated_secretory_ratio: Ciliated to secretory ratio (FT-specific)

    Expected Patterns:
    - Fimbriae: Higher ciliated:secretory ratio (more ciliated cells)
    - Isthmus: Higher smooth muscle content
    - Epithelial:Stromal ratio typically 1.5-3.0 in healthy FT
    """

    @property
    def tissue_name(self) -> str:
        return "fallopian_tube"

    @property
    def metric_names(self) -> List[str]:
        return [
            "epithelial_pct",
            "stromal_pct",
            "smooth_muscle_pct",
            "immune_pct",
            "endothelial_pct",
            "ciliated_pct",
            "secretory_pct",
            "epithelial_stromal_ratio",
            "ciliated_secretory_ratio",
        ]

    def compute_metrics(
        self,
        df: pd.DataFrame,
        cell_type_col: str,
        patterns: PatternConfig,
    ) -> List[BiologyMetric]:
        """Compute FT-specific biology metrics."""
        if len(df) == 0:
            return [
                BiologyMetric(
                    name=m, value=np.nan, category="percentage",
                    description="", is_valid=False
                )
                for m in self.metric_names
            ]

        cell_types = df[cell_type_col]
        total = len(df)

        # Count by category using pattern matching
        epithelial_count = count_by_pattern(cell_types, patterns.epithelial)
        stromal_count = count_by_pattern(cell_types, patterns.stromal)
        smooth_muscle_count = count_by_pattern(cell_types, patterns.smooth_muscle)
        immune_count = count_by_pattern(cell_types, patterns.immune)
        endothelial_count = count_by_pattern(cell_types, patterns.endothelial)
        ciliated_count = count_by_pattern(cell_types, patterns.ciliated)
        secretory_count = count_by_pattern(cell_types, patterns.secretory)

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
                name="smooth_muscle_pct",
                value=compute_percentage(smooth_muscle_count, total),
                category="percentage",
                description="Percentage of smooth muscle cells",
                details={"count": smooth_muscle_count, "total": total},
            ),
            BiologyMetric(
                name="immune_pct",
                value=compute_percentage(immune_count, total),
                category="percentage",
                description="Immune cell infiltration percentage",
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
                name="ciliated_pct",
                value=compute_percentage(ciliated_count, total),
                category="percentage",
                description="Percentage of ciliated epithelial cells",
                details={"count": ciliated_count, "total": total},
            ),
            BiologyMetric(
                name="secretory_pct",
                value=compute_percentage(secretory_count, total),
                category="percentage",
                description="Percentage of secretory epithelial cells",
                details={"count": secretory_count, "total": total},
            ),
            BiologyMetric(
                name="epithelial_stromal_ratio",
                value=compute_ratio(epithelial_count, stromal_count),
                category="ratio",
                description="Epithelial to stromal cell ratio",
                details={"epithelial": epithelial_count, "stromal": stromal_count},
            ),
            BiologyMetric(
                name="ciliated_secretory_ratio",
                value=compute_ratio(ciliated_count, secretory_count),
                category="ratio",
                description="Ciliated to secretory cell ratio (FT-specific)",
                details={"ciliated": ciliated_count, "secretory": secretory_count},
            ),
        ]

        return metrics


def validate_ft_expectations(
    biology_by_region: pd.DataFrame,
    region_col: str = "region",
) -> pd.DataFrame:
    """Validate FT biology metrics against expected patterns.

    Expected FT biology patterns:
    - Fimbriae: Higher ciliated:secretory ratio (more ciliated cells)
    - Isthmus: Higher smooth muscle content
    - Epithelial:Stromal ratio typically 1.5-3.0 in healthy FT

    Parameters
    ----------
    biology_by_region : pd.DataFrame
        Biology metrics with region column
    region_col : str
        Column name for region

    Returns
    -------
    pd.DataFrame
        Validation flags with warnings for unexpected patterns
    """
    if region_col not in biology_by_region.columns:
        return pd.DataFrame()

    flags = []

    for _, row in biology_by_region.iterrows():
        region = row[region_col]
        region_lower = str(region).lower()

        # Check fimbriae expectations
        if "fimbri" in region_lower:
            cil_sec = row.get("ciliated_secretory_ratio", np.nan)
            if not np.isnan(cil_sec) and cil_sec < 1.0:
                flags.append({
                    "region": region,
                    "metric": "ciliated_secretory_ratio",
                    "value": cil_sec,
                    "expected": ">1.0",
                    "flag": "LOW_CILIATED_IN_FIMBRIAE",
                    "severity": "warning",
                })

        # Check isthmus expectations
        if "isthm" in region_lower:
            smooth = row.get("smooth_muscle_pct", np.nan)
            if not np.isnan(smooth) and smooth < 5.0:
                flags.append({
                    "region": region,
                    "metric": "smooth_muscle_pct",
                    "value": smooth,
                    "expected": ">5%",
                    "flag": "LOW_SMOOTH_MUSCLE_IN_ISTHMUS",
                    "severity": "warning",
                })

        # Check E:S ratio
        es_ratio = row.get("epithelial_stromal_ratio", np.nan)
        if not np.isnan(es_ratio):
            if es_ratio < 0.5:
                flags.append({
                    "region": region,
                    "metric": "epithelial_stromal_ratio",
                    "value": es_ratio,
                    "expected": "0.5-5.0",
                    "flag": "LOW_EPITHELIAL_STROMAL_RATIO",
                    "severity": "warning",
                })
            elif es_ratio > 10.0:
                flags.append({
                    "region": region,
                    "metric": "epithelial_stromal_ratio",
                    "value": es_ratio,
                    "expected": "0.5-5.0",
                    "flag": "HIGH_EPITHELIAL_STROMAL_RATIO",
                    "severity": "warning",
                })

    return pd.DataFrame(flags)


# =============================================================================
# Auto-register FallopianTubeMetrics on module import
# =============================================================================
register_organ_metrics(
    organ="fallopian_tube",
    metrics_class=FallopianTubeMetrics,
    aliases=["ft", "fallopian", "oviduct"],
)
