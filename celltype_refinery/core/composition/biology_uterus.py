"""Uterus-specific biology metrics.

This module provides uterus-specific metrics:
- Epithelial:Stromal ratio
- Glandular:Stromal ratio (endometrial health indicator)
- Smooth muscle content (myometrium)
- Vascular and lymphatic endothelium percentages
- Immune cell infiltration with macrophage breakdown

Regional expectations:
- Ectocervix: Higher epithelial (squamous)
- Endocervix: Higher glandular epithelium
- Body/Fundus: Balanced glandular + myometrium
- Lower Segment: Higher smooth muscle

Usage:
    from celltype_refinery.core.composition.biology_uterus import UterusMetrics

    engine = CompositionEngine(biology_metrics=UterusMetrics())
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


class UterusMetrics(TissueBiologyMetrics):
    """Uterus-specific biology metrics.

    Computes:
    - epithelial_pct: Percentage of epithelial cells
    - stromal_pct: Percentage of stromal cells
    - smooth_muscle_pct: Percentage of smooth muscle cells (myometrium)
    - immune_pct: Immune cell infiltration percentage
    - endothelial_pct: Percentage of all endothelial cells
    - glandular_pct: Percentage of glandular epithelium
    - luminal_pct: Percentage of luminal epithelium
    - vascular_pct: Percentage of vascular endothelium
    - lymphatic_pct: Percentage of lymphatic endothelium
    - mesenchymal_pct: Percentage of mesenchymal cells (including MSCs)
    - macrophage_pct: Percentage of macrophages
    - epithelial_stromal_ratio: Epithelial to stromal ratio
    - glandular_stromal_ratio: Glandular to stromal ratio
    - vascular_lymphatic_ratio: Vascular to lymphatic endothelium ratio
    """

    @property
    def tissue_name(self) -> str:
        return "uterus"

    @property
    def metric_names(self) -> List[str]:
        return [
            "epithelial_pct",
            "stromal_pct",
            "smooth_muscle_pct",
            "immune_pct",
            "endothelial_pct",
            "glandular_pct",
            "luminal_pct",
            "vascular_pct",
            "lymphatic_pct",
            "mesenchymal_pct",
            "macrophage_pct",
            "epithelial_stromal_ratio",
            "glandular_stromal_ratio",
            "vascular_lymphatic_ratio",
        ]

    def compute_metrics(
        self,
        df: pd.DataFrame,
        cell_type_col: str,
        patterns: PatternConfig,
    ) -> List[BiologyMetric]:
        """Compute uterus-specific biology metrics."""
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

        # Count by category using pattern matching (standard fields)
        epithelial_count = count_by_pattern(cell_types, patterns.epithelial)
        stromal_count = count_by_pattern(cell_types, patterns.stromal)
        smooth_muscle_count = count_by_pattern(cell_types, patterns.smooth_muscle)
        immune_count = count_by_pattern(cell_types, patterns.immune)
        endothelial_count = count_by_pattern(cell_types, patterns.endothelial)

        # Uterus-specific patterns (from custom dict)
        glandular_count = count_by_pattern(cell_types, patterns.custom.get("glandular", []))
        luminal_count = count_by_pattern(cell_types, patterns.custom.get("luminal", []))
        vascular_count = count_by_pattern(cell_types, patterns.custom.get("vascular", []))
        lymphatic_count = count_by_pattern(cell_types, patterns.custom.get("lymphatic", []))
        mesenchymal_count = count_by_pattern(cell_types, patterns.custom.get("mesenchymal", []))
        macrophage_count = count_by_pattern(cell_types, patterns.custom.get("macrophage", []))

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
                description="Percentage of smooth muscle cells (myometrium)",
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
                description="Percentage of all endothelial cells",
                details={"count": endothelial_count, "total": total},
            ),
            BiologyMetric(
                name="glandular_pct",
                value=compute_percentage(glandular_count, total),
                category="percentage",
                description="Percentage of glandular epithelium",
                details={"count": glandular_count, "total": total},
            ),
            BiologyMetric(
                name="luminal_pct",
                value=compute_percentage(luminal_count, total),
                category="percentage",
                description="Percentage of luminal epithelium",
                details={"count": luminal_count, "total": total},
            ),
            BiologyMetric(
                name="vascular_pct",
                value=compute_percentage(vascular_count, total),
                category="percentage",
                description="Percentage of vascular endothelium",
                details={"count": vascular_count, "total": total},
            ),
            BiologyMetric(
                name="lymphatic_pct",
                value=compute_percentage(lymphatic_count, total),
                category="percentage",
                description="Percentage of lymphatic endothelium",
                details={"count": lymphatic_count, "total": total},
            ),
            BiologyMetric(
                name="mesenchymal_pct",
                value=compute_percentage(mesenchymal_count, total),
                category="percentage",
                description="Percentage of mesenchymal cells (including MSCs)",
                details={"count": mesenchymal_count, "total": total},
            ),
            BiologyMetric(
                name="macrophage_pct",
                value=compute_percentage(macrophage_count, total),
                category="percentage",
                description="Percentage of macrophages",
                details={"count": macrophage_count, "total": total},
            ),
            BiologyMetric(
                name="epithelial_stromal_ratio",
                value=compute_ratio(epithelial_count, stromal_count),
                category="ratio",
                description="Epithelial to stromal cell ratio",
                details={"epithelial": epithelial_count, "stromal": stromal_count},
            ),
            BiologyMetric(
                name="glandular_stromal_ratio",
                value=compute_ratio(glandular_count, stromal_count),
                category="ratio",
                description="Glandular epithelium to stromal ratio (endometrial health)",
                details={"glandular": glandular_count, "stromal": stromal_count},
            ),
            BiologyMetric(
                name="vascular_lymphatic_ratio",
                value=compute_ratio(vascular_count, lymphatic_count),
                category="ratio",
                description="Vascular to lymphatic endothelium ratio",
                details={"vascular": vascular_count, "lymphatic": lymphatic_count},
            ),
        ]

        return metrics


def validate_uterus_expectations(
    biology_by_region: pd.DataFrame,
    region_col: str = "region",
) -> pd.DataFrame:
    """Validate uterus biology metrics against expected patterns.

    Expected uterus biology patterns:
    - Ectocervix: Higher epithelial percentage (squamous epithelium)
    - Endocervix: Higher glandular epithelium
    - Body/Fundus: Balanced glandular + myometrium
    - Lower Segment: Higher smooth muscle

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

        # Check ectocervix expectations (should have higher epithelial)
        if "ectocervix" in region_lower:
            epi_pct = row.get("epithelial_pct", np.nan)
            if not np.isnan(epi_pct) and epi_pct < 10.0:
                flags.append({
                    "region": region,
                    "metric": "epithelial_pct",
                    "value": epi_pct,
                    "expected": ">10%",
                    "flag": "LOW_EPITHELIAL_IN_ECTOCERVIX",
                    "severity": "warning",
                })

        # Check endocervix expectations (should have glandular epithelium)
        if "endocervix" in region_lower:
            gland_pct = row.get("glandular_pct", np.nan)
            if not np.isnan(gland_pct) and gland_pct < 1.0:
                flags.append({
                    "region": region,
                    "metric": "glandular_pct",
                    "value": gland_pct,
                    "expected": ">1%",
                    "flag": "LOW_GLANDULAR_IN_ENDOCERVIX",
                    "severity": "note",
                })

        # Check body/fundus expectations
        if "body" in region_lower or "fundus" in region_lower:
            smooth_pct = row.get("smooth_muscle_pct", np.nan)
            # Body should have some myometrium
            if not np.isnan(smooth_pct) and smooth_pct < 0.1:
                flags.append({
                    "region": region,
                    "metric": "smooth_muscle_pct",
                    "value": smooth_pct,
                    "expected": ">0.1%",
                    "flag": "LOW_SMOOTH_MUSCLE_IN_BODY",
                    "severity": "note",
                })

        # Check E:S ratio (tissue architecture)
        es_ratio = row.get("epithelial_stromal_ratio", np.nan)
        if not np.isnan(es_ratio):
            if es_ratio < 0.1:
                flags.append({
                    "region": region,
                    "metric": "epithelial_stromal_ratio",
                    "value": es_ratio,
                    "expected": ">0.1",
                    "flag": "VERY_LOW_EPITHELIAL_STROMAL_RATIO",
                    "severity": "warning",
                })
            elif es_ratio > 20.0:
                flags.append({
                    "region": region,
                    "metric": "epithelial_stromal_ratio",
                    "value": es_ratio,
                    "expected": "<20",
                    "flag": "VERY_HIGH_EPITHELIAL_STROMAL_RATIO",
                    "severity": "warning",
                })

    return pd.DataFrame(flags)


# =============================================================================
# Auto-register UterusMetrics on module import
# =============================================================================
register_organ_metrics(
    organ="uterus",
    metrics_class=UterusMetrics,
    aliases=["ut", "uterine"],
)
