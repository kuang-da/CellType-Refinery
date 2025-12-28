"""Spatial-based flagging rules.

Rules:
- EXTREME_CLUSTERING: Unusually high Moran's I (potential artifact)
- UNEXPECTED_DISPERSION: Expected clustered type is dispersed
- SPATIAL_ISOLATION: Cell type spatially isolated from all others
"""

from typing import List, TYPE_CHECKING

from .base import BaseRule, FlaggedIssue
from .registry import RuleRegistry

if TYPE_CHECKING:
    from ..aggregator import MetricsAggregator


@RuleRegistry.register
class ExtremeClusteringRule(BaseRule):
    """Flag cell types with unusually high Moran's I."""

    rule_id = "EXTREME_CLUSTERING"
    category = "spatial"
    default_severity = "warning"
    description = "Unusually high spatial clustering (potential artifact)"

    def check(self, aggregator: "MetricsAggregator") -> List[FlaggedIssue]:
        issues = []

        moran_threshold = self.get_threshold("moran_threshold", 0.75)

        morans_data = aggregator.get_morans_i_data()
        if morans_data is None:
            return issues

        for cell_type in aggregator.get_cell_types():
            moran_i = aggregator.get_moran_i(cell_type)
            if moran_i is None:
                continue

            # Check for extreme clustering
            if moran_i > moran_threshold:
                # Check if this is expected for the cell type
                expected = self.template.get_expected_moran_range(cell_type)
                if expected and moran_i <= expected[1] * 1.2:
                    # Within expected range (with 20% buffer)
                    continue

                issues.append(
                    self.create_issue(
                        cell_type=cell_type,
                        description=f"Extreme spatial clustering (Moran's I = {moran_i:.3f})",
                        evidence=f"Moran's I: {moran_i:.3f}, Threshold: {moran_threshold:.2f}",
                        recommendation="Check for segmentation artifacts or autofluorescent debris",
                        value=moran_i,
                        expected=f"<{moran_threshold:.2f}",
                    )
                )

        return issues


@RuleRegistry.register
class UnexpectedDispersionRule(BaseRule):
    """Flag expected clustered types that are dispersed."""

    rule_id = "UNEXPECTED_DISPERSION"
    category = "spatial"
    default_severity = "note"
    description = "Expected clustered type is dispersed"

    def check(self, aggregator: "MetricsAggregator") -> List[FlaggedIssue]:
        issues = []

        expected_min_moran = self.get_threshold("expected_min_moran", 0.30)
        actual_threshold = self.get_threshold("actual_threshold", 0.15)

        for cell_type in aggregator.get_cell_types():
            moran_i = aggregator.get_moran_i(cell_type)
            if moran_i is None:
                continue

            # Check expected Moran's I from template
            expected = self.template.get_expected_moran_range(cell_type)
            if expected is None:
                continue

            exp_min, exp_max = expected

            # Flag if expected to cluster but actually dispersed
            if exp_min >= expected_min_moran and moran_i < actual_threshold:
                issues.append(
                    self.create_issue(
                        cell_type=cell_type,
                        description=f"Unexpectedly dispersed (Moran's I = {moran_i:.3f} vs expected {exp_min:.2f}-{exp_max:.2f})",
                        evidence=f"Moran's I: {moran_i:.3f}, Expected: {exp_min:.2f}-{exp_max:.2f}",
                        recommendation="Check if marker specificity is low or cells are misannotated",
                        value=moran_i,
                        expected=f"{exp_min:.2f}-{exp_max:.2f}",
                    )
                )

        return issues


@RuleRegistry.register
class SpatialIsolationRule(BaseRule):
    """Flag cell types spatially isolated from all others."""

    rule_id = "SPATIAL_ISOLATION"
    category = "spatial"
    default_severity = "note"
    description = "Cell type spatially isolated from all others"

    def check(self, aggregator: "MetricsAggregator") -> List[FlaggedIssue]:
        issues = []

        min_pct_for_flag = self.get_threshold("min_pct_for_flag", 1.0)

        enrichment_data = aggregator.get_enrichment_pairs()
        if enrichment_data is None:
            return issues

        for cell_type in aggregator.get_cell_types():
            pct = aggregator.get_global_pct(cell_type)
            if pct is None or pct < min_pct_for_flag:
                continue

            # Get enrichments involving this cell type
            enrichments = aggregator.get_enrichments_for_type(cell_type)
            if enrichments is None:
                continue

            # Check if enriched with any other type
            n_positive_enrichments = sum(1 for z in enrichments.values() if z > 2.0)

            if n_positive_enrichments == 0:
                # Check if only self-enriched
                self_enrichment = enrichments.get(cell_type)
                if self_enrichment and self_enrichment > 0:
                    issues.append(
                        self.create_issue(
                            cell_type=cell_type,
                            description="Spatially isolated (no enrichment with other cell types)",
                            evidence=f"Self-enrichment: {self_enrichment:.1f}, No positive enrichments with others",
                            recommendation="Check if this type is in a distinct spatial compartment",
                            value=n_positive_enrichments,
                            expected=">0 enrichments",
                        )
                    )

        return issues
