"""Proportion-based flagging rules.

Rules:
- PROPORTION_OUTLIER: Cell type outside expected global proportion range
- REGIONAL_ABSENCE: Expected cell type absent from region
- REGIONAL_EXCESS: Cell type exceeds expected regional proportion
- SAMPLE_ARTIFACT: Cell type concentrated in single sample
- DIVERSITY_ANOMALY: Single cell type dominates composition
"""

from typing import List, TYPE_CHECKING

from .base import BaseRule, FlaggedIssue
from .registry import RuleRegistry

if TYPE_CHECKING:
    from ..aggregator import MetricsAggregator


@RuleRegistry.register
class ProportionOutlierRule(BaseRule):
    """Flag cell types outside expected proportion range."""

    rule_id = "PROPORTION_OUTLIER"
    category = "proportion"
    default_severity = "warning"
    description = "Cell type outside expected proportion range"

    def check(self, aggregator: "MetricsAggregator") -> List[FlaggedIssue]:
        issues = []

        low_factor = self.get_threshold("low_factor", 0.5)
        high_factor = self.get_threshold("high_factor", 1.5)

        for cell_type in aggregator.get_cell_types():
            pct = aggregator.get_global_pct(cell_type)
            if pct is None:
                continue

            expected = self.template.get_expected_pct_range(cell_type)
            if expected is None:
                continue

            exp_min, exp_max = expected

            # Check if below expected range
            if pct < exp_min * low_factor:
                issues.append(
                    self.create_issue(
                        cell_type=cell_type,
                        description=f"Below expected proportion ({pct:.2f}% vs {exp_min:.1f}-{exp_max:.1f}%)",
                        evidence=f"Global: {pct:.2f}%, Expected min: {exp_min:.1f}%",
                        recommendation="Check marker map sensitivity or annotation thresholds",
                        value=pct,
                        expected=f"{exp_min:.1f}-{exp_max:.1f}%",
                    )
                )

            # Check if above expected range
            elif pct > exp_max * high_factor:
                issues.append(
                    self.create_issue(
                        cell_type=cell_type,
                        description=f"Above expected proportion ({pct:.2f}% vs {exp_min:.1f}-{exp_max:.1f}%)",
                        evidence=f"Global: {pct:.2f}%, Expected max: {exp_max:.1f}%",
                        recommendation="Check for over-annotation or marker specificity issues",
                        value=pct,
                        expected=f"{exp_min:.1f}-{exp_max:.1f}%",
                    )
                )

        return issues


@RuleRegistry.register
class RegionalAbsenceRule(BaseRule):
    """Flag cell types absent from expected regions."""

    rule_id = "REGIONAL_ABSENCE"
    category = "proportion"
    default_severity = "warning"
    description = "Expected cell type absent from region"

    def check(self, aggregator: "MetricsAggregator") -> List[FlaggedIssue]:
        issues = []

        min_pct = self.get_threshold("min_pct", 0.1)

        for cell_type in aggregator.get_cell_types():
            ct_info = self.template.get_cell_type(cell_type)
            if ct_info is None or not ct_info.expected_regions:
                continue

            for region, expected_range in ct_info.expected_regions.items():
                regional_pct = aggregator.get_regional_pct(cell_type, region)

                # Skip if region not in data
                if regional_pct is None:
                    continue

                # Check if below minimum
                if regional_pct < min_pct:
                    issues.append(
                        self.create_issue(
                            cell_type=cell_type,
                            description=f"Nearly absent from {region} ({regional_pct:.3f}%)",
                            evidence=f"Regional: {regional_pct:.3f}%, Expected: {expected_range[0]:.1f}-{expected_range[1]:.1f}%",
                            recommendation=f"Verify {cell_type} annotation in {region} samples",
                            region=region,
                            value=regional_pct,
                            expected=f"{expected_range[0]:.1f}-{expected_range[1]:.1f}%",
                        )
                    )

        return issues


@RuleRegistry.register
class RegionalExcessRule(BaseRule):
    """Flag cell types exceeding expected regional proportion."""

    rule_id = "REGIONAL_EXCESS"
    category = "proportion"
    default_severity = "warning"
    description = "Cell type exceeds expected regional proportion"

    def check(self, aggregator: "MetricsAggregator") -> List[FlaggedIssue]:
        issues = []

        excess_factor = self.get_threshold("excess_factor", 2.0)

        for cell_type in aggregator.get_cell_types():
            ct_info = self.template.get_cell_type(cell_type)
            if ct_info is None or not ct_info.expected_regions:
                continue

            for region, expected_range in ct_info.expected_regions.items():
                regional_pct = aggregator.get_regional_pct(cell_type, region)

                # Skip if region not in data
                if regional_pct is None:
                    continue

                exp_max = expected_range[1]

                # Check if exceeds expected max by factor
                if regional_pct > exp_max * excess_factor:
                    issues.append(
                        self.create_issue(
                            cell_type=cell_type,
                            description=f"Excess in {region} ({regional_pct:.1f}% vs max {exp_max:.1f}%)",
                            evidence=f"Regional: {regional_pct:.1f}%, Expected max: {exp_max:.1f}%",
                            recommendation=f"Check for regional batch effect or annotation bias in {region}",
                            region=region,
                            value=regional_pct,
                            expected=f"<{exp_max * excess_factor:.1f}%",
                        )
                    )

        return issues


@RuleRegistry.register
class SampleArtifactRule(BaseRule):
    """Flag cell types concentrated in single sample."""

    rule_id = "SAMPLE_ARTIFACT"
    category = "proportion"
    default_severity = "critical"
    description = "Cell type concentrated in single sample"

    def check(self, aggregator: "MetricsAggregator") -> List[FlaggedIssue]:
        issues = []

        max_concentration = self.get_threshold("max_concentration", 0.50)
        min_samples = self.get_threshold("min_samples", 2)

        for cell_type in aggregator.get_cell_types():
            sample_dist = aggregator.get_sample_distribution(cell_type)
            if sample_dist is None:
                continue

            n_samples = sample_dist.get("n_samples", 0)
            max_sample_pct = sample_dist.get("max_sample_pct", 0)

            # Skip types in few samples (may be rare but valid)
            if n_samples < min_samples:
                continue

            # Check if concentrated in single sample
            if max_sample_pct > max_concentration:
                issues.append(
                    self.create_issue(
                        cell_type=cell_type,
                        description=f"Concentrated in single sample ({max_sample_pct:.1%} of type)",
                        evidence=f"Max sample: {max_sample_pct:.1%}, Present in {n_samples} samples",
                        recommendation="Investigate for sample-specific artifact or batch effect",
                        value=max_sample_pct * 100,
                        expected=f"<{max_concentration * 100:.0f}%",
                        severity="critical",
                    )
                )

        return issues


@RuleRegistry.register
class DiversityAnomalyRule(BaseRule):
    """Flag when single cell type dominates composition."""

    rule_id = "DIVERSITY_ANOMALY"
    category = "proportion"
    default_severity = "warning"
    description = "Single cell type dominates composition"

    def check(self, aggregator: "MetricsAggregator") -> List[FlaggedIssue]:
        issues = []

        max_single_pct = self.get_threshold("max_single_type_pct", 35.0)

        for cell_type in aggregator.get_cell_types():
            pct = aggregator.get_global_pct(cell_type)
            if pct is None:
                continue

            # Skip Unassigned - handled by separate rule
            if cell_type.lower() == "unassigned":
                continue

            if pct > max_single_pct:
                issues.append(
                    self.create_issue(
                        cell_type=cell_type,
                        description=f"Dominates composition ({pct:.1f}% > {max_single_pct:.0f}%)",
                        evidence=f"Global: {pct:.1f}%",
                        recommendation="Verify this is biologically expected for tissue type",
                        value=pct,
                        expected=f"<{max_single_pct:.0f}%",
                    )
                )

        return issues
