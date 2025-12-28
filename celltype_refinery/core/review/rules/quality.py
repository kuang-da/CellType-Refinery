"""Annotation quality flagging rules.

Rules:
- HIGH_UNASSIGNED: High proportion of unassigned cells
- HIGH_HYBRID: High proportion of hybrid cell types
- LOW_CONFIDENCE: Many cells with low confidence assignments
- ORPHAN_QUALITY: Suspicious orphan rescues
- KNOWN_ARTIFACT: Matches known problematic pattern
"""

from typing import List, TYPE_CHECKING

from .base import BaseRule, FlaggedIssue
from .registry import RuleRegistry

if TYPE_CHECKING:
    from ..aggregator import MetricsAggregator


@RuleRegistry.register
class HighUnassignedRule(BaseRule):
    """Flag high proportion of unassigned cells."""

    rule_id = "HIGH_UNASSIGNED"
    category = "quality"
    default_severity = "critical"
    description = "High proportion of unassigned cells"

    def check(self, aggregator: "MetricsAggregator") -> List[FlaggedIssue]:
        issues = []

        max_pct = self.get_threshold(
            "max_pct",
            self.template.get_quality_threshold("max_unassigned_pct", 10.0),
        )
        critical_pct = self.get_threshold("critical_pct", 15.0)

        unassigned_pct = aggregator.get_global_pct("Unassigned")
        if unassigned_pct is None:
            return issues

        if unassigned_pct > critical_pct:
            issues.append(
                self.create_issue(
                    cell_type="Unassigned",
                    description=f"Critical: {unassigned_pct:.1f}% unassigned (target: <{max_pct:.0f}%)",
                    evidence=f"Unassigned: {unassigned_pct:.1f}%, Critical threshold: {critical_pct:.0f}%",
                    recommendation="Review marker map and annotation thresholds; consider orphan rescue",
                    severity="critical",
                    value=unassigned_pct,
                    expected=f"<{max_pct:.0f}%",
                )
            )
        elif unassigned_pct > max_pct:
            issues.append(
                self.create_issue(
                    cell_type="Unassigned",
                    description=f"High unassigned: {unassigned_pct:.1f}% (target: <{max_pct:.0f}%)",
                    evidence=f"Unassigned: {unassigned_pct:.1f}%, Threshold: {max_pct:.0f}%",
                    recommendation="Review marker map and annotation thresholds",
                    severity="warning",
                    value=unassigned_pct,
                    expected=f"<{max_pct:.0f}%",
                )
            )

        return issues


@RuleRegistry.register
class HighHybridRule(BaseRule):
    """Flag high proportion of hybrid cell types."""

    rule_id = "HIGH_HYBRID"
    category = "quality"
    default_severity = "warning"
    description = "High proportion of hybrid cell types"

    def check(self, aggregator: "MetricsAggregator") -> List[FlaggedIssue]:
        issues = []

        max_pct = self.get_threshold(
            "max_pct",
            self.template.get_quality_threshold("max_hybrid_pct", 5.0),
        )
        hybrid_pattern = self.get_threshold("hybrid_pattern", "~")

        # Calculate total hybrid percentage
        total_hybrid_pct = 0
        hybrid_types = []

        for cell_type in aggregator.get_cell_types():
            if hybrid_pattern in cell_type:
                pct = aggregator.get_global_pct(cell_type)
                if pct is not None:
                    total_hybrid_pct += pct
                    hybrid_types.append((cell_type, pct))

        if total_hybrid_pct > max_pct:
            # Sort hybrids by percentage
            hybrid_types.sort(key=lambda x: x[1], reverse=True)
            top_hybrids = [f"{ct}: {p:.1f}%" for ct, p in hybrid_types[:3]]

            issues.append(
                self.create_issue(
                    cell_type="Hybrid types",
                    description=f"High hybrid rate: {total_hybrid_pct:.1f}% (target: <{max_pct:.0f}%)",
                    evidence=f"Total hybrid: {total_hybrid_pct:.1f}%, Top: {', '.join(top_hybrids)}",
                    recommendation="Review marker exclusivity; hybrids may indicate doublets or marker bleed",
                    value=total_hybrid_pct,
                    expected=f"<{max_pct:.0f}%",
                )
            )

        return issues


@RuleRegistry.register
class LowConfidenceRule(BaseRule):
    """Flag many cells with low confidence assignments."""

    rule_id = "LOW_CONFIDENCE"
    category = "quality"
    default_severity = "warning"
    description = "Many cells with low confidence assignments"

    def check(self, aggregator: "MetricsAggregator") -> List[FlaggedIssue]:
        issues = []

        max_low_conf_pct = self.get_threshold(
            "max_low_confidence_pct",
            self.template.get_quality_threshold("max_low_confidence_pct", 20.0),
        )

        confidence_data = aggregator.get_confidence_breakdown()
        if confidence_data is None:
            return issues

        # Calculate low confidence percentage
        total_cells = confidence_data.get("total", 0)
        low_conf_cells = confidence_data.get("low", 0) + confidence_data.get(
            "very_low", 0
        )

        if total_cells > 0:
            low_conf_pct = (low_conf_cells / total_cells) * 100

            if low_conf_pct > max_low_conf_pct:
                issues.append(
                    self.create_issue(
                        cell_type="All types",
                        description=f"High low-confidence rate: {low_conf_pct:.1f}% (target: <{max_low_conf_pct:.0f}%)",
                        evidence=f"Low/very-low confidence: {low_conf_cells:,} cells ({low_conf_pct:.1f}%)",
                        recommendation="Review annotation thresholds or marker map sensitivity",
                        value=low_conf_pct,
                        expected=f"<{max_low_conf_pct:.0f}%",
                    )
                )

        return issues


@RuleRegistry.register
class OrphanQualityRule(BaseRule):
    """Flag suspicious orphan rescues."""

    rule_id = "ORPHAN_QUALITY"
    category = "quality"
    default_severity = "warning"
    description = "Suspicious orphan rescues"

    def check(self, aggregator: "MetricsAggregator") -> List[FlaggedIssue]:
        issues = []

        suspicious_cells_threshold = self.get_threshold("suspicious_cells", 10000)

        orphan_data = aggregator.get_orphan_candidates()
        if orphan_data is None:
            return issues

        # Count suspicious orphans
        suspicious = orphan_data[orphan_data["flag"] == "suspicious_orphan"]

        if len(suspicious) > 0:
            total_suspicious_cells = suspicious["n_cells"].sum()

            if total_suspicious_cells > suspicious_cells_threshold:
                # Get top suspicious subtypes
                top_subtypes = (
                    suspicious.groupby("subtype")["n_cells"]
                    .sum()
                    .sort_values(ascending=False)
                )
                top_list = [
                    f"{st}: {n:,}" for st, n in top_subtypes.head(3).items()
                ]

                issues.append(
                    self.create_issue(
                        cell_type="Orphans",
                        description=f"High suspicious orphan count: {total_suspicious_cells:,} cells",
                        evidence=f"Suspicious orphans: {len(suspicious)} clusters, {total_suspicious_cells:,} cells. Top: {', '.join(top_list)}",
                        recommendation="Review orphan rescue plausibility; may indicate marker map issues",
                        value=total_suspicious_cells,
                        expected=f"<{suspicious_cells_threshold:,}",
                    )
                )

        return issues


@RuleRegistry.register
class KnownArtifactRule(BaseRule):
    """Flag cell types matching known artifacts."""

    rule_id = "KNOWN_ARTIFACT"
    category = "quality"
    default_severity = "warning"
    description = "Matches known problematic pattern"

    def check(self, aggregator: "MetricsAggregator") -> List[FlaggedIssue]:
        issues = []

        known_artifacts = self.template.get_known_artifacts()

        for artifact in known_artifacts:
            name = artifact.get("name", "")
            detection_method = artifact.get("detection_method", "")
            threshold = artifact.get("threshold", 0.5)

            # Check for sample concentration artifact
            if detection_method == "sample_concentration":
                for cell_type in aggregator.get_cell_types():
                    dist = aggregator.get_sample_distribution(cell_type)
                    if dist and dist.get("max_sample_pct", 0) > threshold:
                        issues.append(
                            self.create_issue(
                                cell_type=cell_type,
                                description=f"Matches known artifact: {name}",
                                evidence=artifact.get("description", ""),
                                recommendation=artifact.get("recommendation", ""),
                            )
                        )

            # Check for extreme clustering artifact
            elif detection_method == "extreme_clustering":
                for cell_type in aggregator.get_cell_types():
                    moran_i = aggregator.get_moran_i(cell_type)
                    if moran_i and moran_i > threshold:
                        issues.append(
                            self.create_issue(
                                cell_type=cell_type,
                                description=f"Matches known artifact: {name}",
                                evidence=f"Moran's I: {moran_i:.3f}. {artifact.get('description', '')}",
                                recommendation=artifact.get("recommendation", ""),
                            )
                        )

            # Check for low assignment artifact
            elif detection_method == "low_assignment":
                unassigned_pct = aggregator.get_global_pct("Unassigned")
                if unassigned_pct and unassigned_pct > threshold * 100:
                    issues.append(
                        self.create_issue(
                            cell_type="Unassigned",
                            description=f"Matches known artifact: {name}",
                            evidence=artifact.get("description", ""),
                            recommendation=artifact.get("recommendation", ""),
                        )
                    )

        return issues
