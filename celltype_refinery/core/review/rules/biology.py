"""Tissue-specific biology flagging rules.

These rules are driven by configuration in the tissue template YAML,
allowing tissue-specific biology rules without code changes.

Rules:
- RATIO_VIOLATION: Configurable ratio metrics outside expected range
- REGIONAL_METRIC_LOW: Regional metric below expected threshold
- REGIONAL_GRADIENT_ABSENT: Expected regional variation not observed
"""

from typing import Any, Dict, List, TYPE_CHECKING

import numpy as np

from .base import BaseRule, FlaggedIssue
from .registry import RuleRegistry

if TYPE_CHECKING:
    from ..aggregator import MetricsAggregator


@RuleRegistry.register
class RatioViolationRule(BaseRule):
    """Flag tissue-specific ratio metrics outside expected range.

    Configured via template.biology_rules with format:
    ```yaml
    biology_rules:
      rules:
        - rule_id: "RATIO_VIOLATION"
          metric_name: "epithelial_stromal_ratio"
          description: "E:S ratio check"
          regional_thresholds:
            fimbriae: [0.5, 2.5]
            isthmus: [0.05, 0.3]
    ```
    """

    rule_id = "RATIO_VIOLATION"
    category = "biology"
    default_severity = "warning"
    description = "Ratio metric outside expected range"

    def check(self, aggregator: "MetricsAggregator") -> List[FlaggedIssue]:
        issues = []

        biology_config = self.template.get_biology_rules()
        if not biology_config.enabled:
            return issues

        # Find RATIO_VIOLATION rules in config
        for rule_def in biology_config.rules:
            if rule_def.get("rule_id") != "RATIO_VIOLATION":
                continue

            metric_name = rule_def.get("metric_name", "")
            rule_description = rule_def.get("description", "Ratio check")
            regional_thresholds = rule_def.get("regional_thresholds", {})

            biology_data = aggregator.get_biology_by_region()
            if biology_data is None:
                continue

            for _, row in biology_data.iterrows():
                region = row.get("region", "").lower()
                metric_value = row.get(metric_name)

                if region not in regional_thresholds or metric_value is None:
                    continue

                expected = regional_thresholds[region]
                if isinstance(expected, list) and len(expected) == 2:
                    exp_min, exp_max = expected

                    if metric_value < exp_min or metric_value > exp_max:
                        issues.append(
                            self.create_issue(
                                cell_type=metric_name.replace("_", " ").title(),
                                description=f"{rule_description} in {region} ({metric_value:.2f}) outside expected range",
                                evidence=f"{metric_name}: {metric_value:.2f}, Expected: {exp_min:.2f}-{exp_max:.2f}",
                                recommendation=f"Review annotation in {region}",
                                region=region,
                                value=metric_value,
                                expected=f"{exp_min:.2f}-{exp_max:.2f}",
                            )
                        )

        return issues


@RuleRegistry.register
class RegionalMetricLowRule(BaseRule):
    """Flag regional metrics below expected threshold.

    Configured via template.biology_rules:
    ```yaml
    biology_rules:
      rules:
        - rule_id: "REGIONAL_METRIC_LOW"
          region: "isthmus"
          cell_types: ["Smooth Muscle Cells"]
          min_pct: 35.0
          description: "Isthmus smooth muscle check"
    ```
    """

    rule_id = "REGIONAL_METRIC_LOW"
    category = "biology"
    default_severity = "warning"
    description = "Regional metric below expected threshold"

    def check(self, aggregator: "MetricsAggregator") -> List[FlaggedIssue]:
        issues = []

        biology_config = self.template.get_biology_rules()
        if not biology_config.enabled:
            return issues

        # Find REGIONAL_METRIC_LOW rules in config
        for rule_def in biology_config.rules:
            if rule_def.get("rule_id") != "REGIONAL_METRIC_LOW":
                continue

            region = rule_def.get("region", "")
            cell_types = rule_def.get("cell_types", [])
            min_pct = rule_def.get("min_pct", 0.0)
            rule_description = rule_def.get("description", "Regional metric check")

            # Sum percentages for configured cell types
            total_pct = 0.0
            for ct in cell_types:
                pct = aggregator.get_regional_pct(ct, region)
                if pct is not None:
                    total_pct += pct

            if total_pct < min_pct:
                cell_type_str = ", ".join(cell_types) if len(cell_types) > 1 else cell_types[0]
                issues.append(
                    self.create_issue(
                        cell_type=cell_type_str,
                        description=f"Low {cell_type_str} in {region} ({total_pct:.1f}% < {min_pct:.0f}%)",
                        evidence=f"{region} {cell_type_str}: {total_pct:.1f}%, Expected: >{min_pct:.0f}%",
                        recommendation=rule_description,
                        region=region,
                        value=total_pct,
                        expected=f">{min_pct:.0f}%",
                    )
                )

        return issues


@RuleRegistry.register
class RegionalGradientRule(BaseRule):
    """Flag cell types with no regional variation when expected.

    Configured via template.biology_rules:
    ```yaml
    biology_rules:
      rules:
        - rule_id: "REGIONAL_GRADIENT_ABSENT"
          cell_types: ["Smooth Muscle Cells", "Epithelium"]
          min_cv: 0.3
          description: "Regional gradient check"
    ```
    """

    rule_id = "REGIONAL_GRADIENT_ABSENT"
    category = "biology"
    default_severity = "note"
    description = "Expected regional variation not observed"

    def check(self, aggregator: "MetricsAggregator") -> List[FlaggedIssue]:
        issues = []

        biology_config = self.template.get_biology_rules()
        if not biology_config.enabled:
            return issues

        # Find REGIONAL_GRADIENT_ABSENT rules in config
        for rule_def in biology_config.rules:
            if rule_def.get("rule_id") != "REGIONAL_GRADIENT_ABSENT":
                continue

            cell_types = rule_def.get("cell_types", [])
            min_cv = rule_def.get("min_cv", 0.3)
            rule_description = rule_def.get(
                "description", "This cell type should vary across regions"
            )

            regional_data = aggregator.get_regional_composition()
            if regional_data is None:
                continue

            for cell_type in cell_types:
                # Calculate CV across regions
                regional_pcts = []
                for region in aggregator.get_regions():
                    pct = aggregator.get_regional_pct(cell_type, region)
                    if pct is not None:
                        regional_pcts.append(pct)

                if len(regional_pcts) < 2:
                    continue

                mean_pct = np.mean(regional_pcts)
                std_pct = np.std(regional_pcts)

                if mean_pct > 0:
                    cv = std_pct / mean_pct
                else:
                    continue

                if cv < min_cv:
                    issues.append(
                        self.create_issue(
                            cell_type=cell_type,
                            description=f"No regional gradient (CV = {cv:.2f} < {min_cv:.2f})",
                            evidence=f"Regional CV: {cv:.2f}, Mean: {mean_pct:.1f}%, Std: {std_pct:.2f}%",
                            recommendation=rule_description,
                            value=cv,
                            expected=f"CV > {min_cv:.2f}",
                        )
                    )

        return issues


def create_tissue_specific_biology_rules(
    template_config: Dict[str, Any]
) -> List[Dict[str, Any]]:
    """Helper to create biology rule configurations for a tissue.

    This function can be used to generate YAML configuration for
    tissue-specific biology rules.

    Parameters
    ----------
    template_config : Dict
        Template configuration with expected biology characteristics

    Returns
    -------
    List[Dict]
        List of biology rule definitions for YAML

    Example
    -------
    >>> config = {
    ...     "tissue": "fallopian_tube",
    ...     "biology": {
    ...         "ratios": {
    ...             "epithelial_stromal": {
    ...                 "fimbriae": [0.5, 2.5],
    ...                 "isthmus": [0.05, 0.3],
    ...             }
    ...         },
    ...         "regional_minimums": {
    ...             "isthmus": {"Smooth Muscle Cells": 35.0},
    ...             "fimbriae": {"Epithelium": 12.0},
    ...         },
    ...         "gradient_types": ["Smooth Muscle Cells", "Epithelium"],
    ...     }
    ... }
    >>> rules = create_tissue_specific_biology_rules(config)
    """
    rules = []

    biology = template_config.get("biology", {})

    # Create ratio rules
    for ratio_name, regional_thresholds in biology.get("ratios", {}).items():
        rules.append({
            "rule_id": "RATIO_VIOLATION",
            "metric_name": ratio_name + "_ratio",
            "description": f"{ratio_name.replace('_', ' ').title()} ratio check",
            "regional_thresholds": regional_thresholds,
        })

    # Create regional minimum rules
    for region, type_mins in biology.get("regional_minimums", {}).items():
        for cell_type, min_pct in type_mins.items():
            rules.append({
                "rule_id": "REGIONAL_METRIC_LOW",
                "region": region,
                "cell_types": [cell_type],
                "min_pct": min_pct,
                "description": f"{region} should have high {cell_type}; verify annotation",
            })

    # Create gradient rules
    gradient_types = biology.get("gradient_types", [])
    if gradient_types:
        rules.append({
            "rule_id": "REGIONAL_GRADIENT_ABSENT",
            "cell_types": gradient_types,
            "min_cv": biology.get("gradient_min_cv", 0.3),
            "description": "This cell type should vary across regions",
        })

    return rules
