"""Priority scoring functions for review checklist.

Ported from ft/src/stage_h_coarse_clustering.py lines 1355-1426.
"""

from __future__ import annotations

from typing import Optional

import pandas as pd

from .config import ChecklistConfig


def classify_confidence_band(score: float) -> str:
    """Classify a confidence score into a band.

    Args:
        score: Assignment confidence score.

    Returns:
        Confidence band: "high", "medium", "low", or "very_low".
    """
    if score >= 2.0:
        return "high"
    elif score >= 1.0:
        return "medium"
    elif score >= 0.5:
        return "low"
    else:
        return "very_low"


def compute_review_priority(
    row: pd.Series,
    config: Optional[ChecklistConfig] = None,
) -> int:
    """Compute review priority score for a cluster (higher = higher priority).

    Args:
        row: Row from cluster_assignments DataFrame.
        config: Checklist configuration with priority weights.

    Returns:
        Priority score (0-200 range, higher = review first).
    """
    if config is None:
        config = ChecklistConfig()

    priority = 0

    # Low confidence clusters are high priority
    score = row.get("assigned_score", 0)
    if score < 1.0:
        priority += config.low_conf_weight
    elif score < 1.5:
        priority += config.medium_low_conf_weight

    # Large clusters are higher priority (more impact if wrong)
    n_cells = row.get("n_cells", 0)
    if n_cells > 50000:
        priority += config.large_cluster_weight
    elif n_cells > 10000:
        priority += config.medium_cluster_weight
    elif n_cells > 5000:
        priority += config.small_cluster_weight

    # Ambiguous siblings need review
    stop_reason = row.get("stop_reason", "")
    if stop_reason == "ambiguous_siblings":
        priority += config.ambiguous_siblings_weight

    # Low hierarchical margin (close sibling competition) needs attention
    # Skip if margin_is_infinite (no competition = highest confidence)
    min_margin = row.get("min_margin_along_path", 1e6)
    margin_is_inf = row.get("margin_is_infinite", min_margin >= 1e6)
    if not margin_is_inf and min_margin < 0.3:
        priority += config.low_margin_weight

    # Low coverage may indicate marker issues
    coverage = row.get("coverage", 1.0)
    if coverage < 0.75:
        priority += config.low_coverage_weight
    elif coverage < 0.9:
        priority += config.medium_coverage_weight

    # Needs refinement flag
    if row.get("stopped_before_leaf", False):
        priority += config.stopped_before_leaf_weight

    return priority
