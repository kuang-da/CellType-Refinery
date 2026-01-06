"""Recommendation logic for cluster review actions.

Ported from ft/src/stage_h_coarse_clustering.py lines 1429-1469.
"""

from __future__ import annotations

import pandas as pd


def determine_recommendation(row: pd.Series) -> str:
    """Determine recommended action for a cluster.

    Args:
        row: Row from cluster_assignments DataFrame.

    Returns:
        Recommendation string: "keep", "subcluster", "merge", "manual_override", or "review".
    """
    score = row.get("assigned_score", 0)
    n_cells = row.get("n_cells", 0)
    stop_reason = row.get("stop_reason", "")

    # Check for low hierarchical margin (close sibling competition)
    # margin_is_infinite = True means no competition at any level (highest confidence)
    min_margin = row.get("min_margin_along_path", 1e6)
    margin_is_inf = row.get("margin_is_infinite", min_margin >= 1e6)
    has_low_margin = not margin_is_inf and min_margin < 0.3

    # High confidence, leaf reached - keep as is
    if score >= 2.0 and stop_reason == "leaf_reached":
        return "keep"

    # Large ambiguous cluster - subcluster
    if n_cells > 5000 and (stop_reason in ["ambiguous_siblings", "no_child_passed"] or has_low_margin):
        return "subcluster"

    # Small cluster with low score - consider merge
    if n_cells < 1000 and score < 1.0:
        return "merge"

    # Unassigned - needs manual override
    if stop_reason == "no_root_passed" or row.get("assigned_label", "") == "Unassigned":
        return "manual_override"

    # Medium confidence - just review
    if score < 1.5:
        return "review"

    return "keep"
