"""Orphaned subtype detection and rescue.

An "orphaned subtype" is a cell that:
- Failed the root gate (root_score < threshold)
- Shows strong subtype marker enrichment (subtype_score >= threshold)

This can happen biologically when:
- Lymphatic endothelium has low CD31/VWF but high Podoplanin
- Mature smooth muscle cells lose CD44 but express SMA/ACTA2
- Quiescent fibroblasts have low CD44 but express Caveolin/CD34
"""

from dataclasses import dataclass
from enum import Enum
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd

from .config import OrphanRescueConfig, OrphanRule


class OrphanPlausibility(Enum):
    """Biological plausibility of orphan rescue."""

    HIGH = "high"  # Biologically well-documented
    SUSPICIOUS = "suspicious"  # Unusual, needs expert review
    UNCERTAIN = "uncertain"  # Possible but not well-documented


class OrphanAction(Enum):
    """Action to take for orphan."""

    RESCUE = "rescue"  # Assign subtype label with suffix
    KEEP_UNASSIGNED = "keep_unassigned"  # Keep as Unassigned
    EXPERT_REVIEW = "expert_review"  # Flag for expert review


@dataclass
class OrphanCandidate:
    """A detected orphan candidate.

    Attributes
    ----------
    cluster_id : str
        The cluster ID
    n_cells : int
        Number of cells in cluster
    subtype : str
        The detected subtype
    subtype_score : float
        Score for the subtype
    parent_root : str
        The expected parent root type
    root_score : float
        Score for the parent root
    plausibility : OrphanPlausibility
        Biological plausibility
    action : OrphanAction
        Recommended action
    final_label : Optional[str]
        Final label if rescued (None if not rescued)
    flag : str
        Optional flag for tracking
    """

    cluster_id: str
    n_cells: int
    subtype: str
    subtype_score: float
    parent_root: str
    root_score: float
    plausibility: OrphanPlausibility = OrphanPlausibility.UNCERTAIN
    action: OrphanAction = OrphanAction.EXPERT_REVIEW
    final_label: Optional[str] = None
    flag: str = ""


# Default rescue rules mapping subtype to (parent_root, plausibility, default_action)
# These are biologically-motivated defaults that can be customized via config
DEFAULT_RESCUE_RULES: Dict[str, Tuple[str, OrphanPlausibility, OrphanAction]] = {
    # High plausibility - rescue
    "Lymphatic Endothelium": ("Endothelium", OrphanPlausibility.HIGH, OrphanAction.RESCUE),
    "Smooth Muscle Cells": ("Mesenchymal Cells", OrphanPlausibility.HIGH, OrphanAction.RESCUE),
    "Fibroblasts": ("Mesenchymal Cells", OrphanPlausibility.HIGH, OrphanAction.RESCUE),
    "Vascular Endothelium": ("Endothelium", OrphanPlausibility.HIGH, OrphanAction.RESCUE),
    # Misc subtypes - no parent gate, can be assigned directly
    "Pericytes": ("Misc", OrphanPlausibility.HIGH, OrphanAction.RESCUE),
    "Proliferating Cells": ("Misc", OrphanPlausibility.HIGH, OrphanAction.RESCUE),
    # Suspicious - keep as Unassigned (CD45 should be on all immune cells)
    "Myeloids": ("Immune Cells", OrphanPlausibility.SUSPICIOUS, OrphanAction.KEEP_UNASSIGNED),
    "Lymphoids": ("Immune Cells", OrphanPlausibility.SUSPICIOUS, OrphanAction.KEEP_UNASSIGNED),
    "T Cells": ("Immune Cells", OrphanPlausibility.SUSPICIOUS, OrphanAction.KEEP_UNASSIGNED),
    "B Cells": ("Immune Cells", OrphanPlausibility.SUSPICIOUS, OrphanAction.KEEP_UNASSIGNED),
    "Natural-Killer (NK) Cells": ("Immune Cells", OrphanPlausibility.SUSPICIOUS, OrphanAction.KEEP_UNASSIGNED),
    "Granulocytes": ("Immune Cells", OrphanPlausibility.SUSPICIOUS, OrphanAction.KEEP_UNASSIGNED),
    "Monocytes": ("Immune Cells", OrphanPlausibility.SUSPICIOUS, OrphanAction.KEEP_UNASSIGNED),
    "Neutrophils": ("Immune Cells", OrphanPlausibility.SUSPICIOUS, OrphanAction.KEEP_UNASSIGNED),
    "Dendritic Cells": ("Immune Cells", OrphanPlausibility.SUSPICIOUS, OrphanAction.KEEP_UNASSIGNED),
    "Macrophages": ("Immune Cells", OrphanPlausibility.SUSPICIOUS, OrphanAction.KEEP_UNASSIGNED),
    "Activated T Cells": ("Immune Cells", OrphanPlausibility.SUSPICIOUS, OrphanAction.KEEP_UNASSIGNED),
    # Uncertain - expert review (EMT or dedifferentiation possible)
    "Ciliated Epithelium": ("Epithelium", OrphanPlausibility.UNCERTAIN, OrphanAction.EXPERT_REVIEW),
    "Glandular Epithelium": ("Epithelium", OrphanPlausibility.UNCERTAIN, OrphanAction.EXPERT_REVIEW),
    "Peg Cells": ("Epithelium", OrphanPlausibility.UNCERTAIN, OrphanAction.EXPERT_REVIEW),
}

# Root types that have subtypes
ROOT_TYPES: Set[str] = {"Epithelium", "Endothelium", "Immune Cells", "Mesenchymal Cells", "Misc"}


def detect_orphaned_subtypes(
    marker_scores: pd.DataFrame,
    unassigned_cluster_ids: List[str],
    subtype_threshold: float = 0.8,
    root_fail_threshold: float = 0.5,
    config: Optional[OrphanRescueConfig] = None,
    rescue_rules: Optional[Dict[str, Tuple[str, OrphanPlausibility, OrphanAction]]] = None,
) -> List[OrphanCandidate]:
    """Detect orphaned subtypes in Unassigned clusters.

    Parameters
    ----------
    marker_scores : pd.DataFrame
        DataFrame with columns: cluster_id, label, path, level, score, n_cells
    unassigned_cluster_ids : List[str]
        List of cluster IDs that are currently Unassigned
    subtype_threshold : float
        Minimum subtype score to consider as orphan candidate (default: 0.8)
    root_fail_threshold : float
        Maximum root score to consider as "failed" root gate (default: 0.5)
    config : OrphanRescueConfig, optional
        Optional config with custom thresholds and rules
    rescue_rules : Dict, optional
        Custom rescue rules. If None, uses DEFAULT_RESCUE_RULES.

    Returns
    -------
    List[OrphanCandidate]
        List of detected orphan candidates
    """
    if config is not None:
        subtype_threshold = config.subtype_threshold
        root_fail_threshold = config.root_fail_threshold

    if rescue_rules is None:
        rescue_rules = DEFAULT_RESCUE_RULES

    # Get custom rules from config
    custom_rules = {}
    if config is not None:
        custom_rules = {r.subtype: r for r in config.rules}

    # Filter to unassigned clusters
    unassigned_ids = set(unassigned_cluster_ids)
    unassigned_scores = marker_scores[marker_scores["cluster_id"].isin(unassigned_ids)]

    if len(unassigned_scores) == 0:
        return []

    candidates = []

    for cluster_id in unassigned_ids:
        cluster_scores = unassigned_scores[unassigned_scores["cluster_id"] == cluster_id]
        if len(cluster_scores) == 0:
            continue

        n_cells = cluster_scores["n_cells"].iloc[0]

        # Get root scores
        root_scores = {}
        for root in ROOT_TYPES:
            root_row = cluster_scores[cluster_scores["label"] == root]
            if len(root_row) > 0:
                root_scores[root] = root_row["score"].iloc[0]

        # Check each potential subtype
        best_orphan: Optional[OrphanCandidate] = None
        best_score = -float("inf")

        for _, row in cluster_scores.iterrows():
            subtype = row["label"]

            # Skip root types
            if subtype in ROOT_TYPES:
                continue

            # Skip if subtype not in rescue rules
            if subtype not in rescue_rules:
                continue

            subtype_score = row["score"]
            parent_root, plausibility, default_action = rescue_rules[subtype]
            root_score = root_scores.get(parent_root, 0.0)

            # Check custom rule for this subtype
            if subtype in custom_rules:
                custom_rule = custom_rules[subtype]
                min_score = custom_rule.min_score or subtype_threshold
                action_str = custom_rule.action
                action = OrphanAction(action_str)
            else:
                min_score = subtype_threshold
                action = default_action

            # Check if this is an orphan candidate
            # For Misc subtypes, there's no parent gate to fail
            if parent_root == "Misc":
                # Misc subtypes just need good subtype score
                if subtype_score >= min_score:
                    if subtype_score > best_score:
                        best_score = subtype_score
                        best_orphan = OrphanCandidate(
                            cluster_id=cluster_id,
                            n_cells=n_cells,
                            subtype=subtype,
                            subtype_score=subtype_score,
                            parent_root=parent_root,
                            root_score=0.0,  # No root for Misc
                            plausibility=plausibility,
                            action=action,
                        )
            else:
                # Normal case: failed root gate + good subtype score
                if root_score < root_fail_threshold and subtype_score >= min_score:
                    if subtype_score > best_score:
                        best_score = subtype_score
                        best_orphan = OrphanCandidate(
                            cluster_id=cluster_id,
                            n_cells=n_cells,
                            subtype=subtype,
                            subtype_score=subtype_score,
                            parent_root=parent_root,
                            root_score=root_score,
                            plausibility=plausibility,
                            action=action,
                        )

        if best_orphan is not None:
            candidates.append(best_orphan)

    return candidates


def apply_orphan_rescue(
    candidates: List[OrphanCandidate],
    suffix: str = "(orphan)",
    config: Optional[OrphanRescueConfig] = None,
) -> List[OrphanCandidate]:
    """Apply rescue logic to orphan candidates.

    Parameters
    ----------
    candidates : List[OrphanCandidate]
        List of orphan candidates from detect_orphaned_subtypes
    suffix : str
        Suffix to add to rescued labels (default: "(orphan)")
    config : OrphanRescueConfig, optional
        Optional config with custom suffix per subtype

    Returns
    -------
    List[OrphanCandidate]
        Updated candidates with final_label set for rescued orphans
    """
    # Get custom rules from config
    custom_rules = {}
    if config is not None:
        custom_rules = {r.subtype: r for r in config.rules}

    for candidate in candidates:
        # Determine suffix and flag for this subtype
        if candidate.subtype in custom_rules:
            custom_rule = custom_rules[candidate.subtype]
            subtype_suffix = custom_rule.suffix or suffix
            candidate.flag = custom_rule.flag
        else:
            subtype_suffix = suffix

        # Apply action
        if candidate.action == OrphanAction.RESCUE:
            # Rescue with suffix
            candidate.final_label = f"{candidate.subtype} {subtype_suffix}".strip()
        elif candidate.action == OrphanAction.KEEP_UNASSIGNED:
            # Keep as Unassigned
            candidate.final_label = "Unassigned"
            candidate.flag = candidate.flag or "suspicious_orphan"
        else:
            # Expert review - keep as Unassigned but flag
            candidate.final_label = "Unassigned"
            candidate.flag = candidate.flag or "needs_expert_review"

    return candidates


def summarize_orphans(candidates: List[OrphanCandidate]) -> pd.DataFrame:
    """Create summary DataFrame of orphan candidates.

    Parameters
    ----------
    candidates : List[OrphanCandidate]
        List of orphan candidates

    Returns
    -------
    pd.DataFrame
        Summary with columns: subtype, n_clusters, n_cells, mean_score,
        plausibility, action
    """
    if not candidates:
        return pd.DataFrame(
            columns=[
                "subtype",
                "n_clusters",
                "n_cells",
                "mean_subtype_score",
                "mean_root_score",
                "plausibility",
                "action",
            ]
        )

    data = []
    for c in candidates:
        data.append(
            {
                "cluster_id": c.cluster_id,
                "n_cells": c.n_cells,
                "subtype": c.subtype,
                "subtype_score": c.subtype_score,
                "parent_root": c.parent_root,
                "root_score": c.root_score,
                "plausibility": c.plausibility.value,
                "action": c.action.value,
                "final_label": c.final_label,
                "flag": c.flag,
            }
        )

    df = pd.DataFrame(data)

    # Create summary by subtype
    summary = (
        df.groupby("subtype")
        .agg(
            n_clusters=("cluster_id", "count"),
            n_cells=("n_cells", "sum"),
            mean_subtype_score=("subtype_score", "mean"),
            mean_root_score=("root_score", "mean"),
            plausibility=("plausibility", "first"),
            action=("action", "first"),
        )
        .reset_index()
    )

    return summary.sort_values("n_cells", ascending=False)


def get_orphan_rescue_map(
    candidates: List[OrphanCandidate],
) -> Dict[str, Tuple[str, str]]:
    """Get mapping from cluster_id to (final_label, reason) for rescued orphans.

    Parameters
    ----------
    candidates : List[OrphanCandidate]
        List of orphan candidates (after apply_orphan_rescue)

    Returns
    -------
    Dict[str, Tuple[str, str]]
        Mapping: cluster_id -> (final_label, reason)
    """
    result = {}
    for c in candidates:
        if c.final_label is not None:
            reason = f"Orphan rescue: {c.subtype} (score={c.subtype_score:.2f}, root={c.root_score:.2f})"
            if c.flag:
                reason += f" [flag={c.flag}]"
            result[c.cluster_id] = (c.final_label, reason)
    return result


def compute_all_unassigned_scores(
    marker_scores: pd.DataFrame,
    unassigned_cluster_ids: List[str],
) -> pd.DataFrame:
    """Compute root and subtype scores for ALL Unassigned clusters.

    This function is used for visualization purposes, to show all clusters
    in the orphan detection scatter plot (not just orphan candidates).

    Parameters
    ----------
    marker_scores : pd.DataFrame
        DataFrame with columns: cluster_id, label, path, level, score, n_cells
    unassigned_cluster_ids : List[str]
        List of cluster IDs that are currently Unassigned

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: cluster_id, n_cells, root_score, subtype_score,
        best_root, best_subtype
    """
    unassigned_ids = set(str(c) for c in unassigned_cluster_ids)
    unassigned_scores = marker_scores[marker_scores["cluster_id"].astype(str).isin(unassigned_ids)]

    if len(unassigned_scores) == 0:
        return pd.DataFrame(columns=[
            "cluster_id", "n_cells", "root_score", "subtype_score",
            "best_root", "best_subtype"
        ])

    results = []

    for cluster_id in unassigned_ids:
        cluster_scores = unassigned_scores[unassigned_scores["cluster_id"].astype(str) == cluster_id]
        if len(cluster_scores) == 0:
            continue

        n_cells = cluster_scores["n_cells"].iloc[0]

        # Find best root score
        best_root_score = 0.0
        best_root_label = ""
        for root in ROOT_TYPES:
            root_row = cluster_scores[cluster_scores["label"] == root]
            if len(root_row) > 0:
                score = root_row["score"].iloc[0]
                if score > best_root_score:
                    best_root_score = score
                    best_root_label = root

        # Find best subtype score (excluding root types)
        best_subtype_score = 0.0
        best_subtype_label = ""
        for _, row in cluster_scores.iterrows():
            label = row["label"]
            if label in ROOT_TYPES:
                continue
            score = row["score"]
            if score > best_subtype_score:
                best_subtype_score = score
                best_subtype_label = label

        results.append({
            "cluster_id": cluster_id,
            "n_cells": n_cells,
            "root_score": best_root_score,
            "subtype_score": best_subtype_score,
            "best_root": best_root_label,
            "best_subtype": best_subtype_label,
        })

    return pd.DataFrame(results)
