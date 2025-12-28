"""Hierarchical gating for cell-type annotation.

This module implements hierarchical top-down assignment with sibling-only
comparison. It replaces flat max-score assignment with a hierarchical descent
that uses parent markers as gates and compares only siblings at each level.

Tissue-specific gating parameters (root_hard_requirements, root_veto_markers)
should be provided via configuration files rather than hardcoded defaults.
"""

from __future__ import annotations

import json
import logging
from collections import defaultdict
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy import sparse

try:
    import scanpy as sc
except ImportError:
    sc = None

from .marker_loading import MarkerSet


# Default gating parameters by hierarchy level
# NOTE: root_hard_requirements and root_veto_markers are empty by default.
# These are tissue-specific and should be provided via configs/tissues/*.yaml
DEFAULT_GATING_PARAMS = {
    "min_coverage": {0: 0.5, 1: 0.4, 2: 0.3, 3: 0.3},
    "min_pos_frac": {0: 0.3, 1: 0.20, 2: 0.15, 3: 0.15},
    "min_enrichment": {0: 0.0, 1: -0.5, 2: -1.0, 3: -1.0},
    "min_gap": {0: 0.5, 1: 0.3, 2: 0.2, 3: 0.2},
    "min_frac_markers_on": {0: 0.4, 1: 0.3, 2: 0.25, 3: 0.25},
    "anti_penalty_hard_gate": 1.0,
    "per_cell_gap_threshold": 0.1,
    "anti_agg": "top2mean",
    "root_gap_threshold": 0.25,
    # Tissue-specific: provide via config
    "root_hard_requirements": {},
    "root_veto_markers": {},
}


def _sort_cluster_key(cluster_id: str) -> Tuple[int, int]:
    """Sort key for cluster IDs that handles both plain integers and subclusters.

    Examples:
        "3" -> (3, -1)
        "1:0" -> (1, 0)
        "21:5" -> (21, 5)
    """
    cluster_str = str(cluster_id)
    if ":" in cluster_str:
        parent, sub = cluster_str.split(":", 1)
        try:
            return (int(parent), int(sub))
        except ValueError:
            return (int(1e9), 0)
    else:
        try:
            return (int(cluster_str), -1)
        except ValueError:
            return (int(1e9), 0)


def _build_hierarchy_tree(marker_scores: pd.DataFrame) -> Dict[str, List[str]]:
    """Build parent→children mapping from marker_scores path column."""
    hierarchy: Dict[str, set] = defaultdict(set)

    for path in marker_scores["path"].unique():
        parts = [p.strip() for p in str(path).split(" / ")]
        for i in range(len(parts) - 1):
            parent = parts[i]
            child = parts[i + 1]
            hierarchy[parent].add(child)

    return {k: sorted(list(v)) for k, v in hierarchy.items()}


def _get_root_labels(hierarchy: Dict[str, List[str]], all_labels: set) -> List[str]:
    """Find root labels (labels that are never children of other labels)."""
    all_children = set()
    for children in hierarchy.values():
        all_children.update(children)

    roots = [label for label in all_labels if label not in all_children]
    return sorted(roots)


def _has_anti_conflict(row: pd.Series, params: dict) -> bool:
    """Hard gate: reject node if anti-marker penalty exceeds threshold."""
    threshold = params.get("anti_penalty_hard_gate", 1.0)
    anti_penalty = row.get("anti_penalty", 0.0)
    return anti_penalty > threshold


def _precompute_marker_stats(
    adata: "sc.AnnData",
    cluster_key: str,
    markers: List[str],
    positive_quantile: float = 0.75,
    layer: str = "X",
    logger: Optional[logging.Logger] = None,
) -> Dict[str, Dict[str, Dict[str, float]]]:
    """Pre-compute per-marker statistics for all clusters.

    Used for root-level hard requirements and veto checks.

    Returns:
        Nested dict: {cluster_id: {marker_name: {"pos_frac": float, "enrichment": float}}}
    """
    if layer in adata.layers:
        matrix = adata.layers[layer]
    else:
        matrix = adata.X
    if sparse.issparse(matrix):
        matrix = matrix.toarray()

    var_index = {name: idx for idx, name in enumerate(adata.var_names.astype(str))}
    obs_clusters = adata.obs[cluster_key].astype(str)

    # Compute global stats
    marker_global_stats = {}
    for marker in markers:
        if marker not in var_index:
            continue
        idx = var_index[marker]
        marker_global_stats[marker] = {
            "median": float(np.nanmedian(matrix[:, idx])),
            "std": float(np.nanstd(matrix[:, idx])),
            "threshold": float(np.nanquantile(matrix[:, idx], positive_quantile)),
        }
        if marker_global_stats[marker]["std"] == 0:
            marker_global_stats[marker]["std"] = 1e-6

    # Compute per-cluster stats
    result: Dict[str, Dict[str, Dict[str, float]]] = {}
    clusters = sorted(obs_clusters.unique())

    for cluster in clusters:
        mask = obs_clusters == cluster
        cluster_matrix = matrix[mask.to_numpy(), :]
        result[cluster] = {}

        for marker in markers:
            if marker not in var_index:
                result[cluster][marker] = {"pos_frac": 0.0, "enrichment": 0.0}
                continue

            idx = var_index[marker]
            gs = marker_global_stats[marker]
            cluster_vals = cluster_matrix[:, idx]

            cluster_median = float(np.nanmedian(cluster_vals))
            pos_frac = float(np.mean(cluster_vals >= gs["threshold"]))
            enrichment = (cluster_median - gs["median"]) / gs["std"]

            result[cluster][marker] = {
                "pos_frac": pos_frac,
                "enrichment": float(enrichment),
            }

    if logger:
        logger.debug(
            "Pre-computed stats for %d markers across %d clusters",
            len(markers),
            len(clusters),
        )

    return result


def _per_cell_sibling_vote(
    adata: "sc.AnnData",
    cluster_id: str,
    candidates: List[str],
    marker_scores: pd.DataFrame,
    cluster_key: str = "cluster_lvl0",
    params: Optional[dict] = None,
) -> Dict[str, Any]:
    """Compute per-cell voting among ambiguous sibling candidates.

    For each cell in the cluster, compute a score for each candidate
    using their resolved markers, then vote for the best candidate.

    NOTE: Voting uses Z-scored values from adata.X (the PCA input matrix).
    """
    if params is None:
        params = DEFAULT_GATING_PARAMS

    gap_threshold = params.get("per_cell_gap_threshold", 0.1)

    # Get cells in this cluster
    cluster_mask = adata.obs[cluster_key].astype(str) == str(cluster_id)
    n_cells = int(cluster_mask.sum())

    empty_result = {
        "composition": {"Uncertain": 1.0},
        "candidate_markers": {},
        "expr_stats": {},
        "margin_stats": {"mean_gap": 0.0, "median_gap": 0.0, "pct_uncertain": 100.0},
        "n_cells": n_cells,
        "matrix": "X_zscore",
    }

    if n_cells == 0 or adata.X is None:
        return empty_result

    X = adata.X
    var_names = list(adata.var_names)
    var_index = {name: i for i, name in enumerate(var_names)}

    # Get markers for each candidate
    candidate_markers = {}
    for candidate in candidates:
        rows = marker_scores[
            (marker_scores["cluster_id"].astype(str) == str(cluster_id))
            & (marker_scores["label"] == candidate)
        ]
        if len(rows) == 0:
            continue
        row = rows.iloc[0]
        markers_str = row.get("resolved_markers", "")
        if pd.isna(markers_str) or not markers_str:
            continue
        markers = [m.strip() for m in str(markers_str).split(";") if m.strip()]
        valid_markers = [m for m in markers if m in var_index]
        if valid_markers:
            candidate_markers[candidate] = valid_markers

    if not candidate_markers:
        empty_result["candidate_markers"] = candidate_markers
        return empty_result

    # Initialize vote counts
    votes = {c: 0 for c in candidate_markers.keys()}
    votes["Uncertain"] = 0

    # Get cluster cell indices
    cluster_indices = np.where(cluster_mask)[0]

    if sparse.issparse(X):
        X_dense = X[cluster_indices, :].toarray()
    else:
        X_dense = X[cluster_indices, :]

    per_cell_expr = {c: [] for c in candidate_markers.keys()}
    all_margins = []

    # For each cell, vote
    for cell_idx in range(len(cluster_indices)):
        cell_scores = {}

        for candidate, markers in candidate_markers.items():
            marker_indices = [var_index[m] for m in markers]
            cell_expr = X_dense[cell_idx, marker_indices]
            mean_expr = float(np.nanmean(cell_expr))
            cell_scores[candidate] = mean_expr
            per_cell_expr[candidate].append(mean_expr)

        if not cell_scores:
            votes["Uncertain"] += 1
            continue

        sorted_scores = sorted(cell_scores.items(), key=lambda x: x[1], reverse=True)

        if len(sorted_scores) >= 2:
            gap = sorted_scores[0][1] - sorted_scores[1][1]
            all_margins.append(gap)
            if gap < gap_threshold:
                votes["Uncertain"] += 1
            else:
                votes[sorted_scores[0][0]] += 1
        elif len(sorted_scores) == 1:
            votes[sorted_scores[0][0]] += 1
            all_margins.append(float("inf"))
        else:
            votes["Uncertain"] += 1

    # Convert to fractions
    total = sum(votes.values())
    if total == 0:
        composition = {"Uncertain": 1.0}
    else:
        composition = {k: round(v / total, 3) for k, v in votes.items() if v > 0}

    # Expression stats per candidate
    expr_stats = {}
    for candidate, expr_list in per_cell_expr.items():
        if expr_list:
            arr = np.array(expr_list)
            expr_stats[candidate] = {
                "mean": float(np.nanmean(arr)),
                "std": float(np.nanstd(arr)),
            }
        else:
            expr_stats[candidate] = {"mean": 0.0, "std": 0.0}

    # Margin stats
    if all_margins:
        finite_margins = [m for m in all_margins if m != float("inf")]
        if finite_margins:
            margin_arr = np.array(finite_margins)
            pct_uncertain = (votes.get("Uncertain", 0) / total * 100) if total > 0 else 0.0
            margin_stats = {
                "mean_gap": float(np.mean(margin_arr)),
                "median_gap": float(np.median(margin_arr)),
                "pct_uncertain": round(pct_uncertain, 1),
            }
        else:
            margin_stats = {"mean_gap": 0.0, "median_gap": 0.0, "pct_uncertain": 0.0}
    else:
        margin_stats = {"mean_gap": 0.0, "median_gap": 0.0, "pct_uncertain": 100.0}

    return {
        "composition": composition,
        "candidate_markers": candidate_markers,
        "expr_stats": expr_stats,
        "margin_stats": margin_stats,
        "n_cells": n_cells,
        "matrix": "X_zscore",
    }


def assign_labels_hierarchical(
    marker_scores: pd.DataFrame,
    adata: Optional["sc.AnnData"] = None,
    cluster_key: str = "cluster_lvl0",
    layer: str = "X",
    params: Optional[dict] = None,
    logger: Optional[logging.Logger] = None,
    marker_sets: Optional[Sequence[MarkerSet]] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Hierarchical top-down assignment with sibling-only comparison.

    Algorithm:
    1. Build hierarchy tree from path column
    2. For each cluster:
       a. Find best passing root (top-level category)
       b. Descend through hierarchy, comparing only siblings
       c. At each level: check gates, compare siblings, track decisions
       d. Handle ambiguous siblings with per-cell voting (if adata provided)
       e. Stop when no child passes gate or leaf reached

    Args:
        marker_scores: DataFrame with columns: cluster_id, label, path, level,
                      score, coverage, mean_positive_fraction, mean_enrichment, etc.
        adata: Optional AnnData for per-cell voting in ambiguous cases
        cluster_key: Column name for cluster assignments in adata.obs
        layer: Which layer to use for expression data
        params: Gating parameters dict (uses DEFAULT_GATING_PARAMS if None)
        logger: Logger instance
        marker_sets: Optional list of MarkerSet objects for gating overrides

    Returns:
        Tuple of (cluster_annotations, decision_steps) DataFrames
    """
    if params is None:
        params = DEFAULT_GATING_PARAMS

    if logger is None:
        logger = logging.getLogger(__name__)

    if marker_scores.empty:
        return pd.DataFrame(), pd.DataFrame()

    # Build hierarchy tree
    hierarchy = _build_hierarchy_tree(marker_scores)
    all_labels = set(marker_scores["label"].unique())
    roots = _get_root_labels(hierarchy, all_labels)

    # Build lookup from label to MarkerSet for gating overrides
    mset_lookup: Dict[str, MarkerSet] = {}
    if marker_sets:
        for mset in marker_sets:
            mset_lookup[mset.label] = mset

    logger.info(
        "Hierarchical assignment: %d roots, %d parent nodes in hierarchy",
        len(roots),
        len(hierarchy),
    )

    annotations = []
    all_steps = []

    # Pre-compute marker stats for root gates
    root_marker_stats: Dict[str, Dict[str, Dict[str, float]]] = {}
    if adata is not None:
        root_gate_markers = set()
        for root_req in params.get("root_hard_requirements", {}).values():
            if "marker" in root_req:
                root_gate_markers.add(root_req["marker"])
        for root_veto in params.get("root_veto_markers", {}).values():
            root_gate_markers.update(root_veto.get("markers", []))

        if root_gate_markers:
            root_marker_stats = _precompute_marker_stats(
                adata,
                cluster_key=cluster_key,
                markers=list(root_gate_markers),
                positive_quantile=0.75,
                layer=layer,
                logger=logger,
            )

    def get_threshold(param_dict: dict, level: int) -> float:
        if level in param_dict:
            return param_dict[level]
        max_defined = max(param_dict.keys())
        return param_dict[max_defined]

    # Process each cluster using groupby (O(n) instead of O(n²) filtering)
    # Pre-group marker_scores by cluster_id for efficient iteration
    grouped_scores = marker_scores.groupby("cluster_id", sort=False)
    n_clusters = grouped_scores.ngroups
    logger.info("Processing %d clusters...", n_clusters)

    for cluster_id, cluster_scores in grouped_scores:
        cluster_labels = set(cluster_scores["label"].values)

        first_row = cluster_scores.iloc[0] if len(cluster_scores) > 0 else None
        cluster_n_cells = int(first_row.get("n_cells", 0)) if first_row is not None else 0

        # Pre-build label -> row lookup for O(1) access instead of O(n) filtering
        label_to_row = {row["label"]: row for _, row in cluster_scores.iterrows()}

        def get_row(label: str) -> Optional[pd.Series]:
            return label_to_row.get(label)

        def passes_gate(
            row: pd.Series, level: int, cluster_id_for_root: str = ""
        ) -> Tuple[bool, Optional[str]]:
            label = row.get("label", "")
            mset = mset_lookup.get(label)
            overrides = mset.gating_overrides if mset and mset.gating_overrides else {}

            # Root hard requirements (level 0 only)
            if level == 0 and label in params.get("root_hard_requirements", {}):
                req = params["root_hard_requirements"][label]
                cluster_stats = root_marker_stats.get(str(cluster_id_for_root), {})
                marker_stats = cluster_stats.get(
                    req["marker"], {"pos_frac": 0.0, "enrichment": -999.0}
                )

                if marker_stats["pos_frac"] < req.get("min_pos_frac", 0.0):
                    return False, f"root_hard_req_pos_frac_{req['marker']}"
                if marker_stats["enrichment"] < req.get("min_enrichment", -999.0):
                    return False, f"root_hard_req_enrichment_{req['marker']}"

            # Root veto markers (level 0 only)
            if level == 0 and label in params.get("root_veto_markers", {}):
                veto = params["root_veto_markers"][label]
                cluster_stats = root_marker_stats.get(str(cluster_id_for_root), {})
                for marker in veto.get("markers", []):
                    marker_stats = cluster_stats.get(marker, {"pos_frac": 0.0})
                    if marker_stats["pos_frac"] > veto.get("max_pos_frac", 1.0):
                        return False, f"root_veto_{marker}"

            # Coverage check
            min_cov = overrides.get(
                "min_coverage", get_threshold(params.get("min_coverage", {}), level)
            )
            if row.get("coverage", 0) < min_cov:
                return False, f"low_coverage(<{min_cov:.2f})"

            # Positive fraction check with alternative frac_markers_on gate
            min_pos = overrides.get(
                "min_pos_frac", get_threshold(params.get("min_pos_frac", {}), level)
            )
            min_frac_on = overrides.get(
                "min_frac_markers_on",
                get_threshold(params.get("min_frac_markers_on", {0: 0.3, 1: 0.25}), level),
            )

            pos_frac_ok = row.get("mean_positive_fraction", 0) >= min_pos
            frac_on_ok = row.get("frac_markers_on", 0) >= min_frac_on

            if not pos_frac_ok and not frac_on_ok:
                return False, f"low_pos_frac(<{min_pos:.2f})"

            # Enrichment check
            min_enrich = overrides.get(
                "min_enrichment", get_threshold(params.get("min_enrichment", {}), level)
            )
            if row.get("mean_enrichment", 0) < min_enrich:
                return False, f"low_enrichment(<{min_enrich:.2f})"

            # Stricter anti-marker check if required
            if overrides.get("require_anti_markers_off", False):
                if row.get("mean_anti_positive", 0) > 0.2:
                    return False, "anti_markers_not_off"

            # Hard anti-marker gate
            if _has_anti_conflict(row, params):
                return False, "anti_marker_conflict"

            return True, None

        # Find best passing root
        root_candidates = []
        for root in roots:
            if root not in cluster_labels:
                continue
            row = get_row(root)
            if row is None:
                continue

            passed, reason = passes_gate(row, level=0, cluster_id_for_root=cluster_id)
            all_steps.append({
                "cluster_id": cluster_id,
                "step_idx": 0,
                "parent_label": "ROOT",
                "child_label": root,
                "child_passed_gate": passed,
                "child_score": row["score"],
                "child_coverage": row.get("coverage", np.nan),
                "child_pos_frac": row.get("mean_positive_fraction", np.nan),
                "child_frac_markers_on": row.get("frac_markers_on", np.nan),
                "child_enrichment": row.get("mean_enrichment", np.nan),
                "child_anti_penalty": row.get("anti_penalty", np.nan),
                "selected": False,
                "margin_to_runner_up": np.nan,
                "fail_reason": reason,
            })
            if passed:
                root_candidates.append((root, row))

        if not root_candidates:
            # No root passed - collect fail reasons for each root
            # Extract fail reasons from all_steps for this cluster
            root_fail_reasons = {}
            for step in all_steps:
                if (
                    step["cluster_id"] == cluster_id
                    and step["parent_label"] == "ROOT"
                    and step.get("fail_reason")
                ):
                    root_fail_reasons[step["child_label"]] = step["fail_reason"]

            # Format as semicolon-separated list for CSV readability
            # e.g., "Immune Cells:root_veto_E-cadherin;Epithelium:low_pos_frac(<0.30)"
            root_fail_summary = ";".join(
                f"{root}:{reason}" for root, reason in sorted(root_fail_reasons.items())
            )

            annotations.append({
                "cluster_id": cluster_id,
                "assigned_label": "Unassigned",
                "assigned_path": "",
                "assigned_level": -1,
                "assigned_score": 0.0,
                "root_label": "Unassigned",
                "confidence": 0.0,
                "min_margin_along_path": 0.0,
                "margin_is_infinite": False,
                "stop_reason": "no_root_passed",
                "stopped_before_leaf": True,
                "composition": None,
                "decision_trace": "[]",
                "n_cells": cluster_n_cells,
                "coverage": np.nan,
                "resolved_markers": "",
                "is_ambiguous_root": False,
                "ambiguous_root_candidates": "",
                "root_fail_reasons": root_fail_summary,
            })
            continue

        # Sort by score
        root_candidates.sort(key=lambda x: x[1]["score"], reverse=True)
        current_node, current_row = root_candidates[0]
        root_label = current_node

        # Calculate margin to runner-up root
        if len(root_candidates) >= 2:
            root_margin = root_candidates[0][1]["score"] - root_candidates[1][1]["score"]
        else:
            root_margin = float("inf")

        # Root ambiguity detection
        root_gap_thresh = params.get("root_gap_threshold", 0.25)
        if len(root_candidates) >= 2 and root_margin < root_gap_thresh:
            top_two = f"{root_candidates[0][0]}~{root_candidates[1][0]}"

            ambiguous_root_composition = {
                root_candidates[0][0]: 0.5,
                root_candidates[1][0]: 0.5,
                "_evidence_only": True,
                "_gap": float(root_margin),
                "_threshold": float(root_gap_thresh),
            }

            annotations.append({
                "cluster_id": cluster_id,
                "assigned_label": top_two,
                "assigned_path": top_two,
                "assigned_level": 0,
                "assigned_score": root_candidates[0][1]["score"],
                "root_label": top_two,
                "confidence": root_margin,
                "min_margin_along_path": root_margin,
                "margin_is_infinite": False,
                "stop_reason": "ambiguous_root",
                "stopped_before_leaf": True,
                "composition": json.dumps(ambiguous_root_composition),
                "decision_trace": json.dumps({
                    "type": "ambiguous_root",
                    "candidates": [
                        {"label": c[0], "score": float(c[1]["score"])}
                        for c in root_candidates[:3]
                    ],
                }),
                "n_cells": cluster_n_cells,
                "coverage": np.nan,
                "resolved_markers": "",
                "is_ambiguous_root": True,
                "ambiguous_root_candidates": f"{root_candidates[0][0]};{root_candidates[1][0]}",
                "root_fail_reasons": "",  # N/A for ambiguous root
            })
            continue

        # Mark selected root
        for step in all_steps:
            if (
                step["cluster_id"] == cluster_id
                and step["parent_label"] == "ROOT"
                and step["child_label"] == current_node
            ):
                step["selected"] = True
                step["margin_to_runner_up"] = root_margin

        trace = [(current_node, float(current_row["score"]))]
        min_margin = root_margin
        step_idx = 1
        stop_reason = "leaf_reached"
        composition = None

        # Descend through hierarchy
        while current_node in hierarchy:
            children = hierarchy[current_node]
            child_candidates = []

            for child in children:
                if child not in cluster_labels:
                    continue
                row = get_row(child)
                if row is None:
                    continue

                level = int(row.get("level", 1))
                passed, reason = passes_gate(row, level)

                all_steps.append({
                    "cluster_id": cluster_id,
                    "step_idx": step_idx,
                    "parent_label": current_node,
                    "child_label": child,
                    "child_passed_gate": passed,
                    "child_score": row["score"],
                    "child_coverage": row.get("coverage", np.nan),
                    "child_pos_frac": row.get("mean_positive_fraction", np.nan),
                    "child_frac_markers_on": row.get("frac_markers_on", np.nan),
                    "child_enrichment": row.get("mean_enrichment", np.nan),
                    "child_anti_penalty": row.get("anti_penalty", np.nan),
                    "selected": False,
                    "margin_to_runner_up": np.nan,
                    "fail_reason": reason,
                })

                if passed:
                    child_candidates.append((child, row))

            step_idx += 1

            if not child_candidates:
                stop_reason = "no_child_passed"
                break

            child_candidates.sort(key=lambda x: x[1]["score"], reverse=True)
            best_child, best_row = child_candidates[0]

            if len(child_candidates) >= 2:
                gap = best_row["score"] - child_candidates[1][1]["score"]
            else:
                gap = float("inf")

            min_margin = min(min_margin, gap)

            # Check for ambiguous siblings
            level = int(best_row.get("level", 1))
            min_gap_threshold = get_threshold(params.get("min_gap", {}), level)

            if len(child_candidates) >= 2 and gap < min_gap_threshold:
                # Ambiguous siblings - stop at parent
                score_winner = child_candidates[0][0]
                runner_up = child_candidates[1][0]

                vote_evidence = None
                if adata is not None:
                    vote_evidence = _per_cell_sibling_vote(
                        adata,
                        cluster_id,
                        [c[0] for c in child_candidates[:3]],
                        marker_scores,
                        cluster_key=cluster_key,
                        params=params,
                    )

                stop_reason = "ambiguous_siblings"

                if vote_evidence:
                    composition = vote_evidence["composition"].copy()
                    composition["_evidence_only"] = True
                    composition["_top_by_score"] = score_winner
                    composition["_runner_up_by_score"] = runner_up
                    composition["_gap"] = gap
                    composition["_threshold"] = min_gap_threshold
                    composition["_cell_margins"] = vote_evidence["margin_stats"]
                else:
                    composition = {
                        score_winner: 0.5,
                        runner_up: 0.5,
                        "_evidence_only": True,
                        "_top_by_score": score_winner,
                        "_runner_up_by_score": runner_up,
                        "_gap": gap,
                        "_threshold": min_gap_threshold,
                    }

                break

            # Mark selected child
            for step in all_steps:
                if (
                    step["cluster_id"] == cluster_id
                    and step["step_idx"] == step_idx - 1
                    and step["child_label"] == best_child
                ):
                    step["selected"] = True
                    step["margin_to_runner_up"] = gap

            current_node, current_row = best_child, best_row
            trace.append((current_node, float(current_row["score"])))

        # Build annotation record
        n_cells = int(current_row.get("n_cells", 0)) if current_row is not None else 0
        coverage = float(current_row.get("coverage", np.nan)) if current_row is not None else np.nan
        resolved_markers = current_row.get("resolved_markers", "") if current_row is not None else ""

        MARGIN_SENTINEL = 1e6
        margin_is_infinite = min_margin == float("inf")
        effective_margin = MARGIN_SENTINEL if margin_is_infinite else min_margin

        annotations.append({
            "cluster_id": cluster_id,
            "assigned_label": current_node,
            "assigned_path": " / ".join([t[0] for t in trace]),
            "assigned_level": len(trace) - 1,
            "assigned_score": float(current_row["score"]) if current_row is not None else 0.0,
            "root_label": root_label,
            "confidence": effective_margin,
            "min_margin_along_path": effective_margin,
            "margin_is_infinite": margin_is_infinite,
            "stop_reason": stop_reason,
            "stopped_before_leaf": stop_reason in ("ambiguous_siblings", "no_child_passed"),
            "composition": json.dumps(composition) if composition else None,
            "decision_trace": json.dumps(trace),
            "n_cells": n_cells,
            "coverage": coverage,
            "resolved_markers": resolved_markers,
            "is_ambiguous_root": False,
            "ambiguous_root_candidates": "",
            "root_fail_reasons": "",  # N/A for assigned clusters
        })

    annotations_df = pd.DataFrame(annotations)
    steps_df = pd.DataFrame(all_steps)

    logger.info(
        "Hierarchical assignment complete: %d clusters processed",
        len(annotations_df),
    )
    if "stop_reason" in annotations_df.columns:
        stop_counts = annotations_df["stop_reason"].value_counts()
        for reason, count in stop_counts.items():
            logger.info("  Stop reason '%s': %d clusters", reason, count)

    return annotations_df, steps_df
