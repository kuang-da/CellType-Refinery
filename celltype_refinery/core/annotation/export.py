"""Export functions for annotation results.

Provides export utilities for cluster annotations and enhanced annotations
with multi-level lineage tracking.

Key Functions:
- export_enhanced_annotations: Create 44-column enhanced annotations
- export_cluster_annotations_stage_h_format: Create Stage H format (24 columns)
- get_hierarchical_runner_up: Get runner-up label from decision tree
- export_marker_scores: Export marker scores DataFrame to CSV
- export_composition_stats: Export composition statistics
- export_review_summary: Export review summary JSON
- export_workflow_state: Export workflow state YAML
- generate_enhanced_annotations_44col: Create 44-column format from 24-col input
- run_review_exports: High-level orchestrator for review exports
"""

from __future__ import annotations

import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    import scanpy as sc


# =============================================================================
# Lineage Grouping Constants
# =============================================================================

LINEAGE_GROUPS = {
    'Epithelium': [
        'Epithelium', 'Ciliated Epithelium', 'Glandular Epithelium',
        'Secretory Epithelium', 'Peg Cells'
    ],
    'Endothelium': ['Endothelium', 'Vascular Endothelium', 'Lymphatic Endothelium'],
    'Mesenchymal Cells': [
        'Mesenchymal Cells', 'Smooth Muscle Cells', 'Fibroblasts', 'Pericytes'
    ],
    'Immune Cells': [
        'Immune Cells', 'T Cells', 'B Cells', 'Macrophages', 'Dendritic Cells',
        'Myeloids', 'Lymphoids', 'NK Cells', 'Natural-Killer (NK) Cells',
        'Neutrophils', 'Monocytes', 'Granulocytes', 'Activated T Cells'
    ],
    'Unassigned': ['Unassigned'],
}


# =============================================================================
# Column Detection Utilities
# =============================================================================


def detect_cell_type_column(adata: "sc.AnnData") -> str:
    """Detect best cell type column with fallback chain.

    Priority: cell_type_curated > cell_type_lvl1 > cell_type_lvl0 >
              cell_type_auto > assigned_label > cell_type

    Parameters
    ----------
    adata : sc.AnnData
        AnnData object

    Returns
    -------
    str
        Name of detected cell type column

    Raises
    ------
    ValueError
        If no cell type column found
    """
    priority = [
        'cell_type_curated',
        'cell_type_lvl1',
        'cell_type_lvl0',
        'cell_type_auto',
        'assigned_label',
        'cell_type',
    ]
    for col in priority:
        if col in adata.obs.columns:
            return col
    raise ValueError("No cell type column found in adata.obs")


def detect_cluster_column(adata: "sc.AnnData") -> str:
    """Detect best cluster column with fallback chain.

    Priority: cluster_id > cluster_lvl1 > cluster_lvl0 > leiden

    Parameters
    ----------
    adata : sc.AnnData
        AnnData object

    Returns
    -------
    str
        Name of detected cluster column
    """
    priority = ['cluster_id', 'cluster_lvl1', 'cluster_lvl0', 'leiden']
    for col in priority:
        if col in adata.obs.columns:
            return col
    return 'cluster_id'  # Default


# =============================================================================
# Helper Functions
# =============================================================================


def _sort_cluster_key(cluster_id: str) -> Tuple:
    """Sort key for cluster IDs that handles subclusters (e.g., '1:0', '18:2').

    Examples:
        "1" -> (1, -1)
        "1:0" -> (1, 0)
        "18:2" -> (18, 2)
    """
    if ":" in cluster_id:
        parent, sub = cluster_id.split(":", 1)
        try:
            return (int(parent), int(sub))
        except ValueError:
            return (float("inf"), cluster_id)
    else:
        try:
            return (int(cluster_id), -1)
        except ValueError:
            return (float("inf"), cluster_id)


def _classify_confidence_band(score: float) -> str:
    """Classify score into confidence bands.

    Parameters
    ----------
    score : float
        Marker score

    Returns
    -------
    str
        Confidence level: "high", "medium", "low", or "very_low"

    Note
    ----
    Thresholds aligned with reference ft pipeline:
    - score >= 2.0: "high"
    - score >= 1.0: "medium"
    - score >= 0.5: "low"
    - else: "very_low"
    """
    if score >= 2.0:
        return "high"
    elif score >= 1.0:
        return "medium"
    elif score >= 0.5:
        return "low"
    else:
        return "very_low"


# =============================================================================
# Runner-Up Detection
# =============================================================================


def get_hierarchical_runner_up(
    cluster_id: str,
    decision_steps: pd.DataFrame,
    cluster_annotations: pd.DataFrame,
    marker_scores: Optional[pd.DataFrame] = None,
) -> Tuple[str, float, float]:
    """Get hierarchical runner-up: the losing sibling at the tightest decision point.

    Instead of using flat score ranking (which can incorrectly identify the assigned
    label as runner-up), this function traces the actual decision tree to find which
    alternative was closest to winning at each branch.

    Parameters
    ----------
    cluster_id : str
        Cluster ID to find runner-up for
    decision_steps : pd.DataFrame
        Decision trace with columns: cluster_id, level, candidates, winner, scores, margin
        Each row represents one decision point in the hierarchy.
    cluster_annotations : pd.DataFrame
        Cluster annotations with columns: cluster_id, assigned_label, assigned_path
    marker_scores : pd.DataFrame, optional
        Full marker scores for fallback (flat max if no decision data)

    Returns
    -------
    Tuple[str, float, float]
        (runner_up_label, runner_up_score, margin_from_assigned)

    Examples
    --------
    >>> # For cluster 5 with path "Immune/T Cells" and margin 0.15 at T-cell branch:
    >>> runner_up, score, margin = get_hierarchical_runner_up(
    ...     "5", decision_steps, annotations)
    >>> print(f"Runner-up: {runner_up}, margin: {margin:.2f}")
    Runner-up: B Cells, margin: 0.15

    Notes
    -----
    The margin represents how close the decision was at the tightest point.
    A small margin (e.g., <0.3) indicates the assignment was borderline.
    """
    # Default return for missing data
    default_result = ("Unknown", 0.0, float("inf"))

    # Get annotation for this cluster
    ann_row = cluster_annotations[
        cluster_annotations["cluster_id"].astype(str) == str(cluster_id)
    ]
    if ann_row.empty:
        return default_result

    assigned_label = str(ann_row.iloc[0].get("assigned_label", "Unknown"))
    assigned_score = float(ann_row.iloc[0].get("assigned_score", 0))

    # If no decision steps available, fall back to flat marker scores
    if decision_steps is None or decision_steps.empty:
        if marker_scores is not None and not marker_scores.empty:
            cluster_scores = marker_scores[
                marker_scores["cluster_id"].astype(str) == str(cluster_id)
            ]
            if len(cluster_scores) >= 2:
                sorted_scores = cluster_scores.nlargest(2, "score")
                runner_up_row = sorted_scores.iloc[1]
                runner_up_label = str(runner_up_row.get("label", "Unknown"))
                runner_up_score = float(runner_up_row.get("score", 0))
                margin = assigned_score - runner_up_score
                return (runner_up_label, runner_up_score, margin)
        return default_result

    # Get decision steps for this cluster
    cluster_decisions = decision_steps[
        decision_steps["cluster_id"].astype(str) == str(cluster_id)
    ]
    if cluster_decisions.empty:
        return default_result

    # Find the tightest decision (minimum margin)
    min_margin = float("inf")
    tightest_runner_up = "Unknown"
    tightest_score = 0.0

    for _, step in cluster_decisions.iterrows():
        # Parse candidates and scores
        candidates = step.get("candidates", "")
        scores_str = step.get("scores", "")
        winner = step.get("winner", "")
        margin = step.get("margin", float("inf"))

        # Skip if no margin data
        if pd.isna(margin) or margin == float("inf"):
            continue

        # Parse candidates list
        if isinstance(candidates, str):
            candidate_list = [c.strip() for c in candidates.split(",") if c.strip()]
        elif isinstance(candidates, (list, tuple)):
            candidate_list = list(candidates)
        else:
            continue

        # Parse scores (if available)
        score_list = []
        if isinstance(scores_str, str):
            try:
                score_list = [float(s.strip()) for s in scores_str.split(",") if s.strip()]
            except ValueError:
                pass
        elif isinstance(scores_str, (list, tuple)):
            score_list = [float(s) for s in scores_str]

        # Find runner-up at this decision point
        if len(candidate_list) >= 2 and len(score_list) >= 2:
            # Create candidate-score pairs excluding winner
            pairs = list(zip(candidate_list, score_list))
            non_winners = [(c, s) for c, s in pairs if c != winner]
            if non_winners:
                # Get highest scoring non-winner
                best_non_winner = max(non_winners, key=lambda x: x[1])
                if margin < min_margin:
                    min_margin = margin
                    tightest_runner_up = best_non_winner[0]
                    tightest_score = best_non_winner[1]

    if tightest_runner_up != "Unknown":
        return (tightest_runner_up, tightest_score, min_margin)

    return default_result


# =============================================================================
# Enhanced Annotations Export
# =============================================================================


def export_enhanced_annotations(
    adata: "sc.AnnData",
    output_path: Path,
    marker_scores: Optional[pd.DataFrame] = None,
    decision_steps: Optional[pd.DataFrame] = None,
    cluster_annotations: Optional[pd.DataFrame] = None,
    max_iteration: int = 5,
    include_runner_up: bool = True,
    include_top_markers: bool = True,
    top_n_markers: int = 3,
    logger: Optional[logging.Logger] = None,
) -> pd.DataFrame:
    """Export enhanced cluster annotations with multi-level lineage tracking.

    Creates a comprehensive CSV with 44 columns tracking:
    - Cluster provenance (origin, iteration created)
    - Cell type assignments at each iteration level
    - Marker scores and confidence at each level
    - Runner-up labels and margins
    - Top expressing markers
    - Regional distribution

    Parameters
    ----------
    adata : sc.AnnData
        AnnData with cluster annotations
    output_path : Path
        Output CSV path
    marker_scores : pd.DataFrame, optional
        Marker scores DataFrame. Falls back to adata.uns if not provided.
    decision_steps : pd.DataFrame, optional
        Hierarchical decision steps for runner-up detection
    cluster_annotations : pd.DataFrame, optional
        Pre-computed cluster annotations. Falls back to adata.uns.
    max_iteration : int
        Maximum iteration level to include (default: 5)
    include_runner_up : bool
        Whether to include runner-up labels (default: True)
    include_top_markers : bool
        Whether to include top marker columns (default: True)
    top_n_markers : int
        Number of top markers to include (default: 3)
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    pd.DataFrame
        Enhanced annotations DataFrame with 44 columns
    """
    logger = logger or logging.getLogger(__name__)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Get marker scores from adata.uns if not provided
    if marker_scores is None:
        if "marker_scores_refined" in adata.uns:
            marker_scores = adata.uns["marker_scores_refined"]
            if isinstance(marker_scores, dict):
                marker_scores = pd.DataFrame(marker_scores)
        else:
            marker_scores = pd.DataFrame()

    # Determine cluster column (use highest available level)
    cluster_col = None
    for level in range(max_iteration, -1, -1):
        col = f"cluster_lvl{level}"
        if col in adata.obs:
            cluster_col = col
            break
    if cluster_col is None:
        cluster_col = "cluster_lvl0"

    # Get unique cluster IDs
    all_cluster_ids = sorted(
        adata.obs[cluster_col].astype(str).unique(),
        key=_sort_cluster_key
    )

    logger.info("Generating enhanced annotations for %d clusters", len(all_cluster_ids))

    # Build score lookup
    score_lookup = {}
    if not marker_scores.empty and "cluster_id" in marker_scores.columns:
        for cid in marker_scores["cluster_id"].unique():
            cid_str = str(cid)
            cluster_scores = marker_scores[marker_scores["cluster_id"].astype(str) == cid_str]
            if not cluster_scores.empty:
                # Get all scores for this cluster
                score_lookup[cid_str] = cluster_scores

    # Build records
    records = []
    total_cells = len(adata)

    for cluster_id in all_cluster_ids:
        cluster_id_str = str(cluster_id)
        mask = adata.obs[cluster_col].astype(str) == cluster_id_str
        n_cells = mask.sum()

        if n_cells == 0:
            continue

        # Determine origin cluster and iteration created
        if "cluster_lvl0" in adata.obs:
            origin_cluster = adata.obs.loc[mask, "cluster_lvl0"].astype(str).iloc[0]
        else:
            origin_cluster = cluster_id_str.split(":")[0] if ":" in cluster_id_str else cluster_id_str

        iteration_created = cluster_id_str.count(":")

        # Build record with base columns
        record = {
            "cluster_id": cluster_id_str,
            "origin_cluster": origin_cluster,
            "iteration_created": iteration_created,
            "n_cells": n_cells,
            "proportion": round(n_cells / total_cells, 6),
        }

        # Add cell type and score columns for each iteration level
        for level in range(max_iteration + 1):
            level_suffix = f"_lvl{level}"
            ct_col = f"cell_type_lvl{level}" if level > 0 else "cell_type_lvl0"
            cluster_lvl_col = f"cluster_lvl{level}"

            # Get cell type for this level
            cell_type = ""
            if ct_col in adata.obs:
                ct_values = adata.obs.loc[mask, ct_col].value_counts()
                if len(ct_values) > 0:
                    cell_type = str(ct_values.index[0])

            record[f"cell_type{level_suffix}"] = cell_type

            # Get score for this level (from marker_scores or adata.uns)
            score = 0.0
            assigned_path = ""
            stop_reason = ""
            confidence = ""
            coverage = 0.0
            is_ambiguous_root = False
            root_label = ""
            root_fail_reasons = ""
            min_margin = None
            decision_trace = ""

            # Look for hierarchical annotation data in adata.uns
            hier_key = f"hierarchical_annotations_lvl{level}"
            if hier_key in adata.uns:
                hier_df = adata.uns[hier_key]
                if isinstance(hier_df, pd.DataFrame) and not hier_df.empty:
                    row = hier_df[hier_df["cluster_id"].astype(str) == cluster_id_str]
                    if not row.empty:
                        row = row.iloc[0]
                        score = float(row.get("score", 0))
                        assigned_path = str(row.get("assigned_path", ""))
                        stop_reason = str(row.get("stop_reason", ""))
                        confidence = str(row.get("confidence_level", ""))
                        coverage = float(row.get("coverage", 0))
                        is_ambiguous_root = bool(row.get("is_ambiguous_root", False))
                        root_label = str(row.get("root_label", ""))
                        root_fail_reasons = str(row.get("root_fail_reasons", ""))
                        min_margin = row.get("min_margin_along_path")
                        decision_trace = str(row.get("decision_trace", ""))

            # Fallback to marker_scores if no hierarchical data
            if score == 0.0 and cluster_id_str in score_lookup:
                cluster_scores = score_lookup[cluster_id_str]
                if cell_type and not cluster_scores.empty:
                    label_scores = cluster_scores[cluster_scores["label"] == cell_type]
                    if not label_scores.empty:
                        score = float(label_scores.iloc[0].get("score", 0))
                        coverage = float(label_scores.iloc[0].get("coverage", 0))

            record[f"score{level_suffix}"] = round(score, 3)
            record[f"assigned_path{level_suffix}"] = assigned_path
            record[f"stop_reason{level_suffix}"] = stop_reason
            record[f"confidence{level_suffix}"] = confidence or _classify_confidence_band(score)
            record[f"coverage{level_suffix}"] = round(coverage, 3)
            record[f"is_ambiguous_root{level_suffix}"] = is_ambiguous_root
            record[f"root_label{level_suffix}"] = root_label
            record[f"root_fail_reasons{level_suffix}"] = root_fail_reasons
            if min_margin is not None and not pd.isna(min_margin):
                record[f"min_margin{level_suffix}"] = round(float(min_margin), 3)
            else:
                record[f"min_margin{level_suffix}"] = None
            record[f"decision_trace{level_suffix}"] = decision_trace

        # Add runner-up information (based on highest level)
        if include_runner_up and cluster_annotations is not None:
            runner_up, runner_up_score, margin = get_hierarchical_runner_up(
                cluster_id_str, decision_steps, cluster_annotations, marker_scores
            )
            record["runner_up_label"] = runner_up
            record["runner_up_score"] = round(runner_up_score, 3)
            record["margin_to_runner_up"] = round(margin, 3) if margin != float("inf") else None

        # Add top markers
        if include_top_markers and cluster_id_str in score_lookup:
            cluster_scores = score_lookup[cluster_id_str]
            if not cluster_scores.empty and "label" in cluster_scores.columns:
                top_scores = cluster_scores.nlargest(top_n_markers, "score")
                for i, (_, row) in enumerate(top_scores.iterrows(), 1):
                    record[f"top_marker_{i}"] = str(row.get("label", ""))
                    record[f"top_marker_{i}_score"] = round(float(row.get("score", 0)), 3)

        # Add regional distribution (if region column exists)
        if "region" in adata.obs:
            region_counts = adata.obs.loc[mask, "region"].value_counts()
            for region, count in region_counts.items():
                record[f"region_{region}_count"] = count
                record[f"region_{region}_pct"] = round(100 * count / n_cells, 1)

        # Add mean enrichment and positive fraction
        if cluster_id_str in score_lookup:
            cluster_scores = score_lookup[cluster_id_str]
            # Find the assigned label for this cluster
            assigned_label = record.get("cell_type_lvl1", record.get("cell_type_lvl0", ""))
            if assigned_label and not cluster_scores.empty:
                label_scores = cluster_scores[cluster_scores["label"] == assigned_label]
                if not label_scores.empty:
                    record["mean_enrichment"] = round(
                        float(label_scores.iloc[0].get("mean_enrichment", 0)), 3
                    )
                    record["mean_positive_fraction"] = round(
                        float(label_scores.iloc[0].get("mean_positive_fraction", 0)), 3
                    )

        records.append(record)

    df = pd.DataFrame.from_records(records)
    df.to_csv(output_path, index=False)

    logger.info("Exported enhanced annotations: %d clusters, %d columns → %s",
                len(df), len(df.columns), output_path)
    return df


# =============================================================================
# Stage H Format Export
# =============================================================================


def export_cluster_annotations_stage_h_format(
    adata: "sc.AnnData",
    enhanced_annotations: pd.DataFrame,
    output_path: Path,
    iteration: int = 1,
    logger: Optional[logging.Logger] = None,
) -> pd.DataFrame:
    """Export cluster annotations in Stage H format with 24 columns.

    Provides backward-compatible output matching the original Stage H format,
    derived from enhanced annotations at a specific iteration level.

    Parameters
    ----------
    adata : sc.AnnData
        AnnData with cluster annotations
    enhanced_annotations : pd.DataFrame
        Enhanced annotations DataFrame from export_enhanced_annotations()
    output_path : Path
        Output CSV path
    iteration : int
        Iteration level to export (default: 1)
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    pd.DataFrame
        Stage H format annotations with 24 columns:
        - cluster_id, origin_cluster, iteration_created, proportion
        - mean_enrichment, mean_positive_fraction, n_cells
        - assigned_label, assigned_path, assigned_level, assigned_score
        - root_label, confidence, min_margin_along_path, margin_is_infinite
        - stop_reason, stopped_before_leaf, composition, decision_trace
        - coverage, resolved_markers, is_ambiguous_root
        - ambiguous_root_candidates, root_fail_reasons, confidence_level
    """
    logger = logger or logging.getLogger(__name__)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if enhanced_annotations.empty:
        logger.warning("Empty enhanced annotations; creating empty Stage H format file")
        df = pd.DataFrame()
        df.to_csv(output_path, index=False)
        return df

    records = []
    for _, row in enhanced_annotations.iterrows():
        cluster_id = str(row.get("cluster_id", ""))

        # Get values from the specified iteration level
        cell_type = row.get(f"cell_type_lvl{iteration}", row.get("cell_type_lvl1", ""))
        score = row.get(f"score_lvl{iteration}", row.get("score_lvl1", 0))
        path = row.get(f"assigned_path_lvl{iteration}", row.get("assigned_path_lvl1", ""))
        level = row.get(f"assigned_level_lvl{iteration}", row.get("assigned_level_lvl1", -1))

        # Derive is_ambiguous_root
        is_ambiguous = row.get(f"is_ambiguous_root_lvl{iteration}", False)
        if not is_ambiguous and "~" in str(cell_type):
            is_ambiguous = True

        # Derive root_label
        root_label = row.get(f"root_label_lvl{iteration}", "")
        if not root_label:
            if "~" in str(cell_type):
                root_label = cell_type
            elif "/" in str(path):
                root_label = str(path).split("/")[0].strip()
            else:
                root_label = cell_type

        records.append({
            # Provenance columns
            "cluster_id": cluster_id,
            "origin_cluster": row.get("origin_cluster", str(cluster_id).split(":")[0]),
            "iteration_created": row.get("iteration_created", str(cluster_id).count(":")),
            "proportion": row.get("proportion", 0),
            "mean_enrichment": row.get("mean_enrichment", 0),
            "mean_positive_fraction": row.get("mean_positive_fraction", 0),

            # Stage H format columns
            "n_cells": row.get("n_cells", 0),
            "assigned_label": cell_type,
            "assigned_path": path if path else cell_type,
            "assigned_level": int(level) if level is not None and level != "" else -1,
            "assigned_score": round(float(score) if score else 0.0, 3),
            "root_label": root_label,
            "confidence": row.get(f"confidence_lvl{iteration}", score),
            "min_margin_along_path": row.get(f"min_margin_lvl{iteration}"),
            "margin_is_infinite": row.get(f"margin_is_infinite_lvl{iteration}", False),
            "stop_reason": row.get(f"stop_reason_lvl{iteration}", ""),
            "stopped_before_leaf": row.get(f"stopped_before_leaf_lvl{iteration}", True),
            "composition": row.get(f"composition_lvl{iteration}", ""),
            "decision_trace": row.get(f"decision_trace_lvl{iteration}", ""),
            "coverage": row.get(f"coverage_lvl{iteration}", 0),
            "resolved_markers": row.get(f"resolved_markers_lvl{iteration}", ""),
            "is_ambiguous_root": is_ambiguous,
            "ambiguous_root_candidates": row.get(f"ambiguous_root_candidates_lvl{iteration}", ""),
            "root_fail_reasons": row.get(f"root_fail_reasons_lvl{iteration}", ""),
            "confidence_level": row.get("confidence_level", _classify_confidence_band(float(score) if score else 0)),
        })

    df = pd.DataFrame.from_records(records)
    df.to_csv(output_path, index=False)

    logger.info("Exported Stage H format annotations: %d clusters → %s", len(df), output_path)
    return df


# =============================================================================
# Marker Scores Export
# =============================================================================


def export_marker_scores(
    adata: "sc.AnnData",
    output_path: Path,
    logger: Optional[logging.Logger] = None,
) -> Path:
    """Export marker_scores_refined from adata.uns to CSV.

    Parameters
    ----------
    adata : sc.AnnData
        AnnData with marker scores in uns["marker_scores_refined"]
    output_path : Path
        Output CSV path
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    Path
        Path to the exported CSV file
    """
    logger = logger or logging.getLogger(__name__)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if "marker_scores_refined" in adata.uns:
        scores = adata.uns["marker_scores_refined"]
        if isinstance(scores, pd.DataFrame):
            scores.to_csv(output_path, index=False)
            logger.info("Exported marker scores: %d rows → %s", len(scores), output_path)
        elif isinstance(scores, (dict, list)):
            df = pd.DataFrame(scores)
            df.to_csv(output_path, index=False)
            logger.info("Exported marker scores: %d rows → %s", len(df), output_path)
        else:
            logger.warning("marker_scores_refined has unexpected type: %s", type(scores))
            pd.DataFrame().to_csv(output_path, index=False)
    else:
        logger.warning("No marker_scores_refined in adata.uns; created empty file")
        pd.DataFrame().to_csv(output_path, index=False)

    return output_path


# =============================================================================
# Cluster Label Mapping Export
# =============================================================================


def export_cluster_label_mapping(
    adata: "sc.AnnData",
    output_path: Path,
    marker_scores: Optional[pd.DataFrame] = None,
    cluster_col: str = "cluster_lvl1",
    label_col: str = "cell_type_lvl1",
    logger: Optional[logging.Logger] = None,
) -> pd.DataFrame:
    """Export cluster to label mapping with scores.

    Creates a simple mapping table with cluster_id, assigned_label, assigned_score,
    and n_cells. Sets assigned_score=0.0 for Unassigned clusters.

    Parameters
    ----------
    adata : sc.AnnData
        AnnData with cluster annotations
    output_path : Path
        Output CSV path
    marker_scores : pd.DataFrame, optional
        Marker scores for score lookup
    cluster_col : str
        Cluster column name (default: "cluster_lvl1")
    label_col : str
        Label column name (default: "cell_type_lvl1")
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    pd.DataFrame
        Cluster label mapping with columns: cluster_id, assigned_label, assigned_score, n_cells
    """
    logger = logger or logging.getLogger(__name__)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Get marker scores from adata.uns if not provided
    if marker_scores is None:
        if "marker_scores_refined" in adata.uns:
            marker_scores = adata.uns["marker_scores_refined"]
            if isinstance(marker_scores, dict):
                marker_scores = pd.DataFrame(marker_scores)
        else:
            marker_scores = pd.DataFrame()

    # Build score lookup (max score per cluster)
    score_lookup = {}
    if not marker_scores.empty and "cluster_id" in marker_scores.columns and "score" in marker_scores.columns:
        best_scores = marker_scores.loc[marker_scores.groupby("cluster_id")["score"].idxmax()]
        for _, row in best_scores.iterrows():
            score_lookup[str(row["cluster_id"])] = float(row["score"])

    # Get cluster-label mapping from adata.obs
    if cluster_col not in adata.obs or label_col not in adata.obs:
        logger.warning("Required columns not found: %s or %s", cluster_col, label_col)
        df = pd.DataFrame()
        df.to_csv(output_path, index=False)
        return df

    cluster_label_df = (
        adata.obs.groupby(cluster_col)
        .agg({label_col: "first"})
        .reset_index()
        .rename(columns={cluster_col: "cluster_id", label_col: "assigned_label"})
    )

    # Add n_cells
    cell_counts = adata.obs[cluster_col].value_counts()
    cluster_label_df["n_cells"] = cluster_label_df["cluster_id"].map(cell_counts)

    # Add assigned_score
    cluster_label_df["assigned_score"] = cluster_label_df["cluster_id"].astype(str).map(
        lambda x: score_lookup.get(x, 0.0)
    )

    # Set score=0.0 for Unassigned clusters (no root gate passed = no confidence)
    # This matches the reference pipeline semantics: assigned_score indicates
    # confidence in the assigned label, not best possible match potential
    unassigned_mask = cluster_label_df["assigned_label"] == "Unassigned"
    cluster_label_df.loc[unassigned_mask, "assigned_score"] = 0.0

    # Round scores
    cluster_label_df["assigned_score"] = cluster_label_df["assigned_score"].round(3)

    # Sort by cluster_id
    cluster_label_df = cluster_label_df.sort_values(
        "cluster_id",
        key=lambda x: x.map(lambda c: _sort_cluster_key(str(c)))
    ).reset_index(drop=True)

    # Reorder columns
    cluster_label_df = cluster_label_df[["cluster_id", "assigned_label", "assigned_score", "n_cells"]]

    cluster_label_df.to_csv(output_path, index=False)
    logger.info("Exported cluster label mapping: %d clusters → %s", len(cluster_label_df), output_path)

    return cluster_label_df


# =============================================================================
# High-Level Export Function
# =============================================================================


def run_annotation_exports(
    adata: "sc.AnnData",
    output_dir: Path,
    marker_scores: Optional[pd.DataFrame] = None,
    decision_steps: Optional[pd.DataFrame] = None,
    cluster_annotations: Optional[pd.DataFrame] = None,
    iteration: int = 1,
    cluster_col: str = "cluster_lvl1",
    label_col: str = "cell_type_lvl1",
    logger: Optional[logging.Logger] = None,
) -> Dict[str, Path]:
    """Run all annotation exports.

    Exports enhanced annotations, Stage H format, marker scores, and cluster mapping.

    Parameters
    ----------
    adata : sc.AnnData
        AnnData with annotations
    output_dir : Path
        Output directory
    marker_scores : pd.DataFrame, optional
        Marker scores DataFrame
    decision_steps : pd.DataFrame, optional
        Decision steps for runner-up detection
    cluster_annotations : pd.DataFrame, optional
        Pre-computed cluster annotations
    iteration : int
        Iteration level for Stage H format (default: 1)
    cluster_col : str
        Cluster column name
    label_col : str
        Label column name
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    Dict[str, Path]
        Mapping from export type to output path
    """
    logger = logger or logging.getLogger(__name__)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    output_paths = {}

    # Export enhanced annotations
    enhanced_path = output_dir / "cluster_annotations_enhanced.csv"
    enhanced_df = export_enhanced_annotations(
        adata=adata,
        output_path=enhanced_path,
        marker_scores=marker_scores,
        decision_steps=decision_steps,
        cluster_annotations=cluster_annotations,
        max_iteration=5,
        logger=logger,
    )
    output_paths["enhanced"] = enhanced_path

    # Export Stage H format
    stage_h_path = output_dir / "cluster_annotations.csv"
    export_cluster_annotations_stage_h_format(
        adata=adata,
        enhanced_annotations=enhanced_df,
        output_path=stage_h_path,
        iteration=iteration,
        logger=logger,
    )
    output_paths["stage_h"] = stage_h_path

    # Export marker scores
    marker_scores_path = output_dir / "marker_scores.csv"
    export_marker_scores(adata, marker_scores_path, logger)
    output_paths["marker_scores"] = marker_scores_path

    # Export cluster label mapping
    mapping_path = output_dir / "cluster_label_mapping.csv"
    export_cluster_label_mapping(
        adata=adata,
        output_path=mapping_path,
        marker_scores=marker_scores,
        cluster_col=cluster_col,
        label_col=label_col,
        logger=logger,
    )
    output_paths["mapping"] = mapping_path

    logger.info("Completed annotation exports: %d files → %s", len(output_paths), output_dir)
    return output_paths


# =============================================================================
# Composition Statistics Export
# =============================================================================


def export_composition_stats(
    adata: "sc.AnnData",
    output_dir: Path,
    cell_type_col: str,
    logger: Optional[logging.Logger] = None,
) -> Dict[str, Path]:
    """Export composition statistics (global, by sample, by region, by donor).

    Parameters
    ----------
    adata : sc.AnnData
        AnnData with cell type annotations
    output_dir : Path
        Output directory
    cell_type_col : str
        Column name containing cell type labels
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    Dict[str, Path]
        Mapping from export type to output path
    """
    logger = logger or logging.getLogger(__name__)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_paths = {}

    sample_col = 'sample_id' if 'sample_id' in adata.obs.columns else None
    region_col = 'region' if 'region' in adata.obs.columns else None
    donor_col = 'donor' if 'donor' in adata.obs.columns else None

    # Global composition (matching reference format: proportion as %, total_cells column)
    logger.info('Generating global composition...')
    total_cells = len(adata)
    comp_global = adata.obs[cell_type_col].value_counts().reset_index()
    comp_global.columns = ['cell_type', 'n_cells']
    comp_global['proportion'] = round(comp_global['n_cells'] / total_cells * 100, 2)
    comp_global['total_cells'] = total_cells
    comp_global = comp_global.sort_values('n_cells', ascending=False)
    global_path = output_dir / 'composition_global_refined.csv'
    comp_global.to_csv(global_path, index=False)
    output_paths['global'] = global_path
    logger.info(f'  Wrote composition_global_refined.csv ({len(comp_global)} cell types)')

    # Composition by sample (matching reference format exactly)
    if sample_col:
        logger.info('Generating composition by sample...')
        comp_sample = adata.obs.groupby([sample_col, cell_type_col]).size().unstack(fill_value=0)

        # Long format with reference column names
        comp_sample_long = comp_sample.stack().reset_index()
        comp_sample_long.columns = ['sample_id', 'cell_type', 'n_cells']
        # Add sample totals
        sample_totals = comp_sample_long.groupby('sample_id')['n_cells'].transform('sum')
        comp_sample_long['proportion'] = round(comp_sample_long['n_cells'] / sample_totals * 100, 2)
        comp_sample_long['sample_id_total'] = sample_totals
        # Sort by sample and n_cells descending
        comp_sample_long = comp_sample_long.sort_values(
            ['sample_id', 'n_cells'], ascending=[True, False]
        )
        sample_path = output_dir / 'composition_by_sample_refined.csv'
        comp_sample_long.to_csv(sample_path, index=False)
        output_paths['by_sample'] = sample_path
        logger.info(f'  Wrote composition_by_sample_refined.csv ({len(comp_sample_long)} rows)')

    # Composition by region (matching reference format exactly)
    if region_col:
        logger.info('Generating composition by region...')
        # Long format for internal calculations
        comp_region_long = adata.obs.groupby([region_col, cell_type_col]).size().reset_index(name='n_cells')
        comp_region_long.columns = ['region', 'cell_type', 'n_cells']
        # Add region totals
        region_totals = comp_region_long.groupby('region')['n_cells'].transform('sum')
        comp_region_long['proportion'] = round(
            comp_region_long['n_cells'] / region_totals * 100, 2
        )
        comp_region_long['region_total'] = region_totals
        # Sort by region and n_cells descending
        comp_region_long = comp_region_long.sort_values(
            ['region', 'n_cells'], ascending=[True, False]
        )
        region_path = output_dir / 'composition_by_region_refined.csv'
        comp_region_long.to_csv(region_path, index=False)
        output_paths['by_region'] = region_path
        logger.info(f'  Wrote composition_by_region_refined.csv ({len(comp_region_long)} rows)')

        # === REFERENCE FORMAT PIVOT TABLES ===
        # Use reference region order (anatomical: proximal to distal)
        region_order = ["isthmus", "ampulla", "infundibulum", "fimbriae"]
        available_regions = comp_region_long["region"].unique().tolist()
        region_order = [r for r in region_order if r in available_regions]
        # Add any regions not in the default order
        for r in available_regions:
            if r not in region_order:
                region_order.append(r)

        # === TABLE 1: Raw Counts (cell types as rows, regions as columns) ===
        counts_pivot = comp_region_long.pivot_table(
            index="cell_type",
            columns="region",
            values="n_cells",
            aggfunc="sum"
        ).fillna(0).astype(int)
        # Reorder columns
        counts_pivot = counts_pivot[[r for r in region_order if r in counts_pivot.columns]]
        # Add Total column
        counts_pivot["Total"] = counts_pivot.sum(axis=1)
        # Sort by total descending
        counts_pivot = counts_pivot.sort_values("Total", ascending=False)
        # Add TOTAL row
        col_totals = counts_pivot.sum(axis=0)
        col_totals.name = "TOTAL"
        counts_pivot = pd.concat([counts_pivot, col_totals.to_frame().T])

        counts_path = output_dir / 'composition_counts_by_region.csv'
        counts_pivot.to_csv(counts_path)
        output_paths['counts_by_region'] = counts_path
        logger.info('  Wrote composition_counts_by_region.csv (reference format)')

        # === TABLE 2: % Across Regions (row sums to 100%) ===
        across_pivot = comp_region_long.pivot_table(
            index="cell_type",
            columns="region",
            values="n_cells",
            aggfunc="sum"
        ).fillna(0)
        across_pivot = across_pivot[[r for r in region_order if r in across_pivot.columns]]
        row_totals = across_pivot.sum(axis=1)
        across_pct = across_pivot.div(row_totals, axis=0) * 100
        # Sort by total cells descending
        across_pct["_total"] = row_totals
        across_pct = across_pct.sort_values("_total", ascending=False)
        across_pct = across_pct.drop("_total", axis=1)
        across_pct = across_pct.round(1)

        across_path = output_dir / 'composition_pct_across_regions.csv'
        across_pct.to_csv(across_path)
        output_paths['pct_across_regions'] = across_path
        logger.info('  Wrote composition_pct_across_regions.csv (rows sum to 100%)')

        # === TABLE 3: % Within Regions (column sums to 100%) ===
        within_pivot = comp_region_long.pivot_table(
            index="cell_type",
            columns="region",
            values="n_cells",
            aggfunc="sum"
        ).fillna(0)
        within_pivot = within_pivot[[r for r in region_order if r in within_pivot.columns]]
        col_totals_within = within_pivot.sum(axis=0)
        within_pct = within_pivot.div(col_totals_within, axis=1) * 100
        # Sort by average proportion descending
        within_pct["_avg"] = within_pct.mean(axis=1)
        within_pct = within_pct.sort_values("_avg", ascending=False)
        within_pct = within_pct.drop("_avg", axis=1)
        within_pct = within_pct.round(2)

        within_path = output_dir / 'composition_pct_within_regions.csv'
        within_pct.to_csv(within_path)
        output_paths['pct_within_regions'] = within_path
        logger.info('  Wrote composition_pct_within_regions.csv (columns sum to 100%)')

    # Composition by donor (matching reference format exactly)
    if donor_col:
        logger.info('Generating composition by donor...')
        comp_donor = adata.obs.groupby([donor_col, cell_type_col]).size().unstack(fill_value=0)
        comp_donor_long = comp_donor.stack().reset_index()
        comp_donor_long.columns = ['donor', 'cell_type', 'n_cells']
        # Add donor totals
        donor_totals = comp_donor_long.groupby('donor')['n_cells'].transform('sum')
        comp_donor_long['proportion'] = round(
            comp_donor_long['n_cells'] / donor_totals * 100, 2
        )
        comp_donor_long['donor_total'] = donor_totals
        # Sort by donor and n_cells descending
        comp_donor_long = comp_donor_long.sort_values(
            ['donor', 'n_cells'], ascending=[True, False]
        )
        donor_path = output_dir / 'composition_by_donor_refined.csv'
        comp_donor_long.to_csv(donor_path, index=False)
        output_paths['by_donor'] = donor_path
        logger.info(f'  Wrote composition_by_donor_refined.csv ({len(comp_donor_long)} rows)')

    return output_paths


# =============================================================================
# Marker Scores Export from AnnData
# =============================================================================


def export_marker_scores_from_adata(
    adata: "sc.AnnData",
    output_dir: Path,
    logger: Optional[logging.Logger] = None,
) -> Optional[pd.DataFrame]:
    """Export marker scores from adata.uns.

    Parameters
    ----------
    adata : sc.AnnData
        AnnData with marker scores in uns
    output_dir : Path
        Output directory
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    pd.DataFrame or None
        Marker scores DataFrame if found, else None
    """
    logger = logger or logging.getLogger(__name__)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info('Extracting marker scores...')

    marker_scores_key = None
    for key in ['marker_scores_refined', 'marker_scores', 'scores']:
        if key in adata.uns:
            marker_scores_key = key
            break

    if not marker_scores_key:
        logger.warning('  No marker scores found in adata.uns')
        return None

    scores_data = adata.uns[marker_scores_key]
    if isinstance(scores_data, pd.DataFrame):
        scores_df = scores_data
    elif isinstance(scores_data, list):
        scores_df = pd.DataFrame(scores_data)
    elif isinstance(scores_data, dict):
        # Convert dict format to DataFrame
        records = []
        for cluster_id, labels in scores_data.items():
            if isinstance(labels, dict):
                for label, score in labels.items():
                    records.append({
                        'cluster_id': cluster_id,
                        'label': label,
                        'score': score,
                    })
        scores_df = pd.DataFrame(records)
    else:
        logger.warning(f'  Unexpected marker scores format: {type(scores_data)}')
        return None

    if len(scores_df) > 0:
        scores_path = output_dir / 'marker_scores.csv'
        scores_df.to_csv(scores_path, index=False)
        logger.info(f'  Wrote marker_scores.csv ({len(scores_df)} records)')
        return scores_df

    return None


# =============================================================================
# Simple Cluster Annotations Export
# =============================================================================


def export_cluster_annotations_simple(
    adata: "sc.AnnData",
    output_dir: Path,
    cell_type_col: str,
    cluster_col: str,
    scores_df: Optional[pd.DataFrame],
    logger: Optional[logging.Logger] = None,
) -> None:
    """Export simple cluster annotations with label mapping.

    Parameters
    ----------
    adata : sc.AnnData
        AnnData with cluster annotations
    output_dir : Path
        Output directory
    cell_type_col : str
        Column name containing cell type labels
    cluster_col : str
        Column name containing cluster IDs
    scores_df : pd.DataFrame, optional
        Marker scores DataFrame
    logger : logging.Logger, optional
        Logger instance
    """
    logger = logger or logging.getLogger(__name__)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info('Generating cluster annotations...')

    if cluster_col not in adata.obs.columns:
        logger.warning(f'  Cluster column {cluster_col} not found')
        return

    # Get unique cluster -> cell type mapping
    cluster_mapping = adata.obs.groupby(cluster_col)[cell_type_col].agg(
        lambda x: x.value_counts().index[0] if len(x) > 0 else 'Unknown'
    ).reset_index()
    cluster_mapping.columns = ['cluster_id', 'assigned_label']

    # Add cell counts
    cluster_counts = adata.obs[cluster_col].value_counts().reset_index()
    cluster_counts.columns = ['cluster_id', 'n_cells']
    cluster_mapping = cluster_mapping.merge(cluster_counts, on='cluster_id')

    # Add scores if available
    # IMPORTANT: Use the ASSIGNED LABEL's score, not the highest-scoring label
    # (Fixed 2025-12-30: was incorrectly using idxmax which gave highest score)
    if (scores_df is not None and
            'cluster_id' in scores_df.columns and
            'score' in scores_df.columns):
        scores_df_copy = scores_df.copy()
        scores_df_copy['cluster_id'] = scores_df_copy['cluster_id'].astype(str)
        scores_df_copy['label'] = scores_df_copy['label'].astype(str)
        cluster_mapping['cluster_id'] = cluster_mapping['cluster_id'].astype(str)

        # Build lookup: (cluster_id, label) -> score
        score_lookup = {}
        for _, row in scores_df_copy.iterrows():
            key = (row['cluster_id'], row['label'])
            score_lookup[key] = row['score']

        # Build fallback: cluster_id -> highest score
        best_idx = scores_df_copy.groupby('cluster_id')['score'].idxmax()
        fallback_lookup = dict(zip(
            scores_df_copy.loc[best_idx, 'cluster_id'],
            scores_df_copy.loc[best_idx, 'score']
        ))

        # Get score for assigned_label, fall back to highest score
        def get_score(row):
            key = (str(row['cluster_id']), str(row['assigned_label']))
            if key in score_lookup:
                return score_lookup[key]
            return fallback_lookup.get(str(row['cluster_id']), 0.0)

        cluster_mapping['assigned_score'] = cluster_mapping.apply(get_score, axis=1)

        # Set score=0.0 for Unassigned clusters
        unassigned_mask = cluster_mapping['assigned_label'] == 'Unassigned'
        cluster_mapping.loc[unassigned_mask, 'assigned_score'] = 0.0

    # Save cluster annotations (simple format)
    mapping_path = output_dir / 'cluster_annotations.csv'
    cluster_mapping.to_csv(mapping_path, index=False)
    logger.info(f'  Wrote cluster_annotations.csv ({len(cluster_mapping)} clusters)')

    # Enhanced annotations with regional distribution
    region_col = 'region' if 'region' in adata.obs.columns else None
    enhanced = cluster_mapping.copy()

    if region_col:
        # Add regional distribution
        region_dist = adata.obs.groupby([cluster_col, region_col]).size().unstack(fill_value=0)
        region_dist_pct = region_dist.div(region_dist.sum(axis=1), axis=0) * 100
        region_cols = [f'pct_{r}' for r in region_dist.columns]
        region_dist_pct.columns = region_cols
        region_dist_pct = region_dist_pct.reset_index()
        region_dist_pct[cluster_col] = region_dist_pct[cluster_col].astype(str)
        enhanced = enhanced.merge(
            region_dist_pct,
            left_on='cluster_id',
            right_on=cluster_col,
            how='left',
        )
        if cluster_col != 'cluster_id':
            enhanced = enhanced.drop(columns=[cluster_col])

    enhanced_path = output_dir / 'cluster_annotations_enhanced.csv'
    enhanced.to_csv(enhanced_path, index=False)
    logger.info('  Wrote cluster_annotations_enhanced.csv')


# =============================================================================
# Review Summary Export
# =============================================================================


def export_review_summary(
    adata: "sc.AnnData",
    output_dir: Path,
    cell_type_col: str,
    cluster_col: str,
    comp_global: pd.DataFrame,
    logger: Optional[logging.Logger] = None,
) -> Path:
    """Export review summary JSON with lineage groupings.

    Parameters
    ----------
    adata : sc.AnnData
        AnnData with cell type annotations
    output_dir : Path
        Output directory
    cell_type_col : str
        Column name containing cell type labels
    cluster_col : str
        Column name containing cluster IDs
    comp_global : pd.DataFrame
        Global composition DataFrame
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    Path
        Path to the exported JSON file
    """
    logger = logger or logging.getLogger(__name__)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info('Generating review summary...')

    # Group composition by lineage
    composition_by_group: Dict[str, Dict[str, int]] = {}
    other_types: Dict[str, int] = {}

    for cell_type, count in adata.obs[cell_type_col].value_counts().items():
        assigned = False
        for group, members in LINEAGE_GROUPS.items():
            # Check if cell_type matches any member (case-insensitive, partial match)
            for member in members:
                if (member.lower() in cell_type.lower() or
                        cell_type.lower() in member.lower()):
                    if group not in composition_by_group:
                        composition_by_group[group] = {}
                    composition_by_group[group][cell_type] = int(count)
                    assigned = True
                    break
            if assigned:
                break

        if not assigned:
            # Check for hybrid types (contain ~)
            if '~' in cell_type:
                parts = cell_type.split('~')
                # Assign to first lineage found
                for part in parts:
                    for group, members in LINEAGE_GROUPS.items():
                        for member in members:
                            if member.lower() in part.lower():
                                if group not in composition_by_group:
                                    composition_by_group[group] = {}
                                composition_by_group[group][cell_type] = int(count)
                                assigned = True
                                break
                        if assigned:
                            break
                    if assigned:
                        break

            if not assigned:
                other_types[cell_type] = int(count)

    if other_types:
        composition_by_group['Other'] = other_types

    review_summary = {
        'generated_at': datetime.now().isoformat(),
        'total_cells': int(adata.n_obs),
        'composition_by_group': composition_by_group,
        'cell_type_column': cell_type_col,
        'cluster_column': cluster_col,
        'n_cell_types': len(comp_global),
        'n_clusters': (
            adata.obs[cluster_col].nunique()
            if cluster_col in adata.obs.columns else 0
        ),
    }

    summary_path = output_dir / 'review_summary.json'
    with open(summary_path, 'w') as f:
        json.dump(review_summary, f, indent=2)
    logger.info('  Wrote review_summary.json')

    return summary_path


# =============================================================================
# Workflow State Export
# =============================================================================


def export_workflow_state(
    output_dir: Path,
    marker_map_path: Optional[str] = None,
    input_path: Optional[str] = None,
    stage_h_dir: Optional[str] = None,
    logger: Optional[logging.Logger] = None,
) -> Path:
    """Export workflow state YAML for tracking.

    Parameters
    ----------
    output_dir : Path
        Output directory
    marker_map_path : str, optional
        Path to marker map
    input_path : str, optional
        Path to input h5ad
    stage_h_dir : str, optional
        Path to Stage H directory
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    Path
        Path to the exported YAML file
    """
    import yaml

    logger = logger or logging.getLogger(__name__)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    workflow_state = {
        'version': '1.0',
        'marker_map_path': marker_map_path or '',
        'input_path': input_path or '',
        'stage_h_dir': stage_h_dir or '',
        'out_dir': str(output_dir),
        'created_at': datetime.now().isoformat(),
        'updated_at': datetime.now().isoformat(),
        'review_complete': True,
        'review_timestamp': datetime.now().isoformat(),
    }

    state_path = output_dir / 'workflow_state.yaml'
    with open(state_path, 'w') as f:
        yaml.dump(workflow_state, f, default_flow_style=False)
    logger.info('  Wrote workflow_state.yaml')

    return state_path


# =============================================================================
# 44-Column Enhanced Annotations Format
# =============================================================================


def _map_col_name_44(col: str) -> str:
    """Map enhanced column name to cluster_ann column name."""
    mapping = {
        'cell_type': 'assigned_label',
        'score': 'assigned_score',
        'assigned_path': 'assigned_path',
        'assigned_level': 'assigned_level',
        'root_label': 'root_label',
        'confidence': 'confidence',
        'min_margin_along_path': 'min_margin_along_path',
        'margin_is_infinite': 'margin_is_infinite',
        'stop_reason': 'stop_reason',
        'stopped_before_leaf': 'stopped_before_leaf',
        'composition': 'composition',
        'decision_trace': 'decision_trace',
        'coverage': 'coverage',
        'resolved_markers': 'resolved_markers',
        'is_ambiguous_root': 'is_ambiguous_root',
        'ambiguous_root_candidates': 'ambiguous_root_candidates',
    }
    return mapping.get(col, col)


def _build_reason_string(row: pd.Series) -> str:
    """Build reason string from row data."""
    label = row.get('assigned_label', '')
    score = row.get('assigned_score', 0)
    margin = row.get('min_margin_along_path', 0)

    # Handle None/NaN values
    if score is None or (isinstance(score, float) and pd.isna(score)):
        score = 0
    if margin is None or (isinstance(margin, float) and pd.isna(margin)):
        margin = 0

    if not label or label == 'Unassigned':
        return ''

    # Check if it's a path (contains ' / ')
    path = row.get('assigned_path', '')
    if ' / ' in str(path):
        parts = str(path).split(' / ')
        if len(parts) >= 2:
            return f"{parts[0]} -> {label} (score={score:.2f}, margin={margin:.2f})"

    return f"{label} (score={score:.2f})"


def _get_runner_up_44(
    cluster_id: str,
    assigned_label: str,
    scores_df: Optional[pd.DataFrame],
) -> Tuple[str, float, float]:
    """Get runner-up label, score, and gap for a cluster.

    Parameters
    ----------
    cluster_id : str
        Cluster ID
    assigned_label : str
        Assigned cell type label
    scores_df : pd.DataFrame, optional
        Marker scores DataFrame

    Returns
    -------
    Tuple[str, float, float]
        (runner_up_label, runner_up_score, gap)
    """
    if scores_df is None or scores_df.empty:
        return '', 0.0, 0.0

    if 'cluster_id' not in scores_df.columns or 'score' not in scores_df.columns:
        return '', 0.0, 0.0

    # Filter to this cluster
    cluster_scores = scores_df[
        scores_df['cluster_id'].astype(str) == str(cluster_id)
    ].copy()
    if cluster_scores.empty:
        return '', 0.0, 0.0

    # Sort by score descending
    cluster_scores = cluster_scores.sort_values('score', ascending=False)

    # Get top two scores
    if len(cluster_scores) < 2:
        return '', 0.0, 0.0

    top_row = cluster_scores.iloc[0]
    second_row = cluster_scores.iloc[1]

    top_label = top_row.get('label', '')
    top_score = float(top_row.get('score', 0))
    second_label = second_row.get('label', '')
    second_score = float(second_row.get('score', 0))

    # Runner-up is the second-best label
    # Gap is the difference between assigned and runner-up
    if str(top_label) == str(assigned_label):
        runner_up_label = second_label
        runner_up_score = second_score
        gap = top_score - second_score
    else:
        # Assigned might not be top score (e.g., gating chose a lower-scoring subtype)
        runner_up_label = top_label
        runner_up_score = top_score
        # Find assigned score
        assigned_rows = cluster_scores[cluster_scores['label'] == assigned_label]
        if not assigned_rows.empty:
            assigned_score = float(assigned_rows.iloc[0].get('score', 0))
            gap = assigned_score - top_score
        else:
            gap = 0.0

    return runner_up_label, runner_up_score, gap


def generate_enhanced_annotations_44col(
    cluster_ann: pd.DataFrame,
    scores_df: Optional[pd.DataFrame],
    stage_h_dir: Optional[str] = None,
    iteration: int = 1,
    logger: Optional[logging.Logger] = None,
) -> pd.DataFrame:
    """Generate 44-column enhanced annotations matching reference format.

    Creates hierarchical annotations with:
    - 6 basic columns: cluster_id, origin_cluster, iteration_created, n_cells, proportion, reason
    - 16 lvl0 columns: Stage H annotations (or empty if unavailable)
    - 16 lvl1 columns: Current iteration annotations
    - 3 quality metrics: mean_enrichment, mean_positive_fraction, confidence_level
    - 3 runner-up columns: runner_up_label, runner_up_score, gap

    Parameters
    ----------
    cluster_ann : pd.DataFrame
        Cluster annotations in 24-column format (from cluster_annotations_subcluster)
    scores_df : pd.DataFrame, optional
        Marker scores for runner-up calculation
    stage_h_dir : str, optional
        Path to Stage H directory for lvl0 data
    iteration : int
        Current iteration number (default: 1)
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    pd.DataFrame
        44-column enhanced annotations DataFrame
    """
    logger = logger or logging.getLogger(__name__)

    # Hierarchical columns that exist per level
    HIER_COLS = [
        'cell_type', 'score', 'assigned_path', 'assigned_level',
        'root_label', 'confidence', 'min_margin_along_path', 'margin_is_infinite',
        'stop_reason', 'stopped_before_leaf', 'composition', 'decision_trace',
        'coverage', 'resolved_markers', 'is_ambiguous_root', 'ambiguous_root_candidates',
    ]

    # Load Stage H annotations for lvl0 data if available
    stage_h_lookup = {}
    if stage_h_dir:
        stage_h_path = Path(stage_h_dir) / 'cluster_annotations.csv'
        if stage_h_path.exists():
            try:
                stage_h_df = pd.read_csv(stage_h_path)
                for _, row in stage_h_df.iterrows():
                    cid = str(row.get('cluster_id', ''))
                    stage_h_lookup[cid] = row.to_dict()
                logger.info(f'  Loaded {len(stage_h_lookup)} Stage H annotations for lvl0')
            except Exception as e:
                logger.warning(f'  Failed to load Stage H annotations: {e}')

    records = []
    for _, row in cluster_ann.iterrows():
        cluster_id = str(row.get('cluster_id', ''))
        origin_cluster = str(row.get(
            'origin_cluster',
            cluster_id.split(':')[0] if ':' in cluster_id else cluster_id
        ))
        iteration_created = int(row.get('iteration_created', cluster_id.count(':')))

        # Build base record
        record = {
            'cluster_id': cluster_id,
            'origin_cluster': origin_cluster,
            'iteration_created': iteration_created,
            'n_cells': row.get('n_cells', 0),
            'proportion': row.get('proportion', 0),
            'reason': _build_reason_string(row),
        }

        # Add lvl0 columns (from Stage H or empty)
        h_data = stage_h_lookup.get(origin_cluster, {})
        for col in HIER_COLS:
            src_col = _map_col_name_44(col)
            val = h_data.get(src_col, '')
            # Handle special cases
            if col == 'cell_type':
                val = h_data.get('assigned_label', '')
            elif col == 'score':
                val = h_data.get('assigned_score', 0.0)
            record[f'{col}_lvl0'] = val if val != '' else (
                0.0 if 'score' in col or 'confidence' in col or 'coverage' in col else ''
            )

        # Add lvl1 columns (from current cluster_ann)
        for col in HIER_COLS:
            src_col = _map_col_name_44(col)
            val = row.get(src_col, '')
            # Handle special cases
            if col == 'cell_type':
                val = row.get('assigned_label', '')
            elif col == 'score':
                val = row.get('assigned_score', 0.0)
            record[f'{col}_lvl{iteration}'] = val if val != '' else (
                0.0 if 'score' in col or 'confidence' in col or 'coverage' in col else ''
            )

        # Add quality metrics
        record['mean_enrichment'] = round(
            float(row.get('mean_enrichment', 0) or 0), 3
        )
        record['mean_positive_fraction'] = round(
            float(row.get('mean_positive_fraction', 0) or 0), 3
        )
        record['confidence_level'] = row.get('confidence_level', '')

        # Add runner-up columns
        runner_up_label, runner_up_score, gap = _get_runner_up_44(
            cluster_id=cluster_id,
            assigned_label=row.get('assigned_label', ''),
            scores_df=scores_df,
        )
        record['runner_up_label'] = runner_up_label
        record['runner_up_score'] = round(runner_up_score, 3)
        record['gap'] = round(gap, 3)

        records.append(record)

    df = pd.DataFrame.from_records(records)

    # Ensure column order matches reference (44 columns)
    ref_col_order = [
        'cluster_id', 'origin_cluster', 'iteration_created', 'n_cells', 'proportion', 'reason',
    ]
    # Add lvl0 columns
    for col in HIER_COLS:
        ref_col_order.append(f'{col}_lvl0')
    # Add lvl1 columns
    for col in HIER_COLS:
        ref_col_order.append(f'{col}_lvl{iteration}')
    # Add quality and runner-up columns
    ref_col_order.extend([
        'mean_enrichment', 'mean_positive_fraction', 'confidence_level',
        'runner_up_label', 'runner_up_score', 'gap',
    ])

    # Reorder columns, only keeping those that exist
    final_cols = [c for c in ref_col_order if c in df.columns]
    df = df[final_cols]

    return df


# =============================================================================
# High-Level Review Exports Orchestrator
# =============================================================================


def run_review_exports(
    adata: "sc.AnnData",
    output_dir: Path,
    cell_type_col: Optional[str] = None,
    cluster_col: Optional[str] = None,
    stage_h_dir: Optional[str] = None,
    marker_map_path: Optional[str] = None,
    input_path: Optional[str] = None,
    logger: Optional[logging.Logger] = None,
) -> Dict[str, Path]:
    """Run all review exports (composition, annotations, summary, workflow state).

    This is the high-level orchestrator that combines all review export functions.

    Parameters
    ----------
    adata : sc.AnnData
        AnnData with annotations
    output_dir : Path
        Output directory
    cell_type_col : str, optional
        Cell type column name (auto-detected if not specified)
    cluster_col : str, optional
        Cluster column name (auto-detected if not specified)
    stage_h_dir : str, optional
        Path to Stage H directory for workflow state
    marker_map_path : str, optional
        Path to marker map for workflow state
    input_path : str, optional
        Path to input h5ad for workflow state
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    Dict[str, Path]
        Mapping from export type to output path
    """
    logger = logger or logging.getLogger(__name__)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    output_paths = {}

    # Detect columns if not specified
    if cell_type_col is None:
        cell_type_col = detect_cell_type_column(adata)
    if cluster_col is None:
        cluster_col = detect_cluster_column(adata)

    logger.info(f'Using cell type column: {cell_type_col}')
    logger.info(f'Using cluster column: {cluster_col}')

    # Export composition statistics
    comp_paths = export_composition_stats(adata, output_dir, cell_type_col, logger)
    output_paths.update(comp_paths)

    # Export marker scores
    scores_df = export_marker_scores_from_adata(adata, output_dir, logger)

    # Export cluster annotations from adata.uns if available (fast path)
    cluster_ann_key = 'cluster_annotations_subcluster'
    if cluster_ann_key in adata.uns:
        logger.info('Generating cluster annotations from stored data...')
        cluster_ann = adata.uns[cluster_ann_key].copy()
        total_cells = len(adata)

        # Add missing columns to match reference 24-column format
        if 'origin_cluster' not in cluster_ann.columns:
            cluster_ann['origin_cluster'] = cluster_ann['cluster_id'].astype(str).apply(
                lambda x: x.split(':')[0] if ':' in x else x
            )
        if 'iteration_created' not in cluster_ann.columns:
            cluster_ann['iteration_created'] = cluster_ann['cluster_id'].astype(str).apply(
                lambda x: x.count(':')
            )
        if 'proportion' not in cluster_ann.columns:
            cluster_ann['proportion'] = round(
                cluster_ann['n_cells'] / total_cells * 100, 2
            )

        # Add mean_enrichment and mean_positive_fraction from marker_scores
        # IMPORTANT: Use the ASSIGNED LABEL's metrics, not the highest-scoring label
        # (Fixed 2025-12-30: was incorrectly using idxmax which gave highest score's metrics)
        if scores_df is not None and not scores_df.empty:
            if 'cluster_id' in scores_df.columns and 'score' in scores_df.columns:
                # Build lookup: (cluster_id, label) -> {mean_enrichment, mean_positive_fraction}
                scores_df_copy = scores_df.copy()
                scores_df_copy['cluster_id'] = scores_df_copy['cluster_id'].astype(str)
                scores_df_copy['label'] = scores_df_copy['label'].astype(str)

                # Create lookup dictionary keyed by (cluster_id, label)
                metrics_lookup = {}
                for _, row in scores_df_copy.iterrows():
                    key = (row['cluster_id'], row['label'])
                    metrics_lookup[key] = {
                        'mean_enrichment': row.get('mean_enrichment', 0),
                        'mean_positive_fraction': row.get('mean_positive_fraction', 0),
                    }

                # Also build fallback: cluster_id -> highest score's metrics
                best_idx = scores_df_copy.groupby('cluster_id')['score'].idxmax()
                fallback_lookup = {}
                for idx in best_idx:
                    row = scores_df_copy.loc[idx]
                    fallback_lookup[row['cluster_id']] = {
                        'mean_enrichment': row.get('mean_enrichment', 0),
                        'mean_positive_fraction': row.get('mean_positive_fraction', 0),
                    }

                # Apply: use assigned_label's metrics, fall back to highest score
                if 'mean_enrichment' not in cluster_ann.columns or 'mean_positive_fraction' not in cluster_ann.columns:
                    me_values = []
                    mpf_values = []
                    for _, row in cluster_ann.iterrows():
                        cid = str(row['cluster_id'])
                        label = str(row.get('assigned_label', ''))
                        key = (cid, label)

                        if key in metrics_lookup:
                            # Use assigned label's metrics (correct behavior)
                            me_values.append(metrics_lookup[key]['mean_enrichment'])
                            mpf_values.append(metrics_lookup[key]['mean_positive_fraction'])
                        elif cid in fallback_lookup:
                            # Fallback to highest score (for Unassigned or missing labels)
                            me_values.append(fallback_lookup[cid]['mean_enrichment'])
                            mpf_values.append(fallback_lookup[cid]['mean_positive_fraction'])
                        else:
                            me_values.append(0)
                            mpf_values.append(0)

                    if 'mean_enrichment' not in cluster_ann.columns:
                        cluster_ann['mean_enrichment'] = [round(v, 3) if pd.notna(v) else 0 for v in me_values]
                    if 'mean_positive_fraction' not in cluster_ann.columns:
                        cluster_ann['mean_positive_fraction'] = [round(v, 3) if pd.notna(v) else 0 for v in mpf_values]

        # Add confidence_level if not present
        if 'confidence_level' not in cluster_ann.columns:
            cluster_ann['confidence_level'] = cluster_ann['assigned_score'].apply(
                _classify_confidence_band
            ).str.lower()

        # Reorder columns to match reference format
        ref_cols = [
            'cluster_id', 'origin_cluster', 'iteration_created', 'proportion',
            'mean_enrichment', 'mean_positive_fraction', 'n_cells', 'assigned_label',
            'assigned_path', 'assigned_level', 'assigned_score', 'root_label',
            'confidence', 'min_margin_along_path', 'margin_is_infinite', 'stop_reason',
            'stopped_before_leaf', 'composition', 'decision_trace', 'coverage',
            'resolved_markers', 'is_ambiguous_root', 'ambiguous_root_candidates',
            'confidence_level',
        ]
        final_cols = [c for c in ref_cols if c in cluster_ann.columns]
        cluster_ann = cluster_ann[final_cols]

        # Save cluster_annotations.csv (24 columns)
        stage_h_path = output_dir / 'cluster_annotations.csv'
        cluster_ann.to_csv(stage_h_path, index=False)
        output_paths['cluster_annotations'] = stage_h_path
        logger.info(
            f'  Wrote cluster_annotations.csv ({len(cluster_ann)} clusters, '
            f'{len(final_cols)} columns)'
        )

        # Generate enhanced version with 44-column format
        enhanced_df = generate_enhanced_annotations_44col(
            cluster_ann=cluster_ann,
            scores_df=scores_df,
            stage_h_dir=stage_h_dir,
            iteration=1,
            logger=logger,
        )
        enhanced_path = output_dir / 'cluster_annotations_enhanced.csv'
        enhanced_df.to_csv(enhanced_path, index=False)
        output_paths['cluster_annotations_enhanced'] = enhanced_path
        logger.info(
            f'  Wrote cluster_annotations_enhanced.csv ({len(enhanced_df.columns)} columns)'
        )
    else:
        # Fallback: use simple export
        logger.info('Generating cluster annotations (simple format)...')
        export_cluster_annotations_simple(
            adata, output_dir, cell_type_col, cluster_col, scores_df, logger
        )
        output_paths['cluster_annotations'] = output_dir / 'cluster_annotations.csv'
        output_paths['cluster_annotations_enhanced'] = (
            output_dir / 'cluster_annotations_enhanced.csv'
        )

    # Export review summary
    comp_global = adata.obs[cell_type_col].value_counts()
    summary_path = export_review_summary(
        adata, output_dir, cell_type_col, cluster_col, comp_global, logger
    )
    output_paths['review_summary'] = summary_path

    # Export workflow state
    state_path = export_workflow_state(
        output_dir,
        marker_map_path=marker_map_path,
        input_path=input_path,
        stage_h_dir=stage_h_dir,
        logger=logger,
    )
    output_paths['workflow_state'] = state_path

    logger.info('')
    logger.info('Review export complete!')
    logger.info(f'  Total cells: {adata.n_obs:,}')
    logger.info(f'  Cell types: {len(comp_global)}')
    logger.info(f'  Output files: {len(output_paths)}')

    return output_paths
