"""Refinement module CLI runner.

Enables running Stage I refinement as:
    python -m celltype_refinery.core.refinement --input <h5ad> --output <dir>

This module provides focus-driven cell-type annotation refinement:
1. Load Stage H output (coarse_clusters.h5ad)
2. Generate diagnostic report with refinement recommendations
3. Optionally execute recommendations with --execute flag

Usage Examples:
    # Diagnostic mode (default): show recommendations without modifying data
    python -m celltype_refinery.core.refinement \
        --input output/fallopian_tube/stage_h/coarse_clusters.h5ad \
        --marker-map data/fallopian_tube/FT_cell_type_markers_v9.json \
        --output output/fallopian_tube/stage_i

    # Execution mode: apply recommendations
    python -m celltype_refinery.core.refinement \
        --input output/fallopian_tube/stage_h/coarse_clusters.h5ad \
        --marker-map data/fallopian_tube/FT_cell_type_markers_v9.json \
        --output output/fallopian_tube/stage_i \
        --execute

    # Focus on specific cell types
    python -m celltype_refinery.core.refinement \
        --input output/fallopian_tube/stage_h/coarse_clusters.h5ad \
        --marker-map data/fallopian_tube/FT_cell_type_markers_v9.json \
        --focus-labels "Immune Cells,Epithelium" \
        --output output/fallopian_tube/stage_i \
        --execute

    # Focus on specific clusters
    python -m celltype_refinery.core.refinement \
        --input output/fallopian_tube/stage_h/coarse_clusters.h5ad \
        --marker-map data/fallopian_tube/FT_cell_type_markers_v9.json \
        --focus-clusters "5,12,21" \
        --output output/fallopian_tube/stage_i \
        --execute
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set

import numpy as np
import pandas as pd

try:
    import scanpy as sc
except ImportError:
    sc = None

# Internal imports
from .plan import RefinePlan, SubclusterOp, RescoreOp
from .engine import RefinementEngine, EngineResult
from .policies import AutoPolicy, AutoPolicyConfig
from .schema import InputAdapter, CanonicalSchema
from .provenance import (
    RefinementProvenance,
    create_provenance,
    store_provenance,
    detect_parent_stage,
)
from .validation import ValidationResult
from .cluster_metrics import compute_cluster_marker_heterogeneity


# =============================================================================
# Configuration
# =============================================================================

DEFAULT_LOG_FILENAME = "stage_i_refine.log"


@dataclass
class StageIRunConfig:
    """Configuration for Stage I refinement run.

    Attributes
    ----------
    min_cells : int
        Minimum cluster size for subclustering
    score_threshold : float
        Maximum score for low-confidence criterion
    subtype_signal_threshold : float
        Minimum child score to trigger parent subclustering
    heterogeneity_gap : float
        Score gap threshold for heterogeneity detection
    subcluster_resolution : float
        Leiden resolution for subclustering
    n_pcs : int
        Number of PCs for subclustering
    neighbors_k : int
        Number of neighbors for subclustering
    n_workers : int
        Parallel workers for subclustering
    label_key_out : str
        Output column name for cell type labels
    de_workers : int
        Parallel workers for DE computation
    de_layer : str
        Data layer for DE analysis
    de_within_parent : bool
        Use within-parent DE strategy
    de_tie_correct : bool
        Apply tie correction for Wilcoxon test
    scoring_workers : int
        Parallel workers for marker scoring
    scoring_batch_size : int
        Batch size for parallel scoring
    """

    min_cells: int = 500
    score_threshold: float = 1.0
    subtype_signal_threshold: float = 1.0
    heterogeneity_gap: float = 0.5
    subcluster_resolution: float = 0.4
    n_pcs: int = 30
    neighbors_k: int = 15
    n_workers: int = 1
    label_key_out: str = "cell_type_lvl1"
    de_workers: int = 8
    de_layer: str = "batchcorr"
    de_within_parent: bool = True
    de_tie_correct: bool = True
    scoring_workers: int = 1
    scoring_batch_size: int = 50
    gating_workers: int = 1


# =============================================================================
# Utility Functions
# =============================================================================


def _setup_logging(
    log_dir: Optional[Path],
    verbose: bool = False,
) -> logging.Logger:
    """Setup logging for Stage I refinement.

    Parameters
    ----------
    log_dir : Path, optional
        Directory for log files
    verbose : bool
        Enable debug logging

    Returns
    -------
    logging.Logger
        Configured logger
    """
    logger = logging.getLogger("stage_i_refine")
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)

    # Clear existing handlers
    logger.handlers.clear()

    # Console handler
    console = logging.StreamHandler(sys.stdout)
    console.setLevel(logging.DEBUG if verbose else logging.INFO)
    console_fmt = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    console.setFormatter(console_fmt)
    logger.addHandler(console)

    # File handler
    if log_dir:
        log_dir.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_dir / DEFAULT_LOG_FILENAME)
        file_handler.setLevel(logging.DEBUG)
        file_fmt = logging.Formatter(
            "%(asctime)s [%(levelname)s] %(name)s: %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
        file_handler.setFormatter(file_fmt)
        logger.addHandler(file_handler)

    return logger


def _parse_comma_list(value: Optional[str]) -> Optional[List[str]]:
    """Parse comma-separated string to list."""
    if not value:
        return None
    return [v.strip() for v in value.split(",") if v.strip()]


def _parse_comma_set(value: Optional[str]) -> Optional[Set[str]]:
    """Parse comma-separated string to set."""
    parsed = _parse_comma_list(value)
    return set(parsed) if parsed else None


def _resolve_focus_clusters(
    cluster_annotations: pd.DataFrame,
    focus_labels: Optional[List[str]],
    focus_clusters: Optional[Set[str]],
    logger: logging.Logger,
) -> Optional[Set[str]]:
    """Resolve focus controls to eligible cluster IDs.

    Parameters
    ----------
    cluster_annotations : pd.DataFrame
        Cluster annotations with cluster_id and assigned_label columns
    focus_labels : List[str], optional
        Labels to focus on (matches exact or prefix)
    focus_clusters : Set[str], optional
        Specific cluster IDs to focus on
    logger : logging.Logger
        Logger for messages

    Returns
    -------
    Set[str] or None
        Set of eligible cluster IDs, or None if no focus controls
    """
    if focus_labels is None and focus_clusters is None:
        return None

    eligible = set()

    # Resolve focus_labels to cluster IDs
    if focus_labels:
        for label in focus_labels:
            # Match exact label or path prefix
            mask = (
                (cluster_annotations["assigned_label"] == label)
                | cluster_annotations["assigned_label"].str.startswith(f"{label} / ")
                | cluster_annotations["assigned_label"].str.startswith(f"{label}/")
            )
            matching = cluster_annotations.loc[mask, "cluster_id"].astype(str).unique()
            eligible.update(matching)
            logger.debug("Focus label '%s' matched %d clusters", label, len(matching))

    # Add focus_clusters directly
    if focus_clusters:
        eligible.update(focus_clusters)
        logger.debug("Focus clusters: %s", focus_clusters)

    # If both specified, take intersection
    if focus_labels and focus_clusters:
        label_eligible = eligible.copy()
        eligible = label_eligible.intersection(focus_clusters)
        logger.info(
            "Focus intersection: %d clusters (labels=%d, clusters=%d)",
            len(eligible),
            len(label_eligible),
            len(focus_clusters),
        )

    return eligible if eligible else None


def _print_diagnostic_table(df: pd.DataFrame, logger: logging.Logger) -> None:
    """Print formatted diagnostic table to logger.

    Parameters
    ----------
    df : pd.DataFrame
        Diagnostic report from AutoPolicy.generate_diagnostic_report()
    logger : logging.Logger
        Logger to print to
    """
    if df.empty:
        logger.info("No clusters to diagnose.")
        return

    logger.info("=" * 130)
    logger.info("DIAGNOSTIC REPORT: Cluster Refinement Recommendations")
    logger.info("=" * 130)

    # Header
    logger.info(
        "| %-8s | %-32s | %8s | %8s | %-12s | %-55s |",
        "Cluster",
        "Label",
        "Score",
        "Cells",
        "Recommend",
        "Reason",
    )
    logger.info("-" * 130)

    # Rows
    for _, row in df.iterrows():
        recommendation = row["recommendation"]
        # Highlight SUBCLUSTER and RELABEL
        if recommendation == "SUBCLUSTER":
            rec_str = f"[{recommendation}]"
        elif recommendation == "RELABEL":
            rec_str = f"({recommendation})"
        else:
            rec_str = recommendation

        logger.info(
            "| %-8s | %-32s | %8.3f | %8d | %-12s | %-55s |",
            str(row["cluster_id"])[:8],
            str(row["assigned_label"])[:32],
            row["assigned_score"],
            row["n_cells"],
            rec_str,
            str(row["recommendation_reason"])[:55],
        )

    # Summary
    n_subcluster = (df["recommendation"] == "SUBCLUSTER").sum()
    n_relabel = (df["recommendation"] == "RELABEL").sum()
    n_skip = (df["recommendation"] == "SKIP").sum()

    logger.info("-" * 130)
    logger.info(
        "Summary: %d SUBCLUSTER, %d RELABEL, %d SKIP", n_subcluster, n_relabel, n_skip
    )
    logger.info("=" * 130)


def _print_review_summary(
    df: pd.DataFrame, total_clusters: int, logger: logging.Logger
) -> None:
    """Print review summary block.

    Parameters
    ----------
    df : pd.DataFrame
        Diagnostic report
    total_clusters : int
        Total number of clusters
    logger : logging.Logger
        Logger to print to
    """
    if df.empty:
        return

    n_eligible = len(df)
    n_subcluster = (df["recommendation"] == "SUBCLUSTER").sum()
    n_relabel = (df["recommendation"] == "RELABEL").sum()
    n_skip = (df["recommendation"] == "SKIP").sum()
    n_action = n_subcluster + n_relabel

    logger.info("")
    logger.info("=" * 72)
    logger.info(
        "REVIEW SUMMARY: %d of %d eligible clusters need refinement",
        n_action,
        n_eligible,
    )
    logger.info("=" * 72)
    logger.info(
        "  Eligible for refinement: %d / %d total clusters", n_eligible, total_clusters
    )
    logger.info("  Recommendations:")
    if n_subcluster > 0:
        logger.info(
            "    - SUBCLUSTER: %d clusters (heterogeneous/ambiguous)", n_subcluster
        )
    if n_relabel > 0:
        logger.info("    - RELABEL: %d clusters (strong child signal)", n_relabel)
    if n_skip > 0:
        logger.info("    - SKIP: %d clusters (confident assignments)", n_skip)
    logger.info("")


def _build_plan_from_diagnostic(
    diagnostic_report: pd.DataFrame,
    force_subcluster_ids: Optional[Set[str]],
    config: StageIRunConfig,
    logger: logging.Logger,
) -> RefinePlan:
    """Build execution plan from diagnostic report recommendations.

    Parameters
    ----------
    diagnostic_report : pd.DataFrame
        Diagnostic report from AutoPolicy.generate_diagnostic_report()
    force_subcluster_ids : Set[str], optional
        Cluster IDs to force subcluster even if recommendation is SKIP
    config : StageIRunConfig
        Configuration with subclustering parameters
    logger : logging.Logger
        Logger for messages

    Returns
    -------
    RefinePlan
        Execution plan with subcluster and relabel operations
    """
    plan = RefinePlan(
        metadata={
            "policy": "DiagnosticDriven",
            "created_at": datetime.now().isoformat(),
        }
    )

    force_subcluster_ids = force_subcluster_ids or set()

    for _, row in diagnostic_report.iterrows():
        cluster_id = str(row["cluster_id"])
        recommendation = row["recommendation"]

        # Check if forced
        is_forced = cluster_id in force_subcluster_ids

        if recommendation == "SUBCLUSTER" or is_forced:
            reason = (
                row["recommendation_reason"]
                if not is_forced
                else f"Forced: {row['recommendation_reason']}"
            )
            plan.add_subcluster(
                cluster_id=cluster_id,
                resolution=config.subcluster_resolution,
                n_pcs=config.n_pcs,
                neighbors_k=config.neighbors_k,
                min_cells=config.min_cells,
                reason=reason,
                source="diagnostic",
            )
            if is_forced:
                logger.info("Force subclustering cluster %s", cluster_id)

        elif recommendation == "RELABEL":
            plan.add_relabel(
                cluster_id=cluster_id,
                new_label=row.get("best_child", "Unknown"),
                old_label=row["assigned_label"],
                confidence_score=row.get("best_child_score", 0.0),
                reason=row["recommendation_reason"],
                source="diagnostic",
            )

    # Add rescore if there are modifications
    if plan.operations:
        subclustered = [
            op.cluster_id for op in plan.get_operations_by_type("subcluster")
        ]
        if subclustered or plan.get_operations_by_type("relabel"):
            plan.add_rescore(
                mode="smart",
                target_clusters=subclustered if subclustered else None,
                reason="Rescore after refinement",
            )

    plan.sort_operations()
    return plan


def _load_stage_h_artifacts(
    input_path: Path,
    cluster_annotations_path: Optional[Path],
    marker_scores_path: Optional[Path],
    adapter: InputAdapter,
    logger: logging.Logger,
) -> tuple[Optional[pd.DataFrame], Optional[pd.DataFrame], Optional[Path], Optional[Path]]:
    """Load cluster annotations and marker scores from Stage H output.

    Parameters
    ----------
    input_path : Path
        Path to input AnnData file
    cluster_annotations_path : Path, optional
        Explicit path to cluster_annotations.csv
    marker_scores_path : Path, optional
        Explicit path to marker_scores.csv
    adapter : InputAdapter
        Adapter for input validation and normalization
    logger : logging.Logger
        Logger for messages

    Returns
    -------
    tuple
        (cluster_annotations_df, marker_scores_df, annotations_path, scores_path)
    """
    input_dir = input_path.parent

    # Cluster annotations
    if cluster_annotations_path and cluster_annotations_path.exists():
        annotations_path = cluster_annotations_path
    else:
        candidates = [
            input_dir / "cluster_annotations.csv",
            input_dir / "subcluster_annotations.csv",
            input_dir.parent / "stage_h" / "cluster_annotations.csv",
        ]
        annotations_path = None
        for candidate in candidates:
            if candidate.exists():
                annotations_path = candidate
                break

    cluster_annotations = None
    if annotations_path and annotations_path.exists():
        logger.info("Loading cluster annotations: %s", annotations_path)
        raw_annotations = pd.read_csv(annotations_path)

        # Validate and adapt to canonical column names
        result = adapter.validate_cluster_annotations(raw_annotations)
        result.log_warnings(logger)

        if result.is_valid:
            cluster_annotations = adapter.adapt_cluster_annotations(
                raw_annotations, validate=False
            )
            logger.info("  Adapted columns: %s", list(result.adapted_columns.items()))
        else:
            logger.error("Cluster annotations validation failed:")
            for error in result.errors:
                logger.error("  %s", error)
    else:
        logger.warning("Could not find cluster annotations file")

    # Marker scores
    if marker_scores_path and marker_scores_path.exists():
        scores_path = marker_scores_path
    else:
        candidates = [
            input_dir / "marker_scores.csv",
            input_dir / "subcluster_scores.csv",
            input_dir.parent / "stage_h" / "marker_scores.csv",
        ]
        scores_path = None
        for candidate in candidates:
            if candidate.exists():
                scores_path = candidate
                break

    marker_scores = None
    if scores_path and scores_path.exists():
        logger.info("Loading marker scores: %s", scores_path)
        raw_scores = pd.read_csv(scores_path)

        result = adapter.validate_marker_scores(raw_scores)
        result.log_warnings(logger)

        if result.is_valid:
            marker_scores = adapter.adapt_marker_scores(raw_scores, validate=False)
            logger.info("  Adapted columns: %s", list(result.adapted_columns.items()))
        else:
            logger.error("Marker scores validation failed:")
            for error in result.errors:
                logger.error("  %s", error)
    else:
        logger.warning(
            "Could not find marker scores file (parent subclustering will be limited)"
        )

    return cluster_annotations, marker_scores, annotations_path, scores_path


def _convert_to_native(obj: Any) -> Any:
    """Convert NumPy types to native Python types for JSON serialization."""
    if isinstance(obj, dict):
        return {k: _convert_to_native(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [_convert_to_native(v) for v in obj]
    elif isinstance(obj, tuple):
        return tuple(_convert_to_native(v) for v in obj)
    elif isinstance(obj, (np.integer, np.int64, np.int32)):
        return int(obj)
    elif isinstance(obj, (np.floating, np.float64, np.float32)):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, np.bool_):
        return bool(obj)
    else:
        return obj


# =============================================================================
# Export Functions
# =============================================================================


def export_plan(
    plan: RefinePlan,
    output_path: Path,
    format: str = "yaml",
) -> None:
    """Export RefinePlan to file.

    Parameters
    ----------
    plan : RefinePlan
        Plan to export
    output_path : Path
        Output file path
    format : str
        Format: "yaml" or "json"
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if format == "yaml":
        plan.to_yaml(output_path)
    elif format == "json":
        plan.to_json(output_path)
    else:
        raise ValueError(f"Unknown format: {format}. Use 'yaml' or 'json'.")


def export_curation_log(
    plan: RefinePlan,
    result: EngineResult,
    provenance: RefinementProvenance,
    output_path: Path,
) -> None:
    """Export curation log as JSON.

    Parameters
    ----------
    plan : RefinePlan
        Executed plan
    result : EngineResult
        Execution result
    provenance : RefinementProvenance
        Provenance information
    output_path : Path
        Output file path
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    log = {
        "timestamp": datetime.now().isoformat(),
        "provenance": provenance.to_dict(),
        "plan": {
            "content_checksum": plan.content_checksum(),
            "summary": plan.summary(),
            "operations": [op.to_dict() for op in plan.operations],
        },
        "result": {
            "success": result.success,
            "n_cells_modified": result.n_cells_modified,
            "operations_executed": result.operations_executed,
            "subclusters_created": result.subclusters_created,
            "execution_time_seconds": result.execution_time_seconds,
            "errors": result.errors,
            "warnings": result.warnings,
        },
    }

    log = _convert_to_native(log)

    with open(output_path, "w") as f:
        json.dump(log, f, indent=2)


def export_mapping_table(
    adata: "sc.AnnData",
    output_path: Path,
    label_key: str = "cell_type_curated",
    reason_key: str = "curation_reason",
) -> pd.DataFrame:
    """Export cluster to cell type mapping table.

    Parameters
    ----------
    adata : sc.AnnData
        AnnData object with annotations
    output_path : Path
        Output file path
    label_key : str
        Column with cell type labels
    reason_key : str
        Column with assignment reasons

    Returns
    -------
    pd.DataFrame
        Mapping table
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    cluster_col = "cluster_lvl1" if "cluster_lvl1" in adata.obs else "cluster_lvl0"

    cols = [cluster_col]
    has_label = label_key in adata.obs
    has_reason = reason_key in adata.obs
    if has_label:
        cols.append(label_key)
    if has_reason:
        cols.append(reason_key)

    obs_subset = adata.obs[cols].copy()
    obs_subset[cluster_col] = obs_subset[cluster_col].astype(str)

    def most_common(s):
        if len(s) == 0:
            return ""
        counts = s.value_counts()
        return counts.index[0] if len(counts) > 0 else ""

    grouped = obs_subset.groupby(cluster_col, sort=True)

    data = []
    for cluster_id, group_df in grouped:
        row = {
            "cluster_id": cluster_id,
            "cell_type": most_common(group_df[label_key]) if has_label else "Unknown",
            "n_cells": len(group_df),
            "reason": most_common(group_df[reason_key]) if has_reason else "",
        }
        data.append(row)

    result = pd.DataFrame(data)
    col_order = ["cluster_id", "cell_type", "n_cells", "reason"]
    result = result[[c for c in col_order if c in result.columns]]
    result.to_csv(output_path, index=False)
    return result


def export_annotations_csv(
    adata: "sc.AnnData",
    output_path: Path,
    columns: Optional[List[str]] = None,
) -> pd.DataFrame:
    """Export cell annotations as CSV.

    Parameters
    ----------
    adata : sc.AnnData
        AnnData object
    output_path : Path
        Output file path
    columns : List[str], optional
        Columns to export

    Returns
    -------
    pd.DataFrame
        Exported annotations
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if columns is None:
        columns = []
        for col in [
            "cluster_lvl0",
            "cluster_lvl1",
            "cell_type_curated",
            "curation_reason",
        ]:
            if col in adata.obs:
                columns.append(col)

    df = adata.obs[columns].copy()
    df.to_csv(output_path)
    return df


def generate_next_iteration_template(
    adata: "sc.AnnData",
    output_path: Path,
    label_key: str = "cell_type_curated",
    current_iteration: int = 1,
) -> str:
    """Generate a template config for the next refinement iteration.

    Optimized implementation using groupby for O(n) instead of O(n × k).

    Parameters
    ----------
    adata : sc.AnnData
        AnnData object with current annotations
    output_path : Path
        Output file path for template
    label_key : str
        Column with current cell type labels
    current_iteration : int
        Current iteration number

    Returns
    -------
    str
        Template content
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    lines = [
        "# Next iteration refinement configuration",
        f"# Generated from iteration {current_iteration}",
        f"# Timestamp: {datetime.now().isoformat()}",
        "",
        "version: 1.0",
        f"iteration: {current_iteration + 1}",
        "",
        "# Current cluster assignments (modify as needed):",
        "overrides:",
    ]

    cluster_col = "cluster_lvl1" if "cluster_lvl1" in adata.obs else "cluster_lvl0"
    has_label = label_key in adata.obs

    # Extract only needed columns for efficiency
    cols = [cluster_col]
    if has_label:
        cols.append(label_key)
    obs_subset = adata.obs[cols].copy()
    obs_subset[cluster_col] = obs_subset[cluster_col].astype(str)

    # Helper to get most common value
    def most_common(s):
        if len(s) == 0:
            return "Unknown"
        counts = s.value_counts()
        return counts.index[0] if len(counts) > 0 else "Unknown"

    # Use groupby for O(n) instead of O(n × k) per-cluster masks
    grouped = obs_subset.groupby(cluster_col, sort=True)

    for cluster_id, group_df in grouped:
        n_cells = len(group_df)
        label = most_common(group_df[label_key]) if has_label else "Unknown"

        lines.append(f"  # Cluster {cluster_id}: {n_cells:,} cells")
        lines.append(f"  - cluster_id: \"{cluster_id}\"")
        lines.append(f"    cell_type: \"{label}\"")
        lines.append(f"    reason: \"\"  # Add your reason here")
        lines.append("")

    content = "\n".join(lines)
    with open(output_path, "w") as f:
        f.write(content)

    return content


# =============================================================================
# Main Runner
# =============================================================================


def run_stage_i(
    input_path: Path,
    marker_map_path: Path,
    output_dir: Path,
    config: Optional[StageIRunConfig] = None,
    focus_labels: Optional[List[str]] = None,
    focus_clusters: Optional[Set[str]] = None,
    force_subcluster: Optional[Set[str]] = None,
    execute: bool = False,
    cluster_annotations_path: Optional[Path] = None,
    marker_scores_path: Optional[Path] = None,
    verbose: bool = False,
    log_dir: Optional[Path] = None,
) -> int:
    """Run Stage I refinement.

    Parameters
    ----------
    input_path : Path
        Path to input h5ad file (Stage H output)
    marker_map_path : Path
        Path to marker map JSON file
    output_dir : Path
        Output directory for refined data
    config : StageIRunConfig, optional
        Configuration for refinement
    focus_labels : List[str], optional
        Labels to focus on
    focus_clusters : Set[str], optional
        Cluster IDs to focus on
    force_subcluster : Set[str], optional
        Cluster IDs to force subcluster
    execute : bool
        Whether to execute recommendations (default: diagnostic only)
    cluster_annotations_path : Path, optional
        Explicit path to cluster_annotations.csv
    marker_scores_path : Path, optional
        Explicit path to marker_scores.csv
    verbose : bool
        Enable verbose logging
    log_dir : Path, optional
        Directory for log files

    Returns
    -------
    int
        Exit code (0 = success, non-zero = error)
    """
    if sc is None:
        print("ERROR: scanpy is required. Install with: pip install scanpy")
        return 1

    config = config or StageIRunConfig()
    log_dir = log_dir or output_dir
    logger = _setup_logging(log_dir, verbose)

    # -------------------------------------------------------------------------
    # Phase 1: Initialization
    # -------------------------------------------------------------------------
    logger.info("")
    logger.info("=" * 75)
    logger.info("PHASE 1/7: INITIALIZATION")
    logger.info("=" * 75)
    logger.info("Stage I Refinement: Focus-Driven Cell-Type Refinement")
    logger.info("Timestamp: %s", datetime.now().isoformat())
    logger.info("Input: %s", input_path)
    logger.info("Output: %s", output_dir)
    logger.info("Marker map: %s", marker_map_path)
    logger.info(
        "Focus labels: %s", ", ".join(focus_labels) if focus_labels else "None"
    )
    logger.info(
        "Focus clusters: %s", ", ".join(focus_clusters) if focus_clusters else "None"
    )
    logger.info("Execute: %s", "yes" if execute else "no (diagnostic only)")

    # -------------------------------------------------------------------------
    # Phase 2: Data Loading
    # -------------------------------------------------------------------------
    logger.info("")
    logger.info("=" * 75)
    logger.info("PHASE 2/7: DATA LOADING & VALIDATION")
    logger.info("=" * 75)

    if not input_path.exists():
        logger.error("Input file not found: %s", input_path)
        return 1

    logger.info("Loading input data: %s", input_path)
    adata = sc.read_h5ad(input_path)
    logger.info("Loaded %d cells, %d features", adata.n_obs, adata.n_vars)

    # Create input adapter
    adapter = InputAdapter(logger=logger)
    adata_result = adapter.validate_adata(adata)
    adata_result.log_warnings(logger)

    if not adata_result.is_valid:
        logger.error("Input validation failed:")
        for error in adata_result.errors:
            logger.error("  %s", error)
        return 1

    logger.info("Input validation passed")

    # Detect parent stage
    parent_stage, iteration = detect_parent_stage(adata)
    logger.info("Parent stage: %s (iteration %d)", parent_stage, iteration)

    # -------------------------------------------------------------------------
    # Phase 3: Load Stage H Artifacts
    # -------------------------------------------------------------------------
    logger.info("")
    logger.info("=" * 75)
    logger.info("PHASE 3/7: STAGE H ARTIFACTS & FOCUS")
    logger.info("=" * 75)

    cluster_annotations, marker_scores, annotations_path, scores_path = (
        _load_stage_h_artifacts(
            input_path=input_path,
            cluster_annotations_path=cluster_annotations_path,
            marker_scores_path=marker_scores_path,
            adapter=adapter,
            logger=logger,
        )
    )

    if cluster_annotations is None:
        logger.error("Could not load cluster annotations")
        return 1

    # Resolve focus controls
    eligible_clusters = _resolve_focus_clusters(
        cluster_annotations=cluster_annotations,
        focus_labels=focus_labels,
        focus_clusters=focus_clusters,
        logger=logger,
    )

    if eligible_clusters is not None:
        logger.info(
            "Focus controls active: %d / %d clusters eligible for refinement",
            len(eligible_clusters),
            cluster_annotations["cluster_id"].nunique(),
        )
        if len(eligible_clusters) == 0:
            logger.warning("No clusters match focus criteria. Nothing to refine.")
            return 0

    # -------------------------------------------------------------------------
    # Phase 4: Policy Configuration
    # -------------------------------------------------------------------------
    logger.info("")
    logger.info("=" * 75)
    logger.info("PHASE 4/7: POLICY CONFIGURATION")
    logger.info("=" * 75)

    auto_policy = AutoPolicy(
        min_cells=config.min_cells,
        score_threshold=config.score_threshold,
        subtype_signal_threshold=config.subtype_signal_threshold,
        heterogeneity_gap=config.heterogeneity_gap,
        subcluster_resolution=config.subcluster_resolution,
        n_pcs=config.n_pcs,
        neighbors_k=config.neighbors_k,
        logger=logger,
    )

    logger.info("AutoPolicy configuration:")
    logger.info("  min_cells: %d", config.min_cells)
    logger.info("  score_threshold: %.2f", config.score_threshold)
    logger.info("  subtype_signal_threshold: %.2f", config.subtype_signal_threshold)
    logger.info("  heterogeneity_gap: %.2f", config.heterogeneity_gap)
    logger.info("  subcluster_resolution: %.2f", config.subcluster_resolution)

    # Compute heterogeneity metrics for weak leaf detection
    cluster_col = "cluster_lvl1" if "cluster_lvl1" in adata.obs else "cluster_lvl0"
    logger.info("")
    logger.info("Computing cluster marker heterogeneity for weak leaf detection...")
    try:
        heterogeneity_metrics = compute_cluster_marker_heterogeneity(
            adata,
            cluster_col=cluster_col,
            layer="batchcorr" if "batchcorr" in adata.layers else None,
            logger=logger,
        )
        auto_policy.set_heterogeneity_metrics(heterogeneity_metrics)
        logger.info(
            "Heterogeneity metrics computed for %d clusters (mean=%.3f, max=%.3f)",
            len(heterogeneity_metrics),
            heterogeneity_metrics["marker_heterogeneity"].mean(),
            heterogeneity_metrics["marker_heterogeneity"].max(),
        )
    except Exception as e:
        logger.warning("Failed to compute heterogeneity metrics: %s", e)
        logger.warning("Weak leaf detection will be disabled.")

    # -------------------------------------------------------------------------
    # Phase 5: Diagnostic Report
    # -------------------------------------------------------------------------
    logger.info("")
    logger.info("=" * 75)
    logger.info("PHASE 5/7: DIAGNOSTIC REPORT")
    logger.info("=" * 75)

    diagnostic_report = auto_policy.generate_diagnostic_report(
        cluster_annotations,
        marker_scores if marker_scores is not None else pd.DataFrame(),
        eligible_clusters=eligible_clusters,
    )

    _print_diagnostic_table(diagnostic_report, logger)
    _print_review_summary(
        diagnostic_report, cluster_annotations["cluster_id"].nunique(), logger
    )

    # Save diagnostic report
    output_dir.mkdir(parents=True, exist_ok=True)
    report_path = output_dir / "diagnostic_report.csv"
    diagnostic_report.to_csv(report_path, index=False)
    logger.info("Diagnostic report saved to: %s", report_path)

    if not execute:
        logger.info("-" * 75)
        logger.info(
            "Diagnostic mode complete. Review the report and re-run with --execute to apply."
        )
        return 0

    # -------------------------------------------------------------------------
    # Phase 6: Execution
    # -------------------------------------------------------------------------
    logger.info("")
    logger.info("=" * 75)
    logger.info("PHASE 6/7: REFINEMENT EXECUTION")
    logger.info("=" * 75)

    # Build plan
    merged_plan = _build_plan_from_diagnostic(
        diagnostic_report,
        force_subcluster_ids=force_subcluster,
        config=config,
        logger=logger,
    )

    logger.info("Execution plan: %s", merged_plan.summary())
    logger.info("Affected clusters: %s", merged_plan.get_affected_clusters())

    # Export plan
    plan_path = output_dir / "execution_plan.yaml"
    export_plan(merged_plan, plan_path, format="yaml")
    logger.info("Plan saved to: %s", plan_path)

    # Import DE runner if using within-parent strategy
    from .engine import default_scanpy_de_within_parent_runner

    # Execute
    engine = RefinementEngine(
        marker_map=marker_map_path,
        label_key_out=config.label_key_out,
        adapter=adapter,
        stage_h_scores_path=scores_path,
        stage_h_annotations_path=annotations_path,
        logger=logger,
        output_dir=output_dir,
        n_workers=config.n_workers,
        # DE parameters
        de_workers=config.de_workers,
        de_layer=config.de_layer,
        de_within_parent=config.de_within_parent,
        de_tie_correct=config.de_tie_correct,
        de_within_parent_runner=default_scanpy_de_within_parent_runner if config.de_within_parent else None,
        # Scoring parameters
        scoring_workers=config.scoring_workers,
        scoring_batch_size=config.scoring_batch_size,
        # Gating parameters
        gating_workers=config.gating_workers,
    )

    adata = engine.execute(adata, merged_plan)
    result = engine.get_result()

    if result is None or not result.success:
        logger.error("Execution failed")
        if result:
            logger.error("Errors: %s", result.errors)
        return 1

    # -------------------------------------------------------------------------
    # Phase 7: Output Export
    # -------------------------------------------------------------------------
    logger.info("")
    logger.info("=" * 75)
    logger.info("PHASE 7/7: OUTPUT EXPORT")
    logger.info("=" * 75)

    # Store provenance
    provenance = create_provenance(
        adata=adata,
        plan=merged_plan,
        result=result,
        label_key_out=config.label_key_out,
        auto_policy_config=auto_policy.config.__dict__,
        manual_config_path=None,
        gpu_used=False,
        input_path=str(input_path),
    )
    store_provenance(adata, provenance)

    # Write outputs
    output_h5ad = output_dir / f"refined_v{provenance.iteration}.h5ad"
    logger.info("Writing refined AnnData: %s", output_h5ad)
    adata.write_h5ad(output_h5ad)

    # Also write as refined_final.h5ad for chaining
    final_path = output_dir / "refined_final.h5ad"
    adata.write_h5ad(final_path)
    logger.info("Written final output: %s", final_path)

    # Export curation log
    log_path = output_dir / "curation_log.json"
    export_curation_log(merged_plan, result, provenance, log_path)
    logger.info("Wrote curation log: %s", log_path)

    # Export mapping table
    mapping_path = output_dir / "cluster_label_mapping.csv"
    export_mapping_table(adata, mapping_path, label_key=config.label_key_out)
    logger.info("Wrote mapping table: %s", mapping_path)

    # Export annotations CSV
    annotations_export_path = output_dir / "final_annotations.csv"
    export_annotations_csv(adata, annotations_export_path)
    logger.info("Wrote annotations: %s", annotations_export_path)

    # Generate next iteration template
    template_path = output_dir / "next_iteration_template.yaml"
    generate_next_iteration_template(
        adata,
        template_path,
        label_key=config.label_key_out,
        current_iteration=provenance.iteration,
    )
    logger.info("Wrote next iteration template: %s", template_path)

    # Summary
    logger.info("")
    logger.info("=" * 75)
    logger.info("REFINEMENT COMPLETE")
    logger.info("=" * 75)
    logger.info("  Iteration: %d", provenance.iteration)
    logger.info("  Cells modified: %d", result.n_cells_modified)
    logger.info("  Operations: %s", result.operations_executed)
    logger.info("  Subclusters created: %s", result.subclusters_created)
    logger.info("  Execution time: %.2f seconds", result.execution_time_seconds)
    logger.info("  Output: %s", output_h5ad)

    return 0


# =============================================================================
# CLI
# =============================================================================


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Stage I Refinement: Focus-driven cell-type annotation refinement.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Diagnostic mode (default): show recommendations
  python -m celltype_refinery.core.refinement \\
      --input stage_h/coarse_clusters.h5ad \\
      --marker-map markers.json \\
      --output stage_i

  # Execute recommendations
  python -m celltype_refinery.core.refinement \\
      --input stage_h/coarse_clusters.h5ad \\
      --marker-map markers.json \\
      --output stage_i \\
      --execute

  # Focus on specific labels
  python -m celltype_refinery.core.refinement \\
      --input stage_h/coarse_clusters.h5ad \\
      --marker-map markers.json \\
      --focus-labels "Immune Cells,Epithelium" \\
      --output stage_i \\
      --execute
        """,
    )

    # Input/Output
    parser.add_argument(
        "--input",
        "-i",
        type=Path,
        required=True,
        help="Input AnnData file (Stage H output).",
    )
    parser.add_argument(
        "--marker-map",
        "-m",
        type=Path,
        required=True,
        help="Path to marker map JSON file.",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        required=True,
        help="Output directory for refined data.",
    )

    # Execution controls
    exec_group = parser.add_argument_group("Execution Controls")
    exec_group.add_argument(
        "--execute",
        action="store_true",
        help="Execute recommendations. Without this flag, only shows diagnostic report.",
    )
    exec_group.add_argument(
        "--force-subcluster",
        type=str,
        default=None,
        help="Force subclustering for these clusters (comma-separated IDs).",
    )

    # Focus controls
    focus_group = parser.add_argument_group("Focus Controls")
    focus_group.add_argument(
        "--focus-labels",
        type=str,
        default=None,
        help="Restrict refinement to clusters with these labels (comma-separated).",
    )
    focus_group.add_argument(
        "--focus-clusters",
        type=str,
        default=None,
        help="Restrict refinement to these cluster IDs (comma-separated).",
    )

    # Policy parameters
    policy_group = parser.add_argument_group("Policy Parameters")
    policy_group.add_argument(
        "--min-cells",
        type=int,
        default=500,
        help="Minimum cluster size for subclustering (default: 500).",
    )
    policy_group.add_argument(
        "--score-threshold",
        type=float,
        default=1.0,
        help="Maximum score for low-confidence criterion (default: 1.0).",
    )
    policy_group.add_argument(
        "--subtype-signal-threshold",
        type=float,
        default=1.0,
        help="Minimum child score to trigger parent subclustering (default: 1.0).",
    )
    policy_group.add_argument(
        "--heterogeneity-gap",
        type=float,
        default=0.5,
        help="Score gap threshold for heterogeneity detection (default: 0.5).",
    )
    policy_group.add_argument(
        "--subcluster-resolution",
        type=float,
        default=0.4,
        help="Leiden resolution for subclustering (default: 0.4).",
    )

    # Execution options
    exec_opts = parser.add_argument_group("Execution Options")
    exec_opts.add_argument(
        "--n-workers",
        type=int,
        default=1,
        help="Parallel workers for subclustering (default: 1).",
    )
    exec_opts.add_argument(
        "--label-key-out",
        type=str,
        default="cell_type_lvl1",
        help="Output column name for cell type labels (default: cell_type_lvl1).",
    )
    exec_opts.add_argument(
        "--scoring-workers",
        type=int,
        default=1,
        help="Parallel workers for marker scoring (default: 1). Set >1 for parallel scoring.",
    )
    exec_opts.add_argument(
        "--scoring-batch-size",
        type=int,
        default=50,
        help="Batch size for parallel scoring (default: 50).",
    )
    exec_opts.add_argument(
        "--gating-workers",
        type=int,
        default=1,
        help="Parallel workers for hierarchical gating (default: 1). Set >1 for parallel assignment.",
    )

    # DE options
    de_opts = parser.add_argument_group("Differential Expression Options")
    de_opts.add_argument(
        "--de-workers",
        type=int,
        default=8,
        help="Parallel workers for DE computation (default: 8).",
    )
    de_opts.add_argument(
        "--de-layer",
        type=str,
        default="batchcorr",
        help="Data layer for DE analysis (default: batchcorr).",
    )
    de_opts.add_argument(
        "--de-within-parent",
        action="store_true",
        default=True,
        help="Use within-parent DE strategy (default: True).",
    )
    de_opts.add_argument(
        "--de-global",
        action="store_true",
        default=False,
        help="Use global DE strategy (overrides --de-within-parent).",
    )
    de_opts.add_argument(
        "--de-tie-correct",
        action="store_true",
        default=True,
        help="Apply tie correction for Wilcoxon test (default: True).",
    )
    de_opts.add_argument(
        "--no-de-tie-correct",
        action="store_true",
        default=False,
        help="Disable tie correction for Wilcoxon test.",
    )

    # Stage H artifacts
    artifacts_group = parser.add_argument_group("Stage H Artifacts (Optional)")
    artifacts_group.add_argument(
        "--cluster-annotations",
        type=Path,
        default=None,
        help="Explicit path to cluster_annotations.csv (auto-detected if not provided).",
    )
    artifacts_group.add_argument(
        "--marker-scores",
        type=Path,
        default=None,
        help="Explicit path to marker_scores.csv (auto-detected if not provided).",
    )

    # Logging
    parser.add_argument(
        "--log-dir",
        type=Path,
        default=None,
        help="Directory for log files (default: output directory).",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable verbose logging.",
    )

    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Main entry point."""
    args = parse_args(argv)

    # Determine DE strategy (--de-global overrides --de-within-parent)
    de_within_parent = args.de_within_parent and not args.de_global
    de_tie_correct = args.de_tie_correct and not args.no_de_tie_correct

    config = StageIRunConfig(
        min_cells=args.min_cells,
        score_threshold=args.score_threshold,
        subtype_signal_threshold=args.subtype_signal_threshold,
        heterogeneity_gap=args.heterogeneity_gap,
        subcluster_resolution=args.subcluster_resolution,
        n_workers=args.n_workers,
        label_key_out=args.label_key_out,
        de_workers=args.de_workers,
        de_layer=args.de_layer,
        de_within_parent=de_within_parent,
        de_tie_correct=de_tie_correct,
        scoring_workers=args.scoring_workers,
        scoring_batch_size=args.scoring_batch_size,
        gating_workers=args.gating_workers,
    )

    return run_stage_i(
        input_path=args.input,
        marker_map_path=args.marker_map,
        output_dir=args.output,
        config=config,
        focus_labels=_parse_comma_list(args.focus_labels),
        focus_clusters=_parse_comma_set(args.focus_clusters),
        force_subcluster=_parse_comma_set(args.force_subcluster),
        execute=args.execute,
        cluster_annotations_path=args.cluster_annotations,
        marker_scores_path=args.marker_scores,
        verbose=args.verbose,
        log_dir=args.log_dir,
    )


if __name__ == "__main__":
    sys.exit(main())
