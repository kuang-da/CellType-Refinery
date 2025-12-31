"""Diagnosis module CLI runner.

Enables running cluster diagnosis as:
    python -m celltype_refinery.core.diagnosis --input <h5ad> --output <dir>

This module generates a diagnostic report identifying which clusters need:
- SUBCLUSTER: Re-cluster to resolve heterogeneity
- RELABEL: Simple relabeling without subclustering
- SKIP: No action needed

Usage Examples:
    # Basic usage - generate diagnostic report from refined h5ad
    python -m celltype_refinery.core.diagnosis \\
        --input output/stage_i/refined_final.h5ad \\
        --output output/stage_i \\
        --min-cells 500 \\
        --score-threshold 1.0

    # With explicit annotations and scores
    python -m celltype_refinery.core.diagnosis \\
        --input output/stage_i/refined_final.h5ad \\
        --cluster-annotations output/stage_i/cluster_annotations.csv \\
        --marker-scores output/stage_i/marker_scores.csv \\
        --output output/stage_i

    # Higher thresholds for more aggressive subclustering
    python -m celltype_refinery.core.diagnosis \\
        --input output/stage_i/refined_final.h5ad \\
        --output output/stage_i \\
        --min-cells 2000 \\
        --score-threshold 1.3

Output:
    Creates diagnostic_report.csv in the output directory with columns:
    - cluster_id, origin_cluster, iteration_created
    - assigned_label, assigned_score, n_cells
    - recommendation (SUBCLUSTER | RELABEL | SKIP)
    - recommendation_reason
    - ideal_recommendation, is_blocked_by_size
    - passes_* columns for each criterion
    - confidence_band, triggered_criteria
"""

from __future__ import annotations

import argparse
import logging
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Optional

import pandas as pd
import yaml

try:
    import scanpy as sc
except ImportError:
    sc = None

from .engine import DiagnosticEngine
from .criteria import CriteriaConfig
from ..refinement.cluster_metrics import compute_cluster_marker_heterogeneity


# =============================================================================
# Helper functions for groups export
# =============================================================================


def derive_root_label(assigned_label: str) -> str:
    """Extract root label from hierarchical assignment.

    Handles hierarchical labels like "Epithelium→Ciliated" by extracting
    the first (root) component. Hybrid labels with "~" are kept as-is.

    Parameters
    ----------
    assigned_label : str
        The assigned cell type label (may be hierarchical or hybrid)

    Returns
    -------
    str
        The root label (first component before "→", or original if no "→")

    Examples
    --------
    >>> derive_root_label("Epithelium→Ciliated")
    'Epithelium'
    >>> derive_root_label("Epithelium")
    'Epithelium'
    >>> derive_root_label("Endo~Mesen")
    'Endo~Mesen'
    >>> derive_root_label("Unassigned")
    'Unassigned'
    """
    if pd.isna(assigned_label) or assigned_label == "Unassigned":
        return "Unassigned"
    # Split on "→" and take first part
    parts = str(assigned_label).split("→")
    return parts[0]


def is_ambiguous_root(assigned_label: str) -> bool:
    """Check if label indicates hybrid/ambiguous root.

    Hybrid labels contain "~" to indicate multiple possible lineages
    (e.g., "Endo~Mesen" for Endothelium~Mesenchymal hybrid).

    Parameters
    ----------
    assigned_label : str
        The assigned cell type label

    Returns
    -------
    bool
        True if label contains "~" (hybrid marker)

    Examples
    --------
    >>> is_ambiguous_root("Endo~Mesen")
    True
    >>> is_ambiguous_root("Epithelium")
    False
    >>> is_ambiguous_root("Epithelium→Ciliated")
    False
    """
    return "~" in str(assigned_label)


def setup_logging(log_level: str = "INFO") -> logging.Logger:
    """Setup logging for CLI.

    Parameters
    ----------
    log_level : str
        Logging level (DEBUG, INFO, WARNING, ERROR)

    Returns
    -------
    logging.Logger
        Configured logger
    """
    logging.basicConfig(
        level=getattr(logging, log_level.upper()),
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    return logging.getLogger(__name__)


def update_workflow_state(
    output_dir: Path,
    diagnostic_summary: Dict[str, int],
    diagnostic_report_path: str,
    logger: Optional[logging.Logger] = None,
) -> None:
    """Update workflow_state.yaml with diagnostic summary.

    This function loads the existing workflow_state.yaml (if present),
    updates the diagnose-related fields, and saves the file. This ensures
    the workflow state has accurate diagnostic_summary counts after the
    diagnose step runs.

    Matches reference behavior in ft/src/workflow/state.py:mark_diagnose_complete()

    Parameters
    ----------
    output_dir : Path
        Output directory containing workflow_state.yaml
    diagnostic_summary : Dict[str, int]
        Summary counts {"SUBCLUSTER": n, "RELABEL": n, "SKIP": n}
    diagnostic_report_path : str
        Path to the diagnostic_report.csv file
    logger : logging.Logger, optional
        Logger instance for output
    """
    state_path = output_dir / "workflow_state.yaml"

    # Load existing state or create minimal structure
    if state_path.exists():
        try:
            with open(state_path) as f:
                state = yaml.safe_load(f) or {}
        except Exception as e:
            if logger:
                logger.warning("Failed to load workflow_state.yaml: %s", e)
            state = {}
    else:
        # Create minimal state structure (matches reference schema)
        state = {
            "version": "1.0",
            "out_dir": str(output_dir),
            "created_at": datetime.now().isoformat(),
        }

    # Update diagnose-related fields (matches reference mark_diagnose_complete)
    state["diagnose_complete"] = True
    state["diagnose_timestamp"] = datetime.now().isoformat()
    state["diagnostic_summary"] = diagnostic_summary
    state["diagnostic_report_path"] = diagnostic_report_path
    state["updated_at"] = datetime.now().isoformat()

    # Save updated state
    try:
        with open(state_path, "w") as f:
            yaml.safe_dump(state, f, default_flow_style=False, sort_keys=False)
        if logger:
            logger.info("Updated workflow_state.yaml with diagnostic summary")
    except Exception as e:
        if logger:
            logger.warning("Failed to update workflow_state.yaml: %s", e)


def parse_args(args=None) -> argparse.Namespace:
    """Parse command line arguments.

    Parameters
    ----------
    args : list, optional
        Command line arguments (uses sys.argv if None)

    Returns
    -------
    argparse.Namespace
        Parsed arguments
    """
    parser = argparse.ArgumentParser(
        description="Generate diagnostic report for cluster refinement",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python -m celltype_refinery.core.diagnosis \\
      --input output/stage_i/refined_final.h5ad \\
      --output output/stage_i

  # With custom thresholds
  python -m celltype_refinery.core.diagnosis \\
      --input output/stage_i/refined_final.h5ad \\
      --output output/stage_i \\
      --min-cells 2000 \\
      --score-threshold 1.3
""",
    )

    # Input/Output
    parser.add_argument(
        "--input", "-i",
        type=Path,
        required=True,
        help="Path to refined h5ad file (AnnData with cluster_lvl0/cluster_lvl1 columns)",
    )
    parser.add_argument(
        "--output", "-o",
        type=Path,
        required=True,
        help="Output directory for diagnostic_report.csv",
    )

    # Optional explicit paths (auto-detected from output dir if not provided)
    parser.add_argument(
        "--cluster-annotations",
        type=Path,
        default=None,
        help="Path to cluster_annotations.csv (default: <output>/cluster_annotations.csv)",
    )
    parser.add_argument(
        "--marker-scores",
        type=Path,
        default=None,
        help="Path to marker_scores.csv (default: <output>/marker_scores.csv)",
    )

    # Criteria configuration
    parser.add_argument(
        "--min-cells",
        type=int,
        default=500,
        help="Minimum cells for subclustering (default: 500)",
    )
    parser.add_argument(
        "--score-threshold",
        type=float,
        default=1.0,
        help="Score threshold for low confidence detection (default: 1.0)",
    )
    parser.add_argument(
        "--subtype-signal-threshold",
        type=float,
        default=1.0,
        help="Threshold for subtype signal detection (default: 1.0)",
    )
    parser.add_argument(
        "--heterogeneity-gap",
        type=float,
        default=0.5,
        help="Gap threshold for heterogeneity detection (default: 0.5)",
    )

    # Output options
    parser.add_argument(
        "--report-name",
        type=str,
        default="diagnostic_report.csv",
        help="Name for output report file (default: diagnostic_report.csv)",
    )

    # Groups export options (for workflow integration)
    parser.add_argument(
        "--export-groups",
        action="store_true",
        help="Export groups_derived.yaml and groups/ directory structure after diagnosis",
    )
    parser.add_argument(
        "--marker-map",
        type=Path,
        default=None,
        help="Path to marker map JSON (for groups_derived.yaml metadata)",
    )

    # Logging
    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
        help="Logging level (default: INFO)",
    )
    parser.add_argument(
        "--quiet", "-q",
        action="store_true",
        help="Suppress non-essential output",
    )

    return parser.parse_args(args)


def main(args=None) -> int:
    """Main CLI entry point.

    Parameters
    ----------
    args : list, optional
        Command line arguments (uses sys.argv if None)

    Returns
    -------
    int
        Exit code (0 for success, non-zero for error)
    """
    parsed = parse_args(args)

    # Setup logging
    log_level = "WARNING" if parsed.quiet else parsed.log_level
    logger = setup_logging(log_level)

    if not parsed.quiet:
        logger.info("=" * 70)
        logger.info("Diagnosis CLI: Generating diagnostic report")
        logger.info("=" * 70)

    # Check scanpy
    if sc is None:
        logger.error("scanpy is required but not installed")
        return 1

    # Validate input
    if not parsed.input.exists():
        logger.error("Input file not found: %s", parsed.input)
        return 1

    # Auto-detect paths if not provided
    output_dir = parsed.output
    output_dir.mkdir(parents=True, exist_ok=True)

    cluster_annotations_path = parsed.cluster_annotations or output_dir / "cluster_annotations.csv"
    marker_scores_path = parsed.marker_scores or output_dir / "marker_scores.csv"

    # Load h5ad
    logger.info("Loading h5ad: %s", parsed.input)
    try:
        adata = sc.read_h5ad(parsed.input)
    except Exception as e:
        logger.error("Failed to load h5ad: %s", e)
        return 1

    logger.info("  Cells: %d", adata.n_obs)

    # Determine cluster column
    cluster_col = "cluster_lvl1" if "cluster_lvl1" in adata.obs else "cluster_lvl0"
    logger.info("  Cluster column: %s", cluster_col)

    # Load annotations and scores
    if not cluster_annotations_path.exists():
        logger.error("Cluster annotations not found: %s", cluster_annotations_path)
        logger.error("Hint: Ensure the refinement step has been run first")
        return 1

    if not marker_scores_path.exists():
        logger.error("Marker scores not found: %s", marker_scores_path)
        logger.error("Hint: Ensure the refinement step has been run first")
        return 1

    try:
        cluster_annotations = pd.read_csv(cluster_annotations_path)
        marker_scores = pd.read_csv(marker_scores_path)
    except Exception as e:
        logger.error("Failed to load CSV files: %s", e)
        return 1

    logger.info("Loaded %d cluster annotations, %d marker score records",
                len(cluster_annotations), len(marker_scores))

    # Validate required columns
    required_annot_cols = ["cluster_id", "assigned_label", "assigned_score", "n_cells"]
    missing_annot = [c for c in required_annot_cols if c not in cluster_annotations.columns]
    if missing_annot:
        logger.error("Missing columns in cluster_annotations: %s", missing_annot)
        return 1

    required_score_cols = ["cluster_id", "label", "score"]
    missing_score = [c for c in required_score_cols if c not in marker_scores.columns]
    if missing_score:
        logger.error("Missing columns in marker_scores: %s", missing_score)
        return 1

    # Derive root_label and is_ambiguous_root if needed (for --export-groups)
    if parsed.export_groups:
        if "root_label" not in cluster_annotations.columns:
            cluster_annotations["root_label"] = cluster_annotations["assigned_label"].apply(
                derive_root_label
            )
            logger.info("Derived root_label column from assigned_label")

        if "is_ambiguous_root" not in cluster_annotations.columns:
            cluster_annotations["is_ambiguous_root"] = cluster_annotations["assigned_label"].apply(
                is_ambiguous_root
            )
            logger.info("Derived is_ambiguous_root column from assigned_label")

    # Create diagnostic engine
    config = CriteriaConfig(
        min_cells=parsed.min_cells,
        score_threshold=parsed.score_threshold,
        subtype_signal_threshold=parsed.subtype_signal_threshold,
        heterogeneity_gap=parsed.heterogeneity_gap,
    )
    engine = DiagnosticEngine(config=config, logger=logger)

    logger.info("Criteria config:")
    logger.info("  min_cells: %d", config.min_cells)
    logger.info("  score_threshold: %.2f", config.score_threshold)
    logger.info("  subtype_signal_threshold: %.2f", config.subtype_signal_threshold)
    logger.info("  heterogeneity_gap: %.2f", config.heterogeneity_gap)

    # Get eligible clusters from h5ad
    eligible_clusters = set(adata.obs[cluster_col].astype(str).unique())
    logger.info("Found %d eligible clusters in h5ad", len(eligible_clusters))

    # Compute marker heterogeneity for weak leaf detection
    # This is critical for identifying clusters that need further subclustering
    logger.info("")
    logger.info("Computing cluster marker heterogeneity for weak leaf detection...")
    heterogeneity_metrics = None
    try:
        heterogeneity_metrics = compute_cluster_marker_heterogeneity(
            adata,
            cluster_col=cluster_col,
            layer="batchcorr" if "batchcorr" in adata.layers else None,
            logger=logger,
        )
        logger.info(
            "Heterogeneity metrics computed for %d clusters (mean=%.3f, max=%.3f)",
            len(heterogeneity_metrics),
            heterogeneity_metrics["marker_heterogeneity"].mean(),
            heterogeneity_metrics["marker_heterogeneity"].max(),
        )
    except Exception as e:
        logger.warning("Failed to compute heterogeneity metrics: %s", e)
        logger.warning("Weak leaf detection will be disabled.")

    # Generate diagnostic report
    logger.info("")
    logger.info("Generating diagnostic report...")
    try:
        diagnostic_report = engine.diagnose(
            cluster_annotations=cluster_annotations,
            marker_scores=marker_scores,
            eligible_clusters=eligible_clusters,
            heterogeneity_metrics=heterogeneity_metrics,
        )
    except Exception as e:
        logger.error("Failed to generate diagnostic report: %s", e)
        return 1

    # Save report
    report_path = output_dir / parsed.report_name
    diagnostic_report.to_csv(report_path, index=False)
    logger.info("Saved diagnostic report: %s (%d clusters)", report_path, len(diagnostic_report))

    # Print summary
    rec_counts = diagnostic_report["recommendation"].value_counts()
    logger.info("")
    logger.info("Recommendation Summary:")
    for rec in ["SUBCLUSTER", "RELABEL", "SKIP"]:
        count = rec_counts.get(rec, 0)
        logger.info("  %s: %d clusters", rec, count)

    logger.info("")
    logger.info("Total clusters: %d", len(diagnostic_report))

    # Check for blocked recommendations
    if "is_blocked_by_size" in diagnostic_report.columns:
        n_blocked = diagnostic_report["is_blocked_by_size"].sum()
        if n_blocked > 0:
            logger.info("Note: %d clusters blocked by min_cells=%d", n_blocked, parsed.min_cells)

    # Update workflow_state.yaml with diagnostic summary
    # This ensures the state file has correct counts after diagnose runs
    # (matches reference behavior in ft/src/workflow/state.py:mark_diagnose_complete)
    diagnostic_summary = {
        "SUBCLUSTER": int(rec_counts.get("SUBCLUSTER", 0)),
        "RELABEL": int(rec_counts.get("RELABEL", 0)),
        "SKIP": int(rec_counts.get("SKIP", 0)),
    }
    update_workflow_state(
        output_dir=output_dir,
        diagnostic_summary=diagnostic_summary,
        diagnostic_report_path=str(report_path),
        logger=logger,
    )

    # Export groups structure if requested (--export-groups)
    # This generates groups_derived.yaml and groups/ directory with cluster_label_mapping.csv
    # Using FRESH diagnostic_report ensures subcluster_ids are correctly populated
    if parsed.export_groups:
        logger.info("")
        logger.info("Exporting groups structure...")

        try:
            from ..annotation.export import export_groups_structure
            from ..refinement.__main__ import export_mapping_table

            # ISSUE-003t fix: Use export_mapping_table to generate cluster_label_mapping
            # with the correct schema: cluster_id, cell_type, n_cells, reason
            # This matches the reference implementation in ft/src/refinement/export.py
            # The reason column comes from curation_reason in adata.obs (most common per cluster)
            cluster_label_mapping = export_mapping_table(
                adata,
                output_path=output_dir / "cluster_label_mapping.csv",
                label_key="cell_type_lvl1",
                reason_key="curation_reason",
            )

            result = export_groups_structure(
                cluster_annotations=cluster_annotations,
                cluster_label_mapping=cluster_label_mapping,
                output_dir=output_dir,
                diagnostic_report=diagnostic_report,  # FRESH data - correct IDs!
                marker_map_path=str(parsed.marker_map) if parsed.marker_map else None,
                iteration=1,
                logger=logger,
            )

            logger.info(
                "Groups exported: %d groups → %s",
                result["n_groups"],
                result["groups_derived_path"],
            )
        except Exception as e:
            logger.error("Failed to export groups structure: %s", e)
            # Non-fatal - continue even if groups export fails
            logger.warning("Continuing without groups export")

    return 0


if __name__ == "__main__":
    sys.exit(main())
