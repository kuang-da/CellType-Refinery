"""Stage J Pool – Lineage pooling for cell-type refinement.

This stage pools unresolved clusters by lineage to enable further subclustering
in subsequent Stage I iterations.

Relabeling Rules:
1. Root cell types (Epithelium, Endothelium, etc.) -> Pool_<RootType>
2. Ambiguous types (~) -> Pool_<CanonicalLabel> (alphabetically ordered)
3. SUBCLUSTER carry-over -> Keep original label for next iteration

Usage Examples:
    # Dry-run: show pooling plan without applying
    python -m celltype_refinery.core.pooling \\
        --input output/stage_i/refined.h5ad \\
        --marker-map data/marker_map.json \\
        --dry-run

    # Execute pooling
    python -m celltype_refinery.core.pooling \\
        --input output/stage_i/refined.h5ad \\
        --marker-map data/marker_map.json \\
        --out output/stage_j
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from datetime import datetime
from pathlib import Path
from typing import Optional

import pandas as pd

try:
    import scanpy as sc
except ImportError:
    sc = None

try:
    import yaml
except ImportError:
    yaml = None

from .engine import PoolingEngine, PoolingResult


def setup_logger(
    name: str,
    log_path: Optional[Path] = None,
    level: int = logging.INFO,
) -> logging.Logger:
    """Setup logger with file and console handlers.

    Parameters
    ----------
    name : str
        Logger name
    log_path : Path, optional
        Path to log file
    level : int
        Logging level

    Returns
    -------
    logging.Logger
        Configured logger
    """
    logger = logging.getLogger(name)
    logger.setLevel(level)

    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(level)
    console_formatter = logging.Formatter(
        fmt="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)

    # File handler
    if log_path is not None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_path)
        file_handler.setLevel(level)
        file_handler.setFormatter(console_formatter)
        logger.addHandler(file_handler)

    return logger


def print_pooling_plan(
    pool_summary: pd.DataFrame,
    logger: logging.Logger,
) -> None:
    """Print formatted pooling plan to logger.

    Parameters
    ----------
    pool_summary : pd.DataFrame
        Pool summary from PoolingEngine.execute()
    logger : logging.Logger
        Logger to print to
    """
    if pool_summary.empty:
        logger.info("No clusters to pool.")
        return

    # Split into pooled and unchanged
    pooled = pool_summary[pool_summary["is_pooled"] == True]
    unchanged = pool_summary[pool_summary["is_pooled"] == False]

    logger.info("=" * 100)
    logger.info("POOLING PLAN: Lineage Consolidation")
    logger.info("=" * 100)

    # Show pools
    if not pooled.empty:
        logger.info("")
        logger.info("POOLS TO CREATE (%d clusters -> %d pools):",
                    len(pooled), pooled["new_cluster_id"].nunique())
        logger.info("-" * 100)

        # Group by pool
        for pool_label in sorted(pooled["pool_label"].unique()):
            pool_rows = pooled[pooled["pool_label"] == pool_label]
            pool_id = pool_rows["new_cluster_id"].iloc[0]
            total_cells = pool_rows["n_cells"].sum()
            n_source = len(pool_rows)

            logger.info("")
            logger.info("  %s (%s): %d source clusters, %d cells",
                        pool_label, pool_id, n_source, total_cells)

            for _, row in pool_rows.iterrows():
                logger.info("    <- %s (%s): %d cells",
                            row["original_label"], row["cluster_id"], row["n_cells"])

    # Show unchanged
    if not unchanged.empty:
        logger.info("")
        logger.info("UNCHANGED CLUSTERS (%d):", len(unchanged))
        logger.info("-" * 100)

        # Group by recommendation
        n_subcluster = (unchanged["reason"].str.contains("SUBCLUSTER")).sum()
        n_subtype = len(unchanged) - n_subcluster

        if n_subcluster > 0:
            logger.info("  SUBCLUSTER carry-over: %d clusters", n_subcluster)
        if n_subtype > 0:
            logger.info("  Subtype retained: %d clusters", n_subtype)

    # Summary
    logger.info("")
    logger.info("=" * 100)
    logger.info("SUMMARY:")
    logger.info("  Pools to create: %d", pooled["new_cluster_id"].nunique() if not pooled.empty else 0)
    logger.info("  Clusters pooled: %d", len(pooled))
    logger.info("  Clusters unchanged: %d", len(unchanged))
    logger.info("  Total cells pooled: %d", pooled["n_cells"].sum() if not pooled.empty else 0)
    logger.info("=" * 100)


def save_workflow_state(
    out_dir: Path,
    input_path: Path,
    marker_map_path: Path,
    result: PoolingResult,
    diagnostic_summary: dict,
) -> None:
    """Save workflow state for Stage I compatibility.

    Parameters
    ----------
    out_dir : Path
        Output directory
    input_path : Path
        Input AnnData path
    marker_map_path : Path
        Marker map path
    result : PoolingResult
        Pooling result
    diagnostic_summary : dict
        Summary of diagnostic recommendations
    """
    if yaml is None:
        return

    now = datetime.now().isoformat()
    state = {
        "version": "1.0",
        "stage": "J",
        "marker_map_path": str(marker_map_path),
        "input_path": str(input_path),
        "stage_h_dir": str(out_dir),
        "out_dir": str(out_dir),
        "created_at": now,
        "updated_at": now,
        "diagnose_complete": True,
        "diagnose_timestamp": now,
        "diagnostic_summary": diagnostic_summary,
        "diagnostic_report_path": str(out_dir / "diagnostic_report.csv"),
        "execute_complete": False,
        "execute_timestamp": None,
        "final_output_path": str(out_dir / "refined_final.h5ad"),
        "pooling_complete": result.success,
        "pooling_summary": {
            "n_pools_created": int(result.n_pools_created),
            "n_cells_pooled": int(result.n_cells_pooled),
            "n_cells_unchanged": int(result.n_cells_unchanged),
            "n_clusters_pooled": int(result.n_clusters_pooled),
            "n_clusters_unchanged": int(result.n_clusters_unchanged),
        },
    }

    with open(out_dir / "workflow_state.yaml", "w") as f:
        yaml.dump(state, f, default_flow_style=False, sort_keys=False)


def export_composition(
    adata: "sc.AnnData",
    out_dir: Path,
    cell_type_col: str,
    logger: logging.Logger,
) -> dict:
    """Export composition statistics.

    Parameters
    ----------
    adata : sc.AnnData
        AnnData object
    out_dir : Path
        Output directory
    cell_type_col : str
        Cell type column name
    logger : logging.Logger
        Logger instance

    Returns
    -------
    dict
        Paths to exported files
    """
    obs = adata.obs.copy()
    obs[cell_type_col] = obs[cell_type_col].astype(str)

    # Global composition
    global_comp = obs[cell_type_col].value_counts().reset_index()
    global_comp.columns = ["cell_type", "n_cells"]
    global_comp["proportion"] = global_comp["n_cells"] / len(obs)

    global_path = out_dir / "composition_global.csv"
    global_comp.to_csv(global_path, index=False)
    logger.info(f"  Saved: {global_path}")

    # Sample composition (if sample_id exists)
    if "sample_id" in obs.columns:
        sample_comp = obs.groupby(["sample_id", cell_type_col]).size().reset_index(name="n_cells")
        sample_path = out_dir / "composition_by_sample.csv"
        sample_comp.to_csv(sample_path, index=False)
        logger.info(f"  Saved: {sample_path}")

    return {"global": str(global_path)}


def generate_groups_derived(
    cluster_annotations: pd.DataFrame,
    adata: "sc.AnnData",
    label_col: str,
) -> dict:
    """Generate groups_derived.yaml for Stage I iterate workflow.

    Groups clusters by root_label for organized processing:
    1. Epithelium
    2. Endothelium
    3. Mesenchymal Cells
    4. Immune Cells
    5. Hybrids (ambiguous root clusters)
    6. Unassigned

    Parameters
    ----------
    cluster_annotations : pd.DataFrame
        Cluster annotations with root_label column
    adata : sc.AnnData
        AnnData object for cell count validation
    label_col : str
        Cell type label column name

    Returns
    -------
    dict
        Groups configuration for Stage I
    """
    # Standard group ordering
    GROUP_ORDER = [
        "Epithelium",
        "Endothelium",
        "Mesenchymal Cells",
        "Immune Cells",
        "Hybrids",
        "Unassigned",
    ]

    groups = {}
    group_idx = 1

    # Handle case where root_label may not exist
    if "root_label" not in cluster_annotations.columns:
        # Fall back to parsing from assigned_label
        cluster_annotations = cluster_annotations.copy()
        cluster_annotations["root_label"] = cluster_annotations["assigned_label"].apply(
            lambda x: x.split("~")[0] if "~" in str(x) else (
                str(x)[5:] if str(x).startswith("Pool_") else str(x)
            )
        )

    # Check for ambiguous root flag
    has_ambiguous_flag = "is_ambiguous_root" in cluster_annotations.columns

    # Categorize clusters
    cluster_groups = {}
    for _, row in cluster_annotations.iterrows():
        cluster_id = str(row["cluster_id"])
        root_label = str(row.get("root_label", ""))
        is_ambiguous = bool(row.get("is_ambiguous_root", False)) if has_ambiguous_flag else ("~" in str(row.get("assigned_label", "")))

        if is_ambiguous:
            group_name = "Hybrids"
        elif root_label in GROUP_ORDER:
            group_name = root_label
        elif root_label.startswith("Pool_"):
            # Extract root from Pool_X
            extracted = root_label[5:]
            if "~" in extracted:
                group_name = "Hybrids"
            elif extracted in GROUP_ORDER:
                group_name = extracted
            else:
                group_name = "Unassigned"
        else:
            group_name = "Unassigned"

        if group_name not in cluster_groups:
            cluster_groups[group_name] = []
        cluster_groups[group_name].append(cluster_id)

    # Build groups in standard order
    for group_name in GROUP_ORDER:
        if group_name in cluster_groups:
            cluster_ids = sorted(cluster_groups[group_name])
            n_cells = 0
            for cid in cluster_ids:
                ann_row = cluster_annotations[cluster_annotations["cluster_id"].astype(str) == cid]
                if len(ann_row) > 0:
                    n_cells += int(ann_row["n_cells"].iloc[0])

            groups[f"g{group_idx}_{group_name.lower().replace(' ', '_')}"] = {
                "name": group_name,
                "cluster_ids": cluster_ids,
                "n_clusters": len(cluster_ids),
                "n_cells": n_cells,
                "subcluster_only": True,  # Stage J pools only need subclustering
            }
            group_idx += 1

    return {
        "version": "1.0",
        "source": "Stage J pooling",
        "groups": groups,
    }


def parse_args(argv: Optional[list] = None) -> argparse.Namespace:
    """Parse command line arguments.

    Parameters
    ----------
    argv : list, optional
        Command line arguments (defaults to sys.argv)

    Returns
    -------
    argparse.Namespace
        Parsed arguments
    """
    parser = argparse.ArgumentParser(
        description="Stage J: Pool clusters by lineage for refined subclustering",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Dry-run: show pooling plan
    python -m celltype_refinery.core.pooling \\
        --input output/stage_i/refined.h5ad \\
        --marker-map data/marker_map.json \\
        --dry-run

    # Execute pooling
    python -m celltype_refinery.core.pooling \\
        --input output/stage_i/refined.h5ad \\
        --marker-map data/marker_map.json \\
        --out output/stage_j
        """,
    )

    # Required inputs
    parser.add_argument(
        "--input", "-i",
        type=Path,
        required=True,
        help="Input AnnData from Stage I (refined.h5ad)",
    )
    parser.add_argument(
        "--marker-map", "-m",
        type=Path,
        required=True,
        help="Marker map JSON for root type extraction",
    )
    parser.add_argument(
        "--out", "-o",
        type=Path,
        default=Path("output/stage_j"),
        help="Output directory (default: output/stage_j)",
    )

    # Optional inputs (auto-detected from input directory)
    parser.add_argument(
        "--diagnostic-report",
        type=Path,
        help="Diagnostic report CSV (auto-detected from input dir if not specified)",
    )
    parser.add_argument(
        "--cluster-annotations",
        type=Path,
        help="Cluster annotations CSV (auto-detected from input dir if not specified)",
    )
    parser.add_argument(
        "--marker-scores",
        type=Path,
        help="Marker scores CSV (auto-detected from input dir if not specified)",
    )
    parser.add_argument(
        "--cluster-annotations-enhanced",
        type=Path,
        help="Enhanced cluster annotations CSV (auto-detected from input dir if not specified)",
    )

    # Execution mode
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show pooling plan without applying changes",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing output",
    )

    # Column names
    parser.add_argument(
        "--cluster-col",
        default="cluster_lvl1",
        help="Cluster column to read/write (default: cluster_lvl1)",
    )
    parser.add_argument(
        "--label-col",
        default="cell_type_lvl1",
        help="Label column to read/write (default: cell_type_lvl1)",
    )

    # Logging
    parser.add_argument(
        "--log-dir",
        type=Path,
        help="Log directory (optional)",
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable verbose logging",
    )

    return parser.parse_args(argv)


def main(argv: Optional[list] = None) -> int:
    """Main entry point for Stage J pooling.

    Parameters
    ----------
    argv : list, optional
        Command line arguments

    Returns
    -------
    int
        Exit code (0 for success, non-zero for error)
    """
    args = parse_args(argv)

    # Setup logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    log_path = None
    if args.log_dir is not None:
        log_path = args.log_dir / "stage_j_pool.log"
    logger = setup_logger("stage_j_pool", log_path=log_path, level=log_level)

    logger.info("=" * 72)
    logger.info("STAGE J: Lineage Pooling for Cell-Type Refinement")
    logger.info("=" * 72)
    logger.info(f"Timestamp: {datetime.now().isoformat()}")
    logger.info(f"Input: {args.input}")
    logger.info(f"Marker map: {args.marker_map}")
    logger.info(f"Output: {args.out}")
    logger.info(f"Dry-run: {args.dry_run}")

    # Validate inputs
    if not args.input.exists():
        logger.error(f"Input file not found: {args.input}")
        return 1

    if not args.marker_map.exists():
        logger.error(f"Marker map not found: {args.marker_map}")
        return 1

    # Auto-detect diagnostic report and cluster annotations
    input_dir = args.input.parent

    diagnostic_report_path = args.diagnostic_report
    if diagnostic_report_path is None:
        diagnostic_report_path = input_dir / "diagnostic_report.csv"
        if not diagnostic_report_path.exists():
            logger.error(f"Diagnostic report not found: {diagnostic_report_path}")
            logger.error("Please specify --diagnostic-report explicitly.")
            return 1

    cluster_annotations_path = args.cluster_annotations
    if cluster_annotations_path is None:
        cluster_annotations_path = input_dir / "cluster_annotations.csv"
        if not cluster_annotations_path.exists():
            # Try cluster_label_mapping.csv as fallback
            cluster_annotations_path = input_dir / "cluster_label_mapping.csv"
            if not cluster_annotations_path.exists():
                logger.error(f"Cluster annotations not found in: {input_dir}")
                logger.error("Please specify --cluster-annotations explicitly.")
                return 1

    # Auto-detect marker_scores.csv
    marker_scores_path = args.marker_scores
    if marker_scores_path is None:
        marker_scores_path = input_dir / "marker_scores.csv"
        if not marker_scores_path.exists():
            logger.warning(f"marker_scores.csv not found in: {input_dir}")
            logger.warning("Stage I smart rescore will not work without this file.")
            marker_scores_path = None

    # Auto-detect cluster_annotations_enhanced.csv
    enhanced_annotations_path = args.cluster_annotations_enhanced
    if enhanced_annotations_path is None:
        enhanced_annotations_path = input_dir / "cluster_annotations_enhanced.csv"
        if not enhanced_annotations_path.exists():
            logger.warning(f"cluster_annotations_enhanced.csv not found in: {input_dir}")
            logger.warning("Stage I full output compatibility may be limited.")
            enhanced_annotations_path = None

    logger.info(f"Diagnostic report: {diagnostic_report_path}")
    logger.info(f"Cluster annotations: {cluster_annotations_path}")
    logger.info(f"Marker scores: {marker_scores_path}")
    logger.info(f"Enhanced annotations: {enhanced_annotations_path}")

    # Load inputs
    logger.info("")
    logger.info("Loading inputs...")

    if sc is None:
        logger.error("scanpy is required but not installed")
        return 1

    adata = sc.read_h5ad(args.input)
    logger.info(f"  Loaded AnnData: {adata.n_obs:,} cells, {adata.n_vars} markers")

    diagnostic_report = pd.read_csv(diagnostic_report_path)
    logger.info(f"  Loaded diagnostic report: {len(diagnostic_report)} clusters")

    cluster_annotations = pd.read_csv(cluster_annotations_path)
    logger.info(f"  Loaded cluster annotations: {len(cluster_annotations)} clusters")

    # Load optional marker scores
    input_marker_scores = None
    if marker_scores_path is not None and marker_scores_path.exists():
        input_marker_scores = pd.read_csv(marker_scores_path)
        logger.info(f"  Loaded marker scores: {len(input_marker_scores)} rows")

    # Load optional enhanced annotations
    input_enhanced_annotations = None
    if enhanced_annotations_path is not None and enhanced_annotations_path.exists():
        input_enhanced_annotations = pd.read_csv(enhanced_annotations_path)
        logger.info(f"  Loaded enhanced annotations: {len(input_enhanced_annotations)} rows")

    # Ensure required columns exist
    required_ann_cols = ["cluster_id", "assigned_label"]
    missing_cols = [c for c in required_ann_cols if c not in cluster_annotations.columns]
    if missing_cols:
        # Try to adapt from cluster_label_mapping format
        if "cell_type" in cluster_annotations.columns:
            cluster_annotations = cluster_annotations.rename(columns={"cell_type": "assigned_label"})
        else:
            logger.error(f"Missing required columns in cluster annotations: {missing_cols}")
            return 1

    if "n_cells" not in cluster_annotations.columns:
        # Compute n_cells from adata
        logger.info("  Computing n_cells from AnnData...")
        cluster_col = args.cluster_col
        if cluster_col in adata.obs.columns:
            cell_counts = adata.obs[cluster_col].astype(str).value_counts().to_dict()
            cluster_annotations["n_cells"] = cluster_annotations["cluster_id"].astype(str).map(cell_counts).fillna(0).astype(int)
        else:
            logger.error(f"Cluster column '{cluster_col}' not found in adata.obs")
            return 1

    # Initialize engine
    engine = PoolingEngine(
        marker_map_path=args.marker_map,
        cluster_col=args.cluster_col,
        label_col=args.label_col,
        logger=logger,
    )

    # Build pooling plan
    logger.info("")
    logger.info("Building pooling plan...")

    # Execute to get the plan
    if args.dry_run:
        # For dry-run, create a copy to avoid modifying the original
        adata_copy = adata.copy()
        result = engine.execute(adata_copy, cluster_annotations, diagnostic_report)
    else:
        result = engine.execute(adata, cluster_annotations, diagnostic_report)

    # Print pooling plan
    print_pooling_plan(result.pool_summary, logger)

    if args.dry_run:
        logger.info("")
        logger.info("DRY-RUN MODE: No changes applied.")
        logger.info("Run without --dry-run to apply pooling.")
        return 0

    # Ensure output directory exists
    if args.out.exists() and not args.force:
        logger.warning(f"Output directory already exists: {args.out}")
        logger.warning("Use --force to overwrite.")
        return 1
    args.out.mkdir(parents=True, exist_ok=True)

    # Generate outputs
    logger.info("")
    logger.info("Generating outputs...")

    # 1. Save refined AnnData
    output_h5ad = args.out / "refined_final.h5ad"
    adata.write_h5ad(output_h5ad)
    logger.info(f"  Saved: {output_h5ad}")

    # 2. Save pooling summary
    pooling_summary_path = args.out / "pooling_summary.csv"
    result.pool_summary.to_csv(pooling_summary_path, index=False)
    logger.info(f"  Saved: {pooling_summary_path}")

    # 3. Generate FULL diagnostic report for next Stage I (27 columns)
    new_diagnostic = engine.generate_full_diagnostic_report(
        adata, diagnostic_report, result.pool_summary,
        input_cluster_annotations=cluster_annotations,
    )
    new_diagnostic_path = args.out / "diagnostic_report.csv"
    new_diagnostic.to_csv(new_diagnostic_path, index=False)
    logger.info(f"  Saved: {new_diagnostic_path} ({len(new_diagnostic.columns)} columns)")

    # 4. Generate FULL cluster annotations (24 columns)
    new_annotations = engine.generate_full_cluster_annotations(
        adata, result.pool_summary,
        input_cluster_annotations=cluster_annotations,
    )
    new_annotations_path = args.out / "cluster_annotations.csv"
    new_annotations.to_csv(new_annotations_path, index=False)
    logger.info(f"  Saved: {new_annotations_path} ({len(new_annotations.columns)} columns)")

    # 4b. Generate enhanced annotations (41 columns) for Stage I compatibility
    new_enhanced = engine.generate_enhanced_annotations(
        adata, result.pool_summary,
        input_enhanced_annotations=input_enhanced_annotations,
    )
    if not new_enhanced.empty:
        enhanced_path = args.out / "cluster_annotations_enhanced.csv"
        new_enhanced.to_csv(enhanced_path, index=False)
        logger.info(f"  Saved: {enhanced_path} ({len(new_enhanced.columns)} columns)")
    else:
        logger.warning("  Skipped: cluster_annotations_enhanced.csv (no data)")

    # 5. Generate marker_scores.csv (for Stage I smart rescore)
    if input_marker_scores is not None:
        new_marker_scores = engine.generate_marker_scores(
            result.pool_summary,
            input_marker_scores,
        )
        marker_scores_output_path = args.out / "marker_scores.csv"
        new_marker_scores.to_csv(marker_scores_output_path, index=False)
        logger.info(f"  Saved: {marker_scores_output_path}")
    else:
        logger.warning("  Skipped: marker_scores.csv (no input file)")

    # 6. Generate cluster label mapping
    obs_subset = adata.obs[[args.cluster_col, args.label_col]].copy()
    obs_subset[args.cluster_col] = obs_subset[args.cluster_col].astype(str)
    obs_subset[args.label_col] = obs_subset[args.label_col].astype(str)

    grouped = obs_subset.groupby(args.cluster_col, observed=True)
    mapping_df = grouped.agg(
        cell_type=(args.label_col, lambda x: x.mode().iloc[0] if len(x) > 0 else "Unknown"),
        n_cells=(args.label_col, "size"),
    ).reset_index().rename(columns={args.cluster_col: "cluster_id"})

    # Sort: pools first
    mapping_df["_sort_key"] = mapping_df["cluster_id"].apply(lambda x: (not x.startswith("P"), x))
    mapping_df = mapping_df.sort_values("_sort_key").drop(columns=["_sort_key"])

    mapping_path = args.out / "cluster_label_mapping.csv"
    mapping_df.to_csv(mapping_path, index=False)
    logger.info(f"  Saved: {mapping_path}")

    # 7. Save workflow state
    rec_counts = new_diagnostic["recommendation"].value_counts().to_dict()
    diagnostic_summary = {
        "total_clusters": len(new_diagnostic),
        "subcluster": rec_counts.get("SUBCLUSTER", 0),
        "relabel": rec_counts.get("RELABEL", 0),
        "skip": rec_counts.get("SKIP", 0),
        "n_pools": result.n_pools_created,
        "cells_pooled": int(result.n_cells_pooled),
        "cells_unchanged": int(result.n_cells_unchanged),
    }
    save_workflow_state(args.out, args.input, args.marker_map, result, diagnostic_summary)
    logger.info(f"  Saved: {args.out / 'workflow_state.yaml'}")

    # 8. Generate groups_derived.yaml for Stage I iterate workflow
    if yaml is not None:
        groups_derived = generate_groups_derived(new_annotations, adata, args.label_col)
        groups_path = args.out / "groups_derived.yaml"
        with open(groups_path, "w") as f:
            yaml.dump(groups_derived, f, default_flow_style=False, sort_keys=False)
        logger.info(f"  Saved: {groups_path}")

    # 9. Generate composition statistics
    logger.info("Generating composition statistics...")
    export_composition(adata, args.out, args.label_col, logger)

    # Final summary
    logger.info("")
    logger.info("=" * 72)
    logger.info("STAGE J COMPLETE")
    logger.info("=" * 72)
    logger.info(f"  Pools created: {result.n_pools_created}")
    logger.info(f"  Cells pooled: {result.n_cells_pooled:,}")
    logger.info(f"  Cells unchanged: {result.n_cells_unchanged:,}")
    logger.info(f"  Output: {args.out}")
    logger.info("")
    logger.info("Next step: Run Stage I on the pooled data:")
    logger.info(f"  python -m celltype_refinery.core.refinement \\")
    logger.info(f"    --input {output_h5ad} \\")
    logger.info(f"    --auto --execute \\")
    logger.info(f"    --out output/stage_i_next")
    logger.info("=" * 72)

    return 0


if __name__ == "__main__":
    sys.exit(main())
