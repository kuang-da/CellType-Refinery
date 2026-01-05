"""CLI entry point for the consolidation module.

Usage:
    python -m celltype_refinery.core.consolidation --help

Example:
    python -m celltype_refinery.core.consolidation \\
        --input out/stage_i/refined_final.h5ad \\
        --diagnostic out/stage_i/diagnostic_report.csv \\
        --marker-scores out/stage_i/marker_scores.csv \\
        --config config/consolidation.yaml \\
        --enable-orphan-rescue \\
        --out out/stage_n
"""

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

try:
    import scanpy as sc
except ImportError:
    print("Error: scanpy is required. Install with: pip install scanpy")
    sys.exit(1)

from celltype_refinery.core.consolidation.config import ConsolidationConfig
from celltype_refinery.core.consolidation.engine import ConsolidationEngine
from celltype_refinery.core.consolidation.export import export_all
from celltype_refinery.core.consolidation.harmonize import HarmonizeConfig
from celltype_refinery.io.logging import get_logger

# Project paths
PROJECT_ROOT = Path(__file__).resolve().parents[4]  # celltype_refinery package root
DEFAULT_LOG_DIR = Path("log")
LOG_FILENAME = "stage_n_consolidation.log"


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Final consolidation of cell-type annotations",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic usage (no orphan rescue, no config)
    python -m celltype_refinery.core.consolidation \\
        --input out/stage_i/refined_final.h5ad \\
        --diagnostic out/stage_i/diagnostic_report.csv \\
        --out out/stage_n

    # With orphan rescue
    python -m celltype_refinery.core.consolidation \\
        --input out/stage_i/refined_final.h5ad \\
        --diagnostic out/stage_i/diagnostic_report.csv \\
        --marker-scores out/stage_i/marker_scores.csv \\
        --enable-orphan-rescue \\
        --out out/stage_n

    # With config file
    python -m celltype_refinery.core.consolidation \\
        --input out/stage_i/refined_final.h5ad \\
        --diagnostic out/stage_i/diagnostic_report.csv \\
        --marker-scores out/stage_i/marker_scores.csv \\
        --config config/consolidation.yaml \\
        --enable-orphan-rescue \\
        --out out/stage_n

    # Preview mapping without writing H5AD
    python -m celltype_refinery.core.consolidation \\
        --input out/stage_i/refined_final.h5ad \\
        --diagnostic out/stage_i/diagnostic_report.csv \\
        --dry-run \\
        --out out/stage_n
        """,
    )

    # Required arguments
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Path to input AnnData (refined_final.h5ad from Stage I)",
    )
    parser.add_argument(
        "--diagnostic",
        type=Path,
        required=True,
        help="Path to diagnostic_report.csv from Stage I",
    )

    # Optional inputs
    parser.add_argument(
        "--marker-scores",
        type=Path,
        help="Path to marker_scores.csv (required for orphan rescue)",
    )
    parser.add_argument(
        "--config",
        type=Path,
        help="Path to YAML configuration file with manual overrides",
    )
    parser.add_argument(
        "--harmonize-config",
        type=Path,
        help="Path to YAML config for two-level label harmonization (fine + broad vocab)",
    )

    # Output
    parser.add_argument(
        "--out",
        type=Path,
        required=True,
        help="Output directory",
    )

    # Options
    parser.add_argument(
        "--enable-orphan-rescue",
        action="store_true",
        help="Enable orphan detection and rescue (requires --marker-scores)",
    )
    parser.add_argument(
        "--enable-iel-rescue",
        action="store_true",
        help="Enable IEL (intraepithelial immune) detection and rescue (requires --marker-scores)",
    )
    parser.add_argument(
        "--orphan-suffix",
        type=str,
        default="(orphan)",
        help="Suffix for rescued orphan labels (default: '(orphan)')",
    )
    parser.add_argument(
        "--output-col",
        type=str,
        default="cell_type_phenocycler",
        help="Output column name for final labels (default: cell_type_phenocycler)",
    )
    parser.add_argument(
        "--cluster-col",
        type=str,
        default="cluster_lvl1",
        help="Cluster column in AnnData (default: cluster_lvl1)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Preview mapping without modifying AnnData",
    )
    parser.add_argument(
        "--no-save-h5ad",
        action="store_true",
        help="Skip saving consolidated H5AD file",
    )
    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
        help="Logging level (default: INFO)",
    )
    parser.add_argument(
        "--log-dir",
        type=Path,
        default=DEFAULT_LOG_DIR,
        help=f"Directory for log files (default: {DEFAULT_LOG_DIR})",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable verbose (DEBUG) logging",
    )

    return parser.parse_args()


def main() -> int:
    """Main entry point."""
    args = parse_args()

    # Configure logging with file handler
    log_level = logging.DEBUG if args.verbose else getattr(logging, args.log_level)
    log_dir = Path(args.log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    logger, actual_log_path = get_logger(
        name="stage_n_consolidation",
        log_path=log_dir / LOG_FILENAME,
        level=log_level,
    )

    # Add console handler for dual logging (file + stdout)
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(log_level)
    console_formatter = logging.Formatter(
        fmt="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)

    print(f"[INFO] Log file: {actual_log_path}")
    logger.info("=== Cell Type Consolidation (Stage N) ===")
    logger.info(f"Input: {args.input}")
    logger.info(f"Diagnostic: {args.diagnostic}")
    logger.info(f"Output: {args.out}")

    # Validate inputs
    if not args.input.exists():
        logger.error(f"Input file not found: {args.input}")
        return 1

    if not args.diagnostic.exists():
        logger.error(f"Diagnostic file not found: {args.diagnostic}")
        return 1

    if args.enable_orphan_rescue and not args.marker_scores:
        logger.error("--marker-scores is required when --enable-orphan-rescue is set")
        return 1

    if args.enable_iel_rescue and not args.marker_scores:
        logger.error("--marker-scores is required when --enable-iel-rescue is set")
        return 1

    if args.marker_scores and not args.marker_scores.exists():
        logger.error(f"Marker scores file not found: {args.marker_scores}")
        return 1

    # Load config if provided
    config = None
    if args.config:
        if not args.config.exists():
            logger.error(f"Config file not found: {args.config}")
            return 1
        logger.info(f"Loading config from {args.config}")
        config = ConsolidationConfig.from_yaml(args.config)
        logger.info(config.summary())

    # Load harmonize config
    harmonize_config = None
    if args.harmonize_config:
        if not args.harmonize_config.exists():
            logger.warning(f"Harmonize config not found: {args.harmonize_config}, using defaults")
            harmonize_config = HarmonizeConfig.default()
        else:
            logger.info(f"Loading harmonize config from {args.harmonize_config}")
            harmonize_config = HarmonizeConfig.from_yaml(args.harmonize_config)
        logger.info(f"  Harmonize version: {harmonize_config.version}")
    else:
        # Use default harmonization
        harmonize_config = HarmonizeConfig.default()
        logger.info("Using default harmonize config")

    # Load inputs
    logger.info("Loading AnnData...")
    adata = sc.read_h5ad(args.input)
    logger.info(f"  Loaded {len(adata):,} cells")

    logger.info("Loading diagnostic report...")
    diagnostic_report = pd.read_csv(args.diagnostic)
    logger.info(f"  Loaded {len(diagnostic_report):,} clusters")

    marker_scores = None
    if args.marker_scores:
        logger.info("Loading marker scores...")
        marker_scores = pd.read_csv(args.marker_scores)
        logger.info(f"  Loaded {len(marker_scores):,} score entries")

    # Create engine
    engine = ConsolidationEngine(
        config=config,
        output_col=args.output_col,
        enable_orphan_rescue=args.enable_orphan_rescue,
        orphan_suffix=args.orphan_suffix,
        enable_iel_rescue=args.enable_iel_rescue,
        harmonize_config=harmonize_config,
    )

    # Dry run mode - just preview mapping
    if args.dry_run:
        logger.info("=== DRY RUN MODE ===")
        mapping_df = engine.get_mapping_table(diagnostic_report, marker_scores)

        # Print summary
        logger.info("\nLabel distribution preview:")
        label_counts = mapping_df.groupby("final_label")["n_cells"].sum().sort_values(ascending=False)
        for label, count in label_counts.head(20).items():
            logger.info(f"  {label}: {count:,}")

        # Save mapping table
        args.out.mkdir(parents=True, exist_ok=True)
        mapping_path = args.out / "mapping_preview.csv"
        mapping_df.to_csv(mapping_path, index=False)
        logger.info(f"\nMapping preview saved to {mapping_path}")

        return 0

    # Execute consolidation
    logger.info("Executing consolidation...")
    result = engine.execute(
        adata,
        diagnostic_report,
        marker_scores=marker_scores,
        cluster_col=args.cluster_col,
    )

    if not result.success:
        logger.error("Consolidation failed!")
        for error in result.errors:
            logger.error(f"  {error}")
        return 1

    # Print summary
    logger.info("\n=== CONSOLIDATION SUMMARY ===")
    logger.info(f"Total cells: {result.n_cells_total:,}")
    logger.info(f"Total clusters: {result.n_clusters_total}")
    logger.info(f"Overrides applied: {result.n_overrides_applied}")
    logger.info(f"Relabels applied: {result.n_relabels_applied}")
    logger.info(f"Orphans rescued: {result.n_orphans_rescued}")
    logger.info(f"Orphans flagged: {result.n_orphans_flagged}")
    logger.info(f"IELs rescued: {result.n_iel_rescued}")
    logger.info(f"Execution time: {result.execution_time_seconds:.2f}s")

    logger.info("\nCategory breakdown:")
    for cat, count in result.category_breakdown.items():
        pct = count / result.n_cells_total * 100 if result.n_cells_total > 0 else 0
        logger.info(f"  {cat}: {count:,} ({pct:.1f}%)")

    logger.info("\nTop 15 cell types:")
    sorted_labels = sorted(result.label_distribution.items(), key=lambda x: -x[1])
    for label, count in sorted_labels[:15]:
        pct = count / result.n_cells_total * 100 if result.n_cells_total > 0 else 0
        logger.info(f"  {label}: {count:,} ({pct:.1f}%)")

    if result.warnings:
        logger.warning("\nWarnings:")
        for warning in result.warnings:
            logger.warning(f"  {warning}")

    # Get mapping table for export
    mapping_df = engine.get_mapping_table(diagnostic_report, marker_scores)

    # Export all outputs
    logger.info("\nExporting outputs...")
    outputs = export_all(
        adata,
        result,
        config,
        mapping_df,
        args.input,
        args.out,
        output_col=args.output_col,
    )

    for name, path in outputs.items():
        logger.info(f"  {name}: {path}")

    # Save H5AD
    if not args.no_save_h5ad:
        h5ad_path = args.out / "consolidated.h5ad"
        logger.info(f"Saving AnnData to {h5ad_path}...")
        adata.write_h5ad(h5ad_path)
        logger.info("Done!")

    return 0


if __name__ == "__main__":
    sys.exit(main())
