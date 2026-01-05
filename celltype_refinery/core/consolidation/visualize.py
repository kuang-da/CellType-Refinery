"""CLI entry point for consolidation visualizations.

Usage:
    python -m celltype_refinery.core.consolidation.visualize --help

    # Generate all visualizations
    python -m celltype_refinery.core.consolidation.visualize \
        --input out/stage_n/consolidated.h5ad \
        --mapping out/stage_n/cluster_label_mapping.csv \
        --orphans out/stage_n/orphan_candidates.csv \
        --compute-umap \
        --out out/stage_n/figures

    # Static figures only
    python -m celltype_refinery.core.consolidation.visualize \
        --input out/stage_n/consolidated.h5ad \
        --mapping out/stage_n/cluster_label_mapping.csv \
        --static-only \
        --out out/stage_n/figures
"""

import argparse
import logging
import sys
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate visualizations for Stage N consolidation results",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Generate all visualizations (static + interactive)
    python -m celltype_refinery.core.consolidation.visualize \\
        --input out/stage_n/consolidated.h5ad \\
        --mapping out/stage_n/cluster_label_mapping.csv \\
        --compute-umap \\
        --out out/stage_n

    # Static figures only (faster)
    python -m celltype_refinery.core.consolidation.visualize \\
        --input out/stage_n/consolidated.h5ad \\
        --mapping out/stage_n/cluster_label_mapping.csv \\
        --static-only \\
        --out out/stage_n

    # With orphan visualization
    python -m celltype_refinery.core.consolidation.visualize \\
        --input out/stage_n/consolidated.h5ad \\
        --mapping out/stage_n/cluster_label_mapping.csv \\
        --orphans out/stage_n/orphan_candidates.csv \\
        --out out/stage_n
        """,
    )

    # Required arguments
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Path to consolidated AnnData (consolidated.h5ad)",
    )
    parser.add_argument(
        "--mapping",
        type=Path,
        required=True,
        help="Path to cluster_label_mapping.csv",
    )

    # Optional inputs
    parser.add_argument(
        "--orphans",
        type=Path,
        help="Path to orphan_candidates.csv (optional)",
    )
    parser.add_argument(
        "--iels",
        type=Path,
        help="Path to iel_candidates.csv (optional, for IEL visualizations)",
    )
    parser.add_argument(
        "--marker-scores",
        type=Path,
        help="Path to marker_scores.csv (optional, for complete orphan scatter plot)",
    )
    parser.add_argument(
        "--diagnostic",
        type=Path,
        help="Path to diagnostic_report.csv (optional, for complete orphan scatter plot)",
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
        "--compute-umap",
        action="store_true",
        help="Compute UMAP embedding if not present in AnnData",
    )
    parser.add_argument(
        "--use-gpu",
        action="store_true",
        default=True,
        help="Use GPU acceleration for UMAP (default: True)",
    )
    parser.add_argument(
        "--no-gpu",
        action="store_true",
        help="Disable GPU acceleration for UMAP",
    )
    parser.add_argument(
        "--static-only",
        action="store_true",
        help="Only generate static PNG figures (skip interactive HTML)",
    )
    parser.add_argument(
        "--interactive-only",
        action="store_true",
        help="Only generate interactive HTML figures (skip static PNG)",
    )
    parser.add_argument(
        "--cell-type-col",
        type=str,
        default="cell_type_phenocycler",
        help="Column with cell type labels (default: cell_type_phenocycler)",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=200,
        help="DPI for static figures (default: 200)",
    )
    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
        help="Logging level (default: INFO)",
    )

    return parser.parse_args()


def main() -> int:
    """Main entry point."""
    args = parse_args()

    # Configure logging
    logging.getLogger().setLevel(getattr(logging, args.log_level))

    logger.info("=== Stage N Visualization ===")
    logger.info(f"Input: {args.input}")
    logger.info(f"Mapping: {args.mapping}")
    logger.info(f"Output: {args.out}")

    # Validate inputs
    if not args.input.exists():
        logger.error(f"Input file not found: {args.input}")
        return 1

    if not args.mapping.exists():
        logger.error(f"Mapping file not found: {args.mapping}")
        return 1

    if args.orphans and not args.orphans.exists():
        logger.warning(f"Orphan file not found: {args.orphans}")
        args.orphans = None

    if args.iels and not args.iels.exists():
        logger.warning(f"IEL file not found: {args.iels}")
        args.iels = None

    # Import heavy dependencies
    try:
        import scanpy as sc
        import pandas as pd
    except ImportError as e:
        logger.error(f"Missing dependency: {e}")
        return 1

    # Load data
    logger.info("Loading AnnData...")
    adata = sc.read_h5ad(args.input)
    logger.info(f"  Loaded {len(adata):,} cells")

    logger.info("Loading mapping table...")
    mapping_df = pd.read_csv(args.mapping)
    logger.info(f"  Loaded {len(mapping_df):,} cluster mappings")

    orphan_df = None
    if args.orphans:
        logger.info("Loading orphan candidates...")
        orphan_df = pd.read_csv(args.orphans)
        logger.info(f"  Loaded {len(orphan_df):,} orphan candidates")

    iel_df = None
    if args.iels:
        logger.info("Loading IEL candidates...")
        iel_df = pd.read_csv(args.iels)
        logger.info(f"  Loaded {len(iel_df):,} IEL candidates")

    marker_scores_df = None
    if args.marker_scores and args.marker_scores.exists():
        logger.info("Loading marker scores...")
        marker_scores_df = pd.read_csv(args.marker_scores)
        logger.info(f"  Loaded marker scores for {marker_scores_df['cluster_id'].nunique()} clusters")

    diagnostic_df = None
    if args.diagnostic and args.diagnostic.exists():
        logger.info("Loading diagnostic report...")
        diagnostic_df = pd.read_csv(args.diagnostic)
        logger.info(f"  Loaded {len(diagnostic_df):,} cluster diagnostics")

    # Import visualization module
    from celltype_refinery.core.consolidation.viz import generate_all_figures

    # Determine GPU usage
    use_gpu = not args.no_gpu

    # Generate visualizations
    logger.info("Generating visualizations...")
    outputs = generate_all_figures(
        adata=adata,
        mapping_df=mapping_df,
        orphan_df=orphan_df,
        iel_df=iel_df,
        marker_scores=marker_scores_df,
        diagnostic_df=diagnostic_df,
        output_dir=args.out,
        compute_umap=args.compute_umap,
        use_gpu=use_gpu,
        static_only=args.static_only,
        interactive_only=args.interactive_only,
        dpi=args.dpi,
    )

    # Print summary
    logger.info("")
    logger.info("=== VISUALIZATION SUMMARY ===")
    logger.info(f"Generated {len(outputs)} outputs:")

    for name, path in sorted(outputs.items()):
        if path:
            logger.info(f"  {name}: {path}")

    logger.info("")
    logger.info("Done!")

    return 0


if __name__ == "__main__":
    sys.exit(main())
