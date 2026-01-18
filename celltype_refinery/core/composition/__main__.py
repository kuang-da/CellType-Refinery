"""CLI entry point for composition analysis.

Usage:
    python -m celltype_refinery.core.composition --help

    # NEW DEFAULT: Process all 4 cell type columns
    python -m celltype_refinery.core.composition \\
        --input consolidated.h5ad \\
        --out composition_output

    # Process specific columns
    python -m celltype_refinery.core.composition \\
        --input consolidated.h5ad \\
        --cell-type-cols cell_type_phenocycler cell_type_broad \\
        --out composition_output

    # Single-column mode (backward compatible)
    python -m celltype_refinery.core.composition \\
        --input consolidated.h5ad \\
        --cell-type-col cell_type_phenocycler \\
        --out composition_output
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import List, Optional

# Default cell type columns to analyze
DEFAULT_CELL_TYPE_COLUMNS = [
    "cell_type_phenocycler",
    "cell_type_phenocycler_detailed",
    "cell_type_multiomics",
    "cell_type_broad",
]

# Default log filename
LOG_FILENAME = "composition_analysis.log"


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Cell-type composition analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # NEW DEFAULT: Process all 4 cell type columns (multi-column mode)
    python -m celltype_refinery.core.composition \\
        --input consolidated.h5ad \\
        --out composition_output

    # Process specific columns (multi-column mode)
    python -m celltype_refinery.core.composition \\
        --input consolidated.h5ad \\
        --cell-type-cols cell_type_phenocycler cell_type_broad \\
        --out composition_output

    # Single-column mode (backward compatible)
    python -m celltype_refinery.core.composition \\
        --input consolidated.h5ad \\
        --cell-type-col cell_type_phenocycler \\
        --out composition_output

    # Quick run (skip enrichment tests)
    python -m celltype_refinery.core.composition \\
        --input consolidated.h5ad \\
        --skip-enrichment \\
        --out composition_quick

    # With interactive dashboard
    python -m celltype_refinery.core.composition \\
        --input consolidated.h5ad \\
        --dashboard \\
        --out composition_output
        """,
    )

    # Required arguments
    parser.add_argument(
        "--input", "-i",
        type=Path,
        required=True,
        help="Path to annotated AnnData (e.g., consolidated.h5ad)",
    )
    parser.add_argument(
        "--out", "-o",
        type=Path,
        required=True,
        help="Output directory",
    )

    # Configuration
    parser.add_argument(
        "--config", "-c",
        type=Path,
        help="Path to YAML configuration file",
    )

    # Cell type column options (mutually exclusive behavior)
    parser.add_argument(
        "--cell-type-col",
        type=str,
        default=None,
        help="Single cell type column (enables single-column mode for backward compat)",
    )
    parser.add_argument(
        "--cell-type-cols",
        type=str,
        nargs="+",
        default=None,
        help=(
            "Cell type columns to analyze (multi-column mode). "
            "Default: cell_type_phenocycler cell_type_phenocycler_detailed "
            "cell_type_multiomics cell_type_broad"
        ),
    )
    parser.add_argument(
        "--single-column-mode",
        action="store_true",
        help="Force single-column mode (uses --cell-type-col or default)",
    )

    # Organ-specific biology metrics
    parser.add_argument(
        "--organ",
        type=str,
        default=None,
        help=(
            "Organ type for organ-specific biology metrics. "
            "Available: 'fallopian_tube' (aliases: ft, fallopian, oviduct), 'generic'. "
            "Default: generic"
        ),
    )

    # Skip options
    parser.add_argument(
        "--skip-enrichment",
        action="store_true",
        help="Skip regional enrichment tests (faster execution)",
    )
    parser.add_argument(
        "--skip-biology",
        action="store_true",
        help="Skip organ biology metrics",
    )
    parser.add_argument(
        "--skip-viz",
        action="store_true",
        help="Skip static visualization generation",
    )

    # Interactive dashboard
    parser.add_argument(
        "--dashboard",
        action="store_true",
        help="Generate interactive HTML dashboard (requires plotly)",
    )
    parser.add_argument(
        "--html-dir",
        type=Path,
        default=None,
        help="Output HTML dashboards to this directory (flat structure). "
             "All HTML files will be self-contained and in one folder for easy sharing. "
             "If not specified, HTML files go into the composition output directory.",
    )

    # Visualization options
    parser.add_argument(
        "--dpi",
        type=int,
        default=200,
        help="DPI for figures (default: 200)",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=15,
        help="Number of top cell types to show in plots (default: 15)",
    )

    # Logging
    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
        help="Logging level (default: INFO)",
    )
    parser.add_argument(
        "--log-dir",
        type=Path,
        default=None,
        help="Directory for log files (default: output directory)",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose (DEBUG) logging",
    )

    return parser.parse_args()


def main() -> int:
    """Main entry point."""
    args = parse_args()

    # Determine log directory
    log_dir = args.log_dir if args.log_dir else args.out
    log_dir = Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)

    # Configure logging with file handler
    log_level = logging.DEBUG if args.verbose else getattr(logging, args.log_level)

    from celltype_refinery.io.logging import get_logger

    logger, actual_log_path = get_logger(
        name="composition_analysis",
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
    logger.info("=" * 60)
    logger.info("  Composition Analysis")
    logger.info("=" * 60)
    logger.info(f"Input: {args.input}")
    logger.info(f"Output: {args.out}")

    # Determine mode: single-column vs multi-column
    single_column_mode = args.single_column_mode or args.cell_type_col is not None
    if single_column_mode:
        cell_type_col = args.cell_type_col or "cell_type_phenocycler"
        cell_type_cols: Optional[List[str]] = None
        logger.info("Mode: single-column")
        logger.info(f"Cell type column: {cell_type_col}")
    else:
        cell_type_col = None
        cell_type_cols = args.cell_type_cols or DEFAULT_CELL_TYPE_COLUMNS
        logger.info("Mode: multi-column")
        logger.info(f"Cell type columns: {cell_type_cols}")

    # Validate input
    if not args.input.exists():
        logger.error(f"Input file not found: {args.input}")
        return 1

    # Import heavy dependencies
    try:
        import scanpy as sc
    except ImportError as e:
        logger.error(f"Missing dependency: {e}")
        logger.error("Install with: pip install scanpy")
        return 1

    from .config import CompositionConfig
    from .engine import CompositionEngine
    from .export import export_all, export_all_multi
    from .biology import get_organ_metrics, list_available_organs as list_biology_organs
    from celltype_refinery.config import get_organ_config, list_available_organs as list_organ_configs

    # Load config
    if args.config and args.config.exists():
        logger.info(f"Loading config from {args.config}")
        config = CompositionConfig.from_yaml(args.config)
    else:
        logger.info("Using default configuration")
        config = CompositionConfig.default()

    # Apply organ-specific region ordering and colors
    if args.organ:
        try:
            organ_config = get_organ_config(args.organ)
            config.apply_organ_config(organ_config)
            logger.info(f"Organ config: {organ_config.organ_name}")
            logger.info(f"  Region order: {organ_config.region_order}")
            logger.info(f"  Region colors: {len(organ_config.region_colors)} colors defined")
        except ValueError as e:
            logger.error(f"Invalid organ: {e}")
            logger.error(f"Available organs: {list_organ_configs()}")
            return 1

    # Apply skip flags to config (user preference: adapt __main__ to config)
    config.skip_enrichment = args.skip_enrichment
    config.skip_biology = args.skip_biology

    # Override visualization settings
    config.visualization.dpi = args.dpi
    config.visualization.top_n_cell_types = args.top_n

    logger.info(f"Config: skip_enrichment={config.skip_enrichment}, skip_biology={config.skip_biology}")

    # Get organ-specific biology metrics
    biology_metrics = None
    if args.organ:
        try:
            biology_metrics = get_organ_metrics(args.organ)
            logger.info(f"Biology metrics: {biology_metrics.tissue_name}")
            logger.info(f"  Available metrics: {biology_metrics.metric_names}")
        except ValueError as e:
            logger.error(f"Invalid organ for biology metrics: {e}")
            logger.error(f"Available organs: {list_biology_organs()}")
            return 1
    else:
        biology_metrics = get_organ_metrics(None)  # Generic
        logger.info("Organ: generic (use --organ for organ-specific metrics)")

    # Load data
    logger.info("")
    logger.info("Loading AnnData...")
    adata = sc.read_h5ad(args.input)
    logger.info(f"  Loaded {len(adata):,} cells, {adata.n_vars} features")

    # Create output directory
    args.out.mkdir(parents=True, exist_ok=True)

    # Execute analysis
    logger.info("")
    logger.info("--- Running Composition Analysis ---")

    engine = CompositionEngine(config=config, biology_metrics=biology_metrics)

    if single_column_mode:
        # Single-column mode (backward compatible)
        result = engine.execute(
            adata,
            cell_type_col=cell_type_col,
        )

        # Print summary
        logger.info("")
        logger.info("=" * 60)
        logger.info("  COMPOSITION SUMMARY")
        logger.info("=" * 60)
        logger.info(f"Cell types: {len(result.composition_global)}")
        logger.info(f"Samples: {result.provenance.get('n_samples', 'N/A')}")
        logger.info("")

        # Diversity summary
        if hasattr(result, 'diversity_summary') and isinstance(result.diversity_summary, dict):
            logger.info("Diversity metrics:")
            logger.info(f"  Mean Shannon entropy: {result.diversity_summary.get('mean_shannon', 0):.3f}")
            logger.info(f"  Mean Simpson index: {result.diversity_summary.get('mean_simpson', 0):.3f}")
            logger.info(f"  Mean evenness: {result.diversity_summary.get('mean_evenness', 0):.3f}")
        elif hasattr(result, 'diversity_summary') and hasattr(result.diversity_summary, 'empty'):
            if not result.diversity_summary.empty:
                logger.info("Diversity metrics: computed (see output files)")

        # Enrichment summary
        if result.enrichment is not None and len(result.enrichment) > 0:
            n_significant = result.enrichment["significant"].sum() if "significant" in result.enrichment.columns else 0
            logger.info(f"Significant regional enrichments: {n_significant}")

        logger.info(f"Execution time: {result.provenance.get('duration_seconds', 0):.2f}s")

        # Export outputs (flat structure)
        logger.info("")
        logger.info("--- Exporting Outputs ---")
        outputs = export_all(result, config, args.input, args.out)

        for name, path in outputs.items():
            logger.info(f"  {name}: {path.name}")

        # Generate visualizations
        if not args.skip_viz:
            logger.info("")
            logger.info("--- Generating Static Visualizations ---")

            try:
                from .viz import generate_all_figures

                fig_dir = args.out / "figures"
                fig_outputs = generate_all_figures(
                    adata=adata,
                    result=result,
                    output_dir=fig_dir,
                    cell_type_col=cell_type_col,
                    config=config,
                )

                for name, path in fig_outputs.items():
                    if path:
                        logger.info(f"  {name}: {path.name}")

            except ImportError as e:
                logger.warning(f"Visualization generation skipped (missing dependency): {e}")
            except Exception as e:
                logger.warning(f"Visualization generation failed: {e}")
                logger.warning("Continuing without visualizations...")

    else:
        # Multi-column mode (new default)
        multi_result = engine.execute_multi(
            adata,
            cell_type_columns=cell_type_cols,
        )

        # Print summary for each column
        logger.info("")
        logger.info("=" * 60)
        logger.info("  MULTI-COLUMN COMPOSITION SUMMARY")
        logger.info("=" * 60)
        logger.info(f"Columns processed: {len(multi_result.columns_analyzed)}")
        logger.info("")

        for col in multi_result.columns_analyzed:
            result = multi_result.results[col]
            logger.info(f"  {col}:")
            logger.info(f"    Cell types: {len(result.composition_global)}")

        logger.info("")

        # Export outputs (subdirectory structure)
        logger.info("--- Exporting Outputs ---")
        all_outputs = export_all_multi(multi_result, config, args.input, args.out)

        for col in multi_result.columns_analyzed:
            logger.info(f"  {col}/ (subdirectory with all outputs)")

        # Generate visualizations
        if not args.skip_viz:
            logger.info("")
            logger.info("--- Generating Static Visualizations ---")

            try:
                from .viz import generate_all_figures_multi

                fig_outputs = generate_all_figures_multi(
                    adata=adata,
                    multi_result=multi_result,
                    output_dir=args.out,
                    config=config,
                )

                for col in multi_result.columns_analyzed:
                    logger.info(f"  {col}/figures/ generated")

            except ImportError as e:
                logger.warning(f"Visualization generation skipped (missing dependency): {e}")
            except Exception as e:
                logger.warning(f"Visualization generation failed: {e}")
                logger.warning("Continuing without visualizations...")

    # Generate interactive dashboard
    if args.dashboard:
        logger.info("")
        logger.info("--- Generating Interactive Dashboard ---")

        try:
            from .viz import generate_dashboard

            # Determine output path based on html_dir
            if args.html_dir:
                logger.info(f"  HTML output directory: {args.html_dir}")
                output_path = None  # Will be determined by html_dir
            else:
                output_path = args.out / "dashboard.html"

            dashboard_path = generate_dashboard(
                output_dir=args.out,
                output_path=output_path,
                html_dir=args.html_dir,
                title="Cell-Type Composition Analysis",
            )

            if dashboard_path:
                logger.info(f"  Dashboard: {dashboard_path.name}")
                if args.html_dir:
                    logger.info(f"  All HTML files saved to: {args.html_dir}")
            else:
                logger.warning("Dashboard generation failed (plotly may not be installed)")

        except ImportError as e:
            logger.warning(f"Dashboard generation skipped: {e}")
            logger.warning("Install plotly with: pip install plotly")
        except Exception as e:
            logger.warning(f"Dashboard generation failed: {e}")

    logger.info("")
    logger.info("=" * 60)
    logger.info("  Composition Analysis Complete")
    logger.info("=" * 60)
    logger.info(f"Outputs saved to: {args.out}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
