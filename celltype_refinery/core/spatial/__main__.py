"""CLI entry point for spatial analysis.

Usage:
    python -m celltype_refinery.core.spatial --input <adata.h5ad> --out <output_dir>
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Iterable, List, Optional

import scanpy as sc

from .config import SpatialConfig
from .engine import SpatialEngine
from .export import export_all

# Default cell type columns for multi-column analysis
DEFAULT_CELL_TYPE_COLUMNS: List[str] = [
    "cell_type_phenocycler",
    "cell_type_phenocycler_detailed",
    "cell_type_multiomics",
    "cell_type_broad",
]

LOG_FILENAME = "stage_l_spatial.log"


def setup_logging(log_dir: Path, verbose: bool = False) -> Path:
    """Configure logging with timestamped log files.

    Args:
        log_dir: Directory for log files.
        verbose: If True, use DEBUG level.

    Returns:
        Actual log file path (with timestamp).
    """
    log_dir = Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)

    level = logging.DEBUG if verbose else logging.INFO
    file_handler = None

    # Try to use CellType-Refinery's get_logger
    try:
        from celltype_refinery.io.logging import get_logger

        logger, actual_log_path = get_logger(
            name="stage_l_spatial",
            log_path=log_dir / LOG_FILENAME,
            level=level,
        )
    except ImportError:
        # Fallback to basic logging setup
        actual_log_path = log_dir / LOG_FILENAME
        logging.basicConfig(
            level=level,
            format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
        logger = logging.getLogger("stage_l_spatial")
        file_handler = logging.FileHandler(actual_log_path)
        file_handler.setLevel(level)
        file_handler.setFormatter(
            logging.Formatter(
                fmt="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
                datefmt="%Y-%m-%d %H:%M:%S",
            )
        )
        logger.addHandler(file_handler)

    # Add console handler for dual logging (file + stdout)
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(level)
    console_formatter = logging.Formatter(
        fmt="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)

    # Configure spatial module logger so engine.py etc. output to console/file
    spatial_logger = logging.getLogger("celltype_refinery.core.spatial")
    spatial_logger.setLevel(level)
    spatial_logger.addHandler(console_handler)
    if file_handler:
        spatial_logger.addHandler(file_handler)

    # Reduce noise from other libraries
    logging.getLogger("scanpy").setLevel(logging.WARNING)
    logging.getLogger("anndata").setLevel(logging.WARNING)

    return actual_log_path


def parse_args(argv: Optional[Iterable[str]] = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Stage L: Spatial Analysis",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Input/output
    parser.add_argument(
        "-i",
        "--input",
        dest="input_path",
        type=str,
        required=True,
        help="Path to annotated AnnData file (Stage N output)",
    )
    parser.add_argument(
        "-o",
        "--out",
        dest="output_dir",
        type=str,
        required=True,
        help="Output directory for results",
    )
    parser.add_argument(
        "-g",
        "--graphs-dir",
        dest="graphs_dir",
        type=str,
        required=True,
        help="Directory containing spatial graphs (NPZ files)",
    )
    parser.add_argument(
        "--config",
        dest="config_path",
        type=str,
        default=None,
        help="Path to YAML configuration file",
    )
    parser.add_argument(
        "--log-dir",
        dest="log_dir",
        type=str,
        default=None,
        help="Directory for log files (default: output_dir)",
    )

    # Cell type columns
    parser.add_argument(
        "--cell-type-col",
        dest="cell_type_col",
        type=str,
        default=None,
        help="Single cell type column to analyze (mutually exclusive with --multi)",
    )
    parser.add_argument(
        "--multi",
        dest="multi_column",
        action="store_true",
        default=False,
        help="Run multi-column analysis (all 4 standard columns)",
    )
    parser.add_argument(
        "--cell-type-cols",
        dest="cell_type_cols",
        type=str,
        nargs="+",
        default=None,
        help="Specific columns for multi-column analysis",
    )

    # Permutation settings
    parser.add_argument(
        "-n",
        "--n-permutations",
        dest="n_permutations",
        type=int,
        default=100,
        help="Number of permutations for null model",
    )
    parser.add_argument(
        "-t",
        "--n-threads",
        dest="n_threads",
        type=int,
        default=64,
        help="Threads for parallel permutations",
    )
    parser.add_argument(
        "--seed",
        dest="seed",
        type=int,
        default=42,
        help="Random seed for reproducibility",
    )

    # FDR correction
    parser.add_argument(
        "--correction",
        dest="correction",
        type=str,
        choices=["fdr_bh", "bonferroni", "holm", "none"],
        default="fdr_bh",
        help="Multiple testing correction method",
    )
    parser.add_argument(
        "--alpha",
        dest="alpha",
        type=float,
        default=0.05,
        help="Significance threshold",
    )

    # Skip options
    parser.add_argument(
        "--skip-enrichment",
        dest="skip_enrichment",
        action="store_true",
        default=False,
        help="Skip neighborhood enrichment analysis",
    )
    parser.add_argument(
        "--skip-interaction",
        dest="skip_interaction",
        action="store_true",
        default=False,
        help="Skip interaction score analysis",
    )
    parser.add_argument(
        "--skip-morans",
        dest="skip_morans",
        action="store_true",
        default=False,
        help="Skip Moran's I analysis",
    )
    parser.add_argument(
        "--skip-diversity",
        dest="skip_diversity",
        action="store_true",
        default=False,
        help="Skip local diversity analysis",
    )

    # Visualization
    parser.add_argument(
        "--skip-viz",
        dest="skip_viz",
        action="store_true",
        default=False,
        help="Skip visualization generation",
    )
    parser.add_argument(
        "--dpi",
        dest="dpi",
        type=int,
        default=200,
        help="Figure resolution (DPI)",
    )

    # Other
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        default=False,
        help="Verbose logging",
    )
    parser.add_argument(
        "--organ",
        type=str,
        default=None,
        help=(
            "Organ type for region ordering in visualizations. "
            "Available: 'fallopian_tube' (aliases: ft), 'uterus', etc. "
            "Default: alphabetical order"
        ),
    )

    return parser.parse_args(argv)


def main(argv: Optional[Iterable[str]] = None) -> int:
    """Main entry point."""
    args = parse_args(argv)

    # Setup logging with timestamped log file
    log_dir = Path(args.log_dir) if args.log_dir else Path(args.output_dir)
    actual_log_path = setup_logging(log_dir, verbose=args.verbose)
    logger = logging.getLogger("stage_l_spatial")

    print(f"[INFO] Log file: {actual_log_path}")
    logger.info("=" * 60)
    logger.info("Stage L: Spatial Analysis")
    logger.info("=" * 60)

    # Load or create config
    if args.config_path and Path(args.config_path).exists():
        config = SpatialConfig.from_yaml(Path(args.config_path))
        logger.info(f"Loaded config from {args.config_path}")
    else:
        config = SpatialConfig.default()

    # Override config with CLI args
    config.permutation.n_permutations = args.n_permutations
    config.permutation.n_threads = args.n_threads
    config.permutation.seed = args.seed
    config.correction.method = args.correction
    config.correction.alpha = args.alpha
    config.visualization.dpi = args.dpi

    # Determine cell type columns
    if args.cell_type_cols:
        cell_type_cols = args.cell_type_cols
        multi_column = True
    elif args.multi_column:
        cell_type_cols = DEFAULT_CELL_TYPE_COLUMNS
        multi_column = True
    elif args.cell_type_col:
        cell_type_cols = [args.cell_type_col]
        multi_column = False
    else:
        # Default: single column
        cell_type_cols = [config.cell_type_col]
        multi_column = False

    logger.info(f"Input: {args.input_path}")
    logger.info(f"Graphs: {args.graphs_dir}")
    logger.info(f"Output: {args.output_dir}")
    logger.info(f"Cell type columns: {cell_type_cols}")
    logger.info(f"Multi-column mode: {multi_column}")
    logger.info(f"Permutations: {config.permutation.n_permutations}")
    logger.info(f"Threads: {config.permutation.n_threads}")
    logger.info(f"FDR correction: {config.correction.method}")

    # Load data
    input_path = Path(args.input_path)
    if not input_path.exists():
        logger.error(f"Input file not found: {input_path}")
        return 1

    logger.info(f"Loading data from {input_path}...")
    adata = sc.read_h5ad(input_path)
    logger.info(f"Loaded {adata.n_obs:,} cells with {adata.n_vars} features")

    # Create engine and execute
    engine = SpatialEngine(config)

    graphs_dir = Path(args.graphs_dir)
    output_dir = Path(args.output_dir)

    if multi_column:
        result = engine.execute_multi(
            adata=adata,
            graphs_dir=graphs_dir,
            cell_type_cols=cell_type_cols,
            skip_enrichment=args.skip_enrichment,
            skip_interaction=args.skip_interaction,
            skip_morans=args.skip_morans,
            skip_diversity=args.skip_diversity,
        )
    else:
        result = engine.execute(
            adata=adata,
            graphs_dir=graphs_dir,
            cell_type_col=cell_type_cols[0],
            skip_enrichment=args.skip_enrichment,
            skip_interaction=args.skip_interaction,
            skip_morans=args.skip_morans,
            skip_diversity=args.skip_diversity,
        )

    # Export results
    logger.info(f"\nExporting results to {output_dir}...")
    export_all(
        result=result,
        output_dir=output_dir,
        config=config,
        input_path=input_path,
        graphs_dir=graphs_dir,
        cell_type_cols=cell_type_cols,
    )

    # Generate visualizations
    if not args.skip_viz:
        logger.info("Generating visualizations...")
        try:
            from .viz import generate_all_visualizations

            if multi_column:
                for col, col_result in result.results.items():
                    col_dir = output_dir / col.replace("/", "_").replace(" ", "_")
                    generate_all_visualizations(col_result, col_dir, config, organ=args.organ)
            else:
                generate_all_visualizations(result, output_dir, config, organ=args.organ)
        except ImportError as e:
            logger.warning(f"Visualization module not available, skipping plots: {e}")
        except Exception as e:
            logger.warning(f"Failed to generate visualizations: {e}")

    # Summary
    logger.info("\n" + "=" * 60)
    logger.info("Spatial Analysis Complete!")
    logger.info(f"Output: {output_dir}")

    if multi_column:
        logger.info(f"Columns processed: {result.columns_processed}")
        logger.info(
            f"Total time: {result.total_execution_time_seconds:.1f}s"
        )
    else:
        logger.info(f"Samples: {result.n_samples}")
        logger.info(f"Cells: {result.n_cells_total:,}")
        logger.info(f"Time: {result.execution_time_seconds:.1f}s")

    logger.info("=" * 60)

    return 0 if result.success else 1


if __name__ == "__main__":
    sys.exit(main())
