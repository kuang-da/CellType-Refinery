"""CLI for annotation module exports.

This is a thin CLI wrapper that delegates to export.py for all business logic.

Usage:
    python -m celltype_refinery.core.annotation export-review \
        --input refined.h5ad \
        --output-dir /path/to/output

    python -m celltype_refinery.core.annotation export-all \
        --input refined.h5ad \
        --output-dir /path/to/output
"""

from __future__ import annotations

import argparse
import logging
import sys
import warnings
from pathlib import Path

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')


def setup_logging(verbose: bool = False) -> logging.Logger:
    """Setup logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(message)s',
    )
    return logging.getLogger(__name__)


def cmd_export_review(args: argparse.Namespace) -> int:
    """Export review artifacts from annotated AnnData."""
    logger = setup_logging(args.verbose)

    try:
        import scanpy as sc
    except ImportError:
        logger.error("scanpy is required. Install with: pip install scanpy")
        return 1

    from .export import (
        detect_cell_type_column,
        detect_cluster_column,
        run_review_exports,
    )

    input_path = Path(args.input)
    output_dir = Path(args.output_dir)

    logger.info(f'Loading AnnData from {input_path}...')
    adata = sc.read_h5ad(input_path)
    logger.info(f'  Loaded {adata.n_obs:,} cells, {adata.n_vars} features')

    # Detect columns
    cell_type_col = args.cell_type_col or detect_cell_type_column(adata)
    cluster_col = args.cluster_col or detect_cluster_column(adata)

    # Run all review exports
    run_review_exports(
        adata=adata,
        output_dir=output_dir,
        cell_type_col=cell_type_col,
        cluster_col=cluster_col,
        stage_h_dir=args.stage_h_dir,
        marker_map_path=args.marker_map,
        input_path=str(input_path),
        logger=logger,
    )

    return 0


def cmd_export_all(args: argparse.Namespace) -> int:
    """Export all annotation artifacts including enhanced format."""
    logger = setup_logging(args.verbose)

    try:
        import scanpy as sc
    except ImportError:
        logger.error("scanpy is required. Install with: pip install scanpy")
        return 1

    import pandas as pd

    from .export import (
        detect_cell_type_column,
        detect_cluster_column,
        run_annotation_exports,
        export_composition_stats,
        export_review_summary,
        export_workflow_state,
    )

    input_path = Path(args.input)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f'Loading AnnData from {input_path}...')
    adata = sc.read_h5ad(input_path)
    logger.info(f'  Loaded {adata.n_obs:,} cells, {adata.n_vars} features')

    # Detect columns
    cell_type_col = args.cell_type_col or detect_cell_type_column(adata)
    cluster_col = args.cluster_col or detect_cluster_column(adata)
    logger.info(f'  Using cell type column: {cell_type_col}')
    logger.info(f'  Using cluster column: {cluster_col}')

    # Get marker scores from adata.uns
    marker_scores = None
    for key in ['marker_scores_refined', 'marker_scores', 'scores']:
        if key in adata.uns:
            scores_data = adata.uns[key]
            if isinstance(scores_data, pd.DataFrame):
                marker_scores = scores_data
            elif isinstance(scores_data, (dict, list)):
                marker_scores = pd.DataFrame(scores_data)
            break

    # Run all annotation exports
    output_paths = run_annotation_exports(
        adata=adata,
        output_dir=output_dir,
        marker_scores=marker_scores,
        decision_steps=None,
        cluster_annotations=None,
        iteration=args.iteration,
        cluster_col=cluster_col,
        label_col=cell_type_col,
        logger=logger,
    )

    # Also export composition and review summary
    export_composition_stats(adata, output_dir, cell_type_col, logger)

    comp_global = adata.obs[cell_type_col].value_counts()
    export_review_summary(adata, output_dir, cell_type_col, cluster_col, comp_global, logger)

    export_workflow_state(
        output_dir,
        marker_map_path=args.marker_map,
        input_path=str(input_path),
        stage_h_dir=args.stage_h_dir,
        logger=logger,
    )

    logger.info('')
    logger.info('Full export complete!')
    logger.info(f'  Total cells: {adata.n_obs:,}')
    logger.info(f'  Cell types: {len(comp_global)}')
    logger.info(f'  Export files: {len(output_paths)}')

    return 0


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description='Annotation module CLI for exports',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest='command', help='Available commands')

    # export-review command
    review_parser = subparsers.add_parser(
        'export-review',
        help='Export review artifacts (composition, cluster mapping, summary)',
    )
    review_parser.add_argument(
        '--input', '-i', required=True,
        help='Input h5ad file path',
    )
    review_parser.add_argument(
        '--output-dir', '-o', required=True,
        help='Output directory for exports',
    )
    review_parser.add_argument(
        '--cell-type-col', '-c',
        help='Cell type column name (auto-detected if not specified)',
    )
    review_parser.add_argument(
        '--cluster-col',
        help='Cluster column name (auto-detected if not specified)',
    )
    review_parser.add_argument(
        '--marker-map', '-m',
        help='Path to marker map (for workflow state)',
    )
    review_parser.add_argument(
        '--stage-h-dir',
        help='Path to Stage H directory (for workflow state)',
    )
    review_parser.add_argument(
        '--verbose', '-v', action='store_true',
        help='Enable verbose logging',
    )
    review_parser.set_defaults(func=cmd_export_review)

    # export-all command
    all_parser = subparsers.add_parser(
        'export-all',
        help='Export all annotation artifacts including enhanced format',
    )
    all_parser.add_argument(
        '--input', '-i', required=True,
        help='Input h5ad file path',
    )
    all_parser.add_argument(
        '--output-dir', '-o', required=True,
        help='Output directory for exports',
    )
    all_parser.add_argument(
        '--cell-type-col', '-c',
        help='Cell type column name (auto-detected if not specified)',
    )
    all_parser.add_argument(
        '--cluster-col',
        help='Cluster column name (auto-detected if not specified)',
    )
    all_parser.add_argument(
        '--marker-map', '-m',
        help='Path to marker map (for workflow state)',
    )
    all_parser.add_argument(
        '--stage-h-dir',
        help='Path to Stage H directory (for workflow state)',
    )
    all_parser.add_argument(
        '--iteration', type=int, default=1,
        help='Iteration level for Stage H format (default: 1)',
    )
    all_parser.add_argument(
        '--verbose', '-v', action='store_true',
        help='Enable verbose logging',
    )
    all_parser.set_defaults(func=cmd_export_all)

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        return 1

    return args.func(args)


if __name__ == '__main__':
    sys.exit(main())
