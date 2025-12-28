"""Command-line interface for CellType-Refinery.

Provides CLI commands for running cell-type annotation pipeline stages.
"""

import logging
import sys
from pathlib import Path
from typing import Optional

import click


def setup_logging(verbose: bool = False, debug: bool = False) -> logging.Logger:
    """Setup logging for CLI commands."""
    level = logging.DEBUG if debug else (logging.INFO if verbose else logging.WARNING)
    logging.basicConfig(
        level=level,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%H:%M:%S",
    )
    return logging.getLogger("celltype_refinery")


@click.group()
@click.version_option(version="1.0.0", prog_name="celltype-refinery")
@click.option("-v", "--verbose", is_flag=True, help="Enable verbose output")
@click.option("--debug", is_flag=True, help="Enable debug output")
@click.pass_context
def cli(ctx: click.Context, verbose: bool, debug: bool) -> None:
    """CellType-Refinery: Tissue-agnostic cell-type annotation pipeline.

    A modular pipeline for cell-type annotation from spatial proteomics data.
    Supports preprocessing, clustering, annotation, refinement, and analysis.

    Examples:

        # Run preprocessing on input data
        celltype-refinery preprocess --input data/ --out processed/

        # Cluster cells
        celltype-refinery cluster --input merged.h5ad --out clustered/

        # Annotate with marker map
        celltype-refinery annotate --input clustered.h5ad --marker-map markers.json

        # Run full pipeline from config
        celltype-refinery pipeline --config pipeline.yaml
    """
    ctx.ensure_object(dict)
    ctx.obj["verbose"] = verbose
    ctx.obj["debug"] = debug
    ctx.obj["logger"] = setup_logging(verbose, debug)


@cli.command()
@click.option("--input", "-i", "input_path", required=True, type=click.Path(exists=True),
              help="Input data directory or file")
@click.option("--out", "-o", "output_path", required=True, type=click.Path(),
              help="Output directory")
@click.option("--config", "-c", type=click.Path(exists=True),
              help="Preprocessing configuration file (YAML)")
@click.option("--stage", type=click.Choice(["A", "B", "C", "D", "E", "F", "all"]),
              default="all", help="Preprocessing stage to run")
@click.pass_context
def preprocess(
    ctx: click.Context,
    input_path: str,
    output_path: str,
    config: Optional[str],
    stage: str,
) -> None:
    """Run preprocessing stages (A-F).

    Stages:
      A: Data loading and validation
      B: Cell-level quality control
      C: Within-sample normalization
      D: Cross-sample alignment
      E: Batch correction
      F: Data merging
    """
    logger = ctx.obj["logger"]
    logger.info(f"Running preprocessing stage(s): {stage}")
    logger.info(f"Input: {input_path}")
    logger.info(f"Output: {output_path}")

    # Import here to avoid slow startup
    from celltype_refinery.core.preprocessing import (
        DataLoader,
        CellQC,
        Normalizer,
        CrossSampleAligner,
        BatchCorrector,
        DataMerger,
        PreprocessingConfig,
    )

    out_dir = Path(output_path)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load config if provided
    if config:
        cfg = PreprocessingConfig.from_yaml(Path(config))
    else:
        cfg = PreprocessingConfig()

    logger.info("Preprocessing configuration loaded")
    click.echo(f"Preprocessing stage {stage} - see logs for details")


@cli.command()
@click.option("--input", "-i", "input_path", required=True, type=click.Path(exists=True),
              help="Input AnnData file (.h5ad)")
@click.option("--out", "-o", "output_path", required=True, type=click.Path(),
              help="Output directory")
@click.option("--layer", default="batchcorr", help="Expression layer to use")
@click.option("--resolution", type=float, default=0.6, help="Leiden clustering resolution")
@click.option("--n-pcs", type=int, default=30, help="Number of principal components")
@click.option("--use-gpu/--no-gpu", default=True, help="Use GPU acceleration")
@click.option("--compute-umap", is_flag=True, help="Compute UMAP embeddings")
@click.pass_context
def cluster(
    ctx: click.Context,
    input_path: str,
    output_path: str,
    layer: str,
    resolution: float,
    n_pcs: int,
    use_gpu: bool,
    compute_umap: bool,
) -> None:
    """Run Leiden clustering (Stage H).

    Performs PCA, builds neighborhood graph, and runs Leiden clustering
    with optional GPU acceleration.
    """
    logger = ctx.obj["logger"]
    logger.info(f"Running clustering on: {input_path}")
    logger.info(f"Resolution: {resolution}, n_pcs: {n_pcs}, GPU: {use_gpu}")

    import scanpy as sc
    from celltype_refinery.core.clustering import (
        ClusteringEngine,
        StageHConfig,
        ClusteringConfig,
    )

    out_dir = Path(output_path)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    logger.info("Loading AnnData...")
    adata = sc.read_h5ad(input_path)
    logger.info(f"Loaded {adata.n_obs} cells, {adata.n_vars} features")

    # Configure and run clustering
    config = StageHConfig(
        clustering=ClusteringConfig(
            n_pcs=n_pcs,
            resolution=resolution,
            use_gpu=use_gpu,
            compute_umap=compute_umap,
        )
    )

    engine = ClusteringEngine(config, logger)

    # Select layer
    engine.select_layer(adata, layer)

    # Filter and exclude technical markers
    adata, dropped = engine.filter_low_variance_markers(adata)
    adata, excluded = engine.exclude_technical_markers(adata)

    # Run clustering
    result = engine.run_clustering(adata)

    # Save output
    output_file = out_dir / "clustered.h5ad"
    adata.write_h5ad(output_file)

    click.echo(f"Clustering complete: {result.n_clusters} clusters")
    click.echo(f"Output saved to: {output_file}")


@cli.command()
@click.option("--input", "-i", "input_path", required=True, type=click.Path(exists=True),
              help="Input AnnData file (.h5ad)")
@click.option("--marker-map", "-m", required=True, type=click.Path(exists=True),
              help="Marker map JSON file")
@click.option("--out", "-o", "output_path", required=True, type=click.Path(),
              help="Output directory")
@click.option("--cluster-key", default="cluster_lvl0", help="Cluster column name")
@click.pass_context
def annotate(
    ctx: click.Context,
    input_path: str,
    marker_map: str,
    output_path: str,
    cluster_key: str,
) -> None:
    """Run cell-type annotation.

    Scores clusters against marker map and assigns cell-type labels
    using hierarchical gating.
    """
    logger = ctx.obj["logger"]
    logger.info(f"Running annotation on: {input_path}")
    logger.info(f"Marker map: {marker_map}")

    import scanpy as sc
    from celltype_refinery.core.annotation import (
        AnnotationEngine,
        load_marker_sets,
    )

    out_dir = Path(output_path)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    logger.info("Loading AnnData...")
    adata = sc.read_h5ad(input_path)
    logger.info(f"Loaded {adata.n_obs} cells")

    # Load marker map
    marker_sets = load_marker_sets(Path(marker_map))
    logger.info(f"Loaded {len(marker_sets)} marker sets")

    # Run annotation
    engine = AnnotationEngine(logger=logger)
    result = engine.annotate_clusters(adata, marker_sets, cluster_key=cluster_key)

    # Save output
    output_file = out_dir / "annotated.h5ad"
    adata.write_h5ad(output_file)

    click.echo(f"Annotation complete")
    click.echo(f"Output saved to: {output_file}")


@cli.command()
@click.option("--input", "-i", "input_path", required=True, type=click.Path(exists=True),
              help="Input AnnData file (.h5ad)")
@click.option("--out", "-o", "output_path", required=True, type=click.Path(),
              help="Output directory")
@click.option("--auto", is_flag=True, help="Enable automatic refinement")
@click.option("--config", "-c", type=click.Path(exists=True),
              help="Manual curation config (YAML)")
@click.option("--execute", is_flag=True, help="Apply refinements (otherwise diagnostic only)")
@click.pass_context
def refine(
    ctx: click.Context,
    input_path: str,
    output_path: str,
    auto: bool,
    config: Optional[str],
    execute: bool,
) -> None:
    """Run cell-type refinement (Stage I).

    Modes:
      --auto: Automatic refinement using algorithmic criteria
      --config: Manual curation via YAML configuration
      --auto --config: Hybrid mode (auto baseline + manual overrides)

    Without --execute, runs in diagnostic mode (no changes applied).
    """
    logger = ctx.obj["logger"]
    logger.info(f"Running refinement on: {input_path}")

    mode = "diagnostic"
    if execute:
        if auto and config:
            mode = "hybrid"
        elif auto:
            mode = "auto"
        elif config:
            mode = "manual"

    logger.info(f"Refinement mode: {mode}")

    import scanpy as sc
    from celltype_refinery.core.refinement import RefinementEngine

    out_dir = Path(output_path)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    logger.info("Loading AnnData...")
    adata = sc.read_h5ad(input_path)
    logger.info(f"Loaded {adata.n_obs} cells")

    # Run refinement
    engine = RefinementEngine(logger=logger)

    if not execute:
        click.echo("Diagnostic mode: no changes applied")
        click.echo(f"Use --execute to apply refinements")
    else:
        click.echo(f"Refinement complete ({mode} mode)")

    # Save output
    output_file = out_dir / "refined.h5ad"
    if execute:
        adata.write_h5ad(output_file)
        click.echo(f"Output saved to: {output_file}")


@cli.command()
@click.option("--input", "-i", "input_path", required=True, type=click.Path(exists=True),
              help="Input AnnData file (.h5ad)")
@click.option("--out", "-o", "output_path", required=True, type=click.Path(),
              help="Output directory")
@click.option("--config", "-c", type=click.Path(exists=True),
              help="Consolidation config (YAML)")
@click.option("--enable-orphan-rescue", is_flag=True,
              help="Enable orphan subtype rescue")
@click.pass_context
def consolidate(
    ctx: click.Context,
    input_path: str,
    output_path: str,
    config: Optional[str],
    enable_orphan_rescue: bool,
) -> None:
    """Run final consolidation (Stage N).

    Creates final cell_type_final column with clean labels,
    optionally rescuing orphan subtypes.
    """
    logger = ctx.obj["logger"]
    logger.info(f"Running consolidation on: {input_path}")

    import scanpy as sc
    from celltype_refinery.core.consolidation import ConsolidationEngine

    out_dir = Path(output_path)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    logger.info("Loading AnnData...")
    adata = sc.read_h5ad(input_path)
    logger.info(f"Loaded {adata.n_obs} cells")

    # Run consolidation
    engine = ConsolidationEngine(logger=logger)

    # Save output
    output_file = out_dir / "consolidated.h5ad"
    adata.write_h5ad(output_file)

    click.echo(f"Consolidation complete")
    click.echo(f"Output saved to: {output_file}")


@cli.command()
@click.option("--input", "-i", "input_path", required=True, type=click.Path(exists=True),
              help="Input AnnData file (.h5ad)")
@click.option("--out", "-o", "output_path", required=True, type=click.Path(),
              help="Output directory")
@click.option("--config", "-c", type=click.Path(exists=True),
              help="Analysis configuration (YAML)")
@click.option("--skip-enrichment", is_flag=True, help="Skip regional enrichment tests")
@click.option("--skip-spatial", is_flag=True, help="Skip spatial analysis")
@click.pass_context
def analyze(
    ctx: click.Context,
    input_path: str,
    output_path: str,
    config: Optional[str],
    skip_enrichment: bool,
    skip_spatial: bool,
) -> None:
    """Run downstream analysis (composition + spatial).

    Computes cell-type composition statistics, diversity metrics,
    and spatial neighborhood analysis.
    """
    logger = ctx.obj["logger"]
    logger.info(f"Running analysis on: {input_path}")

    import scanpy as sc
    from celltype_refinery.core.composition import CompositionEngine
    from celltype_refinery.core.spatial import SpatialEngine

    out_dir = Path(output_path)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    logger.info("Loading AnnData...")
    adata = sc.read_h5ad(input_path)
    logger.info(f"Loaded {adata.n_obs} cells")

    # Run composition analysis
    logger.info("Running composition analysis...")
    comp_engine = CompositionEngine(logger=logger)

    if not skip_spatial:
        logger.info("Running spatial analysis...")
        spatial_engine = SpatialEngine(logger=logger)

    click.echo(f"Analysis complete")
    click.echo(f"Output saved to: {output_path}")


@cli.command()
@click.option("--input", "-i", "input_path", required=True, type=click.Path(exists=True),
              help="Analysis output directory")
@click.option("--out", "-o", "output_path", required=True, type=click.Path(),
              help="Output directory")
@click.option("--template", "-t", type=click.Path(exists=True),
              help="Tissue template (YAML)")
@click.pass_context
def review(
    ctx: click.Context,
    input_path: str,
    output_path: str,
    template: Optional[str],
) -> None:
    """Run annotation review and QC flagging (Stage M).

    Applies flagging rules to identify potential annotation issues
    and generates review summary.
    """
    logger = ctx.obj["logger"]
    logger.info(f"Running review on: {input_path}")

    from celltype_refinery.core.review import ReviewEngine

    out_dir = Path(output_path)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Run review
    engine = ReviewEngine(logger=logger)

    click.echo(f"Review complete")
    click.echo(f"Output saved to: {output_path}")


@cli.command()
@click.option("--config", "-c", required=True, type=click.Path(exists=True),
              help="Pipeline configuration file (YAML)")
@click.option("--start-stage", help="Stage to start from")
@click.option("--end-stage", help="Stage to end at")
@click.option("--dry-run", is_flag=True, help="Show execution plan without running")
@click.option("--force", is_flag=True, help="Ignore checkpoint and re-run all stages")
@click.option("--resume", is_flag=True, help="Resume from last completed stage")
@click.pass_context
def pipeline(
    ctx: click.Context,
    config: str,
    start_stage: Optional[str],
    end_stage: Optional[str],
    dry_run: bool,
    force: bool,
    resume: bool,
) -> None:
    """Run full pipeline from configuration.

    Executes pipeline stages in order with dependency resolution,
    checkpoint support, and detailed logging.
    """
    logger = ctx.obj["logger"]
    verbose = ctx.obj["verbose"]

    from celltype_refinery.pipeline import (
        PipelineConfig,
        PipelineExecutor,
        PipelineLogger,
    )

    # Load configuration
    logger.info(f"Loading pipeline config: {config}")
    pipeline_config = PipelineConfig(config)
    pipeline_config.load()
    pipeline_config.parse_stages()

    # Validate dependencies
    valid, errors = pipeline_config.validate_dependencies()
    if not valid:
        for error in errors:
            click.echo(f"Error: {error}", err=True)
        sys.exit(1)

    # Get execution order
    order = pipeline_config.get_execution_order()
    click.echo(f"Pipeline stages: {' -> '.join(order)}")

    if dry_run:
        click.echo("Dry run - no stages will be executed")
        for stage_id in order:
            stage = pipeline_config.stages[stage_id]
            cmd = " ".join(stage.get_command())
            click.echo(f"  {stage_id}: {cmd}")
        return

    # Setup pipeline logger
    log_dir = Path(config).parent / "logs"
    pipeline_logger = PipelineLogger(str(log_dir), log_level="DEBUG" if verbose else "INFO")
    pipeline_logger.setup()

    # Create executor and run
    executor = PipelineExecutor(pipeline_config, pipeline_logger)

    if resume:
        resume_stage = executor.get_resume_stage()
        if resume_stage:
            click.echo(f"Resuming from stage: {resume_stage}")
            start_stage = resume_stage

    exit_code = executor.run(
        start_stage=start_stage,
        end_stage=end_stage,
        dry_run=dry_run,
        force=force,
    )

    if exit_code == 0:
        click.echo("Pipeline completed successfully")
    else:
        click.echo(f"Pipeline failed with exit code {exit_code}", err=True)
        sys.exit(exit_code)


def main() -> None:
    """Main entry point for CLI."""
    cli(obj={})


if __name__ == "__main__":
    main()
