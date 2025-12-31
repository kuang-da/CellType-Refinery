"""Clustering module CLI runner (Stage H).

Enables running Stage H clustering and annotation as:
    python -m celltype_refinery.core.clustering --input <h5ad> --marker-map <json> --output <dir>

Usage Examples:
    # Basic Stage H: Clustering + Annotation
    python -m celltype_refinery.core.clustering \\
        --input output/fallopian_tube/stage_f/merged_data.h5ad \\
        --marker-map data/fallopian_tube/FT_cell_type_markers_v9.json \\
        --output output/fallopian_tube/stage_h

    # With custom parameters
    python -m celltype_refinery.core.clustering \\
        --input output/stage_f/merged_data.h5ad \\
        --marker-map data/markers.json \\
        --output output/stage_h \\
        --layer batchcorr \\
        --resolution 0.6 \\
        --n-pcs 30

    # Annotation-only mode (skip clustering, reuse existing clusters)
    python -m celltype_refinery.core.clustering \\
        --input output/stage_h/coarse_clusters.h5ad \\
        --marker-map data/markers_v2.json \\
        --output output/stage_h_v2 \\
        --annotation-only

    # Clustering-only mode (skip annotation)
    python -m celltype_refinery.core.clustering \\
        --input output/stage_f/merged_data.h5ad \\
        --output output/stage_h \\
        --clustering-only

Note on GPU Non-Determinism:
    GPU-accelerated Leiden clustering (via RAPIDS/cuGraph) is non-deterministic
    even with a fixed seed. Cluster counts may vary between runs (typically 34-42
    for resolution=0.6 on FT data). Recommended workflow:
    1. Run Stage H once to establish baseline clusters
    2. Use --annotation-only for subsequent iterations (preserves cluster IDs)
    3. Focus on cell-type composition similarity rather than exact cluster counts
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Sequence

import numpy as np
import pandas as pd

# GPU acceleration
try:
    import rapids_singlecell as rsc
    import cupy as cp
    GPU_AVAILABLE = cp.cuda.is_available()
except ImportError:
    GPU_AVAILABLE = False
    rsc = None
    cp = None

# Annotation module
try:
    from ..annotation.engine import AnnotationEngine, AnnotationParams
    ANNOTATION_AVAILABLE = True
except ImportError:
    ANNOTATION_AVAILABLE = False
    AnnotationEngine = None
    AnnotationParams = None


@dataclass
class StageHRunConfig:
    """Configuration for Stage H clustering and annotation.

    Attributes
    ----------
    layer : str
        AnnData layer to use for clustering (default: batchcorr)
    n_pcs : int
        Number of principal components for neighbors
    neighbors_k : int
        k for neighborhood graph
    resolution : float
        Leiden clustering resolution
    scale_clip : float
        Value clipping during scaling
    seed : int
        Random seed for reproducibility
    min_marker_std : float
        Drop markers with std below this value
    use_gpu : bool
        Use GPU acceleration if available
    positive_quantile : float
        Quantile for marker positivity threshold
    de_method : str
        Differential expression method
    de_bonus : float
        DE bonus weight for scoring
    de_top_frac : float
        Fraction of panel for DE top-K
    de_min_k : int
        Minimum K for DE top-K
    de_max_k : int
        Maximum K for DE top-K
    anti_weight : float
        Weight for anti-marker penalty
    use_idf : bool
        Use IDF-based marker weighting
    compute_umap : bool
        Compute UMAP embeddings
    annotation_only : bool
        Skip clustering, use existing clusters
    clustering_only : bool
        Skip annotation, only cluster
    expand_markers : bool
        Output per-marker evidence table
    """

    layer: str = "batchcorr"
    n_pcs: int = 30
    neighbors_k: int = 15
    resolution: float = 0.6
    scale_clip: float = 10.0
    seed: int = 1337
    min_marker_std: float = 1e-3
    use_gpu: bool = True
    positive_quantile: float = 0.75
    de_method: str = "wilcoxon"
    de_bonus: float = 0.5
    de_top_frac: float = 0.2
    de_min_k: int = 3
    de_max_k: int = 12
    anti_weight: float = 0.5
    use_idf: bool = True
    compute_umap: bool = False
    annotation_only: bool = False
    clustering_only: bool = False
    expand_markers: bool = True


def setup_logging(
    verbose: bool = False,
    log_dir: Optional[Path] = None,
    log_filename: str = "stage_h.log",
) -> logging.Logger:
    """Setup logging configuration.

    Parameters
    ----------
    verbose : bool
        Enable verbose (DEBUG) logging
    log_dir : Optional[Path]
        Directory for log file (if None, console only)
    log_filename : str
        Log file name

    Returns
    -------
    logging.Logger
        Configured logger
    """
    level = logging.DEBUG if verbose else logging.INFO
    logger = logging.getLogger("stage_h")
    logger.setLevel(level)
    logger.handlers.clear()

    formatter = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # Console handler
    console = logging.StreamHandler()
    console.setLevel(level)
    console.setFormatter(formatter)
    logger.addHandler(console)

    # File handler (optional)
    if log_dir:
        log_dir = Path(log_dir)
        log_dir.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_dir / log_filename)
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        logger.info(f"Log file: {log_dir / log_filename}")

    return logger


def run_stage_h(
    input_path: Path,
    marker_map_path: Optional[Path],
    output_dir: Path,
    config: Optional[StageHRunConfig] = None,
    verbose: bool = False,
    skip_figures: bool = False,
    dpi: int = 200,
    log_dir: Optional[Path] = None,
) -> "anndata.AnnData":
    """Run Stage H: Coarse clustering and cell-type annotation.

    Performs:
    1. Data preprocessing (scaling, variance filtering)
    2. PCA dimensionality reduction
    3. Neighborhood graph construction
    4. Leiden clustering (GPU-accelerated if available)
    5. Differential expression testing
    6. Marker scoring and hierarchical annotation

    Parameters
    ----------
    input_path : Path
        Path to Stage F merged_data.h5ad
    marker_map_path : Optional[Path]
        Path to marker map JSON file (required unless clustering_only)
    output_dir : Path
        Output directory for Stage H artifacts
    config : Optional[StageHRunConfig]
        Clustering and annotation configuration
    verbose : bool
        Enable verbose logging
    skip_figures : bool
        Skip figure generation
    dpi : int
        Figure resolution
    log_dir : Optional[Path]
        Directory for log file

    Returns
    -------
    anndata.AnnData
        AnnData with cluster assignments and annotations
    """
    try:
        import scanpy as sc
    except ImportError:
        raise ImportError("scanpy is required for Stage H. Install with: pip install scanpy")

    # Setup logging
    logger = setup_logging(verbose, log_dir=log_dir, log_filename="stage_h.log")
    logger.info("Stage H: Coarse Clustering and Annotation")
    logger.info(f"Input: {input_path}")
    if marker_map_path:
        logger.info(f"Marker map: {marker_map_path}")
    logger.info(f"Output: {output_dir}")

    # Setup
    config = config or StageHRunConfig()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir = output_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    # Validate inputs
    if not config.clustering_only and marker_map_path is None:
        raise ValueError("--marker-map is required unless --clustering-only is specified")

    # Log configuration
    logger.info(f"Layer: {config.layer}")
    logger.info(f"n_pcs: {config.n_pcs}, neighbors_k: {config.neighbors_k}, resolution: {config.resolution}")
    logger.info(f"DE method: {config.de_method}, DE bonus: {config.de_bonus}")
    logger.info(f"Annotation only: {config.annotation_only}, Clustering only: {config.clustering_only}")

    # Check GPU availability
    use_gpu = config.use_gpu and GPU_AVAILABLE
    if config.use_gpu and not GPU_AVAILABLE:
        logger.warning("GPU requested but not available. Falling back to CPU.")
    logger.info(f"GPU acceleration: {'enabled' if use_gpu else 'disabled'}")

    # Load AnnData
    logger.info("Loading AnnData...")
    adata = sc.read_h5ad(input_path)
    logger.info(f"Loaded: {adata.n_obs:,} cells, {adata.n_vars} markers")

    # Select layer
    layer = config.layer
    if layer not in adata.layers and layer != "X":
        logger.warning(f"Layer '{layer}' not found, falling back to 'X'")
        layer = "X"

    # Prepare data matrix for clustering
    if layer in adata.layers:
        adata.X = adata.layers[layer].copy()
        logger.info(f"Using layer: {layer}")
    else:
        logger.info("Using X matrix")

    # Load marker map (if needed)
    marker_map = None
    if marker_map_path and not config.clustering_only:
        logger.info("Loading marker map...")
        with open(marker_map_path, "r") as f:
            marker_map = json.load(f)
        logger.info(f"Loaded marker map with {len(marker_map)} entries")

    start_time = time.time()

    if not config.annotation_only:
        # =====================================================================
        # CLUSTERING PHASE
        # =====================================================================
        logger.info("=" * 60)
        logger.info("CLUSTERING PHASE")
        logger.info("=" * 60)

        # Filter low variance markers
        logger.info("Filtering low variance markers...")
        var_mask = np.nanstd(adata.X, axis=0) > config.min_marker_std
        n_dropped = adata.n_vars - var_mask.sum()
        if n_dropped > 0:
            dropped_markers = adata.var_names[~var_mask].tolist()
            logger.info(f"Dropping {n_dropped} low-variance markers: {dropped_markers}")
            adata = adata[:, var_mask].copy()

        # Scale data
        logger.info("Scaling data...")
        sc.pp.scale(adata, max_value=config.scale_clip)

        # PCA
        logger.info(f"Running PCA (n_pcs={config.n_pcs})...")
        sc.tl.pca(adata, n_comps=config.n_pcs, random_state=config.seed)

        # Neighbors
        logger.info(f"Building neighbor graph (k={config.neighbors_k})...")
        if use_gpu and rsc is not None:
            try:
                rsc.pp.neighbors(adata, n_neighbors=config.neighbors_k, n_pcs=config.n_pcs)
            except Exception as e:
                logger.warning(f"GPU neighbors failed: {e}. Falling back to CPU.")
                sc.pp.neighbors(adata, n_neighbors=config.neighbors_k, n_pcs=config.n_pcs)
        else:
            sc.pp.neighbors(adata, n_neighbors=config.neighbors_k, n_pcs=config.n_pcs)

        # Leiden clustering
        logger.info(f"Running Leiden clustering (resolution={config.resolution})...")
        if use_gpu and rsc is not None:
            try:
                rsc.tl.leiden(adata, resolution=config.resolution, random_state=config.seed, key_added="cluster_lvl0")
            except Exception as e:
                logger.warning(f"GPU Leiden failed: {e}. Falling back to CPU.")
                sc.tl.leiden(adata, resolution=config.resolution, random_state=config.seed, key_added="cluster_lvl0")
        else:
            sc.tl.leiden(adata, resolution=config.resolution, random_state=config.seed, key_added="cluster_lvl0")

        n_clusters = adata.obs["cluster_lvl0"].nunique()
        logger.info(f"Found {n_clusters} clusters")

        # UMAP (optional)
        if config.compute_umap:
            logger.info("Computing UMAP...")
            if use_gpu and rsc is not None:
                try:
                    rsc.tl.umap(adata, random_state=config.seed)
                except Exception as e:
                    logger.warning(f"GPU UMAP failed: {e}. Falling back to CPU.")
                    sc.tl.umap(adata, random_state=config.seed)
            else:
                sc.tl.umap(adata, random_state=config.seed)

        # Differential Expression
        logger.info(f"Running differential expression ({config.de_method})...")
        de_layer = layer if layer in adata.layers else None
        sc.tl.rank_genes_groups(
            adata,
            groupby="cluster_lvl0",
            method=config.de_method,
            n_genes=config.de_max_k,
            layer=de_layer,
            use_raw=False,
            tie_correct=True,
            pts=True,
            key_added="de_wilcoxon",
        )
        logger.info("DE testing complete")

    else:
        # Annotation-only mode: verify cluster_lvl0 exists
        if "cluster_lvl0" not in adata.obs:
            raise ValueError(
                "--annotation-only requires 'cluster_lvl0' in adata.obs. "
                "Run Stage H without --annotation-only first."
            )
        logger.info("Annotation-only mode: using existing cluster_lvl0")
        n_clusters = adata.obs["cluster_lvl0"].nunique()
        logger.info(f"Found {n_clusters} existing clusters")

    clustering_time = time.time() - start_time
    logger.info(f"Clustering phase completed in {clustering_time:.1f}s")

    if config.clustering_only:
        # Clustering-only mode: skip annotation
        logger.info("Clustering-only mode: skipping annotation")

        # Save AnnData
        h5ad_path = output_dir / "coarse_clusters.h5ad"
        adata.write_h5ad(h5ad_path)
        logger.info(f"Wrote {h5ad_path}")

        logger.info(f"Stage H completed (clustering only): {n_clusters} clusters")
        return adata

    # =========================================================================
    # ANNOTATION PHASE
    # =========================================================================
    if not ANNOTATION_AVAILABLE:
        raise ImportError(
            "Annotation module not available. "
            "Ensure celltype_refinery.core.annotation is installed."
        )

    logger.info("=" * 60)
    logger.info("ANNOTATION PHASE")
    logger.info("=" * 60)

    annotation_start = time.time()

    # Setup annotation parameters
    annotation_params = AnnotationParams(
        cluster_key="cluster_lvl0",
        label_col="cell_type_auto",
        layer=layer if layer in adata.layers else "X",
        positive_quantile=config.positive_quantile,
        de_key="de_wilcoxon",
        de_bonus=config.de_bonus,
        de_top_frac=config.de_top_frac,
        de_min_k=config.de_min_k,
        de_max_k=config.de_max_k,
        anti_weight=config.anti_weight,
        use_idf=config.use_idf,
        expand_markers=config.expand_markers,
    )

    # Run annotation engine
    annotation_engine = AnnotationEngine(
        marker_map=marker_map_path,
        params=annotation_params,
        logger=logger,
    )

    result = annotation_engine.run(
        adata=adata,
        cluster_key="cluster_lvl0",
        output_dir=output_dir,
    )

    annotation_time = time.time() - annotation_start
    logger.info(f"Annotation phase completed in {annotation_time:.1f}s")

    # Save AnnData with annotations
    h5ad_path = output_dir / "coarse_clusters.h5ad"
    adata.write_h5ad(h5ad_path)
    logger.info(f"Wrote {h5ad_path}")

    # Generate figures
    if not skip_figures:
        logger.info("Generating figures...")
        _generate_stage_h_figures(
            adata=adata,
            cluster_annotations=result.cluster_annotations,
            output_dir=figures_dir,
            dpi=dpi,
            logger=logger,
        )

    total_time = time.time() - start_time
    logger.info(f"Stage H completed in {total_time:.1f}s")
    logger.info(f"  Clusters: {n_clusters}")
    logger.info(f"  Cell types: {result.cluster_annotations['assigned_label'].nunique()}")

    return adata


def _generate_stage_h_figures(
    adata: "anndata.AnnData",
    cluster_annotations: pd.DataFrame,
    output_dir: Path,
    dpi: int = 200,
    logger: Optional[logging.Logger] = None,
) -> List[str]:
    """Generate Stage H diagnostic figures.

    Parameters
    ----------
    adata : anndata.AnnData
        AnnData with clusters and annotations
    cluster_annotations : pd.DataFrame
        Cluster annotation table
    output_dir : Path
        Output directory for figures
    dpi : int
        Figure resolution
    logger : Optional[logging.Logger]
        Logger instance

    Returns
    -------
    List[str]
        List of generated figure paths
    """
    import matplotlib.pyplot as plt

    if logger is None:
        logger = logging.getLogger(__name__)

    generated = []
    output_dir = Path(output_dir)

    # Figure 1: Cluster size distribution
    try:
        fig, ax = plt.subplots(figsize=(10, 6))
        cluster_sizes = adata.obs["cluster_lvl0"].value_counts().sort_index()
        ax.bar(range(len(cluster_sizes)), cluster_sizes.values)
        ax.set_xlabel("Cluster ID")
        ax.set_ylabel("Number of Cells")
        ax.set_title("Stage H: Cluster Size Distribution")
        ax.set_yscale("log")
        fig_path = output_dir / "stage_h_cluster_sizes.png"
        fig.savefig(fig_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        generated.append(str(fig_path))
        logger.info(f"Generated: {fig_path.name}")
    except Exception as e:
        logger.warning(f"Failed to generate cluster size figure: {e}")

    # Figure 2: Annotation confidence distribution
    if "confidence" in cluster_annotations.columns:
        try:
            fig, ax = plt.subplots(figsize=(10, 6))
            confidences = cluster_annotations["confidence"].values
            confidences = confidences[np.isfinite(confidences) & (confidences < 1e5)]
            ax.hist(confidences, bins=50, edgecolor="black", alpha=0.7)
            ax.set_xlabel("Confidence (min margin along path)")
            ax.set_ylabel("Number of Clusters")
            ax.set_title("Stage H: Annotation Confidence Distribution")
            ax.axvline(0.3, color="red", linestyle="--", label="Low confidence threshold")
            ax.legend()
            fig_path = output_dir / "stage_h_confidence_distribution.png"
            fig.savefig(fig_path, dpi=dpi, bbox_inches="tight")
            plt.close(fig)
            generated.append(str(fig_path))
            logger.info(f"Generated: {fig_path.name}")
        except Exception as e:
            logger.warning(f"Failed to generate confidence figure: {e}")

    # Figure 3: Cell type composition
    if "cell_type_auto" in adata.obs.columns:
        try:
            fig, ax = plt.subplots(figsize=(12, 8))
            cell_type_counts = adata.obs["cell_type_auto"].value_counts()
            top_n = min(20, len(cell_type_counts))
            top_types = cell_type_counts.head(top_n)
            ax.barh(range(len(top_types)), top_types.values)
            ax.set_yticks(range(len(top_types)))
            ax.set_yticklabels(top_types.index)
            ax.set_xlabel("Number of Cells")
            ax.set_title(f"Stage H: Top {top_n} Cell Types")
            ax.invert_yaxis()
            fig_path = output_dir / "stage_h_cell_type_composition.png"
            fig.savefig(fig_path, dpi=dpi, bbox_inches="tight")
            plt.close(fig)
            generated.append(str(fig_path))
            logger.info(f"Generated: {fig_path.name}")
        except Exception as e:
            logger.warning(f"Failed to generate cell type composition figure: {e}")

    return generated


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Stage H: Coarse Clustering and Cell-Type Annotation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic clustering + annotation
  python -m celltype_refinery.core.clustering \\
      --input output/stage_f/merged_data.h5ad \\
      --marker-map data/markers.json \\
      --output output/stage_h

  # Annotation-only (reuse existing clusters)
  python -m celltype_refinery.core.clustering \\
      --input output/stage_h/coarse_clusters.h5ad \\
      --marker-map data/markers_v2.json \\
      --output output/stage_h_v2 \\
      --annotation-only

  # Clustering-only (no annotation)
  python -m celltype_refinery.core.clustering \\
      --input output/stage_f/merged_data.h5ad \\
      --output output/stage_h \\
      --clustering-only

Note: GPU Leiden clustering is non-deterministic. Use --annotation-only
for reproducible marker map iterations.
        """,
    )

    # Required arguments
    parser.add_argument(
        "--input", "-i",
        type=Path,
        required=True,
        help="Path to input AnnData (Stage F merged_data.h5ad)",
    )
    parser.add_argument(
        "--output", "-o",
        type=Path,
        required=True,
        help="Output directory for Stage H artifacts",
    )

    # Marker map (required unless clustering-only)
    parser.add_argument(
        "--marker-map", "-m",
        type=Path,
        default=None,
        help="Path to marker map JSON file (REQUIRED unless --clustering-only)",
    )

    # Clustering parameters
    parser.add_argument(
        "--layer",
        type=str,
        default="batchcorr",
        help="AnnData layer to use for clustering (default: batchcorr)",
    )
    parser.add_argument(
        "--n-pcs",
        type=int,
        default=30,
        help="Number of principal components (default: 30)",
    )
    parser.add_argument(
        "--neighbors-k",
        type=int,
        default=15,
        help="k for neighborhood graph (default: 15)",
    )
    parser.add_argument(
        "--resolution",
        type=float,
        default=0.6,
        help="Leiden clustering resolution (default: 0.6)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=1337,
        help="Random seed (default: 1337)",
    )

    # GPU options
    parser.add_argument(
        "--use-gpu",
        action="store_true",
        default=True,
        help="Use GPU acceleration if available (default: True)",
    )
    parser.add_argument(
        "--no-gpu",
        action="store_false",
        dest="use_gpu",
        help="Disable GPU acceleration",
    )

    # Mode options
    parser.add_argument(
        "--annotation-only",
        action="store_true",
        help="Skip clustering, use existing clusters from input h5ad",
    )
    parser.add_argument(
        "--clustering-only",
        action="store_true",
        help="Skip annotation, only perform clustering",
    )
    parser.add_argument(
        "--compute-umap",
        action="store_true",
        help="Compute UMAP embeddings for visualization",
    )

    # Scoring parameters
    parser.add_argument(
        "--de-bonus",
        type=float,
        default=0.5,
        help="DE bonus weight for marker scoring (default: 0.5)",
    )
    parser.add_argument(
        "--anti-weight",
        type=float,
        default=0.5,
        help="Weight for anti-marker penalty (default: 0.5)",
    )

    # Output options
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable verbose logging",
    )
    parser.add_argument(
        "--skip-figures",
        action="store_true",
        help="Skip figure generation",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=200,
        help="Figure resolution (default: 200)",
    )
    parser.add_argument(
        "--log-dir", "-l",
        type=Path,
        default=None,
        help="Directory for log file (default: console only)",
    )

    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None):
    """CLI entry point."""
    args = parse_args(argv)

    # Validate marker-map requirement
    if not args.clustering_only and args.marker_map is None:
        print("Error: --marker-map is required unless --clustering-only is specified", file=sys.stderr)
        sys.exit(1)

    # Build config
    config = StageHRunConfig(
        layer=args.layer,
        n_pcs=args.n_pcs,
        neighbors_k=args.neighbors_k,
        resolution=args.resolution,
        seed=args.seed,
        use_gpu=args.use_gpu,
        compute_umap=args.compute_umap,
        annotation_only=args.annotation_only,
        clustering_only=args.clustering_only,
        de_bonus=args.de_bonus,
        anti_weight=args.anti_weight,
    )

    try:
        run_stage_h(
            input_path=args.input,
            marker_map_path=args.marker_map,
            output_dir=args.output,
            config=config,
            verbose=args.verbose,
            skip_figures=args.skip_figures,
            dpi=args.dpi,
            log_dir=args.log_dir,
        )
    except Exception as e:
        logging.error(f"Stage H failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
