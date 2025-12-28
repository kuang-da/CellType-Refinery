"""Clustering engine for cell population identification (Stage H).

Provides Leiden clustering with optional GPU acceleration, and
targeted subclustering functionality.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple
import logging
import warnings

import numpy as np
from scipy import sparse

from .config import ClusteringConfig, StageHConfig


# GPU acceleration support (optional)
try:
    import rapids_singlecell as rsc
    import cupy as cp
    GPU_AVAILABLE = cp.cuda.is_available()
except ImportError:
    GPU_AVAILABLE = False
    rsc = None
    cp = None


@dataclass
class ClusteringResult:
    """Result from clustering operation.

    Attributes
    ----------
    n_clusters : int
        Number of clusters found
    cluster_key : str
        Key in adata.obs containing cluster assignments
    cluster_sizes : Dict[str, int]
        Map of cluster ID to cell count
    dropped_markers : List[str]
        Markers dropped due to low variance
    excluded_markers : List[str]
        Technical markers moved to metadata
    """

    n_clusters: int = 0
    cluster_key: str = "cluster_lvl0"
    cluster_sizes: Dict[str, int] = field(default_factory=dict)
    dropped_markers: List[str] = field(default_factory=list)
    excluded_markers: List[str] = field(default_factory=list)


class ClusteringEngine:
    """Clustering engine with Leiden algorithm and GPU support.

    Provides the core clustering pipeline: scale → PCA → neighbors → Leiden,
    with optional GPU acceleration via rapids-singlecell.

    Parameters
    ----------
    config : StageHConfig, optional
        Stage H configuration. If None, uses defaults.
    logger : logging.Logger, optional
        Logger instance. If None, creates default logger.

    Example
    -------
    >>> from celltype_refinery.core.clustering import ClusteringEngine, StageHConfig
    >>> config = StageHConfig()
    >>> engine = ClusteringEngine(config)
    >>> result = engine.run_clustering(adata, cluster_key="leiden")
    """

    def __init__(
        self,
        config: Optional[StageHConfig] = None,
        logger: Optional[logging.Logger] = None,
    ):
        self.config = config or StageHConfig()
        self.logger = logger or logging.getLogger(__name__)
        self._check_dependencies()

    def _check_dependencies(self) -> None:
        """Check for required dependencies."""
        try:
            import scanpy
        except ImportError:
            raise RuntimeError(
                "Clustering requires scanpy. Install with: pip install scanpy"
            )

    @property
    def gpu_available(self) -> bool:
        """Check if GPU acceleration is available."""
        return GPU_AVAILABLE

    def subset_adata(
        self,
        adata: Any,  # AnnData
        sample_ids: Optional[Sequence[str]] = None,
        donors: Optional[Sequence[str]] = None,
        regions: Optional[Sequence[str]] = None,
        focus_clusters: Optional[Sequence[str]] = None,
        focus_column: str = "cluster_lvl0",
    ) -> Any:
        """Subset AnnData by sample, donor, region, or cluster filters.

        Parameters
        ----------
        adata : AnnData
            Input AnnData object
        sample_ids : Sequence[str], optional
            Sample IDs to keep
        donors : Sequence[str], optional
            Donor IDs to keep
        regions : Sequence[str], optional
            Region names to keep
        focus_clusters : Sequence[str], optional
            Cluster IDs to keep
        focus_column : str
            Column to use for cluster filtering

        Returns
        -------
        AnnData
            Subset AnnData (copy if filtered, original if no filters)

        Raises
        ------
        ValueError
            If filters remove all cells
        """
        mask = np.ones(adata.n_obs, dtype=bool)

        if sample_ids and "sample_id" in adata.obs:
            self.logger.info("Applying sample_id filter (%d values)", len(sample_ids))
            mask &= adata.obs["sample_id"].astype(str).isin(set(sample_ids)).to_numpy()

        if donors and "donor" in adata.obs:
            self.logger.info("Applying donor filter (%d values)", len(donors))
            mask &= adata.obs["donor"].astype(str).isin(set(donors)).to_numpy()

        if regions and "region" in adata.obs:
            self.logger.info("Applying region filter (%d values)", len(regions))
            mask &= adata.obs["region"].astype(str).isin(set(regions)).to_numpy()

        if focus_clusters:
            if focus_column not in adata.obs:
                self.logger.warning(
                    "Focus column %s not present in obs; ignoring focus_clusters",
                    focus_column,
                )
            else:
                self.logger.info(
                    "Applying focus cluster filter (%d values) on column %s",
                    len(focus_clusters),
                    focus_column,
                )
                mask &= (
                    adata.obs[focus_column]
                    .astype(str)
                    .isin(set(focus_clusters))
                    .to_numpy()
                )

        if mask.sum() == 0:
            raise ValueError("Subset filters removed all cells; relax filter arguments.")

        if mask.all():
            return adata

        self.logger.info(
            "Subsetting AnnData: %d -> %d cells after filters",
            adata.n_obs,
            int(mask.sum()),
        )
        return adata[mask].copy()

    def select_layer(
        self,
        adata: Any,  # AnnData
        layer: Optional[str] = None,
    ) -> str:
        """Select and pin expression layer for clustering.

        Parameters
        ----------
        adata : AnnData
            Input AnnData object
        layer : str, optional
            Layer name to use. Falls back to X if not found.

        Returns
        -------
        str
            Name of the layer being used
        """
        if layer and layer in adata.layers:
            base = adata.layers[layer]
            self.logger.info("Using layer '%s' for clustering", layer)
        else:
            if layer:
                self.logger.warning(
                    "Requested layer '%s' not found; falling back to AnnData.X", layer
                )
            base = adata.X
            layer = "X"

        matrix = base.toarray() if sparse.issparse(base) else np.asarray(base)
        adata.layers["stage_h_input"] = matrix.astype(np.float32, copy=True)
        adata.X = matrix.astype(np.float32, copy=True)
        self.logger.info("Pinned working matrix into AnnData.X (shape=%s)", matrix.shape)
        return layer

    def filter_low_variance_markers(
        self,
        adata: Any,  # AnnData
        min_std: Optional[float] = None,
    ) -> Tuple[Any, List[str]]:
        """Filter out markers with low variance.

        Parameters
        ----------
        adata : AnnData
            Input AnnData object
        min_std : float, optional
            Minimum standard deviation threshold. Uses config default if None.

        Returns
        -------
        Tuple[AnnData, List[str]]
            Filtered AnnData and list of dropped marker names
        """
        if min_std is None:
            min_std = self.config.clustering.min_marker_std

        values = adata.X.toarray() if sparse.issparse(adata.X) else np.asarray(adata.X)
        std = np.std(values, axis=0)
        keep = std >= min_std

        dropped = adata.var_names[np.logical_not(keep)].tolist()
        if dropped:
            self.logger.info(
                "Dropping %d low-variance markers (std < %.4f)", len(dropped), min_std
            )
            adata = adata[:, keep].copy()

        return adata, dropped

    def exclude_technical_markers(
        self,
        adata: Any,  # AnnData
        technical_markers: Optional[Sequence[str]] = None,
    ) -> Tuple[Any, List[str]]:
        """Move technical markers from features to cell metadata.

        Technical markers like DAPI are useful for visualization and QC but
        should not be included in biological analyses like PCA, clustering,
        or differential expression.

        Parameters
        ----------
        adata : AnnData
            Input AnnData object
        technical_markers : Sequence[str], optional
            Marker names to exclude. Uses config default if None.

        Returns
        -------
        Tuple[AnnData, List[str]]
            AnnData with technical markers removed, list of excluded markers
        """
        if technical_markers is None:
            technical_markers = self.config.technical_markers

        excluded = []
        for marker in technical_markers:
            if marker in adata.var_names:
                marker_idx = list(adata.var_names).index(marker)
                if sparse.issparse(adata.X):
                    values = adata.X[:, marker_idx].toarray().flatten()
                else:
                    values = adata.X[:, marker_idx].flatten()

                col_name = f"{marker}_intensity"
                adata.obs[col_name] = values
                excluded.append(marker)
                self.logger.info(
                    "Moved '%s' from features to metadata column '%s' (mean=%.2f, std=%.2f)",
                    marker,
                    col_name,
                    np.nanmean(values),
                    np.nanstd(values),
                )

        if excluded:
            keep_markers = [m for m in adata.var_names if m not in excluded]
            adata = adata[:, keep_markers].copy()
            self.logger.info(
                "Excluded %d technical markers from analysis: %s",
                len(excluded),
                ", ".join(excluded),
            )
            self.logger.info(
                "Feature matrix reduced from %d to %d markers",
                len(adata.var_names) + len(excluded),
                len(adata.var_names),
            )

        return adata, excluded

    def run_clustering(
        self,
        adata: Any,  # AnnData
        cluster_key: str = "cluster_lvl0",
        n_pcs: Optional[int] = None,
        neighbors_k: Optional[int] = None,
        resolution: Optional[float] = None,
        scale_clip: Optional[float] = None,
        random_seed: Optional[int] = None,
        use_gpu: Optional[bool] = None,
        compute_umap: Optional[bool] = None,
    ) -> ClusteringResult:
        """Run the full clustering pipeline.

        Pipeline: scale → PCA → neighbors → Leiden (→ UMAP optional)

        Parameters
        ----------
        adata : AnnData
            Input AnnData object (modified in place)
        cluster_key : str
            Key in adata.obs to store cluster assignments
        n_pcs : int, optional
            Number of principal components. Uses config default if None.
        neighbors_k : int, optional
            k for neighborhood graph. Uses config default if None.
        resolution : float, optional
            Leiden resolution. Uses config default if None.
        scale_clip : float, optional
            Value clipping during scaling. Uses config default if None.
        random_seed : int, optional
            Random seed for reproducibility. Uses config default if None.
        use_gpu : bool, optional
            Use GPU acceleration. Uses config default if None.
        compute_umap : bool, optional
            Compute UMAP embeddings. Uses config default if None.

        Returns
        -------
        ClusteringResult
            Clustering result with cluster statistics
        """
        import scanpy as sc

        # Use config defaults where not specified
        cfg = self.config.clustering
        n_pcs = n_pcs if n_pcs is not None else cfg.n_pcs
        neighbors_k = neighbors_k if neighbors_k is not None else cfg.neighbors_k
        resolution = resolution if resolution is not None else cfg.resolution
        scale_clip = scale_clip if scale_clip is not None else cfg.scale_clip
        random_seed = random_seed if random_seed is not None else cfg.random_seed
        use_gpu = use_gpu if use_gpu is not None else cfg.use_gpu
        compute_umap = compute_umap if compute_umap is not None else cfg.compute_umap

        self.logger.info(
            "Running clustering pipeline: n_pcs=%d, neighbors_k=%d, resolution=%.3f, GPU=%s",
            n_pcs,
            neighbors_k,
            resolution,
            GPU_AVAILABLE and use_gpu,
        )

        # Scale (always CPU - minimal impact)
        sc.pp.scale(adata, zero_center=True, max_value=scale_clip)
        adata.uns["scaled"] = True
        adata.uns["scale_clip"] = scale_clip

        # Adjust n_pcs if needed
        use_pcs = min(n_pcs, max(adata.n_vars - 1, 1), max(adata.n_obs - 1, 1))

        # Choose GPU or CPU path
        if GPU_AVAILABLE and use_gpu:
            self.logger.info("Using GPU acceleration for clustering")
            rsc.get.anndata_to_GPU(adata)
            rsc.pp.pca(adata, n_comps=use_pcs, random_state=random_seed)
            rsc.pp.neighbors(adata, n_neighbors=neighbors_k, n_pcs=use_pcs)
            rsc.tl.leiden(
                adata, resolution=resolution, random_state=random_seed, key_added=cluster_key
            )
            if compute_umap:
                rsc.tl.umap(adata, random_state=random_seed)
            rsc.get.anndata_to_CPU(adata)
        else:
            if not use_gpu and GPU_AVAILABLE:
                self.logger.info("GPU available but disabled by user flag")
            self.logger.info("Using CPU for clustering")
            sc.tl.pca(adata, n_comps=use_pcs, svd_solver="arpack", random_state=random_seed)
            sc.pp.neighbors(adata, n_neighbors=neighbors_k, n_pcs=use_pcs)
            sc.tl.leiden(
                adata,
                resolution=resolution,
                random_state=random_seed,
                key_added=cluster_key,
                flavor="igraph",
                n_iterations=2,
                directed=False,
            )
            if compute_umap:
                sc.tl.umap(adata, random_state=random_seed)

        # Build result
        result = ClusteringResult(cluster_key=cluster_key)
        result.n_clusters = adata.obs[cluster_key].nunique()
        result.cluster_sizes = (
            adata.obs[cluster_key].value_counts().to_dict()
        )

        self.logger.info(
            "Computed Leiden clustering with %d clusters", result.n_clusters
        )

        return result

    def subcluster(
        self,
        adata: Any,  # AnnData
        parent_cluster: str,
        parent_key: str = "cluster_lvl0",
        child_key: str = "cluster_lvl1",
        resolution: Optional[float] = None,
        min_cells: Optional[int] = None,
    ) -> Optional[ClusteringResult]:
        """Subcluster a specific parent cluster.

        Parameters
        ----------
        adata : AnnData
            Input AnnData object (modified in place)
        parent_cluster : str
            Cluster ID to subcluster
        parent_key : str
            Column containing parent cluster assignments
        child_key : str
            Column to store subcluster assignments
        resolution : float, optional
            Leiden resolution for subclustering
        min_cells : int, optional
            Minimum cells required for subclustering

        Returns
        -------
        Optional[ClusteringResult]
            Subclustering result, or None if cluster is too small
        """
        import scanpy as sc

        cfg = self.config.subcluster
        resolution = resolution if resolution is not None else cfg.resolution
        min_cells = min_cells if min_cells is not None else cfg.min_cells

        # Get cells in parent cluster
        parent_mask = adata.obs[parent_key].astype(str) == str(parent_cluster)
        n_cells = parent_mask.sum()

        if n_cells < min_cells:
            self.logger.info(
                "Cluster %s has %d cells (< %d minimum), skipping subclustering",
                parent_cluster,
                n_cells,
                min_cells,
            )
            return None

        self.logger.info(
            "Subclustering cluster %s (%d cells) at resolution %.3f",
            parent_cluster,
            n_cells,
            resolution,
        )

        # Extract and cluster subset
        subset = adata[parent_mask].copy()

        # Use CPU for subclustering (smaller datasets, avoid GPU context issues)
        sc.pp.scale(subset, zero_center=True, max_value=self.config.clustering.scale_clip)
        use_pcs = min(
            self.config.clustering.n_pcs,
            max(subset.n_vars - 1, 1),
            max(subset.n_obs - 1, 1),
        )
        sc.tl.pca(subset, n_comps=use_pcs, svd_solver="arpack")
        sc.pp.neighbors(subset, n_neighbors=self.config.clustering.neighbors_k, n_pcs=use_pcs)
        sc.tl.leiden(
            subset,
            resolution=resolution,
            key_added="subcluster",
            flavor="igraph",
            n_iterations=2,
            directed=False,
        )

        # Build hierarchical IDs and update main adata
        sub_ids = subset.obs["subcluster"].astype(str)
        hierarchical_ids = sub_ids.apply(lambda x: f"{parent_cluster}:{x}")

        # Initialize child_key if not present
        if child_key not in adata.obs:
            adata.obs[child_key] = adata.obs[parent_key].astype(str)

        # Update child_key for subclustered cells
        adata.obs.loc[parent_mask, child_key] = hierarchical_ids.values

        # Build result
        result = ClusteringResult(cluster_key=child_key)
        result.n_clusters = subset.obs["subcluster"].nunique()
        result.cluster_sizes = {
            f"{parent_cluster}:{k}": v
            for k, v in subset.obs["subcluster"].value_counts().to_dict().items()
        }

        self.logger.info(
            "Created %d subclusters for parent %s",
            result.n_clusters,
            parent_cluster,
        )

        return result
