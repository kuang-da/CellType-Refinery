"""
Parallel subclustering for cell-type refinement.

This module provides CPU multiprocessing for subclustering operations.
Instead of processing clusters sequentially, it:
1. Extracts minimal data (numpy arrays) for each cluster
2. Runs subclustering in parallel worker processes
3. Merges results back into the main AnnData

Expected speedup: 3-6x with 4 workers (32 clusters: 224s -> 60s)
"""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import sparse

try:
    import scanpy as sc
except ImportError:
    sc = None

from .plan import SubclusterOp


@dataclass
class SubclusterWorkItem:
    """Minimal data for parallel subclustering.

    This contains only numpy arrays to avoid expensive AnnData serialization.
    """
    cluster_id: str
    cell_indices: np.ndarray      # Index positions in parent adata
    X_subset: np.ndarray          # Dense expression matrix (n_cells, n_markers)
    var_names: List[str]          # Marker names
    resolution: float = 0.4
    n_pcs: int = 30
    neighbors_k: int = 15
    min_cells: int = 100
    min_subcluster_cells: int = 20
    max_subclusters_ratio: float = 0.5
    scale_clip: float = 10.0
    seed: int = 1337


@dataclass
class SubclusterResult:
    """Result from parallel subclustering."""
    cluster_id: str
    cell_indices: np.ndarray
    subcluster_labels: np.ndarray  # Local labels (0, 1, 2, ...)
    n_subclusters: int
    success: bool
    error: Optional[str] = None
    timing_seconds: float = 0.0


def extract_work_items(
    adata: "sc.AnnData",
    subcluster_ops: List[SubclusterOp],
    cluster_col: str = "cluster_lvl0",
    layer: Optional[str] = None,
    logger: Optional[logging.Logger] = None,
) -> List[SubclusterWorkItem]:
    """Extract minimal data for parallel processing.

    Parameters
    ----------
    adata : sc.AnnData
        Source AnnData object
    subcluster_ops : List[SubclusterOp]
        Subclustering operations to prepare
    cluster_col : str
        Column with cluster IDs
    layer : str, optional
        Layer to use for expression data (None = adata.X)
    logger : logging.Logger, optional
        Logger for progress tracking

    Returns
    -------
    List[SubclusterWorkItem]
        Work items ready for parallel processing
    """
    _logger = logger or logging.getLogger(__name__)
    work_items = []
    var_names = list(adata.var_names)

    for op in subcluster_ops:
        cluster_id = str(op.cluster_id)

        # Find cells in this cluster (check both lvl0 and lvl1)
        mask = adata.obs[cluster_col].astype(str) == cluster_id
        if "cluster_lvl1" in adata.obs:
            mask = mask | (adata.obs["cluster_lvl1"].astype(str) == cluster_id)

        cell_indices = np.where(mask)[0]
        n_cells = len(cell_indices)

        if n_cells < op.min_cells:
            _logger.debug(
                "Cluster %s has %d cells (min=%d), skipping work item",
                cluster_id, n_cells, op.min_cells
            )
            continue

        # Extract expression matrix as dense array
        if layer and layer in adata.layers:
            X_subset = adata.layers[layer][cell_indices, :]
        else:
            X_subset = adata.X[cell_indices, :]

        # Convert sparse to dense if needed
        if sparse.issparse(X_subset):
            X_subset = X_subset.toarray()
        else:
            X_subset = np.asarray(X_subset).copy()

        work_items.append(SubclusterWorkItem(
            cluster_id=cluster_id,
            cell_indices=cell_indices,
            X_subset=X_subset,
            var_names=var_names,
            resolution=op.resolution,
            n_pcs=op.n_pcs,
            neighbors_k=op.neighbors_k,
            min_cells=op.min_cells,
            min_subcluster_cells=getattr(op, 'min_subcluster_cells', 20),
            max_subclusters_ratio=getattr(op, 'max_subclusters_ratio', 0.5),
            scale_clip=10.0,
            seed=1337,
        ))

        _logger.debug("Created work item for cluster %s (%d cells)", cluster_id, n_cells)

    return work_items


def worker_subcluster(work_item: SubclusterWorkItem) -> SubclusterResult:
    """Process a single cluster subclustering in a worker process.

    This function is designed to be called via multiprocessing.
    It creates a temporary AnnData, runs clustering, and returns results.

    Parameters
    ----------
    work_item : SubclusterWorkItem
        Work item with cluster data and parameters

    Returns
    -------
    SubclusterResult
        Result with subcluster assignments
    """
    start_time = time.time()

    try:
        import numpy as np
        import pandas as pd
        import scanpy as sc
        from anndata import AnnData

        n_cells = work_item.X_subset.shape[0]

        # Create temporary AnnData for clustering
        adata_temp = AnnData(
            X=work_item.X_subset.copy(),
            var=pd.DataFrame(index=work_item.var_names),
        )

        # Scale data
        sc.pp.scale(adata_temp, max_value=work_item.scale_clip)

        # Adjust n_pcs for small datasets
        n_pcs = min(work_item.n_pcs, adata_temp.n_vars - 1, adata_temp.n_obs - 1)
        if n_pcs < 2:
            # Too few features/cells for meaningful PCA
            return SubclusterResult(
                cluster_id=work_item.cluster_id,
                cell_indices=work_item.cell_indices,
                subcluster_labels=np.zeros(n_cells, dtype=np.int32),
                n_subclusters=1,
                success=True,
                timing_seconds=time.time() - start_time,
            )

        # Run CPU clustering
        _run_cpu_clustering(adata_temp, n_pcs, work_item)

        # Extract labels
        labels = adata_temp.obs["leiden"].astype(int).values
        n_subclusters = len(np.unique(labels))

        # Check for pathological subclustering
        if n_subclusters == 1:
            return SubclusterResult(
                cluster_id=work_item.cluster_id,
                cell_indices=work_item.cell_indices,
                subcluster_labels=labels,
                n_subclusters=1,
                success=True,
                timing_seconds=time.time() - start_time,
            )

        max_allowed = max(2, int(n_cells * work_item.max_subclusters_ratio / work_item.min_subcluster_cells))
        if n_subclusters > max_allowed:
            # Pathological subclustering - return single cluster
            return SubclusterResult(
                cluster_id=work_item.cluster_id,
                cell_indices=work_item.cell_indices,
                subcluster_labels=np.zeros(n_cells, dtype=np.int32),
                n_subclusters=1,
                success=True,
                error=f"Pathological: {n_subclusters} > {max_allowed} max",
                timing_seconds=time.time() - start_time,
            )

        # Merge small subclusters
        labels = _merge_small_subclusters(labels, work_item.min_subcluster_cells)
        n_subclusters = len(np.unique(labels))

        return SubclusterResult(
            cluster_id=work_item.cluster_id,
            cell_indices=work_item.cell_indices,
            subcluster_labels=labels,
            n_subclusters=n_subclusters,
            success=True,
            timing_seconds=time.time() - start_time,
        )

    except Exception as e:
        return SubclusterResult(
            cluster_id=work_item.cluster_id,
            cell_indices=work_item.cell_indices,
            subcluster_labels=np.zeros(len(work_item.cell_indices), dtype=np.int32),
            n_subclusters=0,
            success=False,
            error=str(e),
            timing_seconds=time.time() - start_time,
        )


def _run_cpu_clustering(
    adata_temp: "sc.AnnData",
    n_pcs: int,
    work_item: SubclusterWorkItem,
) -> None:
    """Run clustering on CPU using scanpy."""
    import scanpy as sc

    sc.tl.pca(adata_temp, n_comps=n_pcs, svd_solver="arpack", random_state=work_item.seed)

    neighbors_k = min(work_item.neighbors_k, adata_temp.n_obs - 1)
    sc.pp.neighbors(adata_temp, n_neighbors=neighbors_k, n_pcs=n_pcs)
    sc.tl.leiden(
        adata_temp,
        resolution=work_item.resolution,
        random_state=work_item.seed,
        key_added="leiden",
        flavor="igraph",
        n_iterations=2,
        directed=False,
    )


def _merge_small_subclusters(
    labels: np.ndarray,
    min_cells: int,
) -> np.ndarray:
    """Merge subclusters smaller than min_cells into the largest."""
    unique_labels, counts = np.unique(labels, return_counts=True)

    small_clusters = unique_labels[counts < min_cells]
    if len(small_clusters) == 0:
        return labels

    # Find largest cluster to merge into
    largest_cluster = unique_labels[np.argmax(counts)]

    # Merge small clusters
    merged_labels = labels.copy()
    for small in small_clusters:
        if small != largest_cluster:
            merged_labels[merged_labels == small] = largest_cluster

    # Renumber to be contiguous
    remaining = np.unique(merged_labels)
    renumber_map = {old: new for new, old in enumerate(sorted(remaining))}
    merged_labels = np.array([renumber_map[l] for l in merged_labels], dtype=np.int32)

    return merged_labels


def run_subclusters_parallel(
    adata: "sc.AnnData",
    subcluster_ops: List[SubclusterOp],
    n_workers: int = 4,
    cluster_col: str = "cluster_lvl0",
    output_cluster_col: str = "cluster_lvl1",
    layer: Optional[str] = None,
    logger: Optional[logging.Logger] = None,
) -> Tuple[List[str], List[str]]:
    """Run subclustering operations in parallel.

    Parameters
    ----------
    adata : sc.AnnData
        AnnData to modify (output_cluster_col updated in place)
    subcluster_ops : List[SubclusterOp]
        Subclustering operations to execute
    n_workers : int
        Number of parallel workers (default: 4)
    cluster_col : str
        Column to read cluster IDs from (source)
    output_cluster_col : str
        Column to write hierarchical cluster IDs to (destination)
    layer : str, optional
        Layer to use for expression data
    logger : logging.Logger, optional
        Logger for progress tracking

    Returns
    -------
    Tuple[List[str], List[str]]
        (subclustered_cluster_ids, errors)
    """
    _logger = logger or logging.getLogger(__name__)

    if not subcluster_ops:
        return [], []

    # Extract work items
    _logger.info("Extracting work items for %d subclustering operations...", len(subcluster_ops))
    work_items = extract_work_items(
        adata, subcluster_ops, cluster_col, layer=layer, logger=_logger
    )

    if not work_items:
        _logger.warning("No valid work items after filtering")
        return [], []

    _logger.info(
        "Processing %d clusters with %d workers",
        len(work_items), n_workers
    )

    # Run parallel processing
    start_time = time.time()

    if n_workers == 1:
        # Sequential mode (for debugging or single-core systems)
        results = [worker_subcluster(item) for item in work_items]
    else:
        # Parallel mode with joblib
        try:
            from joblib import Parallel, delayed

            # Use 'loky' backend for better process isolation
            # Use verbose=10 for progress logging
            results = Parallel(n_jobs=n_workers, backend='loky', verbose=10)(
                delayed(worker_subcluster)(item) for item in work_items
            )
        except Exception as e:
            _logger.warning("Parallel execution failed: %s. Falling back to sequential.", e)
            results = [worker_subcluster(item) for item in work_items]

    parallel_time = time.time() - start_time
    _logger.info("Parallel subclustering completed in %.2f seconds", parallel_time)

    # Merge results back into adata
    subclustered_ids, errors = _merge_results(
        adata, results, cluster_col, output_cluster_col, _logger
    )

    return subclustered_ids, errors


def _merge_results(
    adata: "sc.AnnData",
    results: List[SubclusterResult],
    cluster_col: str,
    output_cluster_col: str,
    logger: logging.Logger,
) -> Tuple[List[str], List[str]]:
    """Merge parallel subclustering results into AnnData.

    Parameters
    ----------
    adata : sc.AnnData
        AnnData to update in place
    results : List[SubclusterResult]
        Results from parallel workers
    cluster_col : str
        Source cluster column name (to copy from if output doesn't exist)
    output_cluster_col : str
        Output cluster column name (to write hierarchical labels to)
    logger : logging.Logger
        Logger for progress tracking

    Returns
    -------
    Tuple[List[str], List[str]]
        (subclustered_cluster_ids, errors)
    """
    logger.info("Merging %d results into AnnData (output: %s)...", len(results), output_cluster_col)

    # Initialize output column if needed
    if output_cluster_col not in adata.obs:
        adata.obs[output_cluster_col] = adata.obs[cluster_col].astype(str)
    elif hasattr(adata.obs[output_cluster_col], "cat"):
        adata.obs[output_cluster_col] = adata.obs[output_cluster_col].astype(str)

    subclustered_ids = []
    errors = []
    total_cells_updated = 0

    for result in results:
        if result.success and result.n_subclusters > 1:
            # Create hierarchical labels: "parent:subcluster"
            hierarchical_labels = [
                f"{result.cluster_id}:{label}"
                for label in result.subcluster_labels
            ]

            # Update adata.obs using index positions
            # Use .iloc for positional indexing
            output_col_idx = adata.obs.columns.get_loc(output_cluster_col)
            for i, idx in enumerate(result.cell_indices):
                adata.obs.iloc[idx, output_col_idx] = hierarchical_labels[i]

            subclustered_ids.append(result.cluster_id)
            total_cells_updated += len(result.cell_indices)

            logger.info(
                "  Cluster %s: %d subclusters (%d cells, %.2f sec)",
                result.cluster_id,
                result.n_subclusters,
                len(result.cell_indices),
                result.timing_seconds,
            )
        elif result.success and result.n_subclusters == 1:
            logger.info(
                "  Cluster %s: 1 subcluster (no change, %.2f sec)",
                result.cluster_id,
                result.timing_seconds,
            )
        elif not result.success:
            error_msg = f"Cluster {result.cluster_id}: {result.error}"
            errors.append(error_msg)
            logger.error("  Cluster %s: FAILED - %s", result.cluster_id, result.error)

    logger.info(
        "Merged results: %d clusters subclustered, %d cells updated, %d errors",
        len(subclustered_ids), total_cells_updated, len(errors)
    )

    return subclustered_ids, errors
