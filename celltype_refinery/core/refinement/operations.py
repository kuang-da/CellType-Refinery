"""
Shared operations for cell-type annotation refinement.

This module provides common operations used by the RefinementEngine:
- run_subcluster: Re-cluster a specific cluster with custom parameters
- create_hierarchical_labels: Create "parent:subcluster" style labels
"""

from __future__ import annotations

import logging
import warnings
from typing import List, Optional

import numpy as np
import pandas as pd

try:
    import scanpy as sc
except ImportError:
    sc = None


def run_subcluster(
    adata: "sc.AnnData",
    parent_cluster: str,
    cluster_col: str = "cluster_lvl1",
    resolution: float = 0.4,
    n_pcs: int = 30,
    neighbors_k: int = 15,
    focus_markers: Optional[List[str]] = None,
    min_cells: int = 100,
    min_subcluster_cells: int = 20,
    max_subclusters_ratio: float = 0.5,
    scale_clip: float = 10.0,
    seed: int = 1337,
    logger: Optional[logging.Logger] = None,
) -> int:
    """Subcluster a single parent cluster and create hierarchical labels.

    Creates labels like "3:0", "3:1", "3:2" for subclusters of parent cluster "3".
    Updates the specified cluster column in place.

    Parameters
    ----------
    adata : sc.AnnData
        AnnData object to modify in place
    parent_cluster : str
        ID of parent cluster to subcluster
    cluster_col : str
        Name of cluster column to read from and write to (default: "cluster_lvl1").
        For iterative refinement, use "cluster_lvl2", "cluster_lvl3", etc.
    resolution : float
        Leiden clustering resolution (default: 0.4)
    n_pcs : int
        Number of principal components (default: 30)
    neighbors_k : int
        Number of neighbors for k-NN graph (default: 15)
    focus_markers : List[str], optional
        If provided, only use these markers for clustering
    min_cells : int
        Minimum cells to proceed with subclustering (default: 100)
    min_subcluster_cells : int
        Merge subclusters smaller than this (default: 20)
    max_subclusters_ratio : float
        Prevent pathological fragmentation (default: 0.5)
    scale_clip : float
        Clip scaled values at ±this value (default: 10.0)
    seed : int
        Random seed for reproducibility
    logger : logging.Logger, optional
        Logger for tracking operations

    Returns
    -------
    int
        Number of cells modified (0 if subclustering failed or was skipped)
    """
    _logger = logger or logging.getLogger(__name__)
    parent_cluster = str(parent_cluster)

    # Find cells in parent cluster using the specified cluster column
    if cluster_col not in adata.obs:
        _logger.error("Cluster column '%s' not found in adata.obs", cluster_col)
        return 0

    mask = adata.obs[cluster_col].astype(str) == parent_cluster

    n_cells = mask.sum()
    _logger.info("Subclustering cluster %s: %d cells", parent_cluster, n_cells)

    if n_cells < min_cells:
        _logger.warning(
            "Cluster %s has only %d cells (min=%d); skipping subclustering",
            parent_cluster, n_cells, min_cells,
        )
        return 0

    # Subset to parent cluster
    # Suppress SettingWithCopyWarning from anndata's internal category cleanup
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=pd.errors.SettingWithCopyWarning)
        warnings.filterwarnings("ignore", message=".*SettingWithCopyWarning.*")
        subset = adata[mask].copy()

    # Optional: filter to focus markers
    if focus_markers:
        available_markers = [m for m in focus_markers if m in subset.var_names]
        if len(available_markers) < len(focus_markers):
            missing = set(focus_markers) - set(available_markers)
            _logger.warning("Focus markers not found: %s", missing)
        if available_markers:
            subset = subset[:, available_markers].copy()
            _logger.info("Using %d focus markers: %s", len(available_markers), available_markers)

    # Run clustering pipeline
    try:
        success = _run_clustering_pipeline(
            subset=subset,
            n_pcs=n_pcs,
            neighbors_k=neighbors_k,
            resolution=resolution,
            scale_clip=scale_clip,
            seed=seed,
            logger=_logger,
        )
        if not success:
            _logger.warning("Clustering pipeline failed for cluster %s", parent_cluster)
            return 0
    except Exception as e:
        _logger.error("Clustering failed for cluster %s: %s", parent_cluster, e)
        return 0

    # Get subcluster assignments
    local_clusters = subset.obs["local_leiden"].astype(str)
    n_subclusters = local_clusters.nunique()

    _logger.info("Found %d subclusters in cluster %s", n_subclusters, parent_cluster)

    # Check for pathological subclustering
    if n_subclusters == 1:
        _logger.info("Only 1 subcluster found; keeping original cluster label")
        return 0

    max_allowed = max(2, int(n_cells * max_subclusters_ratio / min_subcluster_cells))
    if n_subclusters > max_allowed:
        _logger.warning(
            "Pathological subclustering: %d subclusters > max allowed %d; aborting",
            n_subclusters, max_allowed,
        )
        return 0

    # Merge small subclusters
    local_clusters = _merge_small_subclusters(
        local_clusters=local_clusters,
        min_cells=min_subcluster_cells,
        logger=_logger,
    )

    # Create hierarchical labels
    hierarchical_labels = create_hierarchical_labels(parent_cluster, local_clusters)

    # Ensure cluster column is string dtype, not Categorical
    # (Categorical can't accept new categories on assignment)
    if hasattr(adata.obs[cluster_col], "cat"):
        adata.obs[cluster_col] = adata.obs[cluster_col].astype(str)

    # Write hierarchical labels to the specified cluster column
    adata.obs.loc[subset.obs_names, cluster_col] = hierarchical_labels.values

    final_n_subclusters = hierarchical_labels.nunique()
    _logger.info(
        "Subclustering complete: cluster %s → %d subclusters (%s), %d cells modified",
        parent_cluster,
        final_n_subclusters,
        sorted(hierarchical_labels.unique()),
        n_cells,
    )

    return n_cells


def _run_clustering_pipeline(
    subset: "sc.AnnData",
    n_pcs: int,
    neighbors_k: int,
    resolution: float,
    scale_clip: float,
    seed: int,
    logger: logging.Logger,
) -> bool:
    """Run the clustering pipeline on a subset using CPU (scanpy).

    Returns True if successful.
    """
    return _run_clustering_scanpy(
        subset=subset,
        n_pcs=n_pcs,
        neighbors_k=neighbors_k,
        resolution=resolution,
        scale_clip=scale_clip,
        seed=seed,
        logger=logger,
    )


def _run_clustering_scanpy(
    subset: "sc.AnnData",
    n_pcs: int,
    neighbors_k: int,
    resolution: float,
    scale_clip: float,
    seed: int,
    logger: logging.Logger,
) -> bool:
    """Fallback clustering using pure scanpy."""
    import scanpy as sc

    # Scale features
    sc.pp.scale(subset, max_value=scale_clip)

    # PCA
    n_pcs_actual = min(n_pcs, subset.n_vars - 1, subset.n_obs - 1)
    if n_pcs_actual < 2:
        logger.warning("Not enough features/cells for PCA")
        return False

    sc.tl.pca(subset, n_comps=n_pcs_actual, random_state=seed)

    # Neighbors
    neighbors_k_actual = min(neighbors_k, subset.n_obs - 1)
    sc.pp.neighbors(subset, n_neighbors=neighbors_k_actual, n_pcs=n_pcs_actual, random_state=seed)

    # Leiden
    sc.tl.leiden(
        subset,
        resolution=resolution,
        key_added="local_leiden",
        random_state=seed,
        flavor="igraph",
        n_iterations=2,
        directed=False,
    )

    return True


def _merge_small_subclusters(
    local_clusters: pd.Series,
    min_cells: int,
    logger: logging.Logger,
) -> pd.Series:
    """Merge subclusters smaller than min_cells into the largest subcluster."""
    cluster_counts = local_clusters.value_counts()
    small_clusters = cluster_counts[cluster_counts < min_cells].index.tolist()

    if not small_clusters:
        return local_clusters

    # Find largest cluster to merge into
    largest_cluster = cluster_counts.idxmax()

    n_merged = 0
    merged_clusters = local_clusters.copy()
    for small in small_clusters:
        if small == largest_cluster:
            continue
        n_cells = (local_clusters == small).sum()
        merged_clusters = merged_clusters.replace(small, largest_cluster)
        n_merged += 1
        logger.info("Merged small subcluster %s (%d cells) into %s", small, n_cells, largest_cluster)

    if n_merged > 0:
        logger.info("Merged %d small subclusters", n_merged)

    # Renumber remaining clusters to be contiguous
    remaining = sorted(merged_clusters.unique())
    renumber_map = {old: str(i) for i, old in enumerate(remaining)}
    merged_clusters = merged_clusters.map(renumber_map)

    return merged_clusters


def create_hierarchical_labels(parent_id: str, subcluster_ids: pd.Series) -> pd.Series:
    """Create hierarchical labels like "3:0", "3:1" from parent and subcluster IDs.

    Parameters
    ----------
    parent_id : str
        Parent cluster ID
    subcluster_ids : pd.Series
        Series of local subcluster IDs (0, 1, 2, ...)

    Returns
    -------
    pd.Series
        Hierarchical labels ("3:0", "3:1", etc.)
    """
    return subcluster_ids.astype(str).apply(lambda x: f"{parent_id}:{x}")


def validate_cluster_exists(
    adata: "sc.AnnData",
    cluster_id: str,
    cluster_col: Optional[str] = None,
) -> bool:
    """Check if a cluster ID exists in the data.

    Parameters
    ----------
    adata : sc.AnnData
        AnnData object
    cluster_id : str
        Cluster ID to check
    cluster_col : str, optional
        Specific cluster column to check. If None, checks the highest
        cluster_lvl column found (most refined view).

    Returns
    -------
    bool
        True if cluster exists in the specified (or detected) cluster column
    """
    if cluster_col is None:
        # Use the highest cluster level available
        iteration = get_current_iteration(adata)
        cluster_col = f"cluster_lvl{iteration}" if iteration > 0 else "cluster_lvl0"

    if cluster_col not in adata.obs:
        return False

    return (adata.obs[cluster_col].astype(str) == cluster_id).any()


def get_cluster_cells(
    adata: "sc.AnnData",
    cluster_id: str,
    cluster_col: Optional[str] = None,
) -> int:
    """Get number of cells in a cluster.

    Parameters
    ----------
    adata : sc.AnnData
        AnnData object
    cluster_id : str
        Cluster ID
    cluster_col : str, optional
        Specific cluster column to check. If None, checks the highest
        cluster_lvl column found (most refined view).

    Returns
    -------
    int
        Number of cells in the cluster
    """
    if cluster_col is None:
        # Use the highest cluster level available
        iteration = get_current_iteration(adata)
        cluster_col = f"cluster_lvl{iteration}" if iteration > 0 else "cluster_lvl0"

    if cluster_col not in adata.obs:
        return 0

    mask = adata.obs[cluster_col].astype(str) == cluster_id
    return mask.sum()


def get_current_iteration(adata: "sc.AnnData") -> int:
    """Determine current iteration from existing cluster_lvl columns.

    Scans for cluster_lvl0, cluster_lvl1, cluster_lvl2, ... columns
    and returns the highest iteration number found.

    Parameters
    ----------
    adata : sc.AnnData
        AnnData object to inspect

    Returns
    -------
    int
        Current iteration number:
        - 0 if only cluster_lvl0 exists (no refinement yet)
        - N if cluster_lvl{N} exists (after N iterations of refinement)

    Example
    -------
    >>> # Only cluster_lvl0
    >>> get_current_iteration(adata)
    0
    >>> # After iteration 1 (cluster_lvl0 + cluster_lvl1)
    >>> get_current_iteration(adata)
    1
    >>> # After iteration 2 (cluster_lvl0 + cluster_lvl1 + cluster_lvl2)
    >>> get_current_iteration(adata)
    2
    """
    iteration = 0
    while f"cluster_lvl{iteration + 1}" in adata.obs:
        iteration += 1
    return iteration


def initialize_iteration_level(
    adata: "sc.AnnData",
    logger: Optional[logging.Logger] = None,
) -> int:
    """Initialize columns for the next iteration level.

    Creates new cluster_lvl{N+1} and cell_type_lvl{N+1} columns by copying
    from the previous level. This preserves full lineage while preparing
    for new subclustering operations.

    Parameters
    ----------
    adata : sc.AnnData
        AnnData object to modify in place
    logger : logging.Logger, optional
        Logger for tracking operations

    Returns
    -------
    int
        The new iteration number (N+1)

    Example
    -------
    >>> # After initial clustering (has cluster_lvl0)
    >>> initialize_iteration_level(adata)
    1  # Created cluster_lvl1, cell_type_lvl1
    >>> # After iteration 1 (has cluster_lvl0, cluster_lvl1)
    >>> initialize_iteration_level(adata)
    2  # Created cluster_lvl2, cell_type_lvl2

    Notes
    -----
    The additive levels design preserves full lineage:
    - cluster_lvl0: Original clusters (never modified)
    - cluster_lvl1: After iteration 1 (immutable after it1)
    - cluster_lvl2: After iteration 2 (immutable after it2)
    - etc.

    Non-subclustered cells get their cluster ID copied forward from
    the previous level, ensuring every cell has a valid ID at every level.
    """
    _logger = logger or logging.getLogger(__name__)

    current = get_current_iteration(adata)
    new_level = current + 1

    _logger.info(
        "Initializing iteration level %d (current: %d)",
        new_level, current
    )

    # Initialize cluster column (copy forward from previous level)
    prev_cluster_col = f"cluster_lvl{current}"
    new_cluster_col = f"cluster_lvl{new_level}"

    if prev_cluster_col not in adata.obs:
        raise ValueError(
            f"Cannot initialize iteration {new_level}: "
            f"previous cluster column '{prev_cluster_col}' not found"
        )

    adata.obs[new_cluster_col] = adata.obs[prev_cluster_col].astype(str).copy()
    _logger.debug(
        "Created %s from %s (%d unique clusters)",
        new_cluster_col, prev_cluster_col,
        adata.obs[new_cluster_col].nunique()
    )

    # Initialize cell_type column
    # For level 0, use cell_type_lvl0_auto (initial output)
    # For level N>0, use cell_type_lvl{N}
    if current == 0:
        prev_label_col = "cell_type_lvl0_auto"
    else:
        prev_label_col = f"cell_type_lvl{current}"

    new_label_col = f"cell_type_lvl{new_level}"

    if prev_label_col in adata.obs:
        adata.obs[new_label_col] = adata.obs[prev_label_col].astype(str).copy()
        _logger.debug(
            "Created %s from %s",
            new_label_col, prev_label_col
        )
    else:
        _logger.warning(
            "Previous label column '%s' not found; %s not created",
            prev_label_col, new_label_col
        )

    # Track iteration in uns
    adata.uns["current_iteration"] = new_level

    _logger.info(
        "Iteration level %d initialized: %s, %s",
        new_level, new_cluster_col, new_label_col
    )

    return new_level
