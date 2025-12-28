"""Mock AnnData generators for testing.

Provides functions to create minimal AnnData objects for testing
without requiring real data.
"""

from typing import Dict, List, Optional

import numpy as np
import pandas as pd


def create_mock_adata(
    n_cells: int = 500,
    n_markers: int = 20,
    n_clusters: int = 5,
    n_samples: int = 2,
    seed: int = 42,
    include_spatial: bool = True,
    include_layers: bool = True,
) -> "AnnData":
    """Create a mock AnnData object for testing.

    Parameters
    ----------
    n_cells : int
        Number of cells
    n_markers : int
        Number of markers/features
    n_clusters : int
        Number of pre-assigned clusters
    n_samples : int
        Number of samples
    seed : int
        Random seed for reproducibility
    include_spatial : bool
        Include spatial coordinates in obsm
    include_layers : bool
        Include aligned/batchcorr layers

    Returns
    -------
    AnnData
        Mock AnnData object suitable for testing
    """
    import anndata as ad

    np.random.seed(seed)

    # Create expression matrix with cluster-specific patterns
    X = np.random.lognormal(mean=0, sigma=1, size=(n_cells, n_markers)).astype(np.float32)

    # Add cluster-specific signal
    cluster_ids = np.repeat(np.arange(n_clusters), n_cells // n_clusters + 1)[:n_cells]
    for i in range(n_clusters):
        mask = cluster_ids == i
        marker_idx = i % n_markers
        X[mask, marker_idx] *= 3  # Boost cluster-specific marker

    # Create obs DataFrame
    sample_ids = np.repeat([f"sample_{i}" for i in range(n_samples)], n_cells // n_samples + 1)[:n_cells]
    donor_ids = np.repeat([f"donor_{i}" for i in range(n_samples)], n_cells // n_samples + 1)[:n_cells]
    region_ids = np.random.choice(["region_A", "region_B", "region_C"], n_cells)

    obs = pd.DataFrame({
        "sample_id": pd.Categorical(sample_ids),
        "donor": pd.Categorical(donor_ids),
        "region": pd.Categorical(region_ids),
        "cluster_lvl0": pd.Categorical([str(c) for c in cluster_ids]),
    })
    obs.index = pd.Index([f"cell_{i}" for i in range(n_cells)], name="cell_id")

    # Create var DataFrame
    marker_names = [f"Marker_{i}" for i in range(n_markers)]
    var = pd.DataFrame(index=marker_names)

    # Create AnnData
    adata = ad.AnnData(X=X, obs=obs, var=var)

    # Add layers
    if include_layers:
        adata.layers["aligned"] = X.copy()
        adata.layers["batchcorr"] = X + np.random.normal(0, 0.1, X.shape).astype(np.float32)

    # Add spatial coordinates
    if include_spatial:
        coords = np.random.uniform(0, 1000, size=(n_cells, 2)).astype(np.float32)
        adata.obsm["spatial"] = coords

    return adata


def create_clustered_adata(
    n_cells: int = 500,
    n_markers: int = 20,
    n_clusters: int = 5,
    seed: int = 42,
) -> "AnnData":
    """Create AnnData with Leiden clustering already computed.

    Returns an AnnData that looks like Stage H output.
    """
    adata = create_mock_adata(
        n_cells=n_cells,
        n_markers=n_markers,
        n_clusters=n_clusters,
        seed=seed,
    )

    # Add clustering metadata
    adata.uns["leiden_resolution"] = 0.6
    adata.uns["n_clusters"] = n_clusters

    return adata


def create_annotated_adata(
    n_cells: int = 500,
    n_markers: int = 20,
    seed: int = 42,
    cell_types: Optional[List[str]] = None,
) -> "AnnData":
    """Create AnnData with cell type annotations.

    Returns an AnnData that looks like post-annotation output.
    """
    if cell_types is None:
        cell_types = [
            "Epithelium",
            "Immune Cells",
            "Mesenchymal Cells",
            "Endothelium",
            "Unassigned",
        ]

    n_clusters = len(cell_types)
    adata = create_mock_adata(
        n_cells=n_cells,
        n_markers=n_markers,
        n_clusters=n_clusters,
        seed=seed,
    )

    # Map clusters to cell types
    cluster_to_type = {str(i): cell_types[i] for i in range(n_clusters)}
    adata.obs["cell_type_lvl0_auto"] = adata.obs["cluster_lvl0"].map(cluster_to_type)
    adata.obs["root_label"] = adata.obs["cell_type_lvl0_auto"]

    # Add annotation scores
    np.random.seed(seed)
    adata.obs["assigned_score"] = np.random.uniform(0.5, 3.0, n_cells)

    return adata


def create_minimal_expression_matrix(
    n_cells: int = 100,
    markers: Optional[List[str]] = None,
    cell_id_col: str = "cell_mask_id",
    seed: int = 42,
) -> pd.DataFrame:
    """Create minimal expression matrix DataFrame.

    Returns a DataFrame suitable for preprocessing stages.
    """
    np.random.seed(seed)

    if markers is None:
        markers = [f"Marker_{i}" for i in range(10)]

    data = {cell_id_col: list(range(n_cells))}
    for marker in markers:
        data[marker] = np.random.lognormal(mean=0, sigma=1, size=n_cells)

    return pd.DataFrame(data)


def create_minimal_cell_metadata(
    n_cells: int = 100,
    cell_id_col: str = "cell_mask_id",
    seed: int = 42,
) -> pd.DataFrame:
    """Create minimal cell metadata DataFrame.

    Returns a DataFrame with cell coordinates and morphology.
    """
    np.random.seed(seed)

    return pd.DataFrame({
        cell_id_col: list(range(n_cells)),
        "x": np.random.uniform(0, 1000, n_cells),
        "y": np.random.uniform(0, 1000, n_cells),
        "cell_area": np.random.uniform(50, 500, n_cells),
        "nucleus_area": np.random.uniform(10, 100, n_cells),
    })
