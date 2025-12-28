"""Data merging and spatial graph construction (Stage F).

Merges aligned/corrected data from all samples into a unified
AnnData object with spatial neighborhood graphs.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .config import MergeConfig


@dataclass
class MergeResult:
    """Result from merging samples.

    Attributes
    ----------
    adata : AnnData
        Merged AnnData object
    n_cells : int
        Total number of cells
    n_samples : int
        Number of samples merged
    markers : List[str]
        Marker names
    layers : List[str]
        Layer names (raw, aligned, batchcorr)
    """

    adata: Any = None  # AnnData
    n_cells: int = 0
    n_samples: int = 0
    markers: List[str] = field(default_factory=list)
    layers: List[str] = field(default_factory=list)


class DataMerger:
    """Data merger with spatial graph construction.

    Merges multiple samples into a single AnnData object with
    multiple expression layers and builds per-sample spatial
    neighborhood graphs.

    Parameters
    ----------
    config : MergeConfig
        Merge configuration

    Example
    -------
    >>> from celltype_refinery.core.preprocessing import DataMerger, MergeConfig
    >>> config = MergeConfig(graph_mode="knn", k=10)
    >>> merger = DataMerger(config)
    >>> result = merger.merge_samples(sample_data, metadata)
    """

    def __init__(self, config: Optional[MergeConfig] = None):
        self.config = config or MergeConfig()
        self._check_dependencies()

    def _check_dependencies(self) -> None:
        """Check for required dependencies."""
        try:
            import anndata
        except ImportError:
            raise RuntimeError(
                "Data merging requires anndata. "
                "Install with: pip install anndata"
            )

    def merge_samples(
        self,
        aligned_matrices: Dict[str, pd.DataFrame],
        cell_metadata: Dict[str, pd.DataFrame],
        sample_metadata: Dict[str, Dict[str, Any]],
        raw_matrices: Optional[Dict[str, pd.DataFrame]] = None,
        corrected_matrices: Optional[Dict[str, pd.DataFrame]] = None,
        cell_id_col: str = "cell_mask_id",
    ) -> MergeResult:
        """Merge samples into a single AnnData object.

        Parameters
        ----------
        aligned_matrices : Dict[str, pd.DataFrame]
            Map of sample_id to aligned intensity matrix
        cell_metadata : Dict[str, pd.DataFrame]
            Map of sample_id to cell metadata (coordinates, morphology)
        sample_metadata : Dict[str, Dict[str, Any]]
            Map of sample_id to sample-level metadata (donor, region, etc.)
        raw_matrices : Dict[str, pd.DataFrame], optional
            Map of sample_id to raw/normalized intensity matrix
        corrected_matrices : Dict[str, pd.DataFrame], optional
            Map of sample_id to batch-corrected intensity matrix
        cell_id_col : str
            Cell ID column name

        Returns
        -------
        MergeResult
            Merged result with AnnData object
        """
        import anndata as ad

        result = MergeResult()

        # Collect all data
        merged_tables: List[pd.DataFrame] = []
        raw_tables: List[pd.DataFrame] = []
        corrected_tables: List[pd.DataFrame] = []
        marker_cols: Optional[List[str]] = None

        sample_ids = sorted(aligned_matrices.keys())
        result.n_samples = len(sample_ids)

        for sample_id in sample_ids:
            aligned = aligned_matrices[sample_id]
            meta = cell_metadata.get(sample_id)
            sample_meta = sample_metadata.get(sample_id, {})

            # Get marker columns from first sample
            if marker_cols is None:
                marker_cols = [c for c in aligned.columns if c != cell_id_col]
                result.markers = marker_cols

            # Merge aligned with cell metadata
            if meta is not None:
                # Align on common cell IDs
                common_ids = sorted(
                    set(aligned[cell_id_col]) & set(meta[cell_id_col])
                )
                aligned = aligned[aligned[cell_id_col].isin(common_ids)].copy()
                meta = meta[meta[cell_id_col].isin(common_ids)].copy()

                merged = aligned.merge(
                    meta, on=cell_id_col, how="left", suffixes=("", "_meta")
                )
            else:
                merged = aligned.copy()

            # Add sample-level metadata
            merged["sample_id"] = sample_id
            for key, value in sample_meta.items():
                merged[key] = value

            merged_tables.append(merged)

            # Collect raw layer
            if raw_matrices and sample_id in raw_matrices:
                raw_tables.append(raw_matrices[sample_id])

            # Collect corrected layer
            if corrected_matrices and sample_id in corrected_matrices:
                corrected_tables.append(corrected_matrices[sample_id])

        if marker_cols is None:
            raise ValueError("No markers found in aligned matrices")

        # Concatenate all samples
        merged_df = pd.concat(merged_tables, ignore_index=True)
        result.n_cells = len(merged_df)

        # Build expression matrix
        expression = merged_df[marker_cols].to_numpy(dtype=np.float32)

        # Build obs dataframe
        obs_cols = [cell_id_col, "sample_id"]
        # Add any sample metadata columns
        sample_meta_cols = set()
        for meta in sample_metadata.values():
            sample_meta_cols.update(meta.keys())
        obs_cols.extend(sorted(sample_meta_cols))

        # Add morphology/coordinate columns
        base_set = set(marker_cols + obs_cols)
        extra_cols = [c for c in merged_df.columns if c not in base_set]
        obs_cols.extend(extra_cols)

        obs_df = merged_df[[c for c in obs_cols if c in merged_df.columns]].copy()
        obs_df.index = pd.Index(
            [f"cell_{i}" for i in range(len(obs_df))], name="cell_id"
        )

        # Build var dataframe
        var_df = pd.DataFrame(index=marker_cols)

        # Create AnnData
        adata = ad.AnnData(X=expression, obs=obs_df, var=var_df)
        result.layers.append("X")

        # Add aligned layer (same as X)
        if self.config.include_aligned_layer:
            adata.layers["aligned"] = expression
            result.layers.append("aligned")

        # Add raw layer
        if self.config.include_raw_layer and raw_tables:
            raw_df = pd.concat(raw_tables, ignore_index=True)
            raw_matrix = self._align_layer_to_obs(
                raw_df, adata.obs[cell_id_col], marker_cols, cell_id_col
            )
            adata.layers["raw"] = raw_matrix
            result.layers.append("raw")

        # Add batch-corrected layer
        if self.config.include_batchcorr_layer and corrected_tables:
            corr_df = pd.concat(corrected_tables, ignore_index=True)
            corr_matrix = self._align_layer_to_obs(
                corr_df, adata.obs[cell_id_col], marker_cols, cell_id_col
            )
            adata.layers["batchcorr"] = corr_matrix
            result.layers.append("batchcorr")

        # Add spatial coordinates to obsm
        coord_cols = self._find_coordinate_cols(adata.obs)
        if coord_cols:
            x_col, y_col = coord_cols
            coords = adata.obs[[x_col, y_col]].to_numpy(dtype=np.float32)
            adata.obsm["spatial"] = coords

        result.adata = adata
        return result

    def _align_layer_to_obs(
        self,
        layer_df: pd.DataFrame,
        cell_ids: pd.Series,
        marker_cols: List[str],
        cell_id_col: str,
    ) -> np.ndarray:
        """Align layer data to obs ordering.

        Parameters
        ----------
        layer_df : pd.DataFrame
            Layer data with cell IDs
        cell_ids : pd.Series
            Cell IDs from obs
        marker_cols : List[str]
            Marker column names
        cell_id_col : str
            Cell ID column name

        Returns
        -------
        np.ndarray
            Aligned layer matrix
        """
        layer_df = layer_df.set_index(cell_id_col)
        aligned = layer_df.reindex(cell_ids)
        return aligned[marker_cols].to_numpy(dtype=np.float32)

    def _find_coordinate_cols(
        self, obs: pd.DataFrame
    ) -> Optional[Tuple[str, str]]:
        """Find coordinate columns in obs.

        Returns
        -------
        Optional[Tuple[str, str]]
            (x_col, y_col) or None
        """
        candidates = [
            ("x", "y"),
            ("patch_x", "patch_y"),
            ("cell_centroid_x", "cell_centroid_y"),
        ]
        for x_col, y_col in candidates:
            if x_col in obs.columns and y_col in obs.columns:
                return (x_col, y_col)
        return None

    def build_spatial_graph(
        self,
        coords: np.ndarray,
    ) -> sp.csr_matrix:
        """Build spatial neighborhood graph.

        Parameters
        ----------
        coords : np.ndarray
            Spatial coordinates (n_cells, 2)

        Returns
        -------
        sp.csr_matrix
            Adjacency matrix
        """
        from sklearn.neighbors import NearestNeighbors, radius_neighbors_graph

        if coords.shape[0] == 0:
            return sp.csr_matrix((0, 0))

        if self.config.graph_mode == "knn":
            max_neighbors = coords.shape[0] - 1
            if max_neighbors < 1:
                return sp.csr_matrix((coords.shape[0], coords.shape[0]))

            n_neighbors = min(max(self.config.k, 1), max_neighbors)
            nn = NearestNeighbors(n_neighbors=n_neighbors, metric="euclidean")
            nn.fit(coords)
            graph = nn.kneighbors_graph(mode="connectivity")
            # Remove self-loops
            graph = graph - sp.diags(graph.diagonal())
            return graph.tocsr()

        elif self.config.graph_mode == "ball":
            radius = max(self.config.radius, 1e-6)
            graph = radius_neighbors_graph(
                coords, radius=radius, mode="connectivity", include_self=False
            )
            return graph.tocsr()

        else:
            raise ValueError(f"Unknown graph mode: {self.config.graph_mode}")

    def build_per_sample_graphs(
        self,
        adata: Any,  # AnnData
        output_dir: Optional[Path] = None,
    ) -> Dict[str, sp.csr_matrix]:
        """Build spatial graphs for each sample.

        Parameters
        ----------
        adata : AnnData
            Merged AnnData object with spatial coordinates
        output_dir : Path, optional
            Directory to save graphs as .npz files

        Returns
        -------
        Dict[str, sp.csr_matrix]
            Map of sample_id to adjacency matrix
        """
        graphs: Dict[str, sp.csr_matrix] = {}

        if "spatial" not in adata.obsm:
            return graphs

        sample_ids = adata.obs["sample_id"].unique()

        for sample_id in sample_ids:
            mask = adata.obs["sample_id"] == sample_id
            coords = adata.obsm["spatial"][mask]

            graph = self.build_spatial_graph(coords)
            graphs[sample_id] = graph

            if output_dir is not None:
                output_dir.mkdir(parents=True, exist_ok=True)
                sp.save_npz(output_dir / f"{sample_id}_neighbors.npz", graph)

        return graphs

    @staticmethod
    def compute_degree_distribution(graph: sp.csr_matrix) -> np.ndarray:
        """Compute node degree distribution.

        Parameters
        ----------
        graph : sp.csr_matrix
            Adjacency matrix

        Returns
        -------
        np.ndarray
            Degree per node
        """
        return np.asarray(graph.sum(axis=1)).ravel()
