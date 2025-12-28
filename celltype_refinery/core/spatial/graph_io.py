"""Spatial graph I/O utilities.

Functions for loading and validating spatial neighbor graphs
stored in NPZ format (CSR sparse matrix components).
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple
import logging

import numpy as np
from scipy import sparse

logger = logging.getLogger(__name__)


def load_spatial_graph(
    graphs_dir: Path,
    sample_id: str,
) -> Optional[sparse.csr_matrix]:
    """Load spatial neighbor graph for a sample.

    Graphs are stored as NPZ files containing CSR sparse matrix components:
    - data: non-zero values (all 1.0 for connectivity)
    - indices: column indices
    - indptr: row pointers
    - shape: matrix dimensions

    Parameters
    ----------
    graphs_dir : Path
        Directory containing graph NPZ files
    sample_id : str
        Sample identifier

    Returns
    -------
    sparse.csr_matrix or None
        Sparse adjacency matrix, or None if file not found
    """
    graphs_dir = Path(graphs_dir)
    graph_file = graphs_dir / f"{sample_id}_neighbors.npz"

    if not graph_file.exists():
        logger.warning(f"Graph file not found: {graph_file}")
        return None

    try:
        data = np.load(graph_file)
        graph = sparse.csr_matrix(
            (data["data"], data["indices"], data["indptr"]),
            shape=data["shape"],
        )
        return graph
    except Exception as e:
        logger.error(f"Failed to load graph for {sample_id}: {e}")
        return None


def get_graph_stats(graph: sparse.csr_matrix) -> Dict[str, float]:
    """Get statistics about a spatial graph.

    Parameters
    ----------
    graph : sparse.csr_matrix
        Sparse adjacency matrix

    Returns
    -------
    Dict[str, float]
        Statistics: n_cells, n_edges, mean_degree, min_degree, max_degree
    """
    n_cells = graph.shape[0]
    n_edges = graph.nnz

    degrees = np.diff(graph.indptr)

    return {
        "n_cells": n_cells,
        "n_edges": n_edges,
        "mean_degree": float(np.mean(degrees)) if n_cells > 0 else 0.0,
        "min_degree": int(np.min(degrees)) if n_cells > 0 else 0,
        "max_degree": int(np.max(degrees)) if n_cells > 0 else 0,
    }


def list_available_graphs(graphs_dir: Path) -> List[str]:
    """List all available sample IDs with spatial graphs.

    Parameters
    ----------
    graphs_dir : Path
        Directory containing graph NPZ files

    Returns
    -------
    List[str]
        List of sample IDs with available graphs
    """
    graphs_dir = Path(graphs_dir)
    if not graphs_dir.exists():
        logger.warning(f"Graphs directory not found: {graphs_dir}")
        return []

    sample_ids = []
    for f in graphs_dir.glob("*_neighbors.npz"):
        sample_id = f.stem.replace("_neighbors", "")
        sample_ids.append(sample_id)

    return sorted(sample_ids)


def validate_graphs_dir(
    graphs_dir: Path,
    expected_samples: Optional[List[str]] = None,
) -> Tuple[List[str], List[str]]:
    """Validate graphs directory and check for missing samples.

    Parameters
    ----------
    graphs_dir : Path
        Directory containing graph NPZ files
    expected_samples : List[str], optional
        List of expected sample IDs to check

    Returns
    -------
    Tuple[List[str], List[str]]
        (available_samples, missing_samples)
    """
    graphs_dir = Path(graphs_dir)

    if not graphs_dir.exists():
        logger.error(f"Graphs directory not found: {graphs_dir}")
        if expected_samples:
            return [], expected_samples
        return [], []

    available = set(list_available_graphs(graphs_dir))

    if expected_samples is None:
        return sorted(available), []

    expected = set(expected_samples)
    missing = expected - available

    return sorted(available & expected), sorted(missing)


def load_graphs_for_samples(
    graphs_dir: Path,
    sample_ids: List[str],
) -> Dict[str, sparse.csr_matrix]:
    """Load spatial graphs for multiple samples.

    Parameters
    ----------
    graphs_dir : Path
        Directory containing graph NPZ files
    sample_ids : List[str]
        Sample IDs to load

    Returns
    -------
    Dict[str, sparse.csr_matrix]
        Dictionary mapping sample_id -> graph
    """
    graphs_dir = Path(graphs_dir)
    graphs = {}

    for sample_id in sample_ids:
        graph = load_spatial_graph(graphs_dir, sample_id)
        if graph is not None:
            graphs[sample_id] = graph
        else:
            logger.warning(f"Skipping sample {sample_id}: graph not found")

    logger.info(f"Loaded {len(graphs)}/{len(sample_ids)} graphs")
    return graphs
