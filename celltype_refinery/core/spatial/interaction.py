"""Cell-type interaction score analysis.

Computes normalized contact probability (log2 fold enrichment)
between cell types based on spatial neighbor relationships.
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Tuple

import numpy as np
import pandas as pd
from scipy import sparse


@dataclass
class InteractionResult:
    """Result from interaction score analysis.

    Attributes
    ----------
    interaction_matrix : pd.DataFrame
        Square matrix of log2 fold enrichment scores
    significant_interactions : pd.DataFrame
        Long-format table of interactions above threshold
    cell_types : List[str]
        Cell types in the matrix
    threshold : float
        Threshold used for significance
    """

    interaction_matrix: pd.DataFrame = field(default_factory=pd.DataFrame)
    significant_interactions: pd.DataFrame = field(default_factory=pd.DataFrame)
    cell_types: List[str] = field(default_factory=list)
    threshold: float = 0.5

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            "n_cell_types": len(self.cell_types),
            "n_significant_interactions": len(self.significant_interactions),
            "threshold": self.threshold,
        }


def compute_neighbor_counts(
    graph: sparse.csr_matrix,
    cell_types: np.ndarray,
    unique_types: List[str],
) -> np.ndarray:
    """Count neighbor pairs by cell type.

    Parameters
    ----------
    graph : sparse.csr_matrix
        Sparse adjacency matrix
    cell_types : np.ndarray
        Cell type labels for each cell
    unique_types : List[str]
        List of unique cell types

    Returns
    -------
    np.ndarray
        Matrix of neighbor pair counts (n_types x n_types)
    """
    n_types = len(unique_types)
    type_to_idx = {t: i for i, t in enumerate(unique_types)}

    type_indices = np.array([type_to_idx.get(ct, -1) for ct in cell_types])
    counts = np.zeros((n_types, n_types), dtype=np.float64)

    coo = graph.tocoo()

    for i, j in zip(coo.row, coo.col):
        ti, tj = type_indices[i], type_indices[j]
        if ti >= 0 and tj >= 0:
            counts[ti, tj] += 1

    return counts


def compute_neighbor_counts_vectorized(
    graph: sparse.csr_matrix,
    cell_types: np.ndarray,
    unique_types: List[str],
) -> np.ndarray:
    """Vectorized version of neighbor counting (faster for large graphs).

    Parameters
    ----------
    graph : sparse.csr_matrix
        Sparse adjacency matrix
    cell_types : np.ndarray
        Cell type labels
    unique_types : List[str]
        List of unique cell types

    Returns
    -------
    np.ndarray
        Matrix of neighbor pair counts
    """
    n_types = len(unique_types)
    type_to_idx = {t: i for i, t in enumerate(unique_types)}

    type_indices = np.array([type_to_idx.get(ct, -1) for ct in cell_types])

    coo = graph.tocoo()
    valid_mask = (type_indices[coo.row] >= 0) & (type_indices[coo.col] >= 0)

    src_types = type_indices[coo.row[valid_mask]]
    dst_types = type_indices[coo.col[valid_mask]]

    counts = np.zeros((n_types, n_types), dtype=np.float64)
    flat_idx = src_types * n_types + dst_types
    np.add.at(counts.ravel(), flat_idx, 1)

    return counts


def compute_interaction_scores(
    graph: sparse.csr_matrix,
    cell_types: np.ndarray,
    unique_types: List[str],
) -> np.ndarray:
    """Compute cell-type interaction scores.

    Interaction score is the log2 fold enrichment of contact probability:
        score = log2(P(A|B) / P(A))

    Parameters
    ----------
    graph : sparse.csr_matrix
        Sparse adjacency matrix
    cell_types : np.ndarray
        Cell type labels
    unique_types : List[str]
        List of unique cell types

    Returns
    -------
    np.ndarray
        Matrix of interaction scores
    """
    n_types = len(unique_types)
    type_to_idx = {t: i for i, t in enumerate(unique_types)}

    type_counts = np.zeros(n_types)
    for ct in cell_types:
        idx = type_to_idx.get(ct, -1)
        if idx >= 0:
            type_counts[idx] += 1

    total_cells = type_counts.sum()
    if total_cells == 0:
        return np.zeros((n_types, n_types))

    type_probs = type_counts / total_cells

    neighbor_counts = compute_neighbor_counts_vectorized(graph, cell_types, unique_types)

    row_totals = neighbor_counts.sum(axis=1, keepdims=True)
    row_totals[row_totals == 0] = 1

    conditional_probs = neighbor_counts / row_totals

    expected_probs = type_probs[np.newaxis, :]
    expected_probs = np.where(expected_probs == 0, 1e-10, expected_probs)

    interaction_scores = np.log2(conditional_probs / expected_probs + 1e-10)

    return interaction_scores


def compute_interactions_for_sample(
    sample_id: str,
    graph: sparse.csr_matrix,
    cell_types: np.ndarray,
    valid_types: List[str],
) -> Tuple[np.ndarray, int, int]:
    """Compute interaction scores for a single sample.

    Parameters
    ----------
    sample_id : str
        Sample identifier
    graph : sparse.csr_matrix
        Sparse adjacency matrix
    cell_types : np.ndarray
        Cell type labels
    valid_types : List[str]
        Cell types to include

    Returns
    -------
    Tuple[np.ndarray, int, int]
        (interaction_matrix, n_cells, n_edges)
    """
    interactions = compute_interaction_scores(graph, cell_types, valid_types)
    return interactions, len(cell_types), graph.nnz


def aggregate_interaction_results(
    sample_interactions: List[np.ndarray],
    valid_types: List[str],
    threshold: float = 0.5,
) -> InteractionResult:
    """Aggregate interaction results across samples.

    Parameters
    ----------
    sample_interactions : List[np.ndarray]
        Interaction matrices from each sample
    valid_types : List[str]
        Cell types (column/row labels)
    threshold : float
        Threshold for significant interactions

    Returns
    -------
    InteractionResult
        Aggregated results
    """
    if not sample_interactions:
        return InteractionResult()

    mean_interactions = np.mean(sample_interactions, axis=0)

    interaction_matrix = pd.DataFrame(
        mean_interactions, index=valid_types, columns=valid_types
    )

    sig_interactions = []
    for i, ct1 in enumerate(valid_types):
        for j, ct2 in enumerate(valid_types):
            score = mean_interactions[i, j]
            if abs(score) > threshold:
                sig_interactions.append({
                    "source": ct1,
                    "target": ct2,
                    "interaction_score": round(score, 4),
                    "direction": "enriched" if score > 0 else "depleted",
                })

    sig_df = pd.DataFrame(sig_interactions)
    if not sig_df.empty:
        sig_df = sig_df.sort_values("interaction_score", key=abs, ascending=False)

    return InteractionResult(
        interaction_matrix=interaction_matrix,
        significant_interactions=sig_df,
        cell_types=valid_types,
        threshold=threshold,
    )
