"""Neighborhood enrichment analysis with permutation tests.

Computes neighborhood enrichment z-scores to identify cell-type pairs
that co-locate significantly more or less than random expectation.
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Tuple

import numpy as np
import pandas as pd
from scipy import sparse

from .interaction import compute_neighbor_counts_vectorized


@dataclass
class NeighborhoodEnrichmentResult:
    """Result from neighborhood enrichment analysis.

    Attributes
    ----------
    z_scores : pd.DataFrame
        Square matrix of z-scores (mean across samples)
    p_values_raw : pd.DataFrame
        Raw p-values (two-tailed permutation test)
    p_values_adj : pd.DataFrame
        FDR-corrected p-values
    enrichment_pairs : pd.DataFrame
        Long-format table with all pairs
    sample_stats : pd.DataFrame
        Per-sample statistics
    cell_types : List[str]
        Cell types in the analysis
    n_permutations : int
        Number of permutations used
    correction_method : str
        Multiple testing correction method
    alpha : float
        Significance threshold
    """

    z_scores: pd.DataFrame = field(default_factory=pd.DataFrame)
    p_values_raw: pd.DataFrame = field(default_factory=pd.DataFrame)
    p_values_adj: pd.DataFrame = field(default_factory=pd.DataFrame)
    enrichment_pairs: pd.DataFrame = field(default_factory=pd.DataFrame)
    sample_stats: pd.DataFrame = field(default_factory=pd.DataFrame)
    cell_types: List[str] = field(default_factory=list)
    n_permutations: int = 100
    correction_method: str = "fdr_bh"
    alpha: float = 0.05

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        n_sig = (
            self.enrichment_pairs["significant"].sum()
            if not self.enrichment_pairs.empty
            else 0
        )
        return {
            "n_cell_types": len(self.cell_types),
            "n_permutations": self.n_permutations,
            "correction_method": self.correction_method,
            "alpha": self.alpha,
            "n_significant_pairs": int(n_sig),
            "n_samples": len(self.sample_stats),
        }


def compute_enrichment_zscores(
    observed: np.ndarray,
    null_counts: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute z-scores and p-values from observed vs null distribution.

    Parameters
    ----------
    observed : np.ndarray
        Observed neighbor counts (n_types x n_types)
    null_counts : np.ndarray
        Null distribution (n_permutations x n_types x n_types)

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        (z_scores, p_values) matrices
    """
    null_mean = null_counts.mean(axis=0)
    null_std = null_counts.std(axis=0)
    null_std = np.where(null_std == 0, 1, null_std)

    z_scores = (observed - null_mean) / null_std

    n_types = observed.shape[0]
    n_perms = null_counts.shape[0]
    p_values = np.zeros((n_types, n_types))

    for i in range(n_types):
        for j in range(n_types):
            extreme = np.sum(
                np.abs(null_counts[:, i, j] - null_mean[i, j])
                >= np.abs(observed[i, j] - null_mean[i, j])
            )
            p_values[i, j] = (extreme + 1) / (n_perms + 1)

    return z_scores, p_values


def _run_permutation(
    graph: sparse.csr_matrix,
    cell_types: np.ndarray,
    unique_types: List[str],
    seed: int,
) -> np.ndarray:
    """Run single permutation for null distribution."""
    rng = np.random.default_rng(seed)
    shuffled = rng.permutation(cell_types)
    return compute_neighbor_counts_vectorized(graph, shuffled, unique_types)


def compute_neighborhood_enrichment_sample(
    graph: sparse.csr_matrix,
    cell_types: np.ndarray,
    unique_types: List[str],
    n_permutations: int = 100,
    n_threads: int = 1,
    seed: int = 42,
) -> Tuple[np.ndarray, np.ndarray, int, int]:
    """Compute neighborhood enrichment for a single sample.

    Parameters
    ----------
    graph : sparse.csr_matrix
        Sparse adjacency matrix
    cell_types : np.ndarray
        Cell type labels
    unique_types : List[str]
        List of unique cell types
    n_permutations : int
        Number of permutations for null model
    n_threads : int
        Threads for parallel permutations
    seed : int
        Random seed

    Returns
    -------
    Tuple[np.ndarray, np.ndarray, int, int]
        (z_scores, p_values, n_cells, n_edges)
    """
    observed = compute_neighbor_counts_vectorized(graph, cell_types, unique_types)

    # Run permutations
    n_types = len(unique_types)
    null_counts = np.zeros((n_permutations, n_types, n_types))

    if n_threads > 1:
        try:
            from joblib import Parallel, delayed
            results = Parallel(n_jobs=n_threads, batch_size=16)(
                delayed(_run_permutation)(graph, cell_types, unique_types, seed + i)
                for i in range(n_permutations)
            )
            null_counts = np.array(results)
        except ImportError:
            for i in range(n_permutations):
                null_counts[i] = _run_permutation(
                    graph, cell_types, unique_types, seed + i
                )
    else:
        for i in range(n_permutations):
            null_counts[i] = _run_permutation(
                graph, cell_types, unique_types, seed + i
            )

    z_scores, p_values = compute_enrichment_zscores(observed, null_counts)

    return z_scores, p_values, len(cell_types), graph.nnz


def apply_fdr_correction(
    p_values: np.ndarray,
    method: str = "fdr_bh",
) -> np.ndarray:
    """Apply multiple testing correction to p-values.

    Parameters
    ----------
    p_values : np.ndarray
        Raw p-values
    method : str
        Correction method: "fdr_bh", "bonferroni", "holm", "none"

    Returns
    -------
    np.ndarray
        Corrected p-values
    """
    original_shape = p_values.shape
    flat_pvals = p_values.ravel()
    n_tests = len(flat_pvals)

    if method == "none":
        return p_values.copy()

    elif method == "bonferroni":
        adjusted = np.clip(flat_pvals * n_tests, 0, 1)

    elif method == "fdr_bh":
        sorted_idx = np.argsort(flat_pvals)
        sorted_pvals = flat_pvals[sorted_idx]

        ranks = np.arange(1, n_tests + 1)
        adjusted_sorted = sorted_pvals * n_tests / ranks

        adjusted_sorted = np.minimum.accumulate(adjusted_sorted[::-1])[::-1]
        adjusted_sorted = np.clip(adjusted_sorted, 0, 1)

        adjusted = np.empty(n_tests)
        adjusted[sorted_idx] = adjusted_sorted

    elif method == "holm":
        sorted_idx = np.argsort(flat_pvals)
        sorted_pvals = flat_pvals[sorted_idx]

        adjusted_sorted = sorted_pvals * (n_tests - np.arange(n_tests))

        adjusted_sorted = np.maximum.accumulate(adjusted_sorted)
        adjusted_sorted = np.clip(adjusted_sorted, 0, 1)

        adjusted = np.empty(n_tests)
        adjusted[sorted_idx] = adjusted_sorted

    else:
        raise ValueError(f"Unknown correction method: {method}")

    return adjusted.reshape(original_shape)


def aggregate_enrichment_results(
    sample_results: List[Tuple[np.ndarray, np.ndarray, int, int]],
    sample_ids: List[str],
    valid_types: List[str],
    correction_method: str = "fdr_bh",
    alpha: float = 0.05,
    n_permutations: int = 100,
) -> NeighborhoodEnrichmentResult:
    """Aggregate neighborhood enrichment results across samples.

    Parameters
    ----------
    sample_results : List[Tuple]
        List of (z_scores, p_values, n_cells, n_edges) per sample
    sample_ids : List[str]
        Sample identifiers
    valid_types : List[str]
        Cell types (row/column labels)
    correction_method : str
        FDR correction method
    alpha : float
        Significance threshold
    n_permutations : int
        Number of permutations used

    Returns
    -------
    NeighborhoodEnrichmentResult
        Aggregated results
    """
    if not sample_results:
        return NeighborhoodEnrichmentResult(
            correction_method=correction_method,
            alpha=alpha,
            n_permutations=n_permutations,
        )

    all_z_scores = [r[0] for r in sample_results]
    all_p_values = [r[1] for r in sample_results]
    sample_stats = [
        {
            "sample_id": sid,
            "n_cells": r[2],
            "n_edges": r[3],
            "mean_z_abs": float(np.abs(r[0]).mean()),
        }
        for sid, r in zip(sample_ids, sample_results)
    ]

    mean_z_scores = np.mean(all_z_scores, axis=0)
    mean_p_values = np.mean(all_p_values, axis=0)

    adj_p_values = apply_fdr_correction(mean_p_values, method=correction_method)

    z_scores_df = pd.DataFrame(mean_z_scores, index=valid_types, columns=valid_types)
    p_values_raw_df = pd.DataFrame(
        mean_p_values, index=valid_types, columns=valid_types
    )
    p_values_adj_df = pd.DataFrame(adj_p_values, index=valid_types, columns=valid_types)
    sample_stats_df = pd.DataFrame(sample_stats)

    pairs = []
    for i, ct1 in enumerate(valid_types):
        for j, ct2 in enumerate(valid_types):
            pairs.append({
                "cell_type_1": ct1,
                "cell_type_2": ct2,
                "z_score": round(mean_z_scores[i, j], 4),
                "p_value_raw": round(mean_p_values[i, j], 6),
                "p_value_adj": round(adj_p_values[i, j], 6),
                "significant": adj_p_values[i, j] < alpha,
            })

    enrichment_pairs = pd.DataFrame(pairs)
    enrichment_pairs = enrichment_pairs.sort_values("z_score", ascending=False)

    return NeighborhoodEnrichmentResult(
        z_scores=z_scores_df,
        p_values_raw=p_values_raw_df,
        p_values_adj=p_values_adj_df,
        enrichment_pairs=enrichment_pairs,
        sample_stats=sample_stats_df,
        cell_types=valid_types,
        n_permutations=n_permutations,
        correction_method=correction_method,
        alpha=alpha,
    )
