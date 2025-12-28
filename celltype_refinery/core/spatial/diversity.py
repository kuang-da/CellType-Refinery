"""Local diversity analysis.

Computes Shannon diversity index for each cell's neighborhood,
measuring the variety of cell types in the local spatial context.
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
from scipy import sparse


@dataclass
class LocalDiversityResult:
    """Result from local diversity analysis.

    Attributes
    ----------
    diversity_by_sample : pd.DataFrame
        Per-sample summary statistics
    diversity_by_region : pd.DataFrame
        Aggregated by region
    mean_diversity : float
        Global mean local diversity
    std_diversity : float
        Global std of local diversity
    n_cells : int
        Total cells analyzed
    """

    diversity_by_sample: pd.DataFrame = field(default_factory=pd.DataFrame)
    diversity_by_region: pd.DataFrame = field(default_factory=pd.DataFrame)
    mean_diversity: float = 0.0
    std_diversity: float = 0.0
    n_cells: int = 0

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            "mean_diversity": round(self.mean_diversity, 4),
            "std_diversity": round(self.std_diversity, 4),
            "n_cells": self.n_cells,
            "n_samples": len(self.diversity_by_sample)
            if not self.diversity_by_sample.empty
            else 0,
        }


def compute_local_diversity(
    graph: sparse.csr_matrix,
    cell_types: np.ndarray,
    unique_types: List[str],
    include_self: bool = True,
) -> np.ndarray:
    """Compute local cell-type diversity for each cell.

    For each cell, computes Shannon entropy of its neighborhood:
        H = -sum(p_i * log(p_i))

    Parameters
    ----------
    graph : sparse.csr_matrix
        Sparse adjacency matrix
    cell_types : np.ndarray
        Cell type labels for each cell
    unique_types : List[str]
        List of unique cell types
    include_self : bool
        Whether to include the cell itself in its neighborhood

    Returns
    -------
    np.ndarray
        Local diversity (Shannon entropy) for each cell
    """
    n_cells = len(cell_types)
    n_types = len(unique_types)
    type_to_idx = {t: i for i, t in enumerate(unique_types)}

    local_diversity = np.zeros(n_cells)

    for i in range(n_cells):
        start, end = graph.indptr[i], graph.indptr[i + 1]
        neighbor_idx = graph.indices[start:end]

        if len(neighbor_idx) == 0 and not include_self:
            local_diversity[i] = 0
            continue

        type_counts = np.zeros(n_types)

        for j in neighbor_idx:
            idx = type_to_idx.get(cell_types[j], -1)
            if idx >= 0:
                type_counts[idx] += 1

        if include_self:
            self_idx = type_to_idx.get(cell_types[i], -1)
            if self_idx >= 0:
                type_counts[self_idx] += 1

        total = type_counts.sum()
        if total > 0:
            probs = type_counts / total
            probs = probs[probs > 0]
            local_diversity[i] = -np.sum(probs * np.log(probs))

    return local_diversity


def compute_diversity_for_sample(
    sample_id: str,
    graph: sparse.csr_matrix,
    cell_types: np.ndarray,
    unique_types: List[str],
    region: Optional[str] = None,
) -> Dict[str, Any]:
    """Compute local diversity statistics for a sample.

    Parameters
    ----------
    sample_id : str
        Sample identifier
    graph : sparse.csr_matrix
        Sparse adjacency matrix
    cell_types : np.ndarray
        Cell type labels
    unique_types : List[str]
        List of unique cell types
    region : str, optional
        Region label for this sample

    Returns
    -------
    Dict
        Summary statistics for the sample
    """
    local_div = compute_local_diversity(graph, cell_types, unique_types)

    result = {
        "sample_id": sample_id,
        "local_diversity_mean": float(np.mean(local_div)),
        "local_diversity_std": float(np.std(local_div)),
        "local_diversity_median": float(np.median(local_div)),
        "local_diversity_min": float(np.min(local_div)),
        "local_diversity_max": float(np.max(local_div)),
        "n_cells": len(cell_types),
    }

    if region is not None:
        result["region"] = region

    return result


def aggregate_diversity_results(
    sample_results: List[Dict[str, Any]],
) -> LocalDiversityResult:
    """Aggregate local diversity results across samples.

    Parameters
    ----------
    sample_results : List[Dict]
        Results from compute_diversity_for_sample()

    Returns
    -------
    LocalDiversityResult
        Aggregated results
    """
    if not sample_results:
        return LocalDiversityResult()

    diversity_by_sample = pd.DataFrame(sample_results)

    total_cells = diversity_by_sample["n_cells"].sum()
    weighted_mean = (
        diversity_by_sample["local_diversity_mean"] * diversity_by_sample["n_cells"]
    ).sum() / total_cells

    diversity_by_region = pd.DataFrame()
    if "region" in diversity_by_sample.columns:
        diversity_by_region = (
            diversity_by_sample.groupby("region", observed=True)
            .agg(
                mean=("local_diversity_mean", "mean"),
                std=("local_diversity_mean", "std"),
                n_samples=("sample_id", "count"),
            )
            .reset_index()
        )

    return LocalDiversityResult(
        diversity_by_sample=diversity_by_sample,
        diversity_by_region=diversity_by_region,
        mean_diversity=weighted_mean,
        std_diversity=float(diversity_by_sample["local_diversity_std"].mean()),
        n_cells=int(total_cells),
    )
