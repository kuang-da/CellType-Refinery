"""Moran's I spatial autocorrelation analysis.

Computes global spatial autocorrelation for each cell type,
measuring whether cells of the same type cluster spatially (+1),
disperse (-1), or are randomly distributed (~0).
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List

import numpy as np
import pandas as pd
from scipy import sparse


@dataclass
class MoransResult:
    """Result from Moran's I analysis.

    Attributes
    ----------
    morans_by_sample : pd.DataFrame
        Per-sample, per-cell-type Moran's I values
    morans_global : pd.DataFrame
        Aggregated by cell type across samples
    cell_types : List[str]
        Cell types analyzed
    n_samples : int
        Number of samples analyzed
    """

    morans_by_sample: pd.DataFrame = field(default_factory=pd.DataFrame)
    morans_global: pd.DataFrame = field(default_factory=pd.DataFrame)
    cell_types: List[str] = field(default_factory=list)
    n_samples: int = 0

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            "n_cell_types": len(self.cell_types),
            "n_samples": self.n_samples,
            "mean_morans_i": float(self.morans_global["morans_i_mean"].mean())
            if not self.morans_global.empty
            else None,
        }


def compute_morans_i(
    graph: sparse.csr_matrix,
    values: np.ndarray,
) -> float:
    """Compute Moran's I spatial autocorrelation.

    Moran's I measures the global spatial autocorrelation of a variable.
    Values close to +1 indicate clustering (similar values near each other),
    close to -1 indicate dispersion (dissimilar values near each other),
    and close to 0 indicate random distribution.

    Formula:
        I = (n/W) * (sum_{ij} w_{ij} * z_i * z_j) / (sum_i z_i^2)

    Parameters
    ----------
    graph : sparse.csr_matrix
        Sparse adjacency matrix
    values : np.ndarray
        Variable values for each cell

    Returns
    -------
    float
        Moran's I statistic, or np.nan if computation fails
    """
    n = len(values)
    if n == 0:
        return np.nan

    mean_val = np.mean(values)
    z = values - mean_val

    W = graph.sum()
    if W == 0:
        return np.nan

    numerator = z @ graph @ z

    denominator = np.sum(z**2)
    if denominator == 0:
        return np.nan

    morans_i = (n / W) * (numerator / denominator)

    return float(morans_i)


def compute_morans_i_for_type(
    graph: sparse.csr_matrix,
    cell_types: np.ndarray,
    target_type: str,
) -> float:
    """Compute Moran's I for a specific cell type.

    Creates a binary indicator (1 for target type, 0 otherwise)
    and computes spatial autocorrelation.

    Parameters
    ----------
    graph : sparse.csr_matrix
        Sparse adjacency matrix
    cell_types : np.ndarray
        Cell type labels for each cell
    target_type : str
        Cell type to analyze

    Returns
    -------
    float
        Moran's I for the target cell type
    """
    indicator = (cell_types == target_type).astype(float)
    return compute_morans_i(graph, indicator)


def compute_morans_for_sample(
    sample_id: str,
    graph: sparse.csr_matrix,
    cell_types: np.ndarray,
    valid_types: List[str],
    min_cells_per_type: int = 50,
) -> List[Dict[str, Any]]:
    """Compute Moran's I for all cell types in a sample.

    Parameters
    ----------
    sample_id : str
        Sample identifier
    graph : sparse.csr_matrix
        Sparse adjacency matrix
    cell_types : np.ndarray
        Cell type labels
    valid_types : List[str]
        Cell types to analyze
    min_cells_per_type : int
        Minimum cells required for a cell type

    Returns
    -------
    List[Dict]
        List of results per cell type
    """
    results = []

    for ct in valid_types:
        n_cells = np.sum(cell_types == ct)

        if n_cells < min_cells_per_type:
            continue

        morans_i = compute_morans_i_for_type(graph, cell_types, ct)

        results.append({
            "sample_id": sample_id,
            "cell_type": ct,
            "morans_i": morans_i,
            "n_cells": int(n_cells),
        })

    return results


def aggregate_morans_results(
    sample_results: List[List[Dict[str, Any]]],
) -> MoransResult:
    """Aggregate Moran's I results across samples.

    Parameters
    ----------
    sample_results : List[List[Dict]]
        Results from compute_morans_for_sample() for each sample

    Returns
    -------
    MoransResult
        Aggregated results
    """
    all_results = []
    for sample_result in sample_results:
        all_results.extend(sample_result)

    if not all_results:
        return MoransResult()

    morans_by_sample = pd.DataFrame(all_results)

    morans_global = (
        morans_by_sample.groupby("cell_type", observed=True)
        .agg(
            morans_i_mean=("morans_i", "mean"),
            morans_i_std=("morans_i", "std"),
            morans_i_median=("morans_i", "median"),
            n_samples=("sample_id", "nunique"),
            total_cells=("n_cells", "sum"),
        )
        .reset_index()
    )

    morans_global = morans_global.sort_values("morans_i_mean", ascending=False)

    return MoransResult(
        morans_by_sample=morans_by_sample,
        morans_global=morans_global,
        cell_types=morans_global["cell_type"].tolist(),
        n_samples=morans_by_sample["sample_id"].nunique(),
    )
