"""Diversity metrics for cell-type composition.

This module computes ecological diversity indices for cell-type compositions:
- Shannon entropy: Measures overall diversity
- Simpson index: Probability that two random cells are different types
- Pielou's evenness: How evenly cells are distributed across types
"""

from typing import Optional

import numpy as np
import pandas as pd


def compute_shannon_entropy(counts: np.ndarray) -> float:
    """Compute Shannon entropy for a distribution.

    H = -sum(p_i * log(p_i))

    Parameters
    ----------
    counts : np.ndarray
        Array of counts (non-negative)

    Returns
    -------
    float
        Shannon entropy in nats (using natural log)
    """
    counts = np.asarray(counts, dtype=float)
    total = counts.sum()
    if total == 0:
        return 0.0

    proportions = counts / total
    # Filter out zeros to avoid log(0)
    proportions = proportions[proportions > 0]

    return -np.sum(proportions * np.log(proportions))


def compute_simpson_index(counts: np.ndarray) -> float:
    """Compute Simpson diversity index.

    D = 1 - sum(p_i^2)

    This gives the probability that two randomly selected cells
    are of different types.

    Parameters
    ----------
    counts : np.ndarray
        Array of counts (non-negative)

    Returns
    -------
    float
        Simpson diversity index (0 to 1)
    """
    counts = np.asarray(counts, dtype=float)
    total = counts.sum()
    if total == 0:
        return 0.0

    proportions = counts / total
    return 1.0 - np.sum(proportions ** 2)


def compute_evenness(counts: np.ndarray) -> float:
    """Compute Pielou's evenness index.

    J = H / H_max = H / log(S)

    Where S is the number of types with non-zero counts.

    Parameters
    ----------
    counts : np.ndarray
        Array of counts (non-negative)

    Returns
    -------
    float
        Pielou's evenness (0 to 1)
    """
    counts = np.asarray(counts, dtype=float)
    n_types = np.sum(counts > 0)

    if n_types <= 1:
        return 1.0  # Perfect evenness with 0 or 1 type

    h = compute_shannon_entropy(counts)
    h_max = np.log(n_types)

    return h / h_max if h_max > 0 else 1.0


def compute_diversity_by_group(
    df: pd.DataFrame,
    group_col: str,
    cell_type_col: str = "cell_type",
    min_cells: int = 100,
) -> pd.DataFrame:
    """Compute diversity metrics for each group.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with cell data
    group_col : str
        Column to group by (e.g., "sample_id", "region")
    cell_type_col : str
        Column with cell type labels
    min_cells : int
        Minimum cells per group for valid calculation

    Returns
    -------
    pd.DataFrame
        DataFrame with columns:
        - group_col: Group identifier
        - n_cells: Total cells in group
        - n_types: Number of unique cell types
        - shannon_entropy: Shannon entropy (nats)
        - simpson_index: Simpson diversity index
        - evenness: Pielou's evenness
    """
    results = []

    for group, group_df in df.groupby(group_col):
        n_cells = len(group_df)

        if n_cells < min_cells:
            results.append({
                group_col: group,
                "n_cells": n_cells,
                "n_types": np.nan,
                "shannon_entropy": np.nan,
                "simpson_index": np.nan,
                "evenness": np.nan,
            })
            continue

        counts = group_df[cell_type_col].value_counts().values

        results.append({
            group_col: group,
            "n_cells": n_cells,
            "n_types": len(counts),
            "shannon_entropy": compute_shannon_entropy(counts),
            "simpson_index": compute_simpson_index(counts),
            "evenness": compute_evenness(counts),
        })

    return pd.DataFrame(results)


def compute_diversity_summary(
    diversity_df: pd.DataFrame,
    group_col: Optional[str] = None,
) -> pd.DataFrame:
    """Compute summary statistics for diversity metrics.

    Parameters
    ----------
    diversity_df : pd.DataFrame
        Output from compute_diversity_by_group
    group_col : str, optional
        If provided, compute summary within each group level

    Returns
    -------
    pd.DataFrame
        Summary statistics for each diversity metric
    """
    metrics = ["shannon_entropy", "simpson_index", "evenness", "n_types"]

    if group_col is not None and group_col in diversity_df.columns:
        # Group-level summary
        summaries = []
        for group, gdf in diversity_df.groupby(group_col):
            for metric in metrics:
                values = gdf[metric].dropna()
                if len(values) > 0:
                    summaries.append({
                        group_col: group,
                        "metric": metric,
                        "mean": values.mean(),
                        "std": values.std(),
                        "median": values.median(),
                        "min": values.min(),
                        "max": values.max(),
                        "n": len(values),
                    })
        return pd.DataFrame(summaries)
    else:
        # Global summary
        summaries = []
        for metric in metrics:
            values = diversity_df[metric].dropna()
            if len(values) > 0:
                summaries.append({
                    "metric": metric,
                    "mean": values.mean(),
                    "std": values.std(),
                    "median": values.median(),
                    "min": values.min(),
                    "max": values.max(),
                    "n": len(values),
                })
        return pd.DataFrame(summaries)


def compute_diversity_from_composition(
    composition_df: pd.DataFrame,
    proportion_cols: Optional[list] = None,
) -> pd.DataFrame:
    """Compute diversity from a composition table (wide format).

    Parameters
    ----------
    composition_df : pd.DataFrame
        Composition table with rows as samples and columns as cell types
    proportion_cols : list, optional
        Columns to use as proportions. If None, uses all numeric columns.

    Returns
    -------
    pd.DataFrame
        Diversity metrics for each row (sample)
    """
    if proportion_cols is None:
        proportion_cols = composition_df.select_dtypes(include=[np.number]).columns.tolist()

    results = []
    for idx, row in composition_df.iterrows():
        proportions = row[proportion_cols].values.astype(float)
        # Convert proportions to pseudo-counts
        counts = proportions * 1000  # Scale to avoid numerical issues

        results.append({
            "index": idx,
            "shannon_entropy": compute_shannon_entropy(counts),
            "simpson_index": compute_simpson_index(counts),
            "evenness": compute_evenness(counts),
            "n_types": np.sum(proportions > 0),
        })

    return pd.DataFrame(results).set_index("index")
