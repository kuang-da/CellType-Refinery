"""Aggregation functions for cell-type composition.

This module provides functions to aggregate cell-type counts and proportions
at different levels: sample, region, donor, and global.
"""

from typing import List, Optional

import numpy as np
import pandas as pd


def compute_composition_by_group(
    df: pd.DataFrame,
    group_col: str,
    cell_type_col: str = "cell_type",
    normalize: bool = True,
) -> pd.DataFrame:
    """Compute cell-type composition for each group.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with cell data
    group_col : str
        Column to group by (e.g., "sample_id", "region")
    cell_type_col : str
        Column with cell type labels
    normalize : bool
        If True, return proportions; if False, return counts

    Returns
    -------
    pd.DataFrame
        Composition table with columns:
        - group_col: Group identifier
        - cell_type: Cell type label
        - count: Number of cells
        - proportion: Proportion of cells (if normalize=True)
    """
    # Count cells per group and cell type
    counts = df.groupby([group_col, cell_type_col], observed=True).size().reset_index(name="count")

    # Rename cell type column to canonical "cell_type" for downstream consistency
    if cell_type_col != "cell_type":
        counts = counts.rename(columns={cell_type_col: "cell_type"})

    if normalize:
        # Compute total per group
        totals = counts.groupby(group_col, observed=True)["count"].transform("sum")
        counts["proportion"] = counts["count"] / totals

    return counts


def compute_composition_by_sample(
    adata,
    cell_type_col: str = "cell_type",
    sample_col: str = "sample_id",
    region_col: Optional[str] = "region",
    donor_col: Optional[str] = "donor",
) -> pd.DataFrame:
    """Compute cell-type composition per sample.

    Parameters
    ----------
    adata : AnnData
        AnnData object with cell annotations
    cell_type_col : str
        Column with cell type labels
    sample_col : str
        Column with sample IDs
    region_col : str, optional
        Column with region labels (for metadata)
    donor_col : str, optional
        Column with donor IDs (for metadata)

    Returns
    -------
    pd.DataFrame
        Per-sample composition with count and proportion columns
    """
    df = adata.obs[[sample_col, cell_type_col]].copy()

    # Add optional metadata
    if region_col and region_col in adata.obs.columns:
        df[region_col] = adata.obs[region_col]
    if donor_col and donor_col in adata.obs.columns:
        df[donor_col] = adata.obs[donor_col]

    composition = compute_composition_by_group(df, sample_col, cell_type_col)

    # Add metadata columns
    if region_col and region_col in df.columns:
        sample_meta = df.groupby(sample_col, observed=True)[region_col].first()
        composition[region_col] = composition[sample_col].map(sample_meta)

    if donor_col and donor_col in df.columns:
        sample_meta = df.groupby(sample_col, observed=True)[donor_col].first()
        composition[donor_col] = composition[sample_col].map(sample_meta)

    return composition


def aggregate_by_region(
    composition_df: pd.DataFrame,
    region_col: str = "region",
    region_order: Optional[List[str]] = None,
) -> pd.DataFrame:
    """Aggregate composition statistics by region.

    Parameters
    ----------
    composition_df : pd.DataFrame
        Per-sample composition from compute_composition_by_sample
    region_col : str
        Column with region labels
    region_order : List[str], optional
        Ordered list of regions for consistent output

    Returns
    -------
    pd.DataFrame
        Regional composition with mean_pct, std_pct, median_pct, min_pct, max_pct, n_samples
    """
    if region_col not in composition_df.columns:
        raise ValueError(f"Column '{region_col}' not found in composition DataFrame")

    # Check if we have 'percentage' column (0-100 scale) or 'proportion' column (0-1 scale)
    pct_col = "percentage" if "percentage" in composition_df.columns else "proportion"
    scale_factor = 1.0 if pct_col == "percentage" else 100.0

    agg = composition_df.groupby([region_col, "cell_type"], observed=True)[pct_col].agg([
        ("mean_pct", "mean"),
        ("std_pct", "std"),
        ("median_pct", "median"),
        ("min_pct", "min"),
        ("max_pct", "max"),
        ("n_samples", "count"),
    ]).reset_index()

    # Scale to percentage if needed
    if scale_factor != 1.0:
        for col in ["mean_pct", "std_pct", "median_pct", "min_pct", "max_pct"]:
            agg[col] = agg[col] * scale_factor

    # Order regions if specified
    if region_order:
        agg[region_col] = pd.Categorical(
            agg[region_col],
            categories=region_order,
            ordered=True
        )
        agg = agg.sort_values([region_col, "cell_type"])

    return agg


def aggregate_by_donor(
    composition_df: pd.DataFrame,
    donor_col: str = "donor",
) -> pd.DataFrame:
    """Aggregate composition statistics by donor.

    Parameters
    ----------
    composition_df : pd.DataFrame
        Per-sample composition from compute_composition_by_sample
    donor_col : str
        Column with donor IDs

    Returns
    -------
    pd.DataFrame
        Per-donor composition with mean, std across samples
    """
    if donor_col not in composition_df.columns:
        raise ValueError(f"Column '{donor_col}' not found in composition DataFrame")

    agg = composition_df.groupby([donor_col, "cell_type"], observed=True)["proportion"].agg([
        ("mean", "mean"),
        ("std", "std"),
        ("n_samples", "count"),
    ]).reset_index()

    return agg


def compute_global_summary(
    composition_df: pd.DataFrame,
    min_presence: int = 1,
) -> pd.DataFrame:
    """Compute global cell-type summary.

    Parameters
    ----------
    composition_df : pd.DataFrame
        Per-sample composition
    min_presence : int
        Minimum number of samples a cell type must appear in

    Returns
    -------
    pd.DataFrame
        Global summary with columns:
        - cell_type
        - total_count
        - n_samples_present
        - mean_pct (percentage)
        - std_pct (percentage)
        - median_pct (percentage)
        - min_pct (percentage)
        - max_pct (percentage)
        - global_pct (percentage)
        - cv (coefficient of variation)
    """
    # Check if we have 'percentage' column (0-100 scale) or 'proportion' column (0-1 scale)
    pct_col = "percentage" if "percentage" in composition_df.columns else "proportion"
    scale_factor = 1.0 if pct_col == "percentage" else 100.0

    # Total counts and statistics
    totals = composition_df.groupby("cell_type", observed=True).agg(
        total_count=("count", "sum"),
        n_samples_present=(pct_col, lambda x: (x > 0).sum()),
        mean_pct=(pct_col, "mean"),
        std_pct=(pct_col, "std"),
        median_pct=(pct_col, "median"),
        min_pct=(pct_col, "min"),
        max_pct=(pct_col, "max"),
    ).reset_index()

    # Scale to percentage if needed
    if scale_factor != 1.0:
        for col in ["mean_pct", "std_pct", "median_pct", "min_pct", "max_pct"]:
            totals[col] = totals[col] * scale_factor

    # Compute global percentage
    grand_total = totals["total_count"].sum()
    totals["global_pct"] = 100.0 * totals["total_count"] / grand_total

    # Compute CV (coefficient of variation)
    totals["cv"] = totals["std_pct"] / totals["mean_pct"].replace(0, np.nan)
    totals["cv"] = totals["cv"].replace([np.inf, -np.inf], np.nan)

    # Filter by minimum presence
    totals = totals[totals["n_samples_present"] >= min_presence]

    # Sort by total count (most abundant first)
    totals = totals.sort_values("total_count", ascending=False)

    return totals


def create_composition_wide(
    composition_df: pd.DataFrame,
    index_col: str = "sample_id",
    value_col: str = "proportion",
    fill_value: float = 0.0,
) -> pd.DataFrame:
    """Create wide-format composition matrix.

    Parameters
    ----------
    composition_df : pd.DataFrame
        Long-format composition from compute_composition_by_sample
    index_col : str
        Column to use as index (rows)
    value_col : str
        Column to use as values (proportion or count)
    fill_value : float
        Value for missing cell types

    Returns
    -------
    pd.DataFrame
        Wide-format matrix (samples x cell types)
    """
    wide = composition_df.pivot_table(
        index=index_col,
        columns="cell_type",
        values=value_col,
        fill_value=fill_value,
        aggfunc="first",
        observed=True,
    )

    return wide


def compute_rare_cell_types(
    global_summary: pd.DataFrame,
    threshold: float = 1.0,
) -> pd.DataFrame:
    """Identify rare cell types.

    Parameters
    ----------
    global_summary : pd.DataFrame
        Output from compute_global_summary
    threshold : float
        Percentage threshold for rare classification (default: 1.0 = 1%)

    Returns
    -------
    pd.DataFrame
        Rare cell types with their statistics
    """
    # Support both old (total_proportion) and new (global_pct) column names
    pct_col = "global_pct" if "global_pct" in global_summary.columns else "total_proportion"
    # If using old column (proportion 0-1), convert threshold to match
    if pct_col == "total_proportion":
        threshold = threshold / 100.0
    rare = global_summary[global_summary[pct_col] < threshold].copy()
    rare["is_rare"] = True
    return rare


def compute_cell_type_co_occurrence(
    composition_wide: pd.DataFrame,
    method: str = "pearson",
) -> pd.DataFrame:
    """Compute co-occurrence correlation between cell types.

    Parameters
    ----------
    composition_wide : pd.DataFrame
        Wide-format composition matrix (samples x cell types)
    method : str
        Correlation method: "pearson", "spearman", or "kendall"

    Returns
    -------
    pd.DataFrame
        Correlation matrix (cell types x cell types)
    """
    return composition_wide.corr(method=method)


def compute_composition_differences(
    composition_df: pd.DataFrame,
    group_col: str,
    group1: str,
    group2: str,
) -> pd.DataFrame:
    """Compute differences in composition between two groups.

    Parameters
    ----------
    composition_df : pd.DataFrame
        Per-sample composition
    group_col : str
        Column to compare groups on (e.g., "region", "donor")
    group1 : str
        First group identifier
    group2 : str
        Second group identifier

    Returns
    -------
    pd.DataFrame
        Differences with columns: cell_type, mean_1, mean_2, diff, abs_diff
    """
    g1 = composition_df[composition_df[group_col] == group1]
    g2 = composition_df[composition_df[group_col] == group2]

    mean1 = g1.groupby("cell_type", observed=True)["proportion"].mean()
    mean2 = g2.groupby("cell_type", observed=True)["proportion"].mean()

    diff_df = pd.DataFrame({
        "cell_type": mean1.index.union(mean2.index),
    })
    diff_df[f"mean_{group1}"] = diff_df["cell_type"].map(mean1).fillna(0)
    diff_df[f"mean_{group2}"] = diff_df["cell_type"].map(mean2).fillna(0)
    diff_df["diff"] = diff_df[f"mean_{group1}"] - diff_df[f"mean_{group2}"]
    diff_df["abs_diff"] = diff_df["diff"].abs()

    return diff_df.sort_values("abs_diff", ascending=False)
