"""Regional enrichment analysis.

This module provides statistical tests to identify cell types that are
significantly enriched or depleted in specific regions.
"""

from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats


def compute_regional_enrichment(
    composition_df: pd.DataFrame,
    region_col: str = "region",
    cell_type_col: str = "cell_type",
    value_col: str = "proportion",
    sample_col: str = "sample_id",
    min_samples_per_region: int = 3,
    correction_method: str = "fdr_bh",
    alpha: float = 0.05,
) -> pd.DataFrame:
    """Compute regional enrichment using Mann-Whitney U test.

    For each cell type in each region, tests whether the proportion
    is significantly different from all other regions combined.

    Parameters
    ----------
    composition_df : pd.DataFrame
        Per-sample composition with region information
    region_col : str
        Column with region labels
    cell_type_col : str
        Column with cell type labels
    value_col : str
        Column with values to test (proportion or count)
    sample_col : str
        Column with sample IDs
    min_samples_per_region : int
        Minimum samples per region for valid test
    correction_method : str
        Multiple testing correction: "bonferroni", "fdr_bh", "holm", or "none"
    alpha : float
        Significance threshold

    Returns
    -------
    pd.DataFrame
        Enrichment results with columns:
        - region
        - cell_type
        - mean_in_region
        - mean_out_region
        - fold_change
        - log2_fold_change
        - statistic
        - p_value
        - p_adjusted
        - significant
        - direction (enriched/depleted)
    """
    results = []
    regions = composition_df[region_col].unique()
    cell_types = composition_df[cell_type_col].unique()

    for region in regions:
        # Get samples in this region
        in_region = composition_df[composition_df[region_col] == region]
        out_region = composition_df[composition_df[region_col] != region]

        n_samples_in = in_region[sample_col].nunique()
        n_samples_out = out_region[sample_col].nunique()

        if n_samples_in < min_samples_per_region or n_samples_out < min_samples_per_region:
            continue

        for cell_type in cell_types:
            # Get values for this cell type
            values_in = in_region[in_region[cell_type_col] == cell_type][value_col].values
            values_out = out_region[out_region[cell_type_col] == cell_type][value_col].values

            # Pad with zeros for missing samples
            n_in_expected = n_samples_in
            n_out_expected = n_samples_out

            if len(values_in) < n_in_expected:
                values_in = np.concatenate([values_in, np.zeros(n_in_expected - len(values_in))])
            if len(values_out) < n_out_expected:
                values_out = np.concatenate([values_out, np.zeros(n_out_expected - len(values_out))])

            mean_in = np.mean(values_in) if len(values_in) > 0 else 0
            mean_out = np.mean(values_out) if len(values_out) > 0 else 0

            # Compute fold change
            if mean_out > 0:
                fold_change = mean_in / mean_out
            elif mean_in > 0:
                fold_change = np.inf
            else:
                fold_change = 1.0

            log2_fc = np.log2(fold_change) if fold_change > 0 and np.isfinite(fold_change) else 0

            # Mann-Whitney U test
            if len(values_in) > 0 and len(values_out) > 0:
                try:
                    stat, p_value = stats.mannwhitneyu(
                        values_in, values_out,
                        alternative="two-sided",
                    )
                except ValueError:
                    stat, p_value = np.nan, 1.0
            else:
                stat, p_value = np.nan, 1.0

            results.append({
                "region": region,
                "cell_type": cell_type,
                "mean_in_region": mean_in,
                "mean_out_region": mean_out,
                "fold_change": fold_change,
                "log2_fold_change": log2_fc,
                "statistic": stat,
                "p_value": p_value,
                "n_samples_in": n_samples_in,
                "n_samples_out": n_samples_out,
            })

    if not results:
        return pd.DataFrame(columns=[
            "region", "cell_type", "mean_in_region", "mean_out_region",
            "fold_change", "log2_fold_change", "statistic", "p_value",
            "p_adjusted", "significant", "direction",
        ])

    df = pd.DataFrame(results)

    # Apply multiple testing correction
    df["p_adjusted"] = _apply_multiple_testing_correction(
        df["p_value"].values,
        method=correction_method,
    )

    df["significant"] = df["p_adjusted"] < alpha
    df["direction"] = np.where(
        df["fold_change"] > 1,
        "enriched",
        np.where(df["fold_change"] < 1, "depleted", "neutral"),
    )

    return df.sort_values(["region", "p_adjusted"])


def _apply_multiple_testing_correction(
    p_values: np.ndarray,
    method: str = "fdr_bh",
) -> np.ndarray:
    """Apply multiple testing correction.

    Parameters
    ----------
    p_values : np.ndarray
        Raw p-values
    method : str
        Correction method: "bonferroni", "fdr_bh", "holm", or "none"

    Returns
    -------
    np.ndarray
        Adjusted p-values
    """
    p_values = np.asarray(p_values, dtype=float)
    n = len(p_values)

    if n == 0:
        return p_values

    # Handle NaN values
    valid_mask = ~np.isnan(p_values)
    valid_p = p_values[valid_mask]

    if len(valid_p) == 0:
        return p_values

    if method == "none":
        return p_values

    elif method == "bonferroni":
        adjusted_valid = np.minimum(valid_p * n, 1.0)

    elif method == "fdr_bh":
        # Benjamini-Hochberg
        n_valid = len(valid_p)
        sorted_idx = np.argsort(valid_p)
        sorted_p = valid_p[sorted_idx]

        adjusted_sorted = np.zeros(n_valid)
        adjusted_sorted[-1] = sorted_p[-1]

        for i in range(n_valid - 2, -1, -1):
            adjusted_sorted[i] = min(
                adjusted_sorted[i + 1],
                sorted_p[i] * n_valid / (i + 1),
            )

        adjusted_valid = np.zeros(n_valid)
        adjusted_valid[sorted_idx] = adjusted_sorted
        adjusted_valid = np.minimum(adjusted_valid, 1.0)

    elif method == "holm":
        # Holm-Bonferroni
        n_valid = len(valid_p)
        sorted_idx = np.argsort(valid_p)
        sorted_p = valid_p[sorted_idx]

        adjusted_sorted = sorted_p * (n_valid - np.arange(n_valid))
        adjusted_sorted = np.minimum.accumulate(adjusted_sorted[::-1])[::-1]

        adjusted_valid = np.zeros(n_valid)
        adjusted_valid[sorted_idx] = adjusted_sorted
        adjusted_valid = np.minimum(adjusted_valid, 1.0)

    else:
        raise ValueError(f"Unknown correction method: {method}")

    # Reconstruct full array with NaNs
    adjusted = np.full(n, np.nan)
    adjusted[valid_mask] = adjusted_valid

    return adjusted


def get_significant_enrichments(
    enrichment_df: pd.DataFrame,
    min_fold_change: float = 1.5,
) -> pd.DataFrame:
    """Filter to significant enrichments.

    Parameters
    ----------
    enrichment_df : pd.DataFrame
        Output from compute_regional_enrichment
    min_fold_change : float
        Minimum fold change to include

    Returns
    -------
    pd.DataFrame
        Significant enrichments only
    """
    mask = (
        enrichment_df["significant"] &
        (enrichment_df["fold_change"] >= min_fold_change)
    )
    return enrichment_df[mask].copy()


def get_significant_depletions(
    enrichment_df: pd.DataFrame,
    max_fold_change: float = 0.67,
) -> pd.DataFrame:
    """Filter to significant depletions.

    Parameters
    ----------
    enrichment_df : pd.DataFrame
        Output from compute_regional_enrichment
    max_fold_change : float
        Maximum fold change to include

    Returns
    -------
    pd.DataFrame
        Significant depletions only
    """
    mask = (
        enrichment_df["significant"] &
        (enrichment_df["fold_change"] <= max_fold_change)
    )
    return enrichment_df[mask].copy()


def summarize_enrichments_by_region(
    enrichment_df: pd.DataFrame,
) -> pd.DataFrame:
    """Summarize enrichment results by region.

    Parameters
    ----------
    enrichment_df : pd.DataFrame
        Output from compute_regional_enrichment

    Returns
    -------
    pd.DataFrame
        Summary with n_enriched, n_depleted per region
    """
    summary = enrichment_df.groupby("region").agg(
        n_tests=("cell_type", "count"),
        n_significant=("significant", "sum"),
        n_enriched=("direction", lambda x: (x == "enriched").sum()),
        n_depleted=("direction", lambda x: (x == "depleted").sum()),
        mean_abs_log2fc=("log2_fold_change", lambda x: np.abs(x).mean()),
    ).reset_index()

    return summary


def summarize_enrichments_by_cell_type(
    enrichment_df: pd.DataFrame,
) -> pd.DataFrame:
    """Summarize enrichment results by cell type.

    Parameters
    ----------
    enrichment_df : pd.DataFrame
        Output from compute_regional_enrichment

    Returns
    -------
    pd.DataFrame
        Summary with regions where each cell type is enriched/depleted
    """
    summary = []

    for cell_type, ct_df in enrichment_df.groupby("cell_type"):
        sig = ct_df[ct_df["significant"]]
        enriched_regions = sig[sig["direction"] == "enriched"]["region"].tolist()
        depleted_regions = sig[sig["direction"] == "depleted"]["region"].tolist()

        summary.append({
            "cell_type": cell_type,
            "n_regions_enriched": len(enriched_regions),
            "n_regions_depleted": len(depleted_regions),
            "enriched_in": ", ".join(enriched_regions) if enriched_regions else "",
            "depleted_in": ", ".join(depleted_regions) if depleted_regions else "",
            "max_fold_change": ct_df["fold_change"].max(),
            "min_fold_change": ct_df["fold_change"].min(),
        })

    return pd.DataFrame(summary).sort_values("n_regions_enriched", ascending=False)


def compute_pairwise_enrichment(
    composition_df: pd.DataFrame,
    region1: str,
    region2: str,
    region_col: str = "region",
    cell_type_col: str = "cell_type",
    value_col: str = "proportion",
    sample_col: str = "sample_id",
) -> pd.DataFrame:
    """Compute enrichment between two specific regions.

    Parameters
    ----------
    composition_df : pd.DataFrame
        Per-sample composition
    region1 : str
        First region
    region2 : str
        Second region
    region_col : str
        Column with region labels
    cell_type_col : str
        Column with cell type labels
    value_col : str
        Column with values
    sample_col : str
        Column with sample IDs

    Returns
    -------
    pd.DataFrame
        Pairwise enrichment results
    """
    r1_df = composition_df[composition_df[region_col] == region1]
    r2_df = composition_df[composition_df[region_col] == region2]

    results = []
    cell_types = composition_df[cell_type_col].unique()

    for cell_type in cell_types:
        v1 = r1_df[r1_df[cell_type_col] == cell_type][value_col].values
        v2 = r2_df[r2_df[cell_type_col] == cell_type][value_col].values

        mean1 = np.mean(v1) if len(v1) > 0 else 0
        mean2 = np.mean(v2) if len(v2) > 0 else 0

        if mean2 > 0:
            fold_change = mean1 / mean2
        elif mean1 > 0:
            fold_change = np.inf
        else:
            fold_change = 1.0

        # Mann-Whitney U test
        if len(v1) > 0 and len(v2) > 0:
            try:
                stat, p_value = stats.mannwhitneyu(v1, v2, alternative="two-sided")
            except ValueError:
                stat, p_value = np.nan, 1.0
        else:
            stat, p_value = np.nan, 1.0

        results.append({
            "cell_type": cell_type,
            f"mean_{region1}": mean1,
            f"mean_{region2}": mean2,
            "fold_change": fold_change,
            "log2_fold_change": np.log2(fold_change) if fold_change > 0 and np.isfinite(fold_change) else 0,
            "p_value": p_value,
        })

    df = pd.DataFrame(results)
    df["p_adjusted"] = _apply_multiple_testing_correction(df["p_value"].values, "fdr_bh")
    df["significant"] = df["p_adjusted"] < 0.05

    return df.sort_values("p_adjusted")
