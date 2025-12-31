"""
Cluster-level metrics for refinement diagnostics.

This module computes cluster-level metrics used by AutoPolicy to make
refinement decisions. These metrics complement the taxonomy-based scoring
by detecting heterogeneity that may not be visible through marker hierarchy.

Key metrics:
- marker_heterogeneity: Measures how "mixed" a cluster's marker profile is
- mixed_marker_rate: Fraction of markers with ambiguous positivity (20-80%)

Example:
    >>> from celltype_refinery.core.refinement import compute_cluster_marker_heterogeneity
    >>> metrics = compute_cluster_marker_heterogeneity(adata, cluster_col="cluster_lvl1")
    >>> print(metrics[["cluster_id", "marker_heterogeneity", "mixed_marker_rate"]])
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Optional, Union

import numpy as np
import pandas as pd

try:
    import scanpy as sc
except ImportError:
    sc = None


def compute_cluster_marker_heterogeneity(
    adata: "sc.AnnData",
    cluster_col: str = "cluster_lvl0",
    layer: str = "batchcorr",
    thresholds: Optional[Dict[str, float]] = None,
    marker_evidence_path: Optional[Path] = None,
    mixed_range: tuple = (0.2, 0.8),
    logger: Optional[logging.Logger] = None,
) -> pd.DataFrame:
    """
    Compute marker heterogeneity metrics per cluster.

    For each cluster, this function computes:
    1. p_m = fraction of cells positive for each marker m
    2. marker_heterogeneity = mean([p_m * (1 - p_m)]) across all markers
       - Range: 0 (all markers are decisive) to 0.25 (all markers at 50%)
    3. mixed_marker_rate = fraction of markers with p_m in [0.2, 0.8]
       - Range: 0 (all markers decisive) to 1 (all markers ambiguous)

    Parameters
    ----------
    adata : sc.AnnData
        AnnData object with expression data.
    cluster_col : str
        Column in adata.obs containing cluster assignments.
    layer : str
        Layer to use for expression values (default: "batchcorr").
    thresholds : Dict[str, float], optional
        Per-marker positivity thresholds. If not provided, will try to load
        from marker_evidence_path or use global 75th percentile.
    marker_evidence_path : Path, optional
        Path to marker_evidence.csv from Stage H. Used to extract thresholds
        if not provided explicitly.
    mixed_range : tuple
        Range (low, high) defining "mixed" positivity. Default: (0.2, 0.8).
    logger : logging.Logger, optional
        Logger for diagnostics.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns:
        - cluster_id: Cluster identifier (string)
        - n_cells: Number of cells in cluster
        - marker_heterogeneity: Mean of p*(1-p) across markers (0-0.25)
        - mixed_marker_rate: Fraction of markers in mixed range (0-1)
        - n_markers_mixed: Count of markers in mixed range
        - n_markers_total: Total markers used
        - heterogeneity_band: Categorical band (low/medium/high)

    Notes
    -----
    - Markers with all-zero expression or missing thresholds are excluded.
    - Higher heterogeneity suggests the cluster may contain mixed populations.
    - This metric is independent of the cell-type taxonomy.

    Examples
    --------
    >>> metrics = compute_cluster_marker_heterogeneity(
    ...     adata,
    ...     cluster_col="cluster_lvl1",
    ...     marker_evidence_path=Path("stage_h/marker_evidence.csv")
    ... )
    >>> # Find clusters with high heterogeneity
    >>> high_het = metrics[metrics["heterogeneity_band"] == "high"]
    """
    log = logger or logging.getLogger(__name__)

    # Validate inputs
    if cluster_col not in adata.obs:
        raise ValueError(f"Cluster column '{cluster_col}' not found in adata.obs")

    # Get expression matrix
    if layer and layer in adata.layers:
        X = adata.layers[layer]
        log.info("Using layer '%s' for heterogeneity computation", layer)
    else:
        X = adata.X
        log.info("Using adata.X for heterogeneity computation")

    # Convert to dense if sparse
    if hasattr(X, "toarray"):
        log.info("Converting sparse matrix to dense...")
        X = X.toarray()

    markers = list(adata.var_names)
    n_markers = len(markers)
    log.info("Computing heterogeneity for %d markers", n_markers)

    # Get or compute thresholds
    if thresholds is None:
        thresholds = _get_marker_thresholds(
            adata, markers, marker_evidence_path, logger=log
        )

    # Get cluster assignments
    clusters = adata.obs[cluster_col].astype(str).values
    unique_clusters = np.unique(clusters)
    n_clusters = len(unique_clusters)
    log.info("Computing metrics for %d clusters", n_clusters)

    # Compute metrics per cluster
    records = []
    mixed_low, mixed_high = mixed_range

    for cluster_id in unique_clusters:
        mask = clusters == cluster_id
        n_cells = mask.sum()

        if n_cells == 0:
            continue

        # Get expression for this cluster
        X_cluster = X[mask, :]

        # Compute positivity fractions per marker
        p_values = []
        n_mixed = 0
        n_valid_markers = 0

        for i, marker in enumerate(markers):
            thr = thresholds.get(marker)
            if thr is None:
                continue

            expr = X_cluster[:, i]
            p_m = (expr > thr).mean()
            p_values.append(p_m)
            n_valid_markers += 1

            # Check if in mixed range
            if mixed_low <= p_m <= mixed_high:
                n_mixed += 1

        if n_valid_markers == 0:
            continue

        p_values = np.array(p_values)

        # Compute heterogeneity: mean of p * (1 - p)
        # This is the variance of a Bernoulli distribution
        heterogeneity = np.mean(p_values * (1 - p_values))

        # Compute mixed marker rate
        mixed_rate = n_mixed / n_valid_markers

        records.append({
            "cluster_id": cluster_id,
            "n_cells": int(n_cells),
            "marker_heterogeneity": float(heterogeneity),
            "mixed_marker_rate": float(mixed_rate),
            "n_markers_mixed": int(n_mixed),
            "n_markers_total": int(n_valid_markers),
        })

    df = pd.DataFrame.from_records(records)

    # Add heterogeneity band
    if not df.empty:
        df["heterogeneity_band"] = df["marker_heterogeneity"].apply(_classify_heterogeneity)

    log.info(
        "Heterogeneity computation complete: mean=%.3f, max=%.3f",
        df["marker_heterogeneity"].mean() if not df.empty else 0,
        df["marker_heterogeneity"].max() if not df.empty else 0,
    )

    return df


def _get_marker_thresholds(
    adata: "sc.AnnData",
    markers: List[str],
    marker_evidence_path: Optional[Path] = None,
    logger: Optional[logging.Logger] = None,
) -> Dict[str, float]:
    """
    Get per-marker positivity thresholds.

    Priority:
    1. Load from marker_evidence.csv (Stage H thresholds)
    2. Fall back to global 75th percentile per marker

    Parameters
    ----------
    adata : sc.AnnData
        AnnData object
    markers : List[str]
        List of marker names
    marker_evidence_path : Path, optional
        Path to marker_evidence.csv
    logger : logging.Logger, optional
        Logger

    Returns
    -------
    Dict[str, float]
        Mapping from marker name to threshold value
    """
    log = logger or logging.getLogger(__name__)
    thresholds = {}

    # Try loading from marker_evidence.csv
    if marker_evidence_path and marker_evidence_path.exists():
        log.info("Loading thresholds from %s", marker_evidence_path)
        try:
            evidence_df = pd.read_csv(marker_evidence_path)
            if "marker" in evidence_df.columns and "threshold_used" in evidence_df.columns:
                # Get unique threshold per marker (they should be consistent)
                for marker in markers:
                    marker_rows = evidence_df[evidence_df["marker"] == marker]
                    if not marker_rows.empty:
                        thr = marker_rows["threshold_used"].iloc[0]
                        if pd.notna(thr):
                            thresholds[marker] = float(thr)

                log.info("Loaded %d thresholds from marker_evidence.csv", len(thresholds))
        except Exception as e:
            log.warning("Failed to load marker_evidence.csv: %s", e)

    # Fall back to global percentiles for missing markers
    missing = [m for m in markers if m not in thresholds]
    if missing:
        log.info("Computing 75th percentile thresholds for %d markers", len(missing))
        # Get expression matrix
        X = adata.layers.get("batchcorr", adata.X)
        if hasattr(X, "toarray"):
            X = X.toarray()

        for i, marker in enumerate(markers):
            if marker in missing:
                values = X[:, i]
                thr = np.percentile(values, 75)
                thresholds[marker] = float(thr)

    return thresholds


def _classify_heterogeneity(value: float) -> str:
    """
    Classify heterogeneity value into bands.

    Bands based on the theoretical maximum of 0.25:
    - low: < 0.10 (< 40% of max)
    - medium: 0.10 - 0.18 (40-72% of max)
    - high: >= 0.18 (>= 72% of max)
    """
    if value < 0.10:
        return "low"
    elif value < 0.18:
        return "medium"
    else:
        return "high"


def merge_heterogeneity_with_annotations(
    cluster_annotations: pd.DataFrame,
    heterogeneity_metrics: pd.DataFrame,
) -> pd.DataFrame:
    """
    Merge heterogeneity metrics into cluster annotations DataFrame.

    Parameters
    ----------
    cluster_annotations : pd.DataFrame
        Cluster annotations from Stage H (must have cluster_id column)
    heterogeneity_metrics : pd.DataFrame
        Output from compute_cluster_marker_heterogeneity()

    Returns
    -------
    pd.DataFrame
        cluster_annotations with added heterogeneity columns
    """
    # Ensure cluster_id is string in both
    annotations = cluster_annotations.copy()
    annotations["cluster_id"] = annotations["cluster_id"].astype(str)

    metrics = heterogeneity_metrics.copy()
    metrics["cluster_id"] = metrics["cluster_id"].astype(str)

    # Merge on cluster_id
    merged = annotations.merge(
        metrics[["cluster_id", "marker_heterogeneity", "mixed_marker_rate",
                 "n_markers_mixed", "heterogeneity_band"]],
        on="cluster_id",
        how="left",
    )

    return merged
