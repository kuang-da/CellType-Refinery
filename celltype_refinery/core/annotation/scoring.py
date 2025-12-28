"""Marker scoring for cell-type annotation.

This module computes marker enrichment scores for clusters against
cell-type marker sets.

Supports parallel scoring via joblib for large cluster counts (1000+).
"""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
from scipy import sparse

try:
    import scanpy as sc
except ImportError:
    sc = None

try:
    from joblib import Parallel, delayed
    JOBLIB_AVAILABLE = True
except ImportError:
    JOBLIB_AVAILABLE = False

from .marker_loading import MarkerSet, canonicalize_marker


@dataclass
class ScoringContext:
    """Pre-computed context shared across parallel workers.

    This avoids recomputing global stats for each cluster.
    """
    matrix: np.ndarray  # Dense expression matrix (n_cells, n_markers)
    obs_clusters: np.ndarray  # Cluster assignments as numpy array
    global_median: np.ndarray
    global_std: np.ndarray
    thresholds: np.ndarray  # Positive thresholds per marker
    var_index: Dict[str, int]  # Marker name -> column index
    var_names: List[str]


def compute_marker_idf(
    marker_sets: List[MarkerSet],
    smoothing: float = 1.0,
    logger: Optional[logging.Logger] = None,
) -> Dict[str, float]:
    """Compute IDF (Inverse Document Frequency) weights for markers.

    IDF weights give higher importance to markers that are specific to
    fewer cell types. A marker used by many cell types gets lower weight.

    Formula: IDF(marker) = log(N / (df(marker) + smoothing))
    where N = number of cell types, df = document frequency

    Args:
        marker_sets: List of MarkerSet objects
        smoothing: Smoothing constant to prevent division by zero
        logger: Optional logger instance

    Returns:
        Dict mapping marker name → IDF weight
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    # Count how many cell types use each marker
    doc_freq: Dict[str, int] = {}
    for mset in marker_sets:
        for marker in mset.resolved_markers:
            doc_freq[marker] = doc_freq.get(marker, 0) + 1

    N = len(marker_sets)
    idf_weights: Dict[str, float] = {}

    for marker, df in doc_freq.items():
        idf = np.log(N / (df + smoothing))
        idf_weights[marker] = float(max(idf, 0.1))  # Floor at 0.1

    logger.debug(
        "Computed IDF weights for %d markers (N=%d cell types)",
        len(idf_weights),
        N,
    )

    return idf_weights


def extract_existing_de_results(
    adata: "sc.AnnData",
    key: str,
    logger: logging.Logger,
) -> Dict[str, List[str]]:
    """Extract differential expression results from AnnData.

    Retrieves pre-computed DE results stored in adata.uns by
    scanpy's rank_genes_groups.

    Args:
        adata: AnnData with DE results in .uns[key]
        key: Key in adata.uns containing DE results
        logger: Logger instance

    Returns:
        Dict mapping cluster_id → list of DE gene names
    """
    if key not in adata.uns:
        logger.warning("DE key '%s' not found in adata.uns", key)
        return {}

    de_results = adata.uns[key]
    names = de_results.get("names")
    if names is None:
        logger.warning("No 'names' in DE results at '%s'", key)
        return {}

    result: Dict[str, List[str]] = {}
    clusters = list(names.dtype.names) if hasattr(names.dtype, "names") else []

    for cluster in clusters:
        try:
            cluster_names = names[cluster]
            result[str(cluster)] = [str(n) for n in cluster_names]
        except (KeyError, IndexError):
            if ":" in str(cluster):
                logger.debug(
                    "No DE results for subcluster %s (expected)", cluster
                )
            else:
                logger.warning("No DE results for cluster %s in %s", cluster, key)
            result[str(cluster)] = []

    logger.info(
        "Extracted existing DE results for %d clusters from adata.uns['%s']",
        len(result),
        key,
    )
    return result


def build_de_rank_lookup(
    adata: "sc.AnnData",
    cluster_key: str,
    de_key: str,
    de_top_frac: float,
    de_min_k: int,
    de_max_k: int,
    logger: logging.Logger,
) -> Dict[str, Dict[str, Any]]:
    """Build per-cluster DE rank maps for upregulated markers only.

    Creates a rank-based lookup that enables continuous DE bonus scoring
    (rank-weighted) instead of binary hit/miss. Only upregulated markers
    (Wilcoxon scores > 0) are included.

    Args:
        adata: AnnData with rank_genes_groups results in .uns[de_key]
        cluster_key: Column name for cluster assignments
        de_key: Key in adata.uns containing DE results
        de_top_frac: Fraction of panel size for target K
        de_min_k: Minimum K
        de_max_k: Maximum K / cap
        logger: Logger instance

    Returns:
        Dict[cluster_id] = {
            "K": int,              # Effective K used
            "K_target": int,       # Calculated K before clamping
            "n_pos": int,          # Number of upregulated markers
            "rank_map": Dict[str, int],  # canonical_marker -> rank (1-indexed)
            "topK_names": List[str],
            "topK_scores": List[float],
        }
    """
    if de_key not in adata.uns:
        logger.warning(
            "DE key '%s' not found in adata.uns, returning empty rank lookup",
            de_key,
        )
        return {}

    M = adata.n_vars  # Panel size
    K_target_raw = int(np.ceil(M * de_top_frac))
    K_target = max(de_min_k, min(de_max_k, K_target_raw))
    logger.info(
        "DE rank lookup: panel_size=%d, de_top_frac=%.2f -> K_target=%d",
        M,
        de_top_frac,
        K_target,
    )

    de_rank_lookup: Dict[str, Dict[str, Any]] = {}
    clusters = sorted(adata.obs[cluster_key].astype(str).unique())

    for cluster in clusters:
        try:
            df = sc.get.rank_genes_groups_df(adata, group=str(cluster), key=de_key)
        except (KeyError, ValueError):
            if ":" in str(cluster):
                logger.debug("No DE results for subcluster %s (expected)", cluster)
            else:
                logger.warning("No DE results for cluster %s", cluster)
            de_rank_lookup[str(cluster)] = {
                "K": 0,
                "K_target": K_target,
                "n_pos": 0,
                "rank_map": {},
                "topK_names": [],
                "topK_scores": [],
            }
            continue

        # Sort by scores descending
        df = df.sort_values("scores", ascending=False)

        # Filter to UPREGULATED only (scores > 0)
        df_up = df[df["scores"] > 0]
        n_pos = len(df_up)

        # Effective K: min(K_target, n_pos)
        K = min(K_target, n_pos)

        rank_map: Dict[str, int] = {}
        topK_names: List[str] = []
        topK_scores: List[float] = []

        if K > 0:
            top_df = df_up.head(K)
            for rank_idx, (_, row) in enumerate(top_df.iterrows(), start=1):
                marker_name = str(row["names"])
                canon = canonicalize_marker(marker_name)
                rank_map[canon] = rank_idx
                topK_names.append(marker_name)
                topK_scores.append(float(row["scores"]))

        de_rank_lookup[str(cluster)] = {
            "K": K,
            "K_target": K_target,
            "n_pos": n_pos,
            "rank_map": rank_map,
            "topK_names": topK_names,
            "topK_scores": topK_scores,
        }

    logger.info("Built DE rank lookup for %d clusters", len(de_rank_lookup))
    return de_rank_lookup


def compute_marker_scores(
    adata: "sc.AnnData",
    marker_sets: Sequence[MarkerSet],
    cluster_key: str,
    positive_quantile: float,
    de_lookup: Dict[str, List[str]],
    de_bonus: float,
    anti_weight: float,
    logger: logging.Logger,
    idf_weights: Optional[Dict[str, float]] = None,
    expand_markers: bool = False,
    anti_agg: str = "top2mean",
    de_rank_lookup: Optional[Dict[str, Dict[str, Any]]] = None,
    marker_doc_freq: Optional[Dict[str, int]] = None,
    de_commonness_alpha: float = 0.5,
    layer: str = "X",
) -> Union[pd.DataFrame, Tuple[pd.DataFrame, pd.DataFrame]]:
    """Score clusters against marker sets with anti-marker penalties.

    The scoring formula is:
        score = mean_enrichment + mean_positive + de_component - anti_penalty

    Where:
        - mean_enrichment: z-score of cluster median vs global median
        - mean_positive: fraction of cells above quantile threshold
        - de_component: rank-weighted DE bonus with commonness penalty
        - anti_penalty: weighted penalty for anti-marker expression

    Args:
        adata: AnnData object with expression data
        marker_sets: List of MarkerSet objects to score against
        cluster_key: Column name in adata.obs with cluster labels
        positive_quantile: Quantile for positive threshold
        de_lookup: Dict mapping cluster → list of DE genes (legacy)
        de_bonus: Maximum bonus for DE overlap
        anti_weight: Weight for anti-marker penalty
        logger: Logger instance
        idf_weights: Optional dict mapping marker → IDF weight
        expand_markers: If True, also return per-marker evidence table
        anti_agg: Aggregation mode for anti-markers: "mean", "top2mean", "max"
        de_rank_lookup: Dict mapping cluster → DE rank info
        marker_doc_freq: Dict mapping marker → document frequency
        de_commonness_alpha: Exponent for commonness penalty
        layer: Which layer to use for expression data

    Returns:
        If expand_markers=False: DataFrame with marker scores per cluster
        If expand_markers=True: Tuple of (marker_scores_df, marker_evidence_df)
    """
    idf_status = "IDF-weighted" if idf_weights else "unweighted"
    logger.info(
        "Scoring clusters against %d marker sets (%s)",
        len(marker_sets),
        idf_status,
    )

    if cluster_key not in adata.obs:
        raise ValueError(f"Missing cluster key '{cluster_key}' in AnnData.obs")

    obs_clusters = adata.obs[cluster_key].astype(str)

    # Get expression matrix
    if layer in adata.layers:
        matrix = adata.layers[layer]
    else:
        matrix = adata.X
    matrix = matrix.toarray() if sparse.issparse(matrix) else np.asarray(matrix)

    global_median = np.nanmedian(matrix, axis=0)
    global_std = np.nanstd(matrix, axis=0)
    global_std[global_std == 0] = 1e-6
    thresholds = np.nanquantile(matrix, positive_quantile, axis=0)
    var_index = {name: idx for idx, name in enumerate(adata.var_names.astype(str))}

    records: List[Dict[str, object]] = []
    marker_evidence_records: List[Dict[str, object]] = [] if expand_markers else None

    clusters = sorted(obs_clusters.unique())
    for cluster in clusters:
        mask = obs_clusters == cluster
        n_cells = int(mask.sum())
        if n_cells == 0:
            continue
        cluster_matrix = matrix[mask.to_numpy(), :]

        for mset in marker_sets:
            idxs = [var_index[name] for name in mset.resolved_markers if name in var_index]
            if not idxs:
                continue

            # Positive marker scoring
            cluster_values = cluster_matrix[:, idxs]
            median_vals = np.nanmedian(cluster_values, axis=0)
            positive_vals = np.mean(cluster_values >= thresholds[idxs], axis=0)
            enrichment_vals = (median_vals - global_median[idxs]) / global_std[idxs]
            mean_median = float(np.nanmean(median_vals))
            coverage = len(idxs) / max(len(mset.markers), 1)

            # Compute marker weights (IDF-based if enabled)
            resolved_markers = [m for m in mset.resolved_markers if m in var_index]
            if idf_weights and resolved_markers:
                raw_weights = np.array([idf_weights.get(m, 1.0) for m in resolved_markers])
                marker_weights = raw_weights / raw_weights.sum() if raw_weights.sum() > 0 else None
                marker_weights_str = ";".join(
                    f"{m}:{w:.3f}" for m, w in zip(resolved_markers, raw_weights)
                )
            else:
                marker_weights = None
                marker_weights_str = ""

            # Compute mean_positive and mean_enrichment
            if marker_weights is not None:
                mean_positive = float(np.average(positive_vals, weights=marker_weights))
                mean_enrichment = float(np.average(enrichment_vals, weights=marker_weights))
            else:
                mean_positive = float(np.nanmean(positive_vals))
                mean_enrichment = float(np.nanmean(enrichment_vals))

            # n_markers_on: count markers with enrichment > 0 AND pos_frac > 0.15
            markers_on_mask = (enrichment_vals > 0) & (positive_vals > 0.15)
            n_markers_on = int(np.sum(markers_on_mask))
            frac_markers_on = n_markers_on / len(idxs) if idxs else 0.0

            # Anti-marker penalty calculation
            anti_idxs = [var_index[name] for name in mset.resolved_anti_markers if name in var_index]
            if anti_idxs and anti_weight > 0:
                anti_cluster_values = cluster_matrix[:, anti_idxs]
                anti_median_vals = np.nanmedian(anti_cluster_values, axis=0)
                anti_positive_vals = np.mean(anti_cluster_values >= thresholds[anti_idxs], axis=0)
                anti_enrichment_vals = (anti_median_vals - global_median[anti_idxs]) / global_std[anti_idxs]

                # Clip enrichment to [0, inf) - only penalize when above global median
                anti_enrichment_clipped = np.clip(anti_enrichment_vals, 0, None)

                # Apply aggregation mode
                n_anti = len(anti_enrichment_clipped)
                if anti_agg == "max":
                    agg_anti_enrichment = float(np.max(anti_enrichment_clipped))
                    agg_anti_positive = float(np.max(anti_positive_vals))
                elif anti_agg == "top2mean" and n_anti >= 2:
                    top2_enrich_idx = np.argsort(anti_enrichment_clipped)[-2:]
                    top2_pos_idx = np.argsort(anti_positive_vals)[-2:]
                    agg_anti_enrichment = float(np.mean(anti_enrichment_clipped[top2_enrich_idx]))
                    agg_anti_positive = float(np.mean(anti_positive_vals[top2_pos_idx]))
                else:
                    agg_anti_enrichment = float(np.nanmean(anti_enrichment_clipped))
                    agg_anti_positive = float(np.nanmean(anti_positive_vals))

                mean_anti_enrichment = agg_anti_enrichment
                mean_anti_positive = agg_anti_positive
                anti_penalty = anti_weight * (agg_anti_enrichment + agg_anti_positive)
            else:
                mean_anti_enrichment = 0.0
                mean_anti_positive = 0.0
                anti_penalty = 0.0
                anti_median_vals = np.array([])
                anti_positive_vals = np.array([])
                anti_enrichment_vals = np.array([])

            # Collect per-marker evidence if requested
            if marker_evidence_records is not None:
                cluster_de_markers = set()
                if cluster in de_lookup:
                    cluster_de_markers = {canonicalize_marker(n) for n in de_lookup[cluster]}

                for i, marker in enumerate(resolved_markers):
                    marker_idx = idxs[i]
                    is_de = canonicalize_marker(marker) in cluster_de_markers
                    is_on = enrichment_vals[i] > 0 and positive_vals[i] > 0.15
                    marker_evidence_records.append({
                        "cluster_id": cluster,
                        "cell_type": mset.label,
                        "path": " / ".join(mset.path),
                        "level": mset.level,
                        "marker": marker,
                        "role": "positive",
                        "enrichment": float(enrichment_vals[i]),
                        "pos_frac": float(positive_vals[i]),
                        "median_val": float(median_vals[i]),
                        "threshold_used": float(thresholds[marker_idx]),
                        "global_median": float(global_median[marker_idx]),
                        "global_std": float(global_std[marker_idx]),
                        "is_de_hit": is_de,
                        "idf_weight": idf_weights.get(marker, 1.0) if idf_weights else 1.0,
                        "is_on": is_on,
                        "n_cells": n_cells,
                    })

                # Add anti-markers
                anti_resolved = [m for m in mset.resolved_anti_markers if m in var_index]
                for j, anti_marker in enumerate(anti_resolved):
                    anti_marker_idx = anti_idxs[j]
                    marker_evidence_records.append({
                        "cluster_id": cluster,
                        "cell_type": mset.label,
                        "path": " / ".join(mset.path),
                        "level": mset.level,
                        "marker": anti_marker,
                        "role": "anti",
                        "enrichment": float(anti_enrichment_vals[j]) if len(anti_enrichment_vals) > j else 0.0,
                        "pos_frac": float(anti_positive_vals[j]) if len(anti_positive_vals) > j else 0.0,
                        "median_val": float(anti_median_vals[j]) if len(anti_median_vals) > j else 0.0,
                        "threshold_used": float(thresholds[anti_marker_idx]),
                        "global_median": float(global_median[anti_marker_idx]),
                        "global_std": float(global_std[anti_marker_idx]),
                        "is_de_hit": False,
                        "idf_weight": 1.0,
                        "is_on": False,
                        "n_cells": n_cells,
                    })

            # DE bonus calculation
            de_score = 0.0
            de_component = 0.0
            de_K = 0
            de_n_pos = 0
            de_hits = 0

            if de_rank_lookup is not None and cluster in de_rank_lookup:
                info = de_rank_lookup[cluster]
                de_K = info["K"]
                de_n_pos = info["n_pos"]
                rank_map = info["rank_map"]

                if de_K > 0 and resolved_markers:
                    contribs = []
                    for marker in resolved_markers:
                        canon = canonicalize_marker(marker)
                        r = rank_map.get(canon)
                        if r is None or r > de_K:
                            contribs.append(0.0)
                        else:
                            de_hits += 1
                            w_rank = (de_K - r + 1) / de_K
                            freq = marker_doc_freq.get(marker, 1) if marker_doc_freq else 1
                            freq = max(freq, 1)
                            pen = freq ** de_commonness_alpha
                            contribs.append(w_rank / pen)

                    de_score = float(np.mean(contribs)) if contribs else 0.0
                    de_component = de_bonus * de_score
            elif cluster in de_lookup:
                cluster_markers = {canonicalize_marker(n) for n in de_lookup[cluster]}
                de_hits = sum(
                    1 for marker in mset.resolved_markers if canonicalize_marker(marker) in cluster_markers
                )
                n_resolved = len(mset.resolved_markers)
                de_score = de_hits / n_resolved if n_resolved > 0 else 0.0
                de_component = de_bonus * de_score

            # Final score
            score = mean_enrichment + mean_positive + de_component - anti_penalty

            records.append({
                "cluster_id": cluster,
                "label": mset.label,
                "path": " / ".join(mset.path),
                "level": mset.level,
                "n_cells": n_cells,
                "markers_expected": len(mset.markers),
                "markers_found": len(idxs),
                "coverage": coverage,
                "mean_median": mean_median,
                "mean_positive_fraction": mean_positive,
                "mean_enrichment": mean_enrichment,
                "n_markers_on": n_markers_on,
                "frac_markers_on": frac_markers_on,
                "de_hits": de_hits,
                "de_score": de_score,
                "de_K": de_K,
                "de_n_pos": de_n_pos,
                "de_component": de_component,
                "anti_markers_expected": len(mset.anti_markers),
                "anti_markers_found": len(anti_idxs),
                "mean_anti_enrichment": mean_anti_enrichment,
                "mean_anti_positive": mean_anti_positive,
                "anti_penalty": anti_penalty,
                "score": score,
                "idf_weighted": bool(idf_weights),
                "marker_weights": marker_weights_str,
                "missing_markers": ";".join(mset.missing_markers),
                "resolved_markers": ";".join(mset.resolved_markers),
                "missing_anti_markers": ";".join(mset.missing_anti_markers),
                "resolved_anti_markers": ";".join(mset.resolved_anti_markers),
            })

    if not records:
        logger.warning("No marker scores were computed")
        if expand_markers:
            return pd.DataFrame(), pd.DataFrame()
        return pd.DataFrame()

    marker_scores_df = pd.DataFrame.from_records(records)

    if expand_markers:
        marker_evidence_df = (
            pd.DataFrame.from_records(marker_evidence_records)
            if marker_evidence_records
            else pd.DataFrame()
        )
        logger.info("Generated marker evidence table with %d rows", len(marker_evidence_df))
        return marker_scores_df, marker_evidence_df

    return marker_scores_df


def _score_cluster_batch(
    cluster_batch: List[str],
    matrix: np.ndarray,
    obs_clusters: np.ndarray,
    global_median: np.ndarray,
    global_std: np.ndarray,
    thresholds: np.ndarray,
    var_index: Dict[str, int],
    marker_sets: Sequence[MarkerSet],
    de_lookup: Dict[str, List[str]],
    de_rank_lookup: Optional[Dict[str, Dict[str, Any]]],
    de_bonus: float,
    anti_weight: float,
    idf_weights: Optional[Dict[str, float]],
    marker_doc_freq: Optional[Dict[str, int]],
    de_commonness_alpha: float,
    anti_agg: str,
    expand_markers: bool,
) -> Tuple[List[Dict], List[Dict]]:
    """Score a batch of clusters (worker function for parallel execution).

    This function is designed to run in a separate process via joblib.
    It receives pre-computed global statistics to avoid redundant computation.

    Returns:
        Tuple of (records_list, marker_evidence_list)
    """
    records = []
    marker_evidence_records = [] if expand_markers else []

    for cluster in cluster_batch:
        mask = obs_clusters == cluster
        n_cells = int(mask.sum())
        if n_cells == 0:
            continue
        cluster_matrix = matrix[mask, :]

        for mset in marker_sets:
            idxs = [var_index[name] for name in mset.resolved_markers if name in var_index]
            if not idxs:
                continue

            # Positive marker scoring
            cluster_values = cluster_matrix[:, idxs]
            median_vals = np.nanmedian(cluster_values, axis=0)
            positive_vals = np.mean(cluster_values >= thresholds[idxs], axis=0)
            enrichment_vals = (median_vals - global_median[idxs]) / global_std[idxs]
            mean_median = float(np.nanmean(median_vals))
            coverage = len(idxs) / max(len(mset.markers), 1)

            # Compute marker weights (IDF-based if enabled)
            resolved_markers = [m for m in mset.resolved_markers if m in var_index]
            if idf_weights and resolved_markers:
                raw_weights = np.array([idf_weights.get(m, 1.0) for m in resolved_markers])
                marker_weights = raw_weights / raw_weights.sum() if raw_weights.sum() > 0 else None
                marker_weights_str = ";".join(
                    f"{m}:{w:.3f}" for m, w in zip(resolved_markers, raw_weights)
                )
            else:
                marker_weights = None
                marker_weights_str = ""

            # Compute mean_positive and mean_enrichment
            if marker_weights is not None:
                mean_positive = float(np.average(positive_vals, weights=marker_weights))
                mean_enrichment = float(np.average(enrichment_vals, weights=marker_weights))
            else:
                mean_positive = float(np.nanmean(positive_vals))
                mean_enrichment = float(np.nanmean(enrichment_vals))

            # n_markers_on: count markers with enrichment > 0 AND pos_frac > 0.15
            markers_on_mask = (enrichment_vals > 0) & (positive_vals > 0.15)
            n_markers_on = int(np.sum(markers_on_mask))
            frac_markers_on = n_markers_on / len(idxs) if idxs else 0.0

            # Anti-marker penalty calculation
            anti_idxs = [var_index[name] for name in mset.resolved_anti_markers if name in var_index]
            if anti_idxs and anti_weight > 0:
                anti_cluster_values = cluster_matrix[:, anti_idxs]
                anti_median_vals = np.nanmedian(anti_cluster_values, axis=0)
                anti_positive_vals = np.mean(anti_cluster_values >= thresholds[anti_idxs], axis=0)
                anti_enrichment_vals = (anti_median_vals - global_median[anti_idxs]) / global_std[anti_idxs]

                # Clip enrichment to [0, inf) - only penalize when above global median
                anti_enrichment_clipped = np.clip(anti_enrichment_vals, 0, None)

                # Apply aggregation mode
                n_anti = len(anti_enrichment_clipped)
                if anti_agg == "max":
                    agg_anti_enrichment = float(np.max(anti_enrichment_clipped))
                    agg_anti_positive = float(np.max(anti_positive_vals))
                elif anti_agg == "top2mean" and n_anti >= 2:
                    top2_enrich_idx = np.argsort(anti_enrichment_clipped)[-2:]
                    top2_pos_idx = np.argsort(anti_positive_vals)[-2:]
                    agg_anti_enrichment = float(np.mean(anti_enrichment_clipped[top2_enrich_idx]))
                    agg_anti_positive = float(np.mean(anti_positive_vals[top2_pos_idx]))
                else:
                    agg_anti_enrichment = float(np.nanmean(anti_enrichment_clipped))
                    agg_anti_positive = float(np.nanmean(anti_positive_vals))

                mean_anti_enrichment = agg_anti_enrichment
                mean_anti_positive = agg_anti_positive
                anti_penalty = anti_weight * (agg_anti_enrichment + agg_anti_positive)
            else:
                mean_anti_enrichment = 0.0
                mean_anti_positive = 0.0
                anti_penalty = 0.0
                anti_median_vals = np.array([])
                anti_positive_vals = np.array([])
                anti_enrichment_vals = np.array([])

            # Collect per-marker evidence if requested
            if expand_markers:
                cluster_de_markers = set()
                if cluster in de_lookup:
                    cluster_de_markers = {canonicalize_marker(n) for n in de_lookup[cluster]}

                for i, marker in enumerate(resolved_markers):
                    marker_idx = idxs[i]
                    is_de = canonicalize_marker(marker) in cluster_de_markers
                    is_on = enrichment_vals[i] > 0 and positive_vals[i] > 0.15
                    marker_evidence_records.append({
                        "cluster_id": cluster,
                        "cell_type": mset.label,
                        "path": " / ".join(mset.path),
                        "level": mset.level,
                        "marker": marker,
                        "role": "positive",
                        "enrichment": float(enrichment_vals[i]),
                        "pos_frac": float(positive_vals[i]),
                        "median_val": float(median_vals[i]),
                        "threshold_used": float(thresholds[idxs[i]]),
                        "global_median": float(global_median[idxs[i]]),
                        "global_std": float(global_std[idxs[i]]),
                        "is_de_hit": is_de,
                        "idf_weight": idf_weights.get(marker, 1.0) if idf_weights else 1.0,
                        "is_on": is_on,
                        "n_cells": n_cells,
                    })

                # Add anti-markers
                anti_resolved = [m for m in mset.resolved_anti_markers if m in var_index]
                for j, anti_marker in enumerate(anti_resolved):
                    marker_evidence_records.append({
                        "cluster_id": cluster,
                        "cell_type": mset.label,
                        "path": " / ".join(mset.path),
                        "level": mset.level,
                        "marker": anti_marker,
                        "role": "anti",
                        "enrichment": float(anti_enrichment_vals[j]) if len(anti_enrichment_vals) > j else 0.0,
                        "pos_frac": float(anti_positive_vals[j]) if len(anti_positive_vals) > j else 0.0,
                        "median_val": float(anti_median_vals[j]) if len(anti_median_vals) > j else 0.0,
                        "threshold_used": float(thresholds[anti_idxs[j]]),
                        "global_median": float(global_median[anti_idxs[j]]),
                        "global_std": float(global_std[anti_idxs[j]]),
                        "is_de_hit": False,
                        "idf_weight": 1.0,
                        "is_on": False,
                        "n_cells": n_cells,
                    })

            # DE bonus calculation
            de_score = 0.0
            de_component = 0.0
            de_K = 0
            de_n_pos = 0
            de_hits = 0

            if de_rank_lookup is not None and cluster in de_rank_lookup:
                info = de_rank_lookup[cluster]
                de_K = info["K"]
                de_n_pos = info["n_pos"]
                rank_map = info["rank_map"]

                if de_K > 0 and resolved_markers:
                    contribs = []
                    for marker in resolved_markers:
                        canon = canonicalize_marker(marker)
                        r = rank_map.get(canon)
                        if r is None or r > de_K:
                            contribs.append(0.0)
                        else:
                            de_hits += 1
                            w_rank = (de_K - r + 1) / de_K
                            freq = marker_doc_freq.get(marker, 1) if marker_doc_freq else 1
                            freq = max(freq, 1)
                            pen = freq ** de_commonness_alpha
                            contribs.append(w_rank / pen)

                    de_score = float(np.mean(contribs)) if contribs else 0.0
                    de_component = de_bonus * de_score
            elif cluster in de_lookup:
                cluster_markers = {canonicalize_marker(n) for n in de_lookup[cluster]}
                de_hits = sum(
                    1 for marker in mset.resolved_markers if canonicalize_marker(marker) in cluster_markers
                )
                n_resolved = len(mset.resolved_markers)
                de_score = de_hits / n_resolved if n_resolved > 0 else 0.0
                de_component = de_bonus * de_score

            # Final score
            score = mean_enrichment + mean_positive + de_component - anti_penalty

            records.append({
                "cluster_id": cluster,
                "label": mset.label,
                "path": " / ".join(mset.path),
                "level": mset.level,
                "n_cells": n_cells,
                "markers_expected": len(mset.markers),
                "markers_found": len(idxs),
                "coverage": coverage,
                "mean_median": mean_median,
                "mean_positive_fraction": mean_positive,
                "mean_enrichment": mean_enrichment,
                "n_markers_on": n_markers_on,
                "frac_markers_on": frac_markers_on,
                "de_hits": de_hits,
                "de_score": de_score,
                "de_K": de_K,
                "de_n_pos": de_n_pos,
                "de_component": de_component,
                "anti_markers_expected": len(mset.anti_markers),
                "anti_markers_found": len(anti_idxs),
                "mean_anti_enrichment": mean_anti_enrichment,
                "mean_anti_positive": mean_anti_positive,
                "anti_penalty": anti_penalty,
                "score": score,
                "idf_weighted": bool(idf_weights),
                "marker_weights": marker_weights_str,
                "missing_markers": ";".join(mset.missing_markers),
                "resolved_markers": ";".join(mset.resolved_markers),
                "missing_anti_markers": ";".join(mset.missing_anti_markers),
                "resolved_anti_markers": ";".join(mset.resolved_anti_markers),
            })

    return records, marker_evidence_records


def compute_marker_scores_parallel(
    adata: "sc.AnnData",
    marker_sets: Sequence[MarkerSet],
    cluster_key: str,
    positive_quantile: float,
    de_lookup: Dict[str, List[str]],
    de_bonus: float,
    anti_weight: float,
    logger: logging.Logger,
    idf_weights: Optional[Dict[str, float]] = None,
    expand_markers: bool = False,
    anti_agg: str = "top2mean",
    de_rank_lookup: Optional[Dict[str, Dict[str, Any]]] = None,
    marker_doc_freq: Optional[Dict[str, int]] = None,
    de_commonness_alpha: float = 0.5,
    layer: str = "X",
    n_workers: int = 1,
    batch_size: int = 50,
) -> Union[pd.DataFrame, Tuple[pd.DataFrame, pd.DataFrame]]:
    """Score clusters against marker sets with parallel execution.

    This is a drop-in replacement for compute_marker_scores that adds
    parallel execution capability for large cluster counts.

    When n_workers=1, falls back to sequential processing (identical to
    compute_marker_scores). When n_workers>1, uses joblib to parallelize
    cluster scoring.

    Args:
        adata: AnnData object with expression data
        marker_sets: List of MarkerSet objects to score against
        cluster_key: Column name in adata.obs with cluster labels
        positive_quantile: Quantile for positive threshold
        de_lookup: Dict mapping cluster → list of DE genes
        de_bonus: Maximum bonus for DE overlap
        anti_weight: Weight for anti-marker penalty
        logger: Logger instance
        idf_weights: Optional dict mapping marker → IDF weight
        expand_markers: If True, also return per-marker evidence table
        anti_agg: Aggregation mode for anti-markers
        de_rank_lookup: Dict mapping cluster → DE rank info
        marker_doc_freq: Dict mapping marker → document frequency
        de_commonness_alpha: Exponent for commonness penalty
        layer: Which layer to use for expression data
        n_workers: Number of parallel workers (1=sequential)
        batch_size: Clusters per worker batch (default: 50)

    Returns:
        If expand_markers=False: DataFrame with marker scores per cluster
        If expand_markers=True: Tuple of (marker_scores_df, marker_evidence_df)
    """
    # Fall back to sequential if n_workers=1 or joblib unavailable
    if n_workers <= 1 or not JOBLIB_AVAILABLE:
        return compute_marker_scores(
            adata=adata,
            marker_sets=marker_sets,
            cluster_key=cluster_key,
            positive_quantile=positive_quantile,
            de_lookup=de_lookup,
            de_bonus=de_bonus,
            anti_weight=anti_weight,
            logger=logger,
            idf_weights=idf_weights,
            expand_markers=expand_markers,
            anti_agg=anti_agg,
            de_rank_lookup=de_rank_lookup,
            marker_doc_freq=marker_doc_freq,
            de_commonness_alpha=de_commonness_alpha,
            layer=layer,
        )

    idf_status = "IDF-weighted" if idf_weights else "unweighted"
    logger.info(
        "Scoring clusters against %d marker sets (%s) with %d workers",
        len(marker_sets),
        idf_status,
        n_workers,
    )

    if cluster_key not in adata.obs:
        raise ValueError(f"Missing cluster key '{cluster_key}' in AnnData.obs")

    start_time = time.time()

    # Pre-compute global statistics (done once, shared with workers)
    obs_clusters = adata.obs[cluster_key].astype(str).values  # Convert to numpy

    # Get expression matrix
    if layer in adata.layers:
        matrix = adata.layers[layer]
    else:
        matrix = adata.X
    matrix = matrix.toarray() if sparse.issparse(matrix) else np.asarray(matrix)

    global_median = np.nanmedian(matrix, axis=0)
    global_std = np.nanstd(matrix, axis=0)
    global_std[global_std == 0] = 1e-6
    thresholds = np.nanquantile(matrix, positive_quantile, axis=0)
    var_index = {name: idx for idx, name in enumerate(adata.var_names.astype(str))}

    prep_time = time.time() - start_time
    logger.info("Pre-computed global stats in %.2f sec", prep_time)

    # Get unique clusters and create batches
    clusters = sorted(np.unique(obs_clusters))
    n_clusters = len(clusters)

    # Create batches
    cluster_batches = [
        clusters[i:i + batch_size]
        for i in range(0, n_clusters, batch_size)
    ]

    logger.info(
        "Scoring %d clusters in %d batches (batch_size=%d)",
        n_clusters, len(cluster_batches), batch_size
    )

    # Run parallel scoring
    score_start = time.time()

    # Use loky backend for process isolation, verbose=5 for moderate logging
    results = Parallel(n_jobs=n_workers, backend='loky', verbose=5)(
        delayed(_score_cluster_batch)(
            cluster_batch=batch,
            matrix=matrix,
            obs_clusters=obs_clusters,
            global_median=global_median,
            global_std=global_std,
            thresholds=thresholds,
            var_index=var_index,
            marker_sets=marker_sets,
            de_lookup=de_lookup,
            de_rank_lookup=de_rank_lookup,
            de_bonus=de_bonus,
            anti_weight=anti_weight,
            idf_weights=idf_weights,
            marker_doc_freq=marker_doc_freq,
            de_commonness_alpha=de_commonness_alpha,
            anti_agg=anti_agg,
            expand_markers=expand_markers,
        )
        for batch in cluster_batches
    )

    score_time = time.time() - score_start
    logger.info("Parallel scoring completed in %.2f sec", score_time)

    # Merge results from all batches
    all_records = []
    all_evidence = [] if expand_markers else None

    for records, evidence in results:
        all_records.extend(records)
        if expand_markers:
            all_evidence.extend(evidence)

    if not all_records:
        logger.warning("No marker scores were computed")
        if expand_markers:
            return pd.DataFrame(), pd.DataFrame()
        return pd.DataFrame()

    marker_scores_df = pd.DataFrame.from_records(all_records)

    total_time = time.time() - start_time
    logger.info(
        "Scoring complete: %d records in %.2f sec (%.1fx speedup estimate)",
        len(marker_scores_df),
        total_time,
        (8 * 60) / max(total_time, 1),  # Compare to ~8 min sequential baseline
    )

    if expand_markers:
        marker_evidence_df = (
            pd.DataFrame.from_records(all_evidence)
            if all_evidence
            else pd.DataFrame()
        )
        logger.info("Generated marker evidence table with %d rows", len(marker_evidence_df))
        return marker_scores_df, marker_evidence_df

    return marker_scores_df
