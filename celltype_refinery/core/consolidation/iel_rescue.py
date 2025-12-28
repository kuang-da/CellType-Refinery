"""Intraepithelial Immune Cell (IEL) Rescue.

This module rescues immune cells that were vetoed by epithelial markers
but are likely legitimate intraepithelial lymphocytes (IELs) or
tissue-resident immune cells.

Biological context:
- Intraepithelial lymphocytes (IELs) are T cells that reside within
  epithelial layers, common in mucosal tissues
- In spatial data, these cells may show epithelial marker signal due to:
  - Physical proximity to epithelial cells
  - Segmentation bleed-through
  - Doublet capture
- These are NOT false positives - they are genuine immune cells in
  an epithelial context

Detection criteria:
1. Cluster was vetoed from Immune Cells root (root_veto_*)
2. Cluster ended up as Unassigned (no other root passed)
3. Strong CD45 signal (pos_frac >= threshold)
4. Epithelial marker present (explains the veto)

Classification:
- If high Lymphoid score: "Intraepithelial Lymphocytes"
- If high Myeloid score: "Tissue-Resident Macrophages"
- Otherwise: "Immune (intraepithelial)"
"""

from dataclasses import dataclass
from enum import Enum
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

try:
    import scanpy as sc
except ImportError:
    sc = None


class IELType(Enum):
    """Classification of intraepithelial immune cells."""

    LYMPHOCYTE = "Intraepithelial Lymphocytes"
    MYELOID = "Tissue-Resident Macrophages"
    GENERIC = "Immune (intraepithelial)"


@dataclass
class IELCandidate:
    """A detected IEL rescue candidate.

    Attributes
    ----------
    cluster_id : str
        The cluster ID
    n_cells : int
        Number of cells in cluster
    cd45_pos_frac : float
        CD45 positive fraction
    veto_marker : str
        Which marker triggered the veto
    veto_marker_pos_frac : float
        Positive fraction of the veto marker
    lymphoid_score : float
        Lymphoid marker score
    myeloid_score : float
        Myeloid marker score
    iel_type : IELType
        Classification type
    final_label : str
        Final rescued label
    reason : str
        Reason for rescue
    """

    cluster_id: str
    n_cells: int
    cd45_pos_frac: float
    veto_marker: str
    veto_marker_pos_frac: float
    lymphoid_score: float = 0.0
    myeloid_score: float = 0.0
    iel_type: IELType = IELType.GENERIC
    final_label: str = ""
    reason: str = ""


@dataclass
class IELRescueConfig:
    """Configuration for IEL rescue.

    Attributes
    ----------
    enabled : bool
        Whether IEL rescue is enabled
    cd45_min_pos_frac : float
        Minimum CD45 positive fraction to consider for rescue
    lymphoid_score_threshold : float
        Minimum lymphoid score to classify as IEL
    myeloid_score_threshold : float
        Minimum myeloid score to classify as tissue-resident macrophage
    suffix : str
        Suffix to add to rescued labels (empty by default)
    cd45_marker : str
        Name of CD45 marker in var_names
    epithelial_markers : List[str]
        List of epithelial marker names to check for veto
    """

    enabled: bool = True
    cd45_min_pos_frac: float = 0.30
    lymphoid_score_threshold: float = 1.0
    myeloid_score_threshold: float = 1.0
    suffix: str = ""
    cd45_marker: str = "CD45"
    epithelial_markers: List[str] = None

    def __post_init__(self):
        if self.epithelial_markers is None:
            self.epithelial_markers = ["E-cadherin", "Pan-Cytokeratin"]


def detect_iel_candidates(
    adata: "sc.AnnData",
    diagnostic_report: pd.DataFrame,
    marker_scores: pd.DataFrame,
    cluster_col: str = "cluster_lvl1",
    layer: Optional[str] = None,
    config: Optional[IELRescueConfig] = None,
) -> List[IELCandidate]:
    """Detect IEL rescue candidates from vetoed Immune clusters.

    Parameters
    ----------
    adata : sc.AnnData
        AnnData object with decision_steps_subcluster in uns
    diagnostic_report : pd.DataFrame
        Diagnostic report with cluster assignments
    marker_scores : pd.DataFrame
        Marker scores for all clusters
    cluster_col : str
        Column name for cluster IDs in adata.obs
    layer : str, optional
        Layer to use for marker expression. If None, uses adata.X
    config : IELRescueConfig, optional
        Configuration for rescue thresholds

    Returns
    -------
    List[IELCandidate]
        List of IEL candidates
    """
    if config is None:
        config = IELRescueConfig()

    if not config.enabled:
        return []

    # Check for decision steps
    if "decision_steps_subcluster" not in adata.uns:
        return []

    steps = pd.DataFrame(adata.uns["decision_steps_subcluster"])

    # Find vetoed clusters
    veto_steps = steps[steps["fail_reason"].str.contains("root_veto", na=False)]
    if len(veto_steps) == 0:
        return []

    vetoed_info = {}
    for _, row in veto_steps.iterrows():
        cid = str(row["cluster_id"])
        reason = row["fail_reason"]
        # Extract veto marker from reason (e.g., "root_veto_Pan-Cytokeratin")
        veto_marker = reason.replace("root_veto_", "")
        vetoed_info[cid] = veto_marker

    # Find Unassigned clusters
    unassigned_ids = set(
        diagnostic_report[diagnostic_report["assigned_label"] == "Unassigned"]["cluster_id"].astype(str)
    )

    # Find vetoed clusters that are Unassigned
    vetoed_unassigned = set(vetoed_info.keys()) & unassigned_ids
    if not vetoed_unassigned:
        return []

    # Get marker indices
    var_names = list(adata.var_names)
    var_index = {name: i for i, name in enumerate(var_names)}

    cd45_idx = var_index.get(config.cd45_marker)
    if cd45_idx is None:
        return []

    # Get epithelial marker indices
    epithelial_indices = {}
    for marker in config.epithelial_markers:
        if marker in var_index:
            epithelial_indices[marker] = var_index[marker]

    # Get expression matrix
    if layer and layer in adata.layers:
        matrix = adata.layers[layer]
    else:
        matrix = adata.X

    # Compute global Q75 threshold for CD45
    cd45_vals = matrix[:, cd45_idx]
    if hasattr(cd45_vals, 'toarray'):
        cd45_vals = cd45_vals.toarray().flatten()
    cd45_threshold = float(np.nanquantile(cd45_vals, 0.75))

    # Similarly for epithelial markers
    epithelial_thresholds = {}
    for marker, idx in epithelial_indices.items():
        vals = matrix[:, idx]
        if hasattr(vals, 'toarray'):
            vals = vals.toarray().flatten()
        epithelial_thresholds[marker] = float(np.nanquantile(vals, 0.75))

    # Get cluster assignments
    clusters = adata.obs[cluster_col].astype(str)

    candidates = []

    for cid in vetoed_unassigned:
        mask = clusters == cid
        n_cells = int(mask.sum())
        if n_cells == 0:
            continue

        cluster_matrix = matrix[mask.to_numpy(), :]
        if hasattr(cluster_matrix, 'toarray'):
            cluster_matrix = cluster_matrix.toarray()

        # Compute CD45 pos_frac
        cd45_cluster = cluster_matrix[:, cd45_idx]
        cd45_pos_frac = float(np.mean(cd45_cluster >= cd45_threshold))

        # Check CD45 threshold
        if cd45_pos_frac < config.cd45_min_pos_frac:
            continue

        # Get veto marker stats
        veto_marker = vetoed_info[cid]
        veto_pos_frac = 0.0
        if veto_marker in epithelial_indices:
            veto_idx = epithelial_indices[veto_marker]
            veto_thresh = epithelial_thresholds[veto_marker]
            veto_cluster = cluster_matrix[:, veto_idx]
            veto_pos_frac = float(np.mean(veto_cluster >= veto_thresh))

        # Get lymphoid and myeloid scores from marker_scores
        cid_scores = marker_scores[marker_scores["cluster_id"].astype(str) == cid]
        lymphoid_score = 0.0
        myeloid_score = 0.0

        lymphoid_row = cid_scores[cid_scores["label"] == "Lymphoids"]
        if len(lymphoid_row) > 0:
            lymphoid_score = float(lymphoid_row["score"].iloc[0])

        myeloid_row = cid_scores[cid_scores["label"] == "Myeloids"]
        if len(myeloid_row) > 0:
            myeloid_score = float(myeloid_row["score"].iloc[0])

        # Classify IEL type
        if lymphoid_score >= config.lymphoid_score_threshold:
            iel_type = IELType.LYMPHOCYTE
        elif myeloid_score >= config.myeloid_score_threshold:
            iel_type = IELType.MYELOID
        else:
            iel_type = IELType.GENERIC

        # Create candidate
        candidate = IELCandidate(
            cluster_id=cid,
            n_cells=n_cells,
            cd45_pos_frac=cd45_pos_frac,
            veto_marker=veto_marker,
            veto_marker_pos_frac=veto_pos_frac,
            lymphoid_score=lymphoid_score,
            myeloid_score=myeloid_score,
            iel_type=iel_type,
        )
        candidates.append(candidate)

    return candidates


def apply_iel_rescue(
    candidates: List[IELCandidate],
    config: Optional[IELRescueConfig] = None,
) -> List[IELCandidate]:
    """Apply rescue logic to IEL candidates.

    Parameters
    ----------
    candidates : List[IELCandidate]
        List of IEL candidates
    config : IELRescueConfig, optional
        Configuration

    Returns
    -------
    List[IELCandidate]
        Updated candidates with final_label and reason
    """
    if config is None:
        config = IELRescueConfig()

    suffix = f" {config.suffix}".rstrip() if config.suffix else ""

    for candidate in candidates:
        base_label = candidate.iel_type.value
        candidate.final_label = f"{base_label}{suffix}"

        candidate.reason = (
            f"IEL rescue: CD45={candidate.cd45_pos_frac:.2f}, "
            f"vetoed by {candidate.veto_marker}={candidate.veto_marker_pos_frac:.2f}, "
            f"lymphoid={candidate.lymphoid_score:.2f}, myeloid={candidate.myeloid_score:.2f}"
        )

    return candidates


def get_iel_rescue_map(
    candidates: List[IELCandidate],
) -> Dict[str, Tuple[str, str]]:
    """Get mapping from cluster_id to (final_label, reason) for IEL rescues.

    Parameters
    ----------
    candidates : List[IELCandidate]
        List of IEL candidates (after apply_iel_rescue)

    Returns
    -------
    Dict[str, Tuple[str, str]]
        Mapping: cluster_id -> (final_label, reason)
    """
    return {c.cluster_id: (c.final_label, c.reason) for c in candidates}


def summarize_iel_candidates(candidates: List[IELCandidate]) -> pd.DataFrame:
    """Create summary DataFrame of IEL candidates.

    Parameters
    ----------
    candidates : List[IELCandidate]
        List of IEL candidates

    Returns
    -------
    pd.DataFrame
        Summary with columns for each attribute
    """
    if not candidates:
        return pd.DataFrame(columns=[
            "cluster_id", "n_cells", "cd45_pos_frac", "veto_marker",
            "veto_marker_pos_frac", "lymphoid_score", "myeloid_score",
            "iel_type", "final_label", "reason"
        ])

    data = []
    for c in candidates:
        data.append({
            "cluster_id": c.cluster_id,
            "n_cells": c.n_cells,
            "cd45_pos_frac": c.cd45_pos_frac,
            "veto_marker": c.veto_marker,
            "veto_marker_pos_frac": c.veto_marker_pos_frac,
            "lymphoid_score": c.lymphoid_score,
            "myeloid_score": c.myeloid_score,
            "iel_type": c.iel_type.value,
            "final_label": c.final_label,
            "reason": c.reason,
        })

    return pd.DataFrame(data).sort_values("n_cells", ascending=False)


def summarize_iel_by_type(candidates: List[IELCandidate]) -> pd.DataFrame:
    """Create summary by IEL type.

    Parameters
    ----------
    candidates : List[IELCandidate]
        List of IEL candidates

    Returns
    -------
    pd.DataFrame
        Summary by IEL type
    """
    if not candidates:
        return pd.DataFrame(columns=["iel_type", "n_clusters", "n_cells"])

    df = summarize_iel_candidates(candidates)
    summary = df.groupby("iel_type").agg(
        n_clusters=("cluster_id", "count"),
        n_cells=("n_cells", "sum"),
        mean_cd45=("cd45_pos_frac", "mean"),
        mean_lymphoid=("lymphoid_score", "mean"),
        mean_myeloid=("myeloid_score", "mean"),
    ).reset_index()

    return summary.sort_values("n_cells", ascending=False)
