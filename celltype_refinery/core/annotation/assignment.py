"""Cell-level annotation assignment.

This module handles mapping cluster-level annotations to individual cells
in the AnnData object.
"""

from __future__ import annotations

import logging
from typing import Optional

import pandas as pd

try:
    import scanpy as sc
except ImportError:
    sc = None


def annotate_obs(
    adata: "sc.AnnData",
    cluster_annotations: pd.DataFrame,
    cluster_key: str,
    label_col: str,
    logger: Optional[logging.Logger] = None,
) -> None:
    """Map cluster-level annotations to cells in adata.obs.

    Creates new columns in adata.obs:
    - {label_col}: Assigned cell type label
    - {label_col}_score: Confidence score
    - {label_col}_root: Root-level category (if available)

    Args:
        adata: AnnData object to modify in place
        cluster_annotations: DataFrame with columns:
            - cluster_id: Cluster identifier
            - assigned_label: Cell type label
            - assigned_score: Confidence score
            - root_label: (optional) Root-level category
        cluster_key: Column in adata.obs with cluster assignments
        label_col: Name for the output label column
        logger: Optional logger instance
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    # Create mappings
    mapping = dict(
        zip(
            cluster_annotations["cluster_id"].astype(str),
            cluster_annotations["assigned_label"],
        )
    )
    score_mapping = dict(
        zip(
            cluster_annotations["cluster_id"].astype(str),
            cluster_annotations["assigned_score"],
        )
    )

    # Apply label mapping
    adata.obs[label_col] = (
        adata.obs[cluster_key].astype(str).map(mapping).fillna("Unknown")
    )

    # Apply score mapping
    adata.obs[f"{label_col}_score"] = (
        adata.obs[cluster_key].astype(str).map(score_mapping).fillna(0.0)
    )

    # Add root_label if available
    if "root_label" in cluster_annotations.columns:
        root_mapping = dict(
            zip(
                cluster_annotations["cluster_id"].astype(str),
                cluster_annotations["root_label"],
            )
        )
        adata.obs[f"{label_col}_root"] = (
            adata.obs[cluster_key].astype(str).map(root_mapping).fillna("Unassigned")
        )

    logger.info(
        "Annotated %d cells with %d unique labels",
        len(adata),
        len(adata.obs[label_col].unique()),
    )
