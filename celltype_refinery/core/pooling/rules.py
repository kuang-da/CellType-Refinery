"""
Pooling rules for Stage J lineage consolidation.

This module provides:
- Root type extraction from marker map
- Ambiguous label canonicalization
- Cluster classification for pooling
- Pool mapping construction
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, Optional, Set, Tuple

import pandas as pd


def load_root_types_from_marker_map(marker_map_path: Path) -> Set[str]:
    """Extract root cell types from marker map JSON.

    Root types are top-level keys in the marker map hierarchy,
    excluding metadata keys (starting with "_").

    Parameters
    ----------
    marker_map_path : Path
        Path to marker map JSON file

    Returns
    -------
    Set[str]
        Set of root type names (e.g., {"Epithelium", "Endothelium", ...})
    """
    marker_map_path = Path(marker_map_path)

    with open(marker_map_path) as f:
        marker_map = json.load(f)

    # Top-level keys excluding metadata (keys starting with "_")
    root_types = {k for k in marker_map.keys() if not k.startswith("_")}

    # Always include "Unassigned"
    root_types.add("Unassigned")

    return root_types


def canonicalize_ambiguous_label(label: str) -> str:
    """Canonicalize ambiguous labels by alphabetical ordering.

    Ambiguous labels use "~" to separate multiple types.
    This function ensures consistent ordering for grouping.

    Parameters
    ----------
    label : str
        Label to canonicalize (e.g., "Epithelium~Endothelium")

    Returns
    -------
    str
        Canonicalized label with components in alphabetical order

    Examples
    --------
    >>> canonicalize_ambiguous_label("Epithelium")
    "Epithelium"
    >>> canonicalize_ambiguous_label("Epithelium~Endothelium")
    "Endothelium~Epithelium"
    >>> canonicalize_ambiguous_label("C~A~B")
    "A~B~C"
    """
    if "~" not in label:
        return label

    parts = label.split("~")
    return "~".join(sorted(parts))


def classify_cluster_for_pooling(
    assigned_label: str,
    root_types: Set[str],
    recommendation: str,
) -> Tuple[Optional[str], str]:
    """Classify a cluster for pooling based on relabeling rules.

    Rules are applied in order:
    1. Root cell types -> Pool_<RootType>
    2. Ambiguous types (~) -> Pool_<CanonicalLabel>
    3. SUBCLUSTER carry-over -> Not pooled (keep label)
    4. Default (subtypes) -> Not pooled

    Parameters
    ----------
    assigned_label : str
        Current cell type label for the cluster
    root_types : Set[str]
        Set of root type names from marker map + "Unassigned"
    recommendation : str
        Recommendation from diagnostic report (SUBCLUSTER, RELABEL, SKIP)

    Returns
    -------
    Tuple[Optional[str], str]
        (pool_label, reason) where pool_label is:
        - "Pool_<Name>" if cluster should be pooled
        - None if cluster should not be pooled
    """
    # Rule 1: Root cell types -> Pool
    if assigned_label in root_types:
        return f"Pool_{assigned_label}", f"Root type pooled: {assigned_label}"

    # Rule 1b: Unassigned -> Pool (case-insensitive check)
    if assigned_label.lower() == "unassigned":
        return "Pool_Unassigned", "Unassigned cells pooled"

    # Rule 2: Ambiguous types (~) -> Pool with canonical name
    if "~" in assigned_label:
        canonical = canonicalize_ambiguous_label(assigned_label)
        return f"Pool_{canonical}", f"Ambiguous type pooled: {assigned_label} -> {canonical}"

    # Rule 3 & default: Not pooled
    # Subtypes (e.g., "Ciliated Epithelium", "Smooth Muscle Cells") are not pooled
    # Clusters marked SUBCLUSTER but not caught above are carried over
    if recommendation == "SUBCLUSTER":
        return None, f"SUBCLUSTER carry-over: {assigned_label}"
    else:
        return None, f"Subtype retained: {assigned_label}"


def build_pool_mapping(
    cluster_annotations: pd.DataFrame,
    diagnostic_report: pd.DataFrame,
    root_types: Set[str],
) -> Dict[str, Tuple[Optional[str], str, str]]:
    """Build mapping of cluster_id -> (pool_label, new_cluster_id, reason).

    Parameters
    ----------
    cluster_annotations : pd.DataFrame
        Cluster annotations from Stage I (must have cluster_id, assigned_label)
    diagnostic_report : pd.DataFrame
        Diagnostic report from Stage I (must have cluster_id, recommendation)
    root_types : Set[str]
        Set of root type names (from marker map + "Unassigned")

    Returns
    -------
    Dict[str, Tuple[Optional[str], str, str]]
        Mapping of cluster_id -> (pool_label, new_cluster_id, reason)
        - pool_label: "Pool_<Name>" or None if not pooled
        - new_cluster_id: "P<n>" for pooled, or original_id if not pooled
        - reason: Human-readable reason for the decision
    """
    pool_mapping: Dict[str, Tuple[Optional[str], str, str]] = {}
    pool_id_counter: Dict[str, int] = {}  # pool_label -> assigned P<n> index

    # Merge annotations with diagnostic report to get recommendations
    merged = cluster_annotations.merge(
        diagnostic_report[["cluster_id", "recommendation"]].astype(str),
        on="cluster_id",
        how="left",
    )
    merged["recommendation"] = merged["recommendation"].fillna("SKIP")

    for _, row in merged.iterrows():
        cluster_id = str(row["cluster_id"])
        assigned_label = str(row["assigned_label"])
        recommendation = str(row["recommendation"])

        pool_label, reason = classify_cluster_for_pooling(
            assigned_label, root_types, recommendation
        )

        if pool_label is not None:
            # Assign pool cluster ID: P0, P1, P2, ...
            # Each unique pool_label gets its own P<n>
            if pool_label not in pool_id_counter:
                pool_id_counter[pool_label] = len(pool_id_counter)
            new_cluster_id = f"P{pool_id_counter[pool_label]}"
            pool_mapping[cluster_id] = (pool_label, new_cluster_id, reason)
        else:
            # Keep original cluster ID
            pool_mapping[cluster_id] = (None, cluster_id, reason)

    return pool_mapping


def get_pool_summary(
    pool_mapping: Dict[str, Tuple[Optional[str], str, str]],
    cluster_annotations: pd.DataFrame,
) -> pd.DataFrame:
    """Generate summary DataFrame of pooling decisions.

    Parameters
    ----------
    pool_mapping : Dict
        Pool mapping from build_pool_mapping()
    cluster_annotations : pd.DataFrame
        Cluster annotations with n_cells column

    Returns
    -------
    pd.DataFrame
        Summary with columns: cluster_id, original_label, pool_label,
        new_cluster_id, n_cells, is_pooled, reason
    """
    rows = []

    # Create lookup for n_cells
    n_cells_lookup = cluster_annotations.set_index(
        cluster_annotations["cluster_id"].astype(str)
    )["n_cells"].to_dict()

    for cluster_id, (pool_label, new_cluster_id, reason) in pool_mapping.items():
        # Get original label from annotations
        mask = cluster_annotations["cluster_id"].astype(str) == cluster_id
        if mask.any():
            original_label = cluster_annotations.loc[mask, "assigned_label"].iloc[0]
        else:
            original_label = "Unknown"

        n_cells = n_cells_lookup.get(cluster_id, 0)

        rows.append(
            {
                "cluster_id": cluster_id,
                "original_label": original_label,
                "pool_label": pool_label if pool_label else "(unchanged)",
                "new_cluster_id": new_cluster_id,
                "n_cells": n_cells,
                "is_pooled": pool_label is not None,
                "reason": reason,
            }
        )

    return pd.DataFrame(rows)
