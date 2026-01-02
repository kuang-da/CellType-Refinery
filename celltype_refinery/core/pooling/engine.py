"""
Pooling engine for Stage J lineage consolidation.

This module provides:
- PoolingEngine: Core class for executing pooling operations
- PoolingResult: Result dataclass from pooling execution
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd

try:
    import scanpy as sc
except ImportError:
    sc = None

from .rules import (
    build_pool_mapping,
    load_root_types_from_marker_map,
)


@dataclass
class PoolingResult:
    """Result from pooling execution."""

    success: bool
    n_pools_created: int = 0
    n_cells_pooled: int = 0
    n_cells_unchanged: int = 0
    n_clusters_pooled: int = 0
    n_clusters_unchanged: int = 0
    pool_summary: pd.DataFrame = field(default_factory=pd.DataFrame)
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)


class PoolingEngine:
    """Engine for executing Stage J pooling operations.

    Pools clusters by lineage based on relabeling rules:
    1. Root cell types -> Pool_<RootType>
    2. Ambiguous types (~) -> Pool_<CanonicalLabel>
    3. SUBCLUSTER carry-over -> Keep original label

    Parameters
    ----------
    marker_map_path : Path
        Path to marker map JSON file for root type extraction
    cluster_col : str
        Column name for cluster IDs (default: cluster_lvl1)
    label_col : str
        Column name for cell type labels (default: cell_type_lvl1)
    logger : logging.Logger, optional
        Logger instance
    """

    def __init__(
        self,
        marker_map_path: Path,
        cluster_col: str = "cluster_lvl1",
        label_col: str = "cell_type_lvl1",
        logger: Optional[logging.Logger] = None,
    ):
        self.marker_map_path = Path(marker_map_path)
        self.cluster_col = cluster_col
        self.label_col = label_col
        self.logger = logger or logging.getLogger(__name__)

        # Load root types from marker map
        self.root_types: Set[str] = load_root_types_from_marker_map(self.marker_map_path)
        self.root_types.add("Unassigned")  # Always include Unassigned
        self.logger.info(f"Root types loaded: {sorted(self.root_types)}")

    def execute(
        self,
        adata: "sc.AnnData",
        cluster_annotations: pd.DataFrame,
        diagnostic_report: pd.DataFrame,
    ) -> PoolingResult:
        """Execute pooling on AnnData.

        Modifies adata in place, updating cluster_col and label_col
        for pooled cells.

        Parameters
        ----------
        adata : sc.AnnData
            AnnData to modify in place
        cluster_annotations : pd.DataFrame
            Cluster annotations from Stage I (must have cluster_id, assigned_label, n_cells)
        diagnostic_report : pd.DataFrame
            Diagnostic report from Stage I (must have cluster_id, recommendation)

        Returns
        -------
        PoolingResult
            Result with execution details and summary
        """
        self.logger.info("=" * 70)
        self.logger.info("[STAGE J] Executing lineage pooling...")
        self.logger.info("=" * 70)

        errors: List[str] = []
        warnings: List[str] = []

        # Validate inputs
        if self.cluster_col not in adata.obs:
            errors.append(f"Column '{self.cluster_col}' not found in adata.obs")
            return PoolingResult(success=False, errors=errors)

        if self.label_col not in adata.obs:
            errors.append(f"Column '{self.label_col}' not found in adata.obs")
            return PoolingResult(success=False, errors=errors)

        # Build pool mapping
        pool_mapping = build_pool_mapping(
            cluster_annotations, diagnostic_report, self.root_types
        )

        # Build cluster_id -> n_cells lookup from annotations
        cluster_cell_counts = dict(
            zip(
                cluster_annotations["cluster_id"].astype(str),
                cluster_annotations["n_cells"].astype(int),
            )
        )

        # Count pools vs unchanged
        pools_created: Dict[str, Dict] = {}
        clusters_pooled = []
        clusters_unchanged = []

        for cid, (pool_label, new_cid, reason) in pool_mapping.items():
            n_cells = cluster_cell_counts.get(cid, 0)
            if pool_label is not None:
                if pool_label not in pools_created:
                    pools_created[pool_label] = {"cluster_ids": [], "new_cluster_id": new_cid, "n_cells": 0}
                pools_created[pool_label]["cluster_ids"].append(cid)
                pools_created[pool_label]["n_cells"] += n_cells
                clusters_pooled.append(cid)
            else:
                clusters_unchanged.append(cid)

        self.logger.info(f"Pools to create: {len(pools_created)}")
        for pool_label, info in sorted(pools_created.items()):
            self.logger.info(f"  {pool_label} ({info['new_cluster_id']}): {len(info['cluster_ids'])} clusters, {info['n_cells']:,} cells")
        self.logger.info(f"Clusters unchanged: {len(clusters_unchanged)}")

        # Apply pooling to AnnData using vectorized operations
        self.logger.info("Applying pooling to AnnData...")

        # Build mapping dictionaries for vectorized replace
        cluster_id_map = {}  # old_cluster_id -> new_cluster_id (for pooled only)
        label_map = {}  # old_cluster_id -> new_label (for pooled only)

        for old_cid, (pool_label, new_cid, reason) in pool_mapping.items():
            if pool_label is not None:
                cluster_id_map[old_cid] = new_cid
                label_map[old_cid] = pool_label

        # Convert to string for mapping
        old_clusters = adata.obs[self.cluster_col].astype(str)
        old_labels = adata.obs[self.label_col].astype(str)

        # Vectorized update using map
        new_clusters = old_clusters.map(cluster_id_map).fillna(old_clusters)
        new_labels = old_clusters.map(label_map).fillna(old_labels)

        adata.obs[self.cluster_col] = new_clusters.astype("category")
        adata.obs[self.label_col] = new_labels.astype("category")

        # Count cells
        n_cells_pooled = int((old_clusters.isin(cluster_id_map.keys())).sum())
        n_cells_unchanged = len(adata) - n_cells_pooled

        self.logger.info(f"Cells pooled: {n_cells_pooled:,}")
        self.logger.info(f"Cells unchanged: {n_cells_unchanged:,}")

        # Build pool summary DataFrame
        summary_rows = []
        for old_cid, (pool_label, new_cid, reason) in pool_mapping.items():
            ann_row = cluster_annotations[
                cluster_annotations["cluster_id"].astype(str) == old_cid
            ]
            n_cells = int(ann_row["n_cells"].iloc[0]) if len(ann_row) > 0 else 0
            original_label = str(ann_row["assigned_label"].iloc[0]) if len(ann_row) > 0 else ""

            summary_rows.append(
                {
                    "cluster_id": old_cid,
                    "original_label": original_label,
                    "pool_label": pool_label if pool_label else "(unchanged)",
                    "new_cluster_id": new_cid,
                    "n_cells": n_cells,
                    "is_pooled": pool_label is not None,
                    "reason": reason,
                }
            )

        pool_summary = pd.DataFrame(summary_rows)

        return PoolingResult(
            success=True,
            n_pools_created=len(pools_created),
            n_cells_pooled=n_cells_pooled,
            n_cells_unchanged=n_cells_unchanged,
            n_clusters_pooled=len(clusters_pooled),
            n_clusters_unchanged=len(clusters_unchanged),
            pool_summary=pool_summary,
            errors=errors,
            warnings=warnings,
        )

    def generate_diagnostic_report(
        self,
        adata: "sc.AnnData",
        original_diagnostic: pd.DataFrame,
        pool_summary: pd.DataFrame,
        input_cluster_annotations: Optional[pd.DataFrame] = None,
    ) -> pd.DataFrame:
        """Generate diagnostic report for next Stage I iteration.

        All pools are marked SUBCLUSTER.
        Non-pooled clusters carry over original recommendation and scores.

        Parameters
        ----------
        adata : sc.AnnData
            AnnData after pooling
        original_diagnostic : pd.DataFrame
            Original diagnostic report from Stage I
        pool_summary : pd.DataFrame
            Pool summary from execute()
        input_cluster_annotations : pd.DataFrame, optional
            Original cluster annotations from Stage I input

        Returns
        -------
        pd.DataFrame
            Diagnostic report compatible with Stage I
        """
        # Use vectorized groupby for speed
        obs_subset = adata.obs[[self.cluster_col, self.label_col]].copy()
        obs_subset[self.cluster_col] = obs_subset[self.cluster_col].astype(str)
        obs_subset[self.label_col] = obs_subset[self.label_col].astype(str)

        grouped = obs_subset.groupby(self.cluster_col, observed=True)
        cluster_stats = grouped.agg(
            assigned_label=(self.label_col, lambda x: x.mode().iloc[0] if len(x) > 0 else "Unknown"),
            n_cells=(self.label_col, "size"),
        ).reset_index().rename(columns={self.cluster_col: "cluster_id"})

        # Build original diagnostic lookup
        orig_diag_lookup = {}
        for _, row in original_diagnostic.iterrows():
            cid = str(row["cluster_id"])
            orig_diag_lookup[cid] = {
                "recommendation": str(row.get("recommendation", "SKIP")),
                "recommendation_reason": str(row.get("recommendation_reason", row.get("reason", ""))),
                "confidence_band": str(row.get("confidence_band", "medium")),
                "assigned_score": float(row.get("assigned_score", 0.0) or 0.0),
            }

        # Build input annotations lookup for score information
        input_ann_lookup = {}
        if input_cluster_annotations is not None and not input_cluster_annotations.empty:
            for _, row in input_cluster_annotations.iterrows():
                cid = str(row.get("cluster_id", ""))
                if cid:
                    input_ann_lookup[cid] = {
                        "assigned_score": float(row.get("assigned_score", 0.0) or 0.0),
                        "n_cells": int(row.get("n_cells", 0) or 0),
                    }

        # Build pool score lookup: compute weighted average for each pool
        pool_score_lookup = self._compute_pool_scores(pool_summary, input_ann_lookup)

        # Add recommendation columns
        rows = []
        for _, row in cluster_stats.iterrows():
            cluster_id = row["cluster_id"]
            if cluster_id.startswith("P"):
                recommendation = "SUBCLUSTER"
                recommendation_reason = "Pooled cluster - subclustering recommended"
                confidence_band = "low"
                assigned_score = pool_score_lookup.get(cluster_id, 0.0)
            else:
                orig = orig_diag_lookup.get(cluster_id, {})
                recommendation = orig.get("recommendation", "SKIP")
                recommendation_reason = orig.get("recommendation_reason", "No previous recommendation")
                confidence_band = orig.get("confidence_band", "medium")
                if cluster_id in input_ann_lookup:
                    assigned_score = input_ann_lookup[cluster_id].get("assigned_score", 0.0)
                else:
                    assigned_score = orig.get("assigned_score", 0.0)

            rows.append({
                "cluster_id": cluster_id,
                "assigned_label": row["assigned_label"],
                "assigned_score": assigned_score,
                "n_cells": row["n_cells"],
                "recommendation": recommendation,
                "recommendation_reason": recommendation_reason,
                "confidence_band": confidence_band,
                "is_pool": cluster_id.startswith("P"),
            })

        result = pd.DataFrame(rows)
        result = self._sort_clusters(result)
        return result

    def generate_cluster_annotations(
        self,
        adata: "sc.AnnData",
        pool_summary: pd.DataFrame,
        input_cluster_annotations: Optional[pd.DataFrame] = None,
    ) -> pd.DataFrame:
        """Generate cluster annotations for Stage I compatibility.

        Parameters
        ----------
        adata : sc.AnnData
            AnnData after pooling
        pool_summary : pd.DataFrame
            Pool summary from execute()
        input_cluster_annotations : pd.DataFrame, optional
            Original cluster annotations from Stage I input

        Returns
        -------
        pd.DataFrame
            Cluster annotations compatible with Stage I
        """
        # Use vectorized groupby for speed
        obs_subset = adata.obs[[self.cluster_col, self.label_col]].copy()
        obs_subset[self.cluster_col] = obs_subset[self.cluster_col].astype(str)
        obs_subset[self.label_col] = obs_subset[self.label_col].astype(str)

        grouped = obs_subset.groupby(self.cluster_col, observed=True)
        cluster_stats = grouped.agg(
            assigned_label=(self.label_col, lambda x: x.mode().iloc[0] if len(x) > 0 else "Unknown"),
            n_cells=(self.label_col, "size"),
        ).reset_index().rename(columns={self.cluster_col: "cluster_id"})

        total_cells = len(adata)

        # Build pool origin lookup
        pool_origins = {}
        if not pool_summary.empty and "new_cluster_id" in pool_summary.columns:
            for new_cid in pool_summary["new_cluster_id"].unique():
                source_rows = pool_summary[pool_summary["new_cluster_id"] == new_cid]
                pool_origins[new_cid] = ";".join(source_rows["cluster_id"].astype(str).tolist())

        # Build input annotations lookup
        input_lookup = {}
        if input_cluster_annotations is not None and not input_cluster_annotations.empty:
            for _, row in input_cluster_annotations.iterrows():
                cid = str(row.get("cluster_id", ""))
                if cid:
                    input_lookup[cid] = row.to_dict()

        # Compute pool scores
        pool_score_lookup = self._compute_pool_scores(
            pool_summary,
            {cid: {"assigned_score": d.get("assigned_score", 0.0), "n_cells": d.get("n_cells", 0)}
             for cid, d in input_lookup.items()}
        )

        # Build rows
        rows = []
        for _, row in cluster_stats.iterrows():
            cluster_id = row["cluster_id"]
            n_cells = row["n_cells"]
            assigned_label = row["assigned_label"]
            proportion = n_cells / total_cells if total_cells > 0 else 0.0

            if cluster_id.startswith("P"):
                # Pool row
                root_label = assigned_label
                if assigned_label.startswith("Pool_"):
                    root_label = assigned_label[5:]
                    if "~" in root_label:
                        root_label = root_label.split("~")[0]

                rows.append({
                    "cluster_id": cluster_id,
                    "origin_cluster": cluster_id,
                    "iteration_created": 0,
                    "proportion": proportion,
                    "n_cells": n_cells,
                    "assigned_label": assigned_label,
                    "assigned_path": assigned_label,
                    "assigned_level": 0,
                    "assigned_score": pool_score_lookup.get(cluster_id, 0.0),
                    "root_label": root_label,
                    "confidence": 0.0,
                    "confidence_level": "pool",
                    "is_ambiguous_root": "~" in assigned_label,
                })
            else:
                # Unchanged row
                src = input_lookup.get(cluster_id, {})
                rows.append({
                    "cluster_id": cluster_id,
                    "origin_cluster": str(src.get("origin_cluster", cluster_id)),
                    "iteration_created": int(src.get("iteration_created", 0) or 0),
                    "proportion": proportion,
                    "n_cells": n_cells,
                    "assigned_label": assigned_label,
                    "assigned_path": str(src.get("assigned_path", assigned_label)),
                    "assigned_level": int(src.get("assigned_level", 0) or 0),
                    "assigned_score": float(src.get("assigned_score", 0.0) or 0.0),
                    "root_label": str(src.get("root_label", "")),
                    "confidence": float(src.get("confidence", 0.0) or 0.0),
                    "confidence_level": str(src.get("confidence_level", "medium")),
                    "is_ambiguous_root": bool(src.get("is_ambiguous_root", False)),
                })

        result = pd.DataFrame(rows)
        result = self._sort_clusters(result)
        return result

    def generate_marker_scores(
        self,
        pool_summary: pd.DataFrame,
        input_marker_scores: pd.DataFrame,
    ) -> pd.DataFrame:
        """Generate marker_scores.csv for Stage I smart rescore compatibility.

        For unchanged clusters: copy rows directly from input marker_scores.
        For pools: rows are NOT included (they need fresh scoring in Stage I).

        Parameters
        ----------
        pool_summary : pd.DataFrame
            Pool summary from execute()
        input_marker_scores : pd.DataFrame
            marker_scores.csv from the input Stage I iteration

        Returns
        -------
        pd.DataFrame
            marker_scores.csv for unchanged clusters only
        """
        if input_marker_scores is None or input_marker_scores.empty:
            self.logger.warning("No input marker_scores provided, returning empty DataFrame")
            return pd.DataFrame()

        # Get list of unchanged clusters (not pooled)
        unchanged_clusters = set()
        if not pool_summary.empty:
            unchanged_rows = pool_summary[pool_summary["is_pooled"] == False]
            unchanged_clusters = set(unchanged_rows["cluster_id"].astype(str).tolist())

        self.logger.info(f"Propagating marker_scores for {len(unchanged_clusters)} unchanged clusters")

        # Filter input marker_scores to only include unchanged clusters
        input_marker_scores = input_marker_scores.copy()
        input_marker_scores["cluster_id"] = input_marker_scores["cluster_id"].astype(str)

        result = input_marker_scores[
            input_marker_scores["cluster_id"].isin(unchanged_clusters)
        ].copy()

        self.logger.info(f"Propagated {len(result)} marker_score rows")

        return result

    def _compute_pool_scores(
        self,
        pool_summary: pd.DataFrame,
        input_ann_lookup: Dict[str, Dict],
    ) -> Dict[str, float]:
        """Compute weighted average scores for each pool.

        Parameters
        ----------
        pool_summary : pd.DataFrame
            Pool summary from execute()
        input_ann_lookup : Dict[str, Dict]
            Lookup of cluster_id -> {assigned_score, n_cells}

        Returns
        -------
        Dict[str, float]
            Mapping of pool cluster_id -> weighted average score
        """
        pool_score_lookup = {}
        if not pool_summary.empty and "new_cluster_id" in pool_summary.columns:
            for new_cid in pool_summary["new_cluster_id"].unique():
                if not str(new_cid).startswith("P"):
                    continue
                source_rows = pool_summary[pool_summary["new_cluster_id"] == new_cid]
                source_ids = source_rows["cluster_id"].astype(str).tolist()

                total_weighted_score = 0.0
                total_cells = 0
                for src_id in source_ids:
                    src_data = input_ann_lookup.get(src_id, {})
                    src_score = float(src_data.get("assigned_score", 0.0) or 0.0)
                    src_cells = int(src_data.get("n_cells", 0) or 0)
                    total_weighted_score += src_score * src_cells
                    total_cells += src_cells

                pool_score = total_weighted_score / total_cells if total_cells > 0 else 0.0
                pool_score_lookup[new_cid] = pool_score

        return pool_score_lookup

    def _sort_clusters(self, df: pd.DataFrame) -> pd.DataFrame:
        """Sort DataFrame with pools first, then by cluster_id."""
        df["_sort_key"] = df["cluster_id"].apply(lambda x: (not str(x).startswith("P"), str(x)))
        result = df.sort_values("_sort_key").drop(columns=["_sort_key"])
        return result

    def generate_full_diagnostic_report(
        self,
        adata: "sc.AnnData",
        original_diagnostic: pd.DataFrame,
        pool_summary: pd.DataFrame,
        input_cluster_annotations: Optional[pd.DataFrame] = None,
    ) -> pd.DataFrame:
        """Generate diagnostic report with ALL Stage I columns (27 columns).

        This is an enhanced version that includes all columns expected by Stage I
        for proper diagnostic continuation.

        Parameters
        ----------
        adata : sc.AnnData
            AnnData after pooling
        original_diagnostic : pd.DataFrame
            Original diagnostic report from Stage I
        pool_summary : pd.DataFrame
            Pool summary from execute()
        input_cluster_annotations : pd.DataFrame, optional
            Original cluster annotations from Stage I input

        Returns
        -------
        pd.DataFrame
            Full diagnostic report compatible with Stage I (27 columns)
        """
        # Use vectorized groupby for speed
        obs_subset = adata.obs[[self.cluster_col, self.label_col]].copy()
        obs_subset[self.cluster_col] = obs_subset[self.cluster_col].astype(str)
        obs_subset[self.label_col] = obs_subset[self.label_col].astype(str)

        grouped = obs_subset.groupby(self.cluster_col, observed=True)
        cluster_stats = grouped.agg(
            assigned_label=(self.label_col, lambda x: x.mode().iloc[0] if len(x) > 0 else "Unknown"),
            n_cells=(self.label_col, "size"),
        ).reset_index().rename(columns={self.cluster_col: "cluster_id"})

        # Build lookup from original diagnostic
        orig_diag_lookup = {}
        if not original_diagnostic.empty:
            original_diagnostic = original_diagnostic.copy()
            original_diagnostic["cluster_id"] = original_diagnostic["cluster_id"].astype(str)
            for _, row in original_diagnostic.iterrows():
                orig_diag_lookup[row["cluster_id"]] = row.to_dict()

        # Build input annotations lookup for score information
        input_ann_lookup = {}
        if input_cluster_annotations is not None and not input_cluster_annotations.empty:
            for _, row in input_cluster_annotations.iterrows():
                cid = str(row.get("cluster_id", ""))
                if cid:
                    input_ann_lookup[cid] = row.to_dict()

        # Build pool score lookup
        pool_score_lookup = self._compute_pool_scores(
            pool_summary,
            {cid: {"assigned_score": d.get("assigned_score", 0.0), "n_cells": d.get("n_cells", 0)}
             for cid, d in input_ann_lookup.items()}
        )

        # Build rows with all Stage I columns
        rows = []
        for _, row in cluster_stats.iterrows():
            cluster_id = row["cluster_id"]
            n_cells = row["n_cells"]

            if cluster_id.startswith("P"):
                # Pool row
                rows.append({
                    "cluster_id": cluster_id,
                    "origin_cluster": cluster_id,
                    "iteration_created": 0,
                    "assigned_label": row["assigned_label"],
                    "assigned_score": pool_score_lookup.get(cluster_id, 0.0),
                    "n_cells": n_cells,
                    "is_parent": False,
                    "children": "",
                    "best_child": "",
                    "best_child_score": 0.0,
                    "second_child": "",
                    "second_child_score": 0.0,
                    "max_child_score": 0.0,
                    "all_child_scores": "",
                    "marker_heterogeneity": 0.0,
                    "passes_low_confidence": True,
                    "passes_subtype_signal": False,
                    "passes_heterogeneous": False,
                    "passes_mixed_population": False,
                    "passes_ambiguous_root": "~" in row["assigned_label"],
                    "passes_weak_leaf": False,
                    "recommendation": "SUBCLUSTER",
                    "recommendation_reason": "Pooled cluster - subclustering recommended",
                    "confidence_band": "pool",
                    "triggered_criteria": '["POOL"]',
                    "is_pool": True,
                })
            else:
                # Unchanged row - copy from original diagnostic if available
                orig = orig_diag_lookup.get(cluster_id, {})
                input_ann = input_ann_lookup.get(cluster_id, {})

                # Get score from input annotations (more reliable) or fallback to diagnostic
                assigned_score = float(input_ann.get("assigned_score", orig.get("assigned_score", 0.0)) or 0.0)

                rows.append({
                    "cluster_id": cluster_id,
                    "origin_cluster": str(orig.get("origin_cluster", cluster_id)),
                    "iteration_created": int(orig.get("iteration_created", 0) or 0),
                    "assigned_label": row["assigned_label"],
                    "assigned_score": assigned_score,
                    "n_cells": n_cells,
                    "is_parent": bool(orig.get("is_parent", False)),
                    "children": str(orig.get("children", "") or ""),
                    "best_child": str(orig.get("best_child", "") or ""),
                    "best_child_score": float(orig.get("best_child_score", 0.0) or 0.0),
                    "second_child": str(orig.get("second_child", "") or ""),
                    "second_child_score": float(orig.get("second_child_score", 0.0) or 0.0),
                    "max_child_score": float(orig.get("max_child_score", 0.0) or 0.0),
                    "all_child_scores": str(orig.get("all_child_scores", "") or ""),
                    "marker_heterogeneity": float(orig.get("marker_heterogeneity", 0.0) or 0.0),
                    "passes_low_confidence": bool(orig.get("passes_low_confidence", False)),
                    "passes_subtype_signal": bool(orig.get("passes_subtype_signal", False)),
                    "passes_heterogeneous": bool(orig.get("passes_heterogeneous", False)),
                    "passes_mixed_population": bool(orig.get("passes_mixed_population", False)),
                    "passes_ambiguous_root": bool(orig.get("passes_ambiguous_root", False)),
                    "passes_weak_leaf": bool(orig.get("passes_weak_leaf", False)),
                    "recommendation": str(orig.get("recommendation", "SKIP")),
                    "recommendation_reason": str(orig.get("recommendation_reason", "")),
                    "confidence_band": str(orig.get("confidence_band", "medium")),
                    "triggered_criteria": str(orig.get("triggered_criteria", '["NONE"]')),
                    "is_pool": False,
                })

        result = pd.DataFrame(rows)
        result = self._sort_clusters(result)
        return result

    def generate_full_cluster_annotations(
        self,
        adata: "sc.AnnData",
        pool_summary: pd.DataFrame,
        input_cluster_annotations: Optional[pd.DataFrame] = None,
    ) -> pd.DataFrame:
        """Generate cluster_annotations.csv with ALL Stage I columns (24 columns).

        This is an enhanced version that includes all columns expected by Stage I
        for proper continuation.

        Parameters
        ----------
        adata : sc.AnnData
            AnnData after pooling
        pool_summary : pd.DataFrame
            Pool summary from execute()
        input_cluster_annotations : pd.DataFrame, optional
            Original cluster annotations from Stage I input

        Returns
        -------
        pd.DataFrame
            Full cluster annotations compatible with Stage I (24 columns)
        """
        # Use vectorized groupby for speed
        obs_subset = adata.obs[[self.cluster_col, self.label_col]].copy()
        obs_subset[self.cluster_col] = obs_subset[self.cluster_col].astype(str)
        obs_subset[self.label_col] = obs_subset[self.label_col].astype(str)

        grouped = obs_subset.groupby(self.cluster_col, observed=True)
        cluster_stats = grouped.agg(
            assigned_label=(self.label_col, lambda x: x.mode().iloc[0] if len(x) > 0 else "Unknown"),
            n_cells=(self.label_col, "size"),
        ).reset_index().rename(columns={self.cluster_col: "cluster_id"})

        total_cells = len(adata)

        # Build pool origin lookup
        pool_origins = {}
        if not pool_summary.empty and "new_cluster_id" in pool_summary.columns:
            for new_cid in pool_summary["new_cluster_id"].unique():
                source_rows = pool_summary[pool_summary["new_cluster_id"] == new_cid]
                pool_origins[new_cid] = ";".join(source_rows["cluster_id"].astype(str).tolist())

        # Build input annotations lookup
        input_lookup = {}
        if input_cluster_annotations is not None and not input_cluster_annotations.empty:
            input_cluster_annotations = input_cluster_annotations.copy()
            input_cluster_annotations["cluster_id"] = input_cluster_annotations["cluster_id"].astype(str)
            for _, row in input_cluster_annotations.iterrows():
                input_lookup[row["cluster_id"]] = row.to_dict()

        # Build pool score lookup
        pool_score_lookup = self._compute_pool_scores(
            pool_summary,
            {cid: {"assigned_score": d.get("assigned_score", 0.0), "n_cells": d.get("n_cells", 0)}
             for cid, d in input_lookup.items()}
        )

        # Build rows
        rows = []
        for _, row in cluster_stats.iterrows():
            cluster_id = row["cluster_id"]
            n_cells = row["n_cells"]
            assigned_label = row["assigned_label"]
            proportion = n_cells / total_cells if total_cells > 0 else 0.0

            if cluster_id.startswith("P"):
                # Pool row
                root_label = assigned_label
                if assigned_label.startswith("Pool_"):
                    root_label = assigned_label[5:]
                    if "~" in root_label:
                        root_label = root_label.split("~")[0]

                rows.append({
                    "cluster_id": cluster_id,
                    "origin_cluster": cluster_id,
                    "iteration_created": 0,
                    "proportion": proportion,
                    "mean_enrichment": 0.0,
                    "mean_positive_fraction": 0.0,
                    "n_cells": n_cells,
                    "assigned_label": assigned_label,
                    "assigned_path": assigned_label,
                    "assigned_level": 0,
                    "assigned_score": pool_score_lookup.get(cluster_id, 0.0),
                    "root_label": root_label,
                    "confidence": 0.0,
                    "min_margin_along_path": 0.0,
                    "margin_is_infinite": False,
                    "stop_reason": "pool",
                    "stopped_before_leaf": True,
                    "composition": "",
                    "decision_trace": "",
                    "coverage": 0.0,
                    "resolved_markers": "",
                    "is_ambiguous_root": "~" in assigned_label,
                    "ambiguous_root_candidates": "",
                    "confidence_level": "pool",
                })
            else:
                # Unchanged row - copy from input annotations
                src = input_lookup.get(cluster_id, {})
                rows.append({
                    "cluster_id": cluster_id,
                    "origin_cluster": str(src.get("origin_cluster", cluster_id)),
                    "iteration_created": int(src.get("iteration_created", 0) or 0),
                    "proportion": proportion,
                    "mean_enrichment": float(src.get("mean_enrichment", 0.0) or 0.0),
                    "mean_positive_fraction": float(src.get("mean_positive_fraction", 0.0) or 0.0),
                    "n_cells": n_cells,
                    "assigned_label": assigned_label,
                    "assigned_path": str(src.get("assigned_path", assigned_label)),
                    "assigned_level": int(src.get("assigned_level", 0) or 0),
                    "assigned_score": float(src.get("assigned_score", 0.0) or 0.0),
                    "root_label": str(src.get("root_label", "")),
                    "confidence": float(src.get("confidence", 0.0) or 0.0),
                    "min_margin_along_path": float(src.get("min_margin_along_path", 0.0) or 0.0),
                    "margin_is_infinite": bool(src.get("margin_is_infinite", False)),
                    "stop_reason": str(src.get("stop_reason", "")),
                    "stopped_before_leaf": bool(src.get("stopped_before_leaf", False)),
                    "composition": str(src.get("composition", "") or ""),
                    "decision_trace": str(src.get("decision_trace", "") or ""),
                    "coverage": float(src.get("coverage", 0.0) or 0.0),
                    "resolved_markers": str(src.get("resolved_markers", "") or ""),
                    "is_ambiguous_root": bool(src.get("is_ambiguous_root", False)),
                    "ambiguous_root_candidates": str(src.get("ambiguous_root_candidates", "") or ""),
                    "confidence_level": str(src.get("confidence_level", "medium")),
                })

        result = pd.DataFrame(rows)
        result = self._sort_clusters(result)
        return result

    def generate_enhanced_annotations(
        self,
        adata: "sc.AnnData",
        pool_summary: pd.DataFrame,
        input_enhanced_annotations: Optional[pd.DataFrame] = None,
    ) -> pd.DataFrame:
        """Generate cluster_annotations_enhanced.csv for Stage I compatibility (41 columns).

        For unchanged clusters: copy rows directly from input enhanced annotations.
        For pools: create new rows with pool-specific information.

        Parameters
        ----------
        adata : sc.AnnData
            AnnData after pooling
        pool_summary : pd.DataFrame
            Pool summary from execute()
        input_enhanced_annotations : pd.DataFrame, optional
            cluster_annotations_enhanced.csv from input Stage I iteration

        Returns
        -------
        pd.DataFrame
            Enhanced annotations compatible with Stage I (41 columns)
        """
        # Use vectorized groupby for speed
        obs_subset = adata.obs[[self.cluster_col, self.label_col]].copy()
        obs_subset[self.cluster_col] = obs_subset[self.cluster_col].astype(str)
        obs_subset[self.label_col] = obs_subset[self.label_col].astype(str)

        grouped = obs_subset.groupby(self.cluster_col, observed=True)
        cluster_stats = grouped.agg(
            n_cells=(self.label_col, "size"),
        ).reset_index().rename(columns={self.cluster_col: "cluster_id"})

        # Get unchanged clusters
        unchanged_clusters = set()
        if not pool_summary.empty:
            unchanged_rows = pool_summary[pool_summary["is_pooled"] == False]
            unchanged_clusters = set(unchanged_rows["cluster_id"].astype(str).tolist())

        # If we have input enhanced annotations, use them for unchanged clusters
        # NOTE: We copy rows as-is from input (matching reference behavior in
        # ft/src/pooling/engine.py:574-578). The proportion values come from
        # Stage I's export_enhanced_annotations() which computes percentages.
        rows = []
        if input_enhanced_annotations is not None and not input_enhanced_annotations.empty:
            input_enhanced_annotations = input_enhanced_annotations.copy()
            input_enhanced_annotations["cluster_id"] = input_enhanced_annotations["cluster_id"].astype(str)

            # Copy unchanged cluster rows
            for _, row in input_enhanced_annotations.iterrows():
                cid = row["cluster_id"]
                if cid in unchanged_clusters:
                    rows.append(row.to_dict())

        # Add pool rows with default values
        pool_clusters = set(cluster_stats["cluster_id"]) - unchanged_clusters
        for cid in sorted(pool_clusters):
            if not cid.startswith("P"):
                continue
            n_cells = int(cluster_stats[cluster_stats["cluster_id"] == cid]["n_cells"].iloc[0])
            # Get label from adata
            cell_type = adata.obs[adata.obs[self.cluster_col].astype(str) == cid][self.label_col].iloc[0] if n_cells > 0 else "Unknown"

            # Extract root label from pool label
            root_label = cell_type
            if cell_type.startswith("Pool_"):
                root_label = cell_type[5:]
                if "~" in root_label:
                    root_label = root_label.split("~")[0]

            rows.append({
                "cluster_id": cid,
                "origin_cluster": cid,
                "iteration_created": 0,
                "n_cells": n_cells,
                # NOTE: Pool proportion uses raw fraction (n_cells / len(adata)),
                # matching reference behavior in ft/src/pooling/engine.py:601.
                # This differs from non-Pool clusters which use percentages
                # (from Stage I export_enhanced_annotations()).
                "proportion": n_cells / len(adata),
                "reason": "Pool created by Stage J",
                "cell_type_lvl0": root_label,
                "score_lvl0": 0.0,
                "assigned_path_lvl0": root_label,
                "assigned_level_lvl0": 0,
                "root_label_lvl0": root_label,
                "confidence_lvl0": 0.0,
                "min_margin_along_path_lvl0": 0.0,
                "margin_is_infinite_lvl0": False,
                "stop_reason_lvl0": "pool",
                "stopped_before_leaf_lvl0": True,
                "composition_lvl0": "",
                "decision_trace_lvl0": "",
                "coverage_lvl0": 0.0,
                "resolved_markers_lvl0": "",
                "is_ambiguous_root_lvl0": "~" in cell_type,
                "ambiguous_root_candidates_lvl0": "",
                "cell_type_lvl1": cell_type,
                "score_lvl1": 0.0,
                "assigned_path_lvl1": cell_type,
                "assigned_level_lvl1": 0,
                "root_label_lvl1": root_label,
                "confidence_lvl1": 0.0,
                "min_margin_along_path_lvl1": 0.0,
                "margin_is_infinite_lvl1": False,
                "stop_reason_lvl1": "pool",
                "stopped_before_leaf_lvl1": True,
                "composition_lvl1": "",
                "decision_trace_lvl1": "",
                "coverage_lvl1": 0.0,
                "resolved_markers_lvl1": "",
                "is_ambiguous_root_lvl1": "~" in cell_type,
                "ambiguous_root_candidates_lvl1": "",
                "mean_enrichment": 0.0,
                "mean_positive_fraction": 0.0,
                "confidence_level": "pool",
                "runner_up_label": "",
                "runner_up_score": 0.0,
                "gap": 0.0,
            })

        if not rows:
            return pd.DataFrame()

        result = pd.DataFrame(rows)
        result = self._sort_clusters(result)
        return result
