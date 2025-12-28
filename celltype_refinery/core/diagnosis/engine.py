"""Diagnostic engine for cluster refinement decisions.

This module provides the DiagnosticEngine class for generating diagnostic reports
that recommend SUBCLUSTER, RELABEL, or SKIP actions for each cluster.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Set

import pandas as pd

from .criteria import CriteriaConfig, CriteriaEvaluator, build_hierarchy_from_scores


class DiagnosticEngine:
    """Generates diagnostic reports for cluster refinement.

    This engine evaluates clusters against refinement criteria and generates
    recommendations for each cluster: SUBCLUSTER, RELABEL, or SKIP.

    Example:
        >>> engine = DiagnosticEngine(config=CriteriaConfig(score_threshold=1.0))
        >>> report = engine.diagnose(cluster_annotations, marker_scores)
        >>> print(report["recommendation"].value_counts())
        SUBCLUSTER    25
        SKIP          10
        RELABEL        5
    """

    def __init__(
        self,
        config: Optional[CriteriaConfig] = None,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize diagnostic engine.

        Parameters
        ----------
        config : CriteriaConfig, optional
            Configuration for diagnostic criteria. Uses defaults if not provided.
        logger : logging.Logger, optional
            Logger for diagnostic messages.
        """
        self.config = config or CriteriaConfig()
        self.logger = logger or logging.getLogger(__name__)
        self._evaluator = CriteriaEvaluator(config=self.config, logger=self.logger)
        self._hierarchy: Dict[str, List[str]] = {}

    def diagnose(
        self,
        cluster_annotations: pd.DataFrame,
        marker_scores: pd.DataFrame,
        eligible_clusters: Optional[Set[str]] = None,
        heterogeneity_metrics: Optional[pd.DataFrame] = None,
    ) -> pd.DataFrame:
        """Generate diagnostic report with recommendations for each cluster.

        Parameters
        ----------
        cluster_annotations : pd.DataFrame
            DataFrame with columns: cluster_id, assigned_label, assigned_score, n_cells.
            Output from clustering/annotation stage or refinement.
        marker_scores : pd.DataFrame
            Full marker scores for all clusters.
            Columns: cluster_id, label, score, path.
        eligible_clusters : Set[str], optional
            If provided, only evaluate these clusters. Others get SKIP.
        heterogeneity_metrics : pd.DataFrame, optional
            Pre-computed marker heterogeneity metrics for weak leaf detection.
            Columns: cluster_id, marker_heterogeneity.

        Returns
        -------
        pd.DataFrame
            Diagnostic report with columns:
            - cluster_id, assigned_label, assigned_score, n_cells
            - is_parent, children, best_child, best_child_score, etc.
            - passes_* columns for each criterion
            - recommendation (SUBCLUSTER | RELABEL | SKIP)
            - recommendation_reason
            - confidence_band, triggered_criteria
        """
        # Build hierarchy from marker scores
        self._hierarchy = build_hierarchy_from_scores(marker_scores)

        # Set heterogeneity metrics if provided
        if heterogeneity_metrics is not None:
            self._evaluator.set_heterogeneity_metrics(heterogeneity_metrics)

        results = []
        for _, row in cluster_annotations.iterrows():
            cluster_id = str(row["cluster_id"])

            # Skip ineligible clusters if filter provided
            if eligible_clusters is not None and cluster_id not in eligible_clusters:
                results.append(self._skip_result(row, "Not in eligible clusters"))
                continue

            # Evaluate criteria
            result = self._evaluator.evaluate(row, marker_scores, self._hierarchy)
            results.append(result)

        report = pd.DataFrame(results)

        # Add lineage columns for iterative refinement tracking
        report = self._add_lineage_columns(report, cluster_annotations)

        # Log summary
        n_sub = (report["recommendation"] == "SUBCLUSTER").sum()
        n_rel = (report["recommendation"] == "RELABEL").sum()
        n_skip = (report["recommendation"] == "SKIP").sum()
        self.logger.info("Diagnostic Summary:")
        self.logger.info("  SUBCLUSTER: %d clusters", n_sub)
        self.logger.info("  RELABEL: %d clusters", n_rel)
        self.logger.info("  SKIP: %d clusters", n_skip)

        return report

    def _add_lineage_columns(
        self,
        report: pd.DataFrame,
        cluster_annotations: pd.DataFrame,
    ) -> pd.DataFrame:
        """Add lineage tracking columns to diagnostic report.

        Parameters
        ----------
        report : pd.DataFrame
            Diagnostic report DataFrame
        cluster_annotations : pd.DataFrame
            Cluster annotations (may contain origin_cluster column)

        Returns
        -------
        pd.DataFrame
            Report with added columns:
            - origin_cluster: Original cluster ID (e.g., "3" for "3:0:1")
            - iteration_created: Iteration when cluster was created (count of ":")
        """
        # Compute origin_cluster from cluster_id format
        # "3" → "3", "3:0" → "3", "3:0:1" → "3"
        def get_origin(cid: str) -> str:
            return str(cid).split(":")[0] if ":" in str(cid) else str(cid)

        # Compute iteration_created from cluster_id format
        # "3" → 0, "3:0" → 1, "3:0:1" → 2
        def get_iteration_created(cid: str) -> int:
            return str(cid).count(":")

        # Add columns to report
        report["origin_cluster"] = report["cluster_id"].apply(get_origin)
        report["iteration_created"] = report["cluster_id"].apply(get_iteration_created)

        # Try to get origin_cluster from cluster_annotations if available (more accurate)
        if "origin_cluster" in cluster_annotations.columns:
            origin_lookup = dict(zip(
                cluster_annotations["cluster_id"].astype(str),
                cluster_annotations["origin_cluster"].astype(str)
            ))
            report["origin_cluster"] = report["cluster_id"].apply(
                lambda cid: origin_lookup.get(str(cid), get_origin(str(cid)))
            )

        # Reorder columns to put lineage columns after cluster_id
        cols = report.columns.tolist()
        # Move origin_cluster and iteration_created after cluster_id
        for col in ["iteration_created", "origin_cluster"]:
            if col in cols:
                cols.remove(col)
                idx = cols.index("cluster_id") + 1
                cols.insert(idx, col)
        report = report[cols]

        return report

    def _skip_result(self, row: pd.Series, reason: str) -> Dict[str, Any]:
        """Generate a SKIP result for a cluster."""
        return {
            "cluster_id": str(row["cluster_id"]),
            "assigned_label": row["assigned_label"],
            "assigned_score": float(row.get("assigned_score", 0.0) or 0.0),
            "n_cells": int(row["n_cells"]),
            "is_parent": False,
            "children": "",
            "best_child": "",
            "best_child_score": None,
            "second_child": "",
            "second_child_score": None,
            "max_child_score": None,
            "all_child_scores": "",
            "marker_heterogeneity": None,
            "passes_low_confidence": False,
            "passes_subtype_signal": False,
            "passes_heterogeneous": False,
            "passes_mixed_population": False,
            "passes_ambiguous_root": False,
            "passes_weak_leaf": False,
            "recommendation": "SKIP",
            "recommendation_reason": reason,
            "ideal_recommendation": "SKIP",
            "is_blocked_by_size": False,
            "confidence_band": "",
            "triggered_criteria": '["EXCLUDED"]',
        }

    def get_subcluster_ids(self, report: pd.DataFrame) -> List[str]:
        """Get cluster IDs recommended for subclustering.

        Parameters
        ----------
        report : pd.DataFrame
            Diagnostic report from diagnose().

        Returns
        -------
        List[str]
            Cluster IDs with SUBCLUSTER recommendation.
        """
        return report[report["recommendation"] == "SUBCLUSTER"]["cluster_id"].tolist()

    def get_relabel_ids(self, report: pd.DataFrame) -> List[str]:
        """Get cluster IDs recommended for relabeling.

        Parameters
        ----------
        report : pd.DataFrame
            Diagnostic report from diagnose().

        Returns
        -------
        List[str]
            Cluster IDs with RELABEL recommendation.
        """
        return report[report["recommendation"] == "RELABEL"]["cluster_id"].tolist()

    def derive_groups(
        self,
        report: pd.DataFrame,
        cluster_annotations: pd.DataFrame,
        group_order: Optional[List[str]] = None,
    ) -> Dict[str, Any]:
        """Derive execution groups from diagnostic report.

        Groups clusters by their root_label for focused execution.

        Parameters
        ----------
        report : pd.DataFrame
            Diagnostic report from diagnose().
        cluster_annotations : pd.DataFrame
            Cluster annotations with root_label column.
        group_order : List[str], optional
            Custom ordering for groups. Defaults to standard biological order.

        Returns
        -------
        Dict[str, Any]
            Group configuration with structure:
            {
                "groups": [
                    {
                        "name": "Epithelium",
                        "focus_label": "Epithelium",
                        "cluster_ids": ["8", "9", "10"],
                        "subcluster_ids": ["8", "9"],
                        "n_clusters": 3,
                        "n_subcluster": 2,
                    },
                    ...
                ]
            }
        """
        # Merge report with annotations to get root_label
        merged = report.merge(
            cluster_annotations[["cluster_id", "root_label"]].astype(str),
            on="cluster_id",
            how="left",
        )

        # Handle missing root_label
        if "root_label" not in merged.columns:
            merged["root_label"] = merged["assigned_label"]

        # Fill NaN root_label with "Unassigned"
        merged["root_label"] = merged["root_label"].fillna("Unassigned")

        # Default group ordering (biological order)
        if group_order is None:
            group_order = [
                "Epithelium",
                "Endothelium",
                "Mesenchymal Cells",
                "Immune Cells",
                "Hybrids",
                "Unassigned",
            ]

        # Group by root_label
        groups = []
        root_labels = merged["root_label"].unique()

        # Process in standard order first
        for label in group_order:
            if label not in root_labels:
                continue

            group_data = merged[merged["root_label"] == label]
            cluster_ids = group_data["cluster_id"].tolist()
            subcluster_ids = group_data[group_data["recommendation"] == "SUBCLUSTER"][
                "cluster_id"
            ].tolist()

            groups.append(
                {
                    "name": label,
                    "focus_label": label,
                    "cluster_ids": cluster_ids,
                    "subcluster_ids": subcluster_ids,
                    "n_clusters": len(cluster_ids),
                    "n_subcluster": len(subcluster_ids),
                }
            )

        # Add any remaining labels not in standard order
        for label in sorted(root_labels):
            if label in group_order:
                continue

            group_data = merged[merged["root_label"] == label]
            cluster_ids = group_data["cluster_id"].tolist()
            subcluster_ids = group_data[group_data["recommendation"] == "SUBCLUSTER"][
                "cluster_id"
            ].tolist()

            groups.append(
                {
                    "name": label,
                    "focus_label": label,
                    "cluster_ids": cluster_ids,
                    "subcluster_ids": subcluster_ids,
                    "n_clusters": len(cluster_ids),
                    "n_subcluster": len(subcluster_ids),
                }
            )

        return {"groups": groups}

    def save_report(self, report: pd.DataFrame, path: Path) -> None:
        """Save diagnostic report to CSV.

        Parameters
        ----------
        report : pd.DataFrame
            Diagnostic report from diagnose().
        path : Path
            Output path for CSV file.
        """
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        report.to_csv(path, index=False)
        self.logger.info("Saved diagnostic report: %s", path)
