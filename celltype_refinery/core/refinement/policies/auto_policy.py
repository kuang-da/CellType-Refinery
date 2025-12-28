"""
AutoPolicy: Automatic criteria-based candidate selection for cell-type refinement.

This policy implements the automatic selection logic:
1. Low confidence: score < threshold → SUBCLUSTER
2. Homogeneous parent: gap >= threshold → RELABEL (instant, no clustering)
3. Heterogeneous parent: gap < threshold → SUBCLUSTER
4. Mixed population: high parent score but ALL children have low/negative scores → SUBCLUSTER
"""

from __future__ import annotations

import json
import logging
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import pandas as pd
import yaml

from ..plan import (
    RefinePlan,
    SubclusterOp,
    RelabelOp,
    RescoreOp,
)


@dataclass
class AutoPolicyConfig:
    """Configuration for AutoPolicy candidate selection."""
    min_cells: int = 500
    score_threshold: float = 1.0
    subtype_signal_threshold: float = 1.0
    heterogeneity_gap: float = 0.5
    subcluster_resolution: float = 0.4
    enable_parent_subclustering: bool = True
    n_pcs: int = 30
    neighbors_k: int = 15
    min_subcluster_cells: int = 20
    add_rescore: bool = True
    rescore_mode: str = "smart"
    # Mixed population detection: high-score parent with all low-score children
    mixed_population_parent_threshold: float = 1.0  # Parent must have score >= this
    mixed_population_max_child_threshold: float = 1.0  # All children must have score < this
    # Weak leaf detection: leaf nodes with borderline scores and high heterogeneity
    weak_leaf_score_threshold: float = 1.5  # Leaf with score < this is "weak"
    weak_leaf_heterogeneity_threshold: float = 0.12  # Heterogeneity > this triggers subcluster
    enable_weak_leaf_detection: bool = True  # Enable weak leaf criterion


class AutoPolicy:
    """Generates RefinePlan using automatic criteria.

    This policy automatically selects clusters for refinement based on:
    1. Low confidence: assigned_score < score_threshold
    2. Parent with subtype signal: parent type with detectable child scores

    For parent types, the action depends on heterogeneity:
    - Homogeneous (one dominant child): RELABEL directly (fast)
    - Heterogeneous (competing children): SUBCLUSTER (expensive but necessary)

    Example:
        >>> policy = AutoPolicy(score_threshold=1.0, heterogeneity_gap=0.5)
        >>> plan = policy.generate_plan(cluster_annotations, marker_scores)
        >>> print(plan.summary())
        {'subcluster': 3, 'relabel': 2, 'rescore': 1}
    """

    def __init__(
        self,
        min_cells: int = 500,
        score_threshold: float = 1.0,
        subtype_signal_threshold: float = 1.0,
        heterogeneity_gap: float = 0.5,
        subcluster_resolution: float = 0.4,
        enable_parent_subclustering: bool = True,
        n_pcs: int = 30,
        neighbors_k: int = 15,
        min_subcluster_cells: int = 20,
        add_rescore: bool = True,
        rescore_mode: str = "smart",
        mixed_population_parent_threshold: float = 1.0,
        mixed_population_max_child_threshold: float = 1.0,
        weak_leaf_score_threshold: float = 1.5,
        weak_leaf_heterogeneity_threshold: float = 0.12,
        enable_weak_leaf_detection: bool = True,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize AutoPolicy with selection parameters.

        Parameters
        ----------
        min_cells : int
            Minimum cluster size for subclustering (default: 500)
        score_threshold : float
            Maximum score for low-confidence criterion (default: 1.0)
        subtype_signal_threshold : float
            Minimum child score to trigger parent subclustering (default: 1.0)
        heterogeneity_gap : float
            Score gap threshold for detecting heterogeneous clusters (default: 0.5).
            If gap >= this value: homogeneous → RELABEL
            If gap < this value: heterogeneous → SUBCLUSTER
        subcluster_resolution : float
            Leiden resolution for subclustering (default: 0.4)
        enable_parent_subclustering : bool
            Enable parent-aware subclustering (default: True)
        n_pcs : int
            Number of principal components for subclustering
        neighbors_k : int
            Number of neighbors for k-NN graph
        min_subcluster_cells : int
            Minimum cells per subcluster
        add_rescore : bool
            Whether to add rescore operation at end of plan
        rescore_mode : str
            Rescoring mode: "smart" or "full"
        mixed_population_parent_threshold : float
            Minimum parent score to trigger mixed population detection (default: 1.0)
        mixed_population_max_child_threshold : float
            Maximum child score for mixed population (default: 1.0).
            If parent >= parent_threshold AND all children < this, it's a mixed population.
        weak_leaf_score_threshold : float
            Maximum score for weak leaf detection (default: 1.5).
            Leaf clusters with score < this are considered "weak" candidates.
        weak_leaf_heterogeneity_threshold : float
            Minimum marker heterogeneity to trigger weak leaf subclustering (default: 0.12).
            Only weak leaves with heterogeneity > this will be subclustered.
        enable_weak_leaf_detection : bool
            Enable weak leaf criterion (default: True)
        logger : logging.Logger, optional
            Logger for tracking operations
        """
        self.config = AutoPolicyConfig(
            min_cells=min_cells,
            score_threshold=score_threshold,
            subtype_signal_threshold=subtype_signal_threshold,
            heterogeneity_gap=heterogeneity_gap,
            subcluster_resolution=subcluster_resolution,
            enable_parent_subclustering=enable_parent_subclustering,
            n_pcs=n_pcs,
            neighbors_k=neighbors_k,
            min_subcluster_cells=min_subcluster_cells,
            add_rescore=add_rescore,
            rescore_mode=rescore_mode,
            mixed_population_parent_threshold=mixed_population_parent_threshold,
            mixed_population_max_child_threshold=mixed_population_max_child_threshold,
            weak_leaf_score_threshold=weak_leaf_score_threshold,
            weak_leaf_heterogeneity_threshold=weak_leaf_heterogeneity_threshold,
            enable_weak_leaf_detection=enable_weak_leaf_detection,
        )
        self.logger = logger or logging.getLogger(__name__)
        self._heterogeneity_metrics: Optional[pd.DataFrame] = None

    def set_heterogeneity_metrics(self, metrics: pd.DataFrame) -> None:
        """Set pre-computed marker heterogeneity metrics for weak leaf detection.

        Parameters
        ----------
        metrics : pd.DataFrame
            DataFrame with columns: cluster_id, marker_heterogeneity, mixed_marker_rate.
            Output from compute_cluster_marker_heterogeneity().

        Notes
        -----
        This must be called before generate_plan() or generate_diagnostic_report()
        for weak_leaf detection to work. If not set, weak_leaf criterion will be
        skipped with a warning.
        """
        if "cluster_id" not in metrics.columns:
            raise ValueError("metrics must have 'cluster_id' column")
        if "marker_heterogeneity" not in metrics.columns:
            raise ValueError("metrics must have 'marker_heterogeneity' column")

        # Ensure cluster_id is string for consistent matching
        metrics = metrics.copy()
        metrics["cluster_id"] = metrics["cluster_id"].astype(str)
        self._heterogeneity_metrics = metrics
        self.logger.info(
            "Heterogeneity metrics set: %d clusters, mean=%.3f, max=%.3f",
            len(metrics),
            metrics["marker_heterogeneity"].mean(),
            metrics["marker_heterogeneity"].max(),
        )

    def _get_heterogeneity(self, cluster_id: str) -> Optional[float]:
        """Get marker heterogeneity for a cluster.

        Returns None if heterogeneity metrics are not set or cluster not found.
        """
        if self._heterogeneity_metrics is None:
            return None
        row = self._heterogeneity_metrics[
            self._heterogeneity_metrics["cluster_id"] == str(cluster_id)
        ]
        if row.empty:
            return None
        return float(row.iloc[0]["marker_heterogeneity"])

    def generate_plan(
        self,
        cluster_annotations: pd.DataFrame,
        marker_scores: pd.DataFrame,
        hierarchy_config: Optional[Dict] = None,
        eligible_clusters: Optional[Set[str]] = None,
    ) -> RefinePlan:
        """Generate RefinePlan using automatic selection criteria.

        Parameters
        ----------
        cluster_annotations : pd.DataFrame
            Cluster assignments with CANONICAL column names:
            cluster_id, assigned_label, assigned_score, n_cells.
            Use InputAdapter.adapt_cluster_annotations() to normalize first.
        marker_scores : pd.DataFrame
            Full marker scores with CANONICAL column names:
            cluster_id, label, path, score.
            Use InputAdapter.adapt_marker_scores() to normalize first.
        hierarchy_config : Dict, optional
            Cell type hierarchy configuration with category_lookup
        eligible_clusters : Set[str], optional
            Set of cluster IDs eligible for refinement (focus controls).
            If provided, only clusters in this set will be considered.
            If None, all clusters are eligible.

        Returns
        -------
        RefinePlan
            Plan with subcluster and relabel operations

        Raises
        ------
        ValueError
            If required columns are missing from DataFrames

        Notes
        -----
        This method expects DataFrames with canonical column names.
        If your input uses different column names, first adapt them:

            >>> from celltype_refinery.core.refinement.schema import InputAdapter
            >>> adapter = InputAdapter()
            >>> annotations = adapter.adapt_cluster_annotations(raw_annotations)
            >>> scores = adapter.adapt_marker_scores(raw_scores)
            >>> plan = policy.generate_plan(annotations, scores)
        """
        # Validate required columns in cluster_annotations
        required_annotation_cols = ["cluster_id", "assigned_label", "assigned_score", "n_cells"]
        missing_annotation = [c for c in required_annotation_cols if c not in cluster_annotations.columns]
        if missing_annotation:
            raise ValueError(
                f"cluster_annotations missing required columns: {missing_annotation}. "
                f"Available columns: {list(cluster_annotations.columns)}. "
                f"Ensure data was adapted via InputAdapter.adapt_cluster_annotations()"
            )

        # Validate required columns in marker_scores
        required_score_cols = ["cluster_id", "label", "path", "score"]
        missing_score = [c for c in required_score_cols if c not in marker_scores.columns]
        if missing_score:
            raise ValueError(
                f"marker_scores missing required columns: {missing_score}. "
                f"Available columns: {list(marker_scores.columns)}. "
                f"Ensure data was adapted via InputAdapter.adapt_marker_scores()"
            )

        plan = RefinePlan(
            metadata={
                "policy": "AutoPolicy",
                "created_at": datetime.now().isoformat(),
                "config": {
                    "min_cells": self.config.min_cells,
                    "score_threshold": self.config.score_threshold,
                    "subtype_signal_threshold": self.config.subtype_signal_threshold,
                    "heterogeneity_gap": self.config.heterogeneity_gap,
                    "subcluster_resolution": self.config.subcluster_resolution,
                    "enable_parent_subclustering": self.config.enable_parent_subclustering,
                },
            }
        )

        # Build parent-child hierarchy from marker scores
        hierarchy = self._build_hierarchy_from_scores(marker_scores) if self.config.enable_parent_subclustering else {}

        # Select candidates (respecting focus controls)
        candidates = self._select_candidates(
            cluster_annotations, marker_scores, hierarchy, eligible_clusters
        )

        # Track subclustered cluster IDs for rescore
        subclustered_clusters: List[str] = []

        # Generate operations from candidates
        for candidate in candidates:
            if candidate.action == "relabel":
                # Add relabel operation
                plan.add_relabel(
                    cluster_id=candidate.cluster_id,
                    new_label=candidate.best_child_label,
                    old_label=candidate.assigned_label,
                    confidence_score=candidate.max_child_score,
                    reason=f"Homogeneous parent: {candidate.assigned_label} → {candidate.best_child_label} "
                           f"(score={candidate.max_child_score:.3f})",
                    source="auto",
                )
            elif candidate.action == "subcluster":
                # Add subcluster operation
                plan.add_subcluster(
                    cluster_id=candidate.cluster_id,
                    resolution=self.config.subcluster_resolution,
                    n_pcs=self.config.n_pcs,
                    neighbors_k=self.config.neighbors_k,
                    min_cells=self.config.min_subcluster_cells,
                    reason=self._format_subcluster_reason(candidate),
                    source="auto",
                )
                subclustered_clusters.append(candidate.cluster_id)

        # Add rescore operation if requested and there are modifications
        if self.config.add_rescore and (subclustered_clusters or plan.get_operations_by_type("relabel")):
            # When subclustering creates new clusters, we need to run DE tests on them
            # to compute proper DE-boosted scores. Without recompute_de=True, subclusters
            # won't have DE results and scoring will be suboptimal.
            needs_de_recompute = bool(subclustered_clusters)
            plan.add_rescore(
                mode=self.config.rescore_mode,
                recompute_de=needs_de_recompute,
                target_clusters=subclustered_clusters if subclustered_clusters else None,
                reason="Rescore modified clusters" + (" (with DE recomputation)" if needs_de_recompute else ""),
            )

        # Sort operations in deterministic order
        plan.sort_operations()

        # Update metadata with results
        plan.metadata["n_candidates"] = len(candidates)
        plan.metadata["n_subclusters"] = len([c for c in candidates if c.action == "subcluster"])
        plan.metadata["n_relabels"] = len([c for c in candidates if c.action == "relabel"])

        return plan

    def _build_hierarchy_from_scores(self, marker_scores: pd.DataFrame) -> Dict[str, List[str]]:
        """Build parent-child relationships from marker_scores.csv path column.

        The path column contains hierarchical paths like:
        - "Epithelium"
        - "Epithelium / Ciliated Epithelium"
        - "Immune Cells / Myeloids / Granulocytes"

        Returns mapping from parent label to list of direct child labels.
        """
        parents_children: Dict[str, Set[str]] = defaultdict(set)

        if "path" not in marker_scores.columns:
            return {}

        for path in marker_scores["path"].unique():
            parts = [p.strip() for p in str(path).split(" / ")]
            for i in range(len(parts) - 1):
                parent = parts[i]
                child = parts[i + 1]
                parents_children[parent].add(child)

        return {k: sorted(list(v)) for k, v in parents_children.items()}

    def _select_candidates(
        self,
        cluster_annotations: pd.DataFrame,
        marker_scores: pd.DataFrame,
        hierarchy: Dict[str, List[str]],
        eligible_clusters: Optional[Set[str]] = None,
    ) -> List["_Candidate"]:
        """Select clusters for subclustering or relabeling.

        Parameters
        ----------
        cluster_annotations : pd.DataFrame
            Cluster assignments with canonical columns.
        marker_scores : pd.DataFrame
            Marker scores with canonical columns.
        hierarchy : Dict[str, List[str]]
            Parent-child relationships from marker hierarchy.
        eligible_clusters : Set[str], optional
            If provided, only clusters in this set are considered.

        Returns
        -------
        List[_Candidate]
            List of candidate objects with selection metadata.
        """
        candidates: List[_Candidate] = []

        self.logger.info("=" * 75)
        self.logger.info("AutoPolicy candidate selection")
        self.logger.info("=" * 75)

        # Log focus controls status
        if eligible_clusters is not None:
            self.logger.info("Focus controls: %d clusters eligible", len(eligible_clusters))
            # Log table of eligible clusters
            if len(eligible_clusters) > 0:
                eligible_rows = cluster_annotations[
                    cluster_annotations["cluster_id"].astype(str).isin(eligible_clusters)
                ].copy()
                if not eligible_rows.empty:
                    self.logger.info("  | cluster_id | assigned_label                       | assigned_score | n_cells   |")
                    self.logger.info("  |------------|--------------------------------------|----------------|-----------|")
                    for _, row in eligible_rows.iterrows():
                        self.logger.info(
                            "  | %-10s | %-36s | %14.3f | %9d |",
                            str(row["cluster_id"]),
                            str(row["assigned_label"])[:36],
                            row["assigned_score"],
                            row["n_cells"],
                        )
        else:
            self.logger.info("Focus controls: disabled (all clusters eligible)")

        self.logger.info("Criteria:")
        self.logger.info("  1. Low confidence: assigned_score < %.2f AND n_cells >= %d → SUBCLUSTER",
                        self.config.score_threshold, self.config.min_cells)
        if self.config.enable_parent_subclustering:
            self.logger.info("  2. Parent with subtype signal: best_child_score > %.2f",
                           self.config.subtype_signal_threshold)
            self.logger.info("     2a. Homogeneous → RELABEL if: (only 1 passing child) OR")
            self.logger.info("         (second_child < %.2f AND gap >= %.2f)",
                           self.config.subtype_signal_threshold, self.config.heterogeneity_gap)
            self.logger.info("     2b. Heterogeneous → SUBCLUSTER otherwise")
            self.logger.info("  3. Mixed population: assigned_score >= %.2f AND max_child_score < %.2f → SUBCLUSTER",
                           self.config.mixed_population_parent_threshold,
                           self.config.mixed_population_max_child_threshold)
        if self.config.enable_weak_leaf_detection:
            self.logger.info("  4. Weak leaf: is_leaf AND score < %.2f AND heterogeneity > %.2f → SUBCLUSTER",
                           self.config.weak_leaf_score_threshold,
                           self.config.weak_leaf_heterogeneity_threshold)
            if self._heterogeneity_metrics is None:
                self.logger.warning("     (Heterogeneity metrics not set - weak leaf detection disabled)")

        for _, row in cluster_annotations.iterrows():
            cluster_id = str(row["cluster_id"])
            assigned_label = row["assigned_label"]
            assigned_score = row["assigned_score"]
            n_cells = row["n_cells"]

            # Skip if not in eligible set (focus controls)
            if eligible_clusters is not None and cluster_id not in eligible_clusters:
                continue

            # Skip if too few cells
            if n_cells < self.config.min_cells:
                continue

            # Criterion 1: Low confidence → always subcluster
            if assigned_score < self.config.score_threshold:
                candidates.append(_Candidate(
                    cluster_id=cluster_id,
                    assigned_label=assigned_label,
                    assigned_score=assigned_score,
                    n_cells=n_cells,
                    reason="low_confidence",
                    action="subcluster",
                ))
                continue

            # Criterion 2 & 3: Parent-based criteria (if enabled)
            if self.config.enable_parent_subclustering and assigned_label in hierarchy:
                children = hierarchy[assigned_label]

                # Get scores for child labels for this cluster
                cluster_scores = marker_scores[marker_scores["cluster_id"] == row["cluster_id"]]
                child_scores = cluster_scores[cluster_scores["label"].isin(children)].copy()

                if not child_scores.empty and len(child_scores) >= 1:
                    # Compute max child score from ALL children (for criterion 3)
                    max_child_score_all = child_scores["score"].max()
                    best_child_all = child_scores.loc[child_scores["score"].idxmax()]
                    best_child_label_all = best_child_all["label"]

                    # Filter to children with positive scores (basic gate: score > 0)
                    # This ensures we only consider children with some marker evidence
                    passing_children = child_scores[child_scores["score"] > 0].copy()

                    # Criterion 2: Check for subtype signal (requires positive children)
                    criterion_2_triggered = False
                    if not passing_children.empty:
                        # Sort by score descending
                        passing_children = passing_children.sort_values("score", ascending=False)

                        best_child_score = passing_children.iloc[0]["score"]
                        best_child_label = passing_children.iloc[0]["label"]

                        # Get second best if available (only from passing children)
                        if len(passing_children) >= 2:
                            second_child_score = passing_children.iloc[1]["score"]
                            second_child_label = passing_children.iloc[1]["label"]
                        else:
                            second_child_score = -999.0
                            second_child_label = None

                        # Check if best child is above threshold
                        if best_child_score > self.config.subtype_signal_threshold:
                            criterion_2_triggered = True
                            gap = best_child_score - second_child_score

                            # Decide action based on gap
                            # Homogeneous if: (a) only one passing child OR (b) second < threshold AND gap >= heterogeneity_gap
                            if second_child_label is None or (second_child_score < self.config.subtype_signal_threshold and gap >= self.config.heterogeneity_gap):
                                # Homogeneous → RELABEL
                                candidates.append(_Candidate(
                                    cluster_id=cluster_id,
                                    assigned_label=assigned_label,
                                    assigned_score=assigned_score,
                                    n_cells=n_cells,
                                    reason="relabel_to_subtype",
                                    action="relabel",
                                    max_child_score=best_child_score,
                                    best_child_label=best_child_label,
                                    second_child_score=second_child_score if second_child_label else None,
                                    second_child_label=second_child_label,
                                ))
                            else:
                                # Heterogeneous → SUBCLUSTER
                                candidates.append(_Candidate(
                                    cluster_id=cluster_id,
                                    assigned_label=assigned_label,
                                    assigned_score=assigned_score,
                                    n_cells=n_cells,
                                    reason="parent_with_subtype_signal",
                                    action="subcluster",
                                    max_child_score=best_child_score,
                                    best_child_label=best_child_label,
                                    second_child_score=second_child_score,
                                    second_child_label=second_child_label,
                                ))

                    # Criterion 3: Mixed population
                    # High parent score but ALL children have low/negative scores
                    if not criterion_2_triggered:
                        if (assigned_score >= self.config.mixed_population_parent_threshold and
                            max_child_score_all < self.config.mixed_population_max_child_threshold):
                            candidates.append(_Candidate(
                                cluster_id=cluster_id,
                                assigned_label=assigned_label,
                                assigned_score=assigned_score,
                                n_cells=n_cells,
                                reason="mixed_population",
                                action="subcluster",
                                max_child_score=max_child_score_all,
                                best_child_label=best_child_label_all,
                            ))
                            continue

            # Criterion 4: Weak leaf - leaf node with borderline score and high heterogeneity
            # Only applies to clusters that are NOT parents (no children in hierarchy)
            if (self.config.enable_weak_leaf_detection and
                assigned_label not in hierarchy and
                assigned_score < self.config.weak_leaf_score_threshold):
                # Check marker heterogeneity
                het = self._get_heterogeneity(cluster_id)
                if het is not None and het > self.config.weak_leaf_heterogeneity_threshold:
                    candidates.append(_Candidate(
                        cluster_id=cluster_id,
                        assigned_label=assigned_label,
                        assigned_score=assigned_score,
                        n_cells=n_cells,
                        reason="weak_leaf",
                        action="subcluster",
                        marker_heterogeneity=het,
                    ))

        # Log results
        low_conf = [c for c in candidates if c.reason == "low_confidence"]
        relabels = [c for c in candidates if c.action == "relabel"]
        heterogeneous = [c for c in candidates if c.reason == "parent_with_subtype_signal"]
        mixed_pop = [c for c in candidates if c.reason == "mixed_population"]
        weak_leaves = [c for c in candidates if c.reason == "weak_leaf"]

        self.logger.info("-" * 75)
        self.logger.info("Selection results:")
        self.logger.info("  Low confidence (subcluster): %d clusters", len(low_conf))
        self.logger.info("  Homogeneous parent (relabel): %d clusters", len(relabels))
        self.logger.info("  Heterogeneous parent (subcluster): %d clusters", len(heterogeneous))
        self.logger.info("  Mixed population (subcluster): %d clusters", len(mixed_pop))
        self.logger.info("  Weak leaf (subcluster): %d clusters", len(weak_leaves))
        self.logger.info("  Total selected: %d clusters", len(candidates))

        return candidates

    def _format_subcluster_reason(self, candidate: "_Candidate") -> str:
        """Format a human-readable reason for subclustering."""
        if candidate.reason == "low_confidence":
            return f"Low confidence: {candidate.assigned_label} (score={candidate.assigned_score:.3f})"
        elif candidate.reason == "parent_with_subtype_signal":
            return (f"Heterogeneous parent: {candidate.assigned_label} → "
                   f"{candidate.best_child_label} vs {candidate.second_child_label} "
                   f"(scores={candidate.max_child_score:.3f} vs {candidate.second_child_score:.3f})")
        elif candidate.reason == "mixed_population":
            return (f"Mixed population: {candidate.assigned_label} "
                   f"(parent_score={candidate.assigned_score:.3f}, "
                   f"max_child={candidate.max_child_score:.3f} < {self.config.mixed_population_max_child_threshold})")
        elif candidate.reason == "weak_leaf":
            return (f"Weak leaf: {candidate.assigned_label} "
                   f"(score={candidate.assigned_score:.3f}, "
                   f"heterogeneity={candidate.marker_heterogeneity:.3f} > {self.config.weak_leaf_heterogeneity_threshold})")
        else:
            return candidate.reason

    def generate_diagnostic_report(
        self,
        cluster_annotations: pd.DataFrame,
        marker_scores: pd.DataFrame,
        eligible_clusters: Optional[Set[str]] = None,
    ) -> pd.DataFrame:
        """Generate diagnostic table showing all criteria for each cluster.

        Parameters
        ----------
        cluster_annotations : pd.DataFrame
            Cluster assignments with canonical columns.
        marker_scores : pd.DataFrame
            Marker scores with canonical columns.
        eligible_clusters : Set[str], optional
            If provided, only include these clusters in the report.

        Returns
        -------
        pd.DataFrame
            Diagnostic report with columns:
            - cluster_id, assigned_label, assigned_score, n_cells
            - is_parent, children (semicolon-separated list)
            - best_child, best_child_score, max_child_score
            - passes_low_confidence, passes_subtype_signal, passes_mixed_population
            - recommendation (SUBCLUSTER | RELABEL | SKIP)
            - recommendation_reason
        """
        # Build hierarchy from marker scores
        hierarchy = self._build_hierarchy_from_scores(marker_scores) if self.config.enable_parent_subclustering else {}

        results = []
        for _, row in cluster_annotations.iterrows():
            result = self._evaluate_cluster_criteria(row, marker_scores, hierarchy)
            results.append(result)

        df = pd.DataFrame(results)

        # Apply focus controls filter if provided
        if eligible_clusters is not None:
            df = df[df["cluster_id"].astype(str).isin(eligible_clusters)]

        return df

    def _evaluate_cluster_criteria(
        self,
        row: pd.Series,
        marker_scores: pd.DataFrame,
        hierarchy: Dict[str, List[str]],
    ) -> Dict[str, Any]:
        """Evaluate all criteria for a single cluster.

        Returns a dictionary with all diagnostic fields.
        """
        import numpy as np

        cluster_id = str(row["cluster_id"])
        assigned_label = row["assigned_label"]
        assigned_score = row["assigned_score"]
        n_cells = row["n_cells"]

        result = {
            "cluster_id": cluster_id,
            "assigned_label": assigned_label,
            "assigned_score": assigned_score,
            "n_cells": n_cells,
            "is_parent": assigned_label in hierarchy,
            "children": "",
            "best_child": "",
            "best_child_score": None,
            "second_child": "",
            "second_child_score": None,
            "max_child_score": None,
            "all_child_scores": "",
            # Heterogeneity metrics
            "marker_heterogeneity": None,
            # Criteria results
            "passes_low_confidence": False,
            "passes_subtype_signal": False,
            "passes_heterogeneous": False,
            "passes_mixed_population": False,
            "passes_ambiguous_root": False,
            "passes_weak_leaf": False,
            # Recommendation
            "recommendation": "SKIP",
            "recommendation_reason": "",
            # Audit columns
            "confidence_band": "",
            "triggered_criteria": "",
        }

        # Check cell count threshold
        if n_cells < self.config.min_cells:
            result["recommendation_reason"] = f"Too few cells ({n_cells} < {self.config.min_cells})"
            result["triggered_criteria"] = '["TOO_SMALL"]'
            return result

        # Criterion 1: Low confidence
        result["passes_low_confidence"] = assigned_score < self.config.score_threshold

        # Get child data if parent
        if result["is_parent"]:
            children = hierarchy[assigned_label]
            result["children"] = ";".join(children)

            cluster_scores = marker_scores[marker_scores["cluster_id"] == row["cluster_id"]]
            child_scores = cluster_scores[cluster_scores["label"].isin(children)].copy()

            if not child_scores.empty:
                # Sort by score
                child_scores = child_scores.sort_values("score", ascending=False)

                result["max_child_score"] = child_scores["score"].max()
                result["best_child"] = child_scores.iloc[0]["label"]
                result["best_child_score"] = child_scores.iloc[0]["score"]

                if len(child_scores) >= 2:
                    result["second_child"] = child_scores.iloc[1]["label"]
                    result["second_child_score"] = child_scores.iloc[1]["score"]

                result["all_child_scores"] = ";".join(
                    f"{r['label']}={r['score']:.2f}" for _, r in child_scores.iterrows()
                )

                # Criterion 2: Subtype signal (positive children above threshold)
                passing_children = child_scores[child_scores["score"] > 0]
                if not passing_children.empty:
                    best_passing_score = passing_children.iloc[0]["score"]
                    if best_passing_score > self.config.subtype_signal_threshold:
                        result["passes_subtype_signal"] = True

                        # Check heterogeneity
                        # Heterogeneous if: (2+ passing children) AND (second >= threshold OR gap < heterogeneity_gap)
                        if len(passing_children) >= 2:
                            second_passing_score = passing_children.iloc[1]["score"]
                            gap = best_passing_score - second_passing_score
                            if second_passing_score >= self.config.subtype_signal_threshold or gap < self.config.heterogeneity_gap:
                                result["passes_heterogeneous"] = True

                # Criterion 3: Mixed population
                if (assigned_score >= self.config.mixed_population_parent_threshold and
                    result["max_child_score"] is not None and
                    result["max_child_score"] < self.config.mixed_population_max_child_threshold):
                    result["passes_mixed_population"] = True

        # Criterion 4: Ambiguous root (from input or detected via ~ in label)
        is_ambiguous = row.get("is_ambiguous_root", False)
        if not is_ambiguous and "~" in str(row.get("assigned_label", "")):
            is_ambiguous = True
        if is_ambiguous:
            result["passes_ambiguous_root"] = True

        # Criterion 5: Weak leaf - leaf with borderline score and high heterogeneity
        if self.config.enable_weak_leaf_detection and not result["is_parent"]:
            het = self._get_heterogeneity(cluster_id)
            if het is not None:
                result["marker_heterogeneity"] = het
                if (assigned_score < self.config.weak_leaf_score_threshold and
                    het > self.config.weak_leaf_heterogeneity_threshold):
                    result["passes_weak_leaf"] = True

        # Generate recommendation
        result["recommendation"], result["recommendation_reason"] = self._get_recommendation(result)

        # Populate audit columns
        result["confidence_band"] = self._classify_confidence_band(assigned_score)
        result["triggered_criteria"] = self._get_reason_codes(result)

        return result

    def _get_recommendation(self, result: Dict[str, Any]) -> Tuple[str, str]:
        """Get recommendation based on evaluated criteria."""
        # Priority order: low_confidence > ambiguous_root > heterogeneous > mixed_population > weak_leaf > subtype_signal (relabel)

        if result["passes_low_confidence"]:
            return "SUBCLUSTER", f"Low confidence: score={result['assigned_score']:.3f} < {self.config.score_threshold}"

        if result.get("passes_ambiguous_root", False):
            return "SUBCLUSTER", "Ambiguous root: multiple competing cell types at root level"

        if result["passes_heterogeneous"]:
            return "SUBCLUSTER", (
                f"Heterogeneous: {result['best_child']} vs {result['second_child']} "
                f"({result['best_child_score']:.2f} vs {result['second_child_score']:.2f})"
            )

        if result["passes_mixed_population"]:
            return "SUBCLUSTER", (
                f"Mixed population: max_child={result['max_child_score']:.2f} < {self.config.mixed_population_max_child_threshold}"
            )

        if result.get("passes_weak_leaf", False):
            het = result.get("marker_heterogeneity", 0)
            return "SUBCLUSTER", (
                f"Weak leaf: score={result['assigned_score']:.2f}, heterogeneity={het:.3f} > {self.config.weak_leaf_heterogeneity_threshold}"
            )

        if result["passes_subtype_signal"]:
            # Include runner-up info and gap threshold for clearer reasoning
            second_child = result.get("second_child")
            second_score = result.get("second_child_score")
            if second_child and second_score is not None:
                gap = result['best_child_score'] - second_score
                return "RELABEL", (
                    f"Clear subtype: → {result['best_child']} "
                    f"(score={result['best_child_score']:.2f} vs runner-up {second_child}={second_score:.2f}, "
                    f"gap={gap:.2f} ≥ {self.config.heterogeneity_gap})"
                )
            else:
                return "RELABEL", f"Clear subtype: → {result['best_child']} (score={result['best_child_score']:.2f}, no runner-up)"

        return "SKIP", "No refinement criteria met"

    def _classify_confidence_band(self, score: float) -> str:
        """Classify score into confidence bands.

        Bands:
        - very_low: score < 0.5
        - low: 0.5 <= score < 1.0
        - medium: 1.0 <= score < 2.0
        - high: score >= 2.0
        """
        if score < 0.5:
            return "very_low"
        elif score < 1.0:
            return "low"
        elif score < 2.0:
            return "medium"
        else:
            return "high"

    def _get_reason_codes(self, result: Dict[str, Any]) -> str:
        """Get JSON-formatted list of triggered criteria.

        Returns a JSON string with all criteria that passed evaluation.
        """
        codes = []

        if result.get("passes_low_confidence", False):
            codes.append("LOW_CONFIDENCE")
        if result.get("passes_ambiguous_root", False):
            codes.append("AMBIGUOUS_ROOT")
        if result.get("passes_subtype_signal", False):
            codes.append("SUBTYPE_SIGNAL")
        if result.get("passes_heterogeneous", False):
            codes.append("HETEROGENEOUS_CHILDREN")
        if result.get("passes_mixed_population", False):
            codes.append("MIXED_POPULATION")
        if result.get("passes_weak_leaf", False):
            codes.append("WEAK_LEAF")

        # Add constraint codes
        if result.get("n_cells", 0) < self.config.min_cells:
            codes.append("TOO_SMALL")

        if not codes:
            codes.append("NONE")

        return json.dumps(codes)

    @classmethod
    def from_yaml(cls, config_path: Path, logger: Optional[logging.Logger] = None) -> "AutoPolicy":
        """Create AutoPolicy from YAML configuration file.

        YAML format:
        ```yaml
        auto_policy:
          min_cells: 500
          score_threshold: 1.0
          subtype_signal_threshold: 0.5
          heterogeneity_gap: 0.3
          subcluster_resolution: 0.4
          enable_parent_subclustering: true
          mixed_population_parent_threshold: 1.0
          mixed_population_max_child_threshold: 0.3
        ```
        """
        with open(config_path, "r") as f:
            config = yaml.safe_load(f)

        auto_config = config.get("auto_policy", config)
        return cls(
            min_cells=auto_config.get("min_cells", 500),
            score_threshold=auto_config.get("score_threshold", 1.0),
            subtype_signal_threshold=auto_config.get("subtype_signal_threshold", 0.5),
            heterogeneity_gap=auto_config.get("heterogeneity_gap", 0.3),
            subcluster_resolution=auto_config.get("subcluster_resolution", 0.4),
            enable_parent_subclustering=auto_config.get("enable_parent_subclustering", True),
            n_pcs=auto_config.get("n_pcs", 30),
            neighbors_k=auto_config.get("neighbors_k", 15),
            min_subcluster_cells=auto_config.get("min_subcluster_cells", 20),
            add_rescore=auto_config.get("add_rescore", True),
            rescore_mode=auto_config.get("rescore_mode", "smart"),
            mixed_population_parent_threshold=auto_config.get("mixed_population_parent_threshold", 1.0),
            mixed_population_max_child_threshold=auto_config.get("mixed_population_max_child_threshold", 0.3),
            logger=logger,
        )


@dataclass
class _Candidate:
    """Internal dataclass for candidate tracking."""
    cluster_id: str
    assigned_label: str
    assigned_score: float
    n_cells: int
    reason: str  # "low_confidence", "parent_with_subtype_signal", "relabel_to_subtype", "mixed_population", "weak_leaf"
    action: str = "subcluster"  # "subcluster" or "relabel"
    max_child_score: Optional[float] = None
    best_child_label: Optional[str] = None
    second_child_score: Optional[float] = None
    second_child_label: Optional[str] = None
    marker_heterogeneity: Optional[float] = None  # For weak_leaf criterion
