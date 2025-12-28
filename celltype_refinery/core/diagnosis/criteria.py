"""Diagnostic criteria configuration and evaluation for cluster refinement.

This module provides the core criteria evaluation logic used for diagnostic
reports that recommend SUBCLUSTER, RELABEL, or SKIP actions.
"""

from __future__ import annotations

import json
import logging
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Set, Tuple

import pandas as pd


@dataclass
class CriteriaConfig:
    """Configuration for diagnostic criteria thresholds.

    These thresholds control when clusters are recommended for
    SUBCLUSTER, RELABEL, or SKIP.
    """

    # Core thresholds
    min_cells: int = 500
    score_threshold: float = 1.0
    subtype_signal_threshold: float = 1.0
    heterogeneity_gap: float = 0.5

    # Subclustering parameters
    subcluster_resolution: float = 0.4
    n_pcs: int = 30
    neighbors_k: int = 15
    min_subcluster_cells: int = 20

    # Parent subclustering
    enable_parent_subclustering: bool = True

    # Mixed population detection
    mixed_population_parent_threshold: float = 1.0
    mixed_population_max_child_threshold: float = 1.0

    # Weak leaf detection
    weak_leaf_score_threshold: float = 1.5
    weak_leaf_heterogeneity_threshold: float = 0.12
    enable_weak_leaf_detection: bool = True


class CriteriaEvaluator:
    """Evaluates refinement criteria for clusters.

    Example:
        >>> config = CriteriaConfig(score_threshold=1.0)
        >>> evaluator = CriteriaEvaluator(config)
        >>> result = evaluator.evaluate(row, marker_scores, hierarchy)
        >>> print(result["recommendation"])
        'SUBCLUSTER'
    """

    def __init__(
        self,
        config: Optional[CriteriaConfig] = None,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize evaluator with configuration.

        Parameters
        ----------
        config : CriteriaConfig, optional
            Thresholds and parameters. Uses defaults if not provided.
        logger : logging.Logger, optional
            Logger for debug messages.
        """
        self.config = config or CriteriaConfig()
        self.logger = logger or logging.getLogger(__name__)
        self._heterogeneity_metrics: Optional[pd.DataFrame] = None

    def set_heterogeneity_metrics(self, metrics: pd.DataFrame) -> None:
        """Set pre-computed marker heterogeneity metrics.

        Parameters
        ----------
        metrics : pd.DataFrame
            DataFrame with columns: cluster_id, marker_heterogeneity.
            Output from compute_cluster_marker_heterogeneity().
        """
        if "cluster_id" not in metrics.columns:
            raise ValueError("metrics must have 'cluster_id' column")
        if "marker_heterogeneity" not in metrics.columns:
            raise ValueError("metrics must have 'marker_heterogeneity' column")

        self._heterogeneity_metrics = metrics.copy()
        self._heterogeneity_metrics["cluster_id"] = (
            self._heterogeneity_metrics["cluster_id"].astype(str)
        )

    def evaluate(
        self,
        row: pd.Series,
        marker_scores: pd.DataFrame,
        hierarchy: Dict[str, List[str]],
    ) -> Dict[str, Any]:
        """Evaluate all criteria for a single cluster.

        Parameters
        ----------
        row : pd.Series
            Row from cluster_annotations DataFrame.
        marker_scores : pd.DataFrame
            Full marker scores for all clusters.
        hierarchy : Dict[str, List[str]]
            Mapping from parent label to list of child labels.

        Returns
        -------
        Dict[str, Any]
            Evaluation result with all diagnostic fields.
        """
        cluster_id = str(row["cluster_id"])
        assigned_label = row["assigned_label"]
        assigned_score = float(row.get("assigned_score", 0.0) or 0.0)
        n_cells = int(row["n_cells"])

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
            # Recommendation (final, after constraints applied)
            "recommendation": "SKIP",
            "recommendation_reason": "",
            # Ideal recommendation (before n_cells constraint)
            "ideal_recommendation": "SKIP",
            # Constraint flags
            "is_blocked_by_size": False,
            # Audit columns
            "confidence_band": "",
            "triggered_criteria": "",
        }

        # NOTE: n_cells check moved to END of evaluation (after _get_recommendation)
        # This allows us to compute ideal_recommendation first, then apply constraints

        # Criterion 1: Low confidence
        result["passes_low_confidence"] = assigned_score < self.config.score_threshold

        # Get child data if parent
        if result["is_parent"]:
            self._evaluate_parent_criteria(result, row, marker_scores, hierarchy)

        # Criterion 4: Ambiguous root (from input or detected via ~ in label)
        is_ambiguous = row.get("is_ambiguous_root", False)
        if not is_ambiguous and "~" in str(row.get("assigned_label", "")):
            is_ambiguous = True
        if is_ambiguous:
            result["passes_ambiguous_root"] = True

        # Criterion 5: Weak leaf
        if self.config.enable_weak_leaf_detection and not result["is_parent"]:
            self._evaluate_weak_leaf(result, cluster_id, assigned_score)

        # Generate IDEAL recommendation (before applying n_cells constraint)
        ideal_rec, ideal_reason = self._get_recommendation(result)
        result["ideal_recommendation"] = ideal_rec

        # Apply n_cells constraint as FINAL GATE
        # Only SUBCLUSTER is constrained by n_cells (RELABEL and SKIP don't need many cells)
        if ideal_rec == "SUBCLUSTER" and n_cells < self.config.min_cells:
            result["is_blocked_by_size"] = True
            result["recommendation"] = "SKIP"
            result["recommendation_reason"] = (
                f"Would {ideal_rec} ({ideal_reason}) but too few cells "
                f"({n_cells} < {self.config.min_cells})"
            )
        else:
            result["recommendation"] = ideal_rec
            result["recommendation_reason"] = ideal_reason

        # Populate audit columns
        result["confidence_band"] = self._classify_confidence_band(assigned_score)
        result["triggered_criteria"] = self._get_reason_codes(result)

        return result

    def _evaluate_parent_criteria(
        self,
        result: Dict[str, Any],
        row: pd.Series,
        marker_scores: pd.DataFrame,
        hierarchy: Dict[str, List[str]],
    ) -> None:
        """Evaluate parent-specific criteria (subtype signal, heterogeneity, mixed population)."""
        assigned_label = result["assigned_label"]
        children = hierarchy[assigned_label]
        result["children"] = ";".join(children)

        cluster_scores = marker_scores[marker_scores["cluster_id"] == row["cluster_id"]]
        child_scores = cluster_scores[cluster_scores["label"].isin(children)].copy()

        if child_scores.empty:
            return

        # Sort by score
        child_scores = child_scores.sort_values("score", ascending=False)

        result["max_child_score"] = float(child_scores["score"].max())
        result["best_child"] = child_scores.iloc[0]["label"]
        result["best_child_score"] = float(child_scores.iloc[0]["score"])

        if len(child_scores) >= 2:
            result["second_child"] = child_scores.iloc[1]["label"]
            result["second_child_score"] = float(child_scores.iloc[1]["score"])

        result["all_child_scores"] = ";".join(
            f"{r['label']}={r['score']:.2f}" for _, r in child_scores.iterrows()
        )

        # Criterion 2: Subtype signal (positive children above threshold)
        passing_children = child_scores[child_scores["score"] > 0]
        if not passing_children.empty:
            best_passing_score = float(passing_children.iloc[0]["score"])
            if best_passing_score > self.config.subtype_signal_threshold:
                result["passes_subtype_signal"] = True

                # Check heterogeneity
                if len(passing_children) >= 2:
                    second_passing_score = float(passing_children.iloc[1]["score"])
                    gap = best_passing_score - second_passing_score
                    if (
                        second_passing_score >= self.config.subtype_signal_threshold
                        or gap < self.config.heterogeneity_gap
                    ):
                        result["passes_heterogeneous"] = True

        # Criterion 3: Mixed population
        if (
            result["assigned_score"] >= self.config.mixed_population_parent_threshold
            and result["max_child_score"] is not None
            and result["max_child_score"] < self.config.mixed_population_max_child_threshold
        ):
            result["passes_mixed_population"] = True

    def _evaluate_weak_leaf(
        self,
        result: Dict[str, Any],
        cluster_id: str,
        assigned_score: float,
    ) -> None:
        """Evaluate weak leaf criterion (borderline score + high heterogeneity)."""
        het = self._get_heterogeneity(cluster_id)
        if het is not None:
            result["marker_heterogeneity"] = het
            if (
                assigned_score < self.config.weak_leaf_score_threshold
                and het > self.config.weak_leaf_heterogeneity_threshold
            ):
                result["passes_weak_leaf"] = True

    def _get_heterogeneity(self, cluster_id: str) -> Optional[float]:
        """Get marker heterogeneity for a cluster."""
        if self._heterogeneity_metrics is None:
            return None
        row = self._heterogeneity_metrics[
            self._heterogeneity_metrics["cluster_id"] == str(cluster_id)
        ]
        if row.empty:
            return None
        return float(row.iloc[0]["marker_heterogeneity"])

    def _get_recommendation(self, result: Dict[str, Any]) -> Tuple[str, str]:
        """Get recommendation based on evaluated criteria.

        Priority order:
        1. low_confidence → SUBCLUSTER
        2. ambiguous_root → SUBCLUSTER
        3. heterogeneous → SUBCLUSTER
        4. mixed_population → SUBCLUSTER
        5. weak_leaf → SUBCLUSTER
        6. subtype_signal → RELABEL
        7. none → SKIP
        """
        cfg = self.config

        if result["passes_low_confidence"]:
            return (
                "SUBCLUSTER",
                f"Low confidence: score={result['assigned_score']:.3f} < {cfg.score_threshold}",
            )

        if result.get("passes_ambiguous_root", False):
            return (
                "SUBCLUSTER",
                "Ambiguous root: multiple competing cell types at root level",
            )

        if result["passes_heterogeneous"]:
            return (
                "SUBCLUSTER",
                f"Heterogeneous: {result['best_child']} vs {result['second_child']} "
                f"({result['best_child_score']:.2f} vs {result['second_child_score']:.2f})",
            )

        if result["passes_mixed_population"]:
            return (
                "SUBCLUSTER",
                f"Mixed population: max_child={result['max_child_score']:.2f} "
                f"< {cfg.mixed_population_max_child_threshold}",
            )

        if result.get("passes_weak_leaf", False):
            het = result.get("marker_heterogeneity", 0)
            return (
                "SUBCLUSTER",
                f"Weak leaf: score={result['assigned_score']:.2f}, "
                f"heterogeneity={het:.3f} > {cfg.weak_leaf_heterogeneity_threshold}",
            )

        if result["passes_subtype_signal"]:
            second_child = result.get("second_child")
            second_score = result.get("second_child_score")
            if second_child and second_score is not None:
                gap = result["best_child_score"] - second_score
                return (
                    "RELABEL",
                    f"Clear subtype: → {result['best_child']} "
                    f"(score={result['best_child_score']:.2f} vs {second_child}={second_score:.2f}, "
                    f"gap={gap:.2f} >= {cfg.heterogeneity_gap})",
                )
            else:
                return (
                    "RELABEL",
                    f"Clear subtype: → {result['best_child']} "
                    f"(score={result['best_child_score']:.2f}, no runner-up)",
                )

        return "SKIP", "No refinement criteria met"

    def _classify_confidence_band(self, score: float) -> str:
        """Classify score into confidence bands."""
        if score < 0.5:
            return "very_low"
        elif score < 1.0:
            return "low"
        elif score < 2.0:
            return "medium"
        else:
            return "high"

    def _get_reason_codes(self, result: Dict[str, Any]) -> str:
        """Get JSON-formatted list of triggered criteria and blocking constraints.

        Format: ["CRITERION1", "CRITERION2", "BLOCKED_CONSTRAINT"]
        - Criteria codes: what triggered the ideal recommendation
        - BLOCKED_ prefix: constraints that prevented the ideal recommendation

        Examples:
            - ["NONE"] - no criteria triggered, SKIP
            - ["LOW_CONFIDENCE"] - low score triggered SUBCLUSTER
            - ["LOW_CONFIDENCE", "BLOCKED_TOO_SMALL"] - would SUBCLUSTER but blocked by size
            - ["SUBTYPE_SIGNAL"] - clear subtype triggered RELABEL
        """
        codes = []

        # Add triggered criteria (what would cause action)
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

        # Add blocking constraints (what prevented ideal action)
        if result.get("is_blocked_by_size", False):
            codes.append("BLOCKED_TOO_SMALL")

        if not codes:
            codes.append("NONE")

        return json.dumps(codes)


def build_hierarchy_from_scores(marker_scores: pd.DataFrame) -> Dict[str, List[str]]:
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

    return {k: list(v) for k, v in parents_children.items()}
