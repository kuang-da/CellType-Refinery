"""Audit engine for orchestrating audit card generation."""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Optional

import pandas as pd

from .card import generate_all_audit_cards
from .config import AuditConfig


@dataclass
class AuditResult:
    """Result from audit card generation.

    Attributes:
        html_path: Path to the generated HTML file.
        n_cards: Number of audit cards generated.
        n_high_confidence: Number of clusters with high confidence.
        n_needs_review: Number of clusters needing review.
        provenance: Provenance information for reproducibility.
    """

    html_path: Path
    n_cards: int
    n_high_confidence: int
    n_needs_review: int
    provenance: Dict[str, Any] = field(default_factory=dict)


class AuditEngine:
    """Engine for generating audit cards.

    This engine orchestrates the generation of HTML audit cards
    for reviewing cell-type annotations.

    Example:
        >>> engine = AuditEngine(AuditConfig(), logger)
        >>> result = engine.execute(
        ...     cluster_annotations=annotations_df,
        ...     decision_steps=steps_df,
        ...     marker_scores=scores_df,
        ...     marker_evidence=evidence_df,
        ...     output_dir=Path("output/"),
        ... )
        >>> print(f"Generated {result.n_cards} audit cards")
    """

    def __init__(
        self,
        config: Optional[AuditConfig] = None,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize audit engine.

        Args:
            config: Audit configuration (uses defaults if None).
            logger: Logger instance.
        """
        self.config = config or AuditConfig()
        self.logger = logger or logging.getLogger(__name__)

    def execute(
        self,
        cluster_annotations: pd.DataFrame,
        decision_steps: Optional[pd.DataFrame],
        marker_scores: pd.DataFrame,
        marker_evidence: Optional[pd.DataFrame],
        output_dir: Path,
    ) -> AuditResult:
        """Generate audit cards.

        Args:
            cluster_annotations: DataFrame with cluster assignments.
            decision_steps: DataFrame with decision steps (can be None).
            marker_scores: DataFrame with marker scores.
            marker_evidence: DataFrame with per-marker evidence (optional).
            output_dir: Directory to save the HTML file.

        Returns:
            AuditResult with generation metadata.
        """
        self.logger.info("Generating audit cards...")

        # Generate the HTML file
        html_path = generate_all_audit_cards(
            cluster_annotations=cluster_annotations,
            decision_steps=decision_steps if self.config.include_decision_steps else None,
            marker_scores=marker_scores,
            marker_evidence=marker_evidence if self.config.include_marker_evidence else None,
            output_dir=output_dir,
            logger=self.logger,
            stage_name=self.config.stage_name,
        )

        # Compute summary stats
        n_cards = len(cluster_annotations)
        n_high_confidence = 0
        n_needs_review = 0

        if "assigned_score" in cluster_annotations.columns:
            n_high_confidence = len(
                cluster_annotations[
                    cluster_annotations["assigned_score"] >= self.config.high_confidence_threshold
                ]
            )

        if "stopped_before_leaf" in cluster_annotations.columns:
            n_needs_review = len(
                cluster_annotations[cluster_annotations["stopped_before_leaf"] == True]
            )

        # Build provenance
        provenance = {
            "generated_at": datetime.now().isoformat(),
            "config": {
                "stage_name": self.config.stage_name,
                "include_marker_evidence": self.config.include_marker_evidence,
                "include_decision_steps": self.config.include_decision_steps,
                "high_confidence_threshold": self.config.high_confidence_threshold,
            },
            "input_shapes": {
                "cluster_annotations": cluster_annotations.shape,
                "decision_steps": decision_steps.shape if decision_steps is not None else None,
                "marker_scores": marker_scores.shape,
                "marker_evidence": marker_evidence.shape if marker_evidence is not None else None,
            },
            "summary": {
                "n_cards": n_cards,
                "n_high_confidence": n_high_confidence,
                "n_needs_review": n_needs_review,
            },
        }

        self.logger.info(
            "Generated %d audit cards (%d high confidence, %d needs review)",
            n_cards,
            n_high_confidence,
            n_needs_review,
        )

        return AuditResult(
            html_path=html_path,
            n_cards=n_cards,
            n_high_confidence=n_high_confidence,
            n_needs_review=n_needs_review,
            provenance=provenance,
        )
