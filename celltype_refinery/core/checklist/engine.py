"""Checklist engine for orchestrating review checklist generation."""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Optional

import pandas as pd

from .config import ChecklistConfig
from .export import export_checklist_md, export_priority_table
from .priority import classify_confidence_band, compute_review_priority
from .recommendation import determine_recommendation


@dataclass
class ChecklistResult:
    """Result from checklist generation.

    Attributes:
        md_path: Path to the generated markdown checklist.
        csv_path: Path to the generated priority table CSV.
        n_clusters: Total number of clusters.
        n_high_priority: Number of high priority clusters.
        n_medium_priority: Number of medium priority clusters.
        n_low_priority: Number of low priority clusters.
        priority_table: DataFrame with priority scores.
        provenance: Provenance information for reproducibility.
    """

    md_path: Path
    csv_path: Path
    n_clusters: int
    n_high_priority: int
    n_medium_priority: int
    n_low_priority: int
    priority_table: pd.DataFrame
    provenance: Dict[str, Any] = field(default_factory=dict)


class ChecklistEngine:
    """Engine for generating review checklists.

    This engine orchestrates the computation of review priorities
    and generation of checklists for reviewing cell-type annotations.

    Example:
        >>> engine = ChecklistEngine(ChecklistConfig(), logger)
        >>> result = engine.execute(
        ...     cluster_annotations=annotations_df,
        ...     output_dir=Path("output/"),
        ... )
        >>> print(f"Generated checklist with {result.n_high_priority} high priority items")
    """

    def __init__(
        self,
        config: Optional[ChecklistConfig] = None,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize checklist engine.

        Args:
            config: Checklist configuration (uses defaults if None).
            logger: Logger instance.
        """
        self.config = config or ChecklistConfig()
        self.logger = logger or logging.getLogger(__name__)

    def execute(
        self,
        cluster_annotations: pd.DataFrame,
        output_dir: Path,
    ) -> ChecklistResult:
        """Generate review checklist.

        Args:
            cluster_annotations: DataFrame with cluster assignments.
            output_dir: Directory to save output files.

        Returns:
            ChecklistResult with generation metadata.
        """
        self.logger.info("Generating review checklist...")

        # Create priority table
        priority_table = self._build_priority_table(cluster_annotations)

        # Count by priority band
        n_high = len(priority_table[priority_table["priority"] >= self.config.high_priority_threshold])
        n_medium = len(priority_table[
            (priority_table["priority"] >= self.config.medium_priority_threshold) &
            (priority_table["priority"] < self.config.high_priority_threshold)
        ])
        n_low = len(priority_table[priority_table["priority"] < self.config.medium_priority_threshold])

        # Export priority table CSV
        csv_path = export_priority_table(
            priority_table=priority_table,
            output_path=output_dir / "priority_table.csv",
            logger=self.logger,
        )

        # Export markdown checklist
        md_path = export_checklist_md(
            priority_table=priority_table,
            output_path=output_dir / "REVIEW_CHECKLIST.md",
            logger=self.logger,
            config=self.config,
        )

        # Build provenance
        provenance = {
            "generated_at": datetime.now().isoformat(),
            "config": {
                "stage_name": self.config.stage_name,
                "high_priority_threshold": self.config.high_priority_threshold,
                "medium_priority_threshold": self.config.medium_priority_threshold,
            },
            "summary": {
                "n_clusters": len(priority_table),
                "n_high_priority": n_high,
                "n_medium_priority": n_medium,
                "n_low_priority": n_low,
            },
        }

        self.logger.info(
            "Generated checklist: %d high, %d medium, %d low priority",
            n_high,
            n_medium,
            n_low,
        )

        return ChecklistResult(
            md_path=md_path,
            csv_path=csv_path,
            n_clusters=len(priority_table),
            n_high_priority=n_high,
            n_medium_priority=n_medium,
            n_low_priority=n_low,
            priority_table=priority_table,
            provenance=provenance,
        )

    def _build_priority_table(self, cluster_annotations: pd.DataFrame) -> pd.DataFrame:
        """Build priority table with scores and recommendations.

        Args:
            cluster_annotations: DataFrame with cluster assignments.

        Returns:
            DataFrame with priority scores and recommendations.
        """
        if cluster_annotations.empty:
            return pd.DataFrame()

        # Select base columns
        base_columns = ["cluster_id", "assigned_label", "assigned_score", "n_cells", "stop_reason"]
        available_columns = [c for c in base_columns if c in cluster_annotations.columns]

        priority_table = cluster_annotations[available_columns].copy()

        # Add confidence band
        if "assigned_score" in priority_table.columns:
            priority_table["confidence_band"] = priority_table["assigned_score"].apply(classify_confidence_band)
        else:
            priority_table["confidence_band"] = "unknown"

        # Compute priority score
        priority_table["priority"] = cluster_annotations.apply(
            lambda row: compute_review_priority(row, self.config),
            axis=1,
        )

        # Compute recommendation
        priority_table["recommendation"] = cluster_annotations.apply(
            determine_recommendation,
            axis=1,
        )

        # Sort by priority descending
        priority_table = priority_table.sort_values("priority", ascending=False)

        return priority_table
