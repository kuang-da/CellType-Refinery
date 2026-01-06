"""Export functions for review checklist.

Ported from ft/src/stage_h_coarse_clustering.py lines 1591-1693.
"""

from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path

import pandas as pd

from .config import ChecklistConfig


def export_priority_table(
    priority_table: pd.DataFrame,
    output_path: Path,
    logger: logging.Logger,
) -> Path:
    """Export priority table to CSV.

    Args:
        priority_table: DataFrame with priority scores.
        output_path: Path to save the CSV file.
        logger: Logger instance.

    Returns:
        Path to the generated CSV file.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    priority_table.to_csv(output_path, index=False)
    logger.info("Wrote priority table: %s", output_path)
    return output_path


def export_checklist_md(
    priority_table: pd.DataFrame,
    output_path: Path,
    logger: logging.Logger,
    config: ChecklistConfig = None,
) -> Path:
    """Generate markdown checklist prioritized by review need.

    Args:
        priority_table: DataFrame with priority scores and recommendations.
        output_path: Path to save the markdown file.
        logger: Logger instance.
        config: Checklist configuration.

    Returns:
        Path to the generated checklist file.
    """
    if config is None:
        config = ChecklistConfig()

    if priority_table.empty:
        logger.warning("Empty priority table; generating empty checklist.")
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(f"# {config.stage_name} Review Checklist\n\nNo clusters to review.\n")
        return output_path

    # Sort by priority_score descending (highest priority first)
    if "priority" in priority_table.columns:
        sorted_df = priority_table.sort_values("priority", ascending=False)
    else:
        sorted_df = priority_table.copy()

    # Group by priority bands
    high_priority = sorted_df[sorted_df.get("priority", 0) >= config.high_priority_threshold] if "priority" in sorted_df.columns else pd.DataFrame()
    medium_priority = sorted_df[
        (sorted_df.get("priority", 0) >= config.medium_priority_threshold) &
        (sorted_df.get("priority", 0) < config.high_priority_threshold)
    ] if "priority" in sorted_df.columns else pd.DataFrame()
    low_priority = sorted_df[sorted_df.get("priority", 0) < config.medium_priority_threshold] if "priority" in sorted_df.columns else sorted_df

    lines = [
        f"# {config.stage_name} Review Checklist",
        "",
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        f"**Total clusters**: {len(sorted_df)}",
        f"**High priority**: {len(high_priority)}",
        f"**Medium priority**: {len(medium_priority)}",
        f"**Low priority**: {len(low_priority)}",
        "",
    ]

    def format_cluster_item(row):
        """Format a single cluster as a checklist item."""
        cluster_id = row.get("cluster_id", "?")
        label = row.get("assigned_label", "Unknown")
        score = row.get("assigned_score", 0)
        n_cells = row.get("n_cells", 0)
        confidence = row.get("confidence_band", "?")
        action = row.get("recommendation", "review")
        stop_reason = row.get("stop_reason", "")

        item = [
            f"- [ ] **Cluster {cluster_id}**: {label} (score={score:.2f}, cells={n_cells:,})",
            f"  - Confidence: {confidence}",
            f"  - Action: {action}",
        ]
        if stop_reason:
            item.append(f"  - Stop reason: {stop_reason}")
        return "\n".join(item)

    # High priority section
    if not high_priority.empty:
        lines.append(f"## High Priority (Score >= {config.high_priority_threshold})")
        lines.append("")
        lines.append("These clusters require immediate attention due to low confidence, large size, or ambiguous assignment.")
        lines.append("")
        for _, row in high_priority.iterrows():
            lines.append(format_cluster_item(row))
            lines.append("")

    # Medium priority section
    if not medium_priority.empty:
        lines.append(f"## Medium Priority (Score {config.medium_priority_threshold}-{config.high_priority_threshold - 1})")
        lines.append("")
        lines.append("These clusters have moderate concerns and should be reviewed.")
        lines.append("")
        for _, row in medium_priority.iterrows():
            lines.append(format_cluster_item(row))
            lines.append("")

    # Low priority section
    if not low_priority.empty:
        lines.append(f"## Low Priority (Score < {config.medium_priority_threshold})")
        lines.append("")
        lines.append("These clusters have high confidence and can likely be auto-approved.")
        lines.append("")
        for _, row in low_priority.head(config.low_priority_max_items).iterrows():
            lines.append(format_cluster_item(row))
            lines.append("")
        if len(low_priority) > config.low_priority_max_items:
            lines.append(f"... and {len(low_priority) - config.low_priority_max_items} more low-priority clusters.")
            lines.append("")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines))
    logger.info("Generated review checklist with %d clusters to %s", len(sorted_df), output_path)
    return output_path
