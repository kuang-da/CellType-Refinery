"""Orphan rescue visualization functions.

Provides:
- Rescue summary bar chart (rescued orphans by subtype)
- Orphan decision pie chart (rescue vs flagged vs unassigned)
- Root vs Subtype score scatter plot
"""

from pathlib import Path
from typing import Dict, Optional, Union
import logging

import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


def plot_orphan_rescue_summary(
    orphan_df: pd.DataFrame,
    output_path: Union[str, Path],
    dpi: int = 200,
) -> Path:
    """Plot bar chart of rescued orphans by subtype.

    Args:
        orphan_df: DataFrame with orphan candidates
            Expected columns: subtype, action, n_cells
        output_path: Path to save figure
        dpi: Figure resolution

    Returns:
        Path to saved figure
    """
    import matplotlib.pyplot as plt
    from .style import set_publication_style, save_figure, ORPHAN_ACTION_COLORS

    set_publication_style()

    if orphan_df is None or len(orphan_df) == 0:
        logger.warning("No orphan data to visualize")
        return None

    # Filter to rescued orphans
    action_col = "action" if "action" in orphan_df.columns else "rescue_action"
    if action_col not in orphan_df.columns:
        logger.warning(f"Action column not found in orphan_df")
        return None

    rescued = orphan_df[orphan_df[action_col] == "rescue"].copy()

    if len(rescued) == 0:
        logger.warning("No rescued orphans to visualize")
        return None

    # Aggregate by subtype
    n_cells_col = "n_cells" if "n_cells" in rescued.columns else None
    if n_cells_col:
        counts = rescued.groupby("subtype")[n_cells_col].sum().sort_values(ascending=True)
    else:
        counts = rescued["subtype"].value_counts().sort_values(ascending=True)

    # Create horizontal bar chart
    fig, ax = plt.subplots(figsize=(10, max(6, len(counts) * 0.4)))

    colors = []
    for subtype in counts.index:
        # Color by parent lineage
        if "Endothelium" in subtype or "Lymphatic" in subtype:
            colors.append("#e74c3c")
        elif "Muscle" in subtype or "Fibroblast" in subtype or "Mesenchymal" in subtype:
            colors.append("#3498db")
        elif "Pericyte" in subtype or "Prolif" in subtype:
            colors.append("#95a5a6")
        else:
            colors.append("#27ae60")

    bars = ax.barh(range(len(counts)), counts.values, color=colors)

    # Add value labels
    for bar, count in zip(bars, counts.values):
        if count >= 1000:
            label = f"{count/1000:.1f}K"
        else:
            label = str(count)
        ax.text(
            bar.get_width() + max(counts) * 0.02,
            bar.get_y() + bar.get_height() / 2,
            label,
            va="center",
            fontsize=9,
        )

    ax.set_yticks(range(len(counts)))
    ax.set_yticklabels(counts.index)
    ax.set_xlabel("Number of Cells")
    ax.set_title(f"Rescued Orphans by Subtype (n={counts.sum():,} cells)")
    ax.set_xlim(0, max(counts) * 1.2)

    plt.tight_layout()
    return save_figure(fig, output_path, dpi=dpi)


def plot_orphan_decision_pie(
    orphan_df: pd.DataFrame,
    output_path: Union[str, Path],
    dpi: int = 200,
) -> Path:
    """Plot pie chart of orphan fate decisions.

    Args:
        orphan_df: DataFrame with orphan candidates
            Expected columns: action, n_cells
        output_path: Path to save figure
        dpi: Figure resolution

    Returns:
        Path to saved figure
    """
    import matplotlib.pyplot as plt
    from .style import set_publication_style, save_figure, ORPHAN_ACTION_COLORS

    set_publication_style()

    if orphan_df is None or len(orphan_df) == 0:
        logger.warning("No orphan data to visualize")
        return None

    # Get action column
    action_col = "action" if "action" in orphan_df.columns else "rescue_action"
    if action_col not in orphan_df.columns:
        logger.warning(f"Action column not found")
        return None

    # Aggregate by action
    n_cells_col = "n_cells" if "n_cells" in orphan_df.columns else None
    if n_cells_col:
        action_counts = orphan_df.groupby(action_col)[n_cells_col].sum()
    else:
        action_counts = orphan_df[action_col].value_counts()

    # Map action names to display names
    action_labels = {
        "rescue": "Rescued",
        "keep_unassigned": "Kept Unassigned",
        "expert_review": "Flagged for Review",
    }

    labels = [action_labels.get(a, a.replace("_", " ").title()) for a in action_counts.index]
    colors = [ORPHAN_ACTION_COLORS.get(a, "#bdc3c7") for a in action_counts.index]

    # Create pie chart
    fig, ax = plt.subplots(figsize=(8, 8))

    wedges, texts, autotexts = ax.pie(
        action_counts.values,
        labels=labels,
        colors=colors,
        autopct=lambda pct: f"{pct:.1f}%\n({int(pct/100*action_counts.sum()):,})",
        startangle=90,
        explode=[0.02] * len(action_counts),
        textprops={"fontsize": 10},
    )

    # Make percentage text bold
    for autotext in autotexts:
        autotext.set_fontweight("bold")

    ax.set_title(f"Orphan Fate Decisions\n(n={action_counts.sum():,} clusters)")

    plt.tight_layout()
    return save_figure(fig, output_path, dpi=dpi)


def plot_orphan_score_scatter(
    orphan_df: pd.DataFrame,
    output_path: Union[str, Path],
    all_unassigned_scores: Optional[pd.DataFrame] = None,
    dpi: int = 200,
) -> Path:
    """Plot scatter of root score vs subtype score for orphans.

    Args:
        orphan_df: DataFrame with orphan candidates
            Expected columns: root_score, subtype_score, action
        output_path: Path to save figure
        all_unassigned_scores: Optional DataFrame with ALL Unassigned cluster scores
            (including those that passed root gate or had low subtype score).
            If provided, non-orphan clusters shown in gray background.
            Expected columns: cluster_id, root_score, subtype_score, n_cells
        dpi: Figure resolution

    Returns:
        Path to saved figure
    """
    import matplotlib.pyplot as plt
    from .style import set_publication_style, save_figure, ORPHAN_ACTION_COLORS

    set_publication_style()

    if orphan_df is None or len(orphan_df) == 0:
        logger.warning("No orphan data to visualize")
        return None

    # Check required columns
    required = ["root_score", "subtype_score"]
    for col in required:
        if col not in orphan_df.columns:
            logger.warning(f"Column '{col}' not found in orphan_df")
            return None

    action_col = "action" if "action" in orphan_df.columns else "rescue_action"
    subtype_col = "subtype" if "subtype" in orphan_df.columns else None

    # Create scatter plot
    fig, ax = plt.subplots(figsize=(12, 9))

    # First, plot non-orphan clusters in gray background (if provided)
    if all_unassigned_scores is not None and len(all_unassigned_scores) > 0:
        # Get cluster IDs that are orphan candidates
        orphan_cluster_ids = set(orphan_df["cluster_id"].astype(str)) if "cluster_id" in orphan_df.columns else set()

        # Filter to non-orphan clusters
        non_orphan_mask = ~all_unassigned_scores["cluster_id"].astype(str).isin(orphan_cluster_ids)
        non_orphan_df = all_unassigned_scores[non_orphan_mask]

        if len(non_orphan_df) > 0:
            # Size by n_cells if available
            if "n_cells" in non_orphan_df.columns:
                sizes = np.clip(non_orphan_df["n_cells"] / 100, 10, 150)
            else:
                sizes = 30

            ax.scatter(
                non_orphan_df["root_score"],
                non_orphan_df["subtype_score"],
                c="#7f8c8d",  # Darker gray
                s=sizes,
                alpha=0.5,
                label=f"Other Unassigned ({len(non_orphan_df)})",
                edgecolors="white",
                linewidth=0.3,
                zorder=1,  # Behind orphan candidates
            )

    # Plot orphan candidates by action (on top)
    actions = orphan_df[action_col].unique() if action_col in orphan_df.columns else ["unknown"]

    for action in actions:
        if action_col in orphan_df.columns:
            mask = orphan_df[action_col] == action
            subset = orphan_df[mask]
        else:
            subset = orphan_df

        color = ORPHAN_ACTION_COLORS.get(action, "#bdc3c7")
        label = action.replace("_", " ").title()

        # Size by n_cells if available
        if "n_cells" in subset.columns:
            sizes = np.clip(subset["n_cells"] / 100, 10, 200)
        else:
            sizes = 50

        ax.scatter(
            subset["root_score"],
            subset["subtype_score"],
            c=color,
            s=sizes,
            alpha=0.6,
            label=label,
            edgecolors="white",
            linewidth=0.5,
            zorder=2,  # On top of gray dots
        )

        # Add subtype labels for clarity (only for larger clusters)
        if subtype_col and "n_cells" in subset.columns:
            for _, row in subset.iterrows():
                if row["n_cells"] > 1000:  # Only label significant clusters
                    # Truncate long subtype names
                    subtype_short = row[subtype_col][:12] + "..." if len(row[subtype_col]) > 12 else row[subtype_col]
                    ax.annotate(
                        subtype_short,
                        (row["root_score"], row["subtype_score"]),
                        fontsize=7,
                        alpha=0.7,
                        xytext=(3, 3),
                        textcoords="offset points",
                        zorder=3,
                    )

    # Compute axis limits (include all_unassigned_scores if provided)
    x_max = orphan_df["root_score"].max()
    y_max = orphan_df["subtype_score"].max()
    if all_unassigned_scores is not None and len(all_unassigned_scores) > 0:
        x_max = max(x_max, all_unassigned_scores["root_score"].max())
        y_max = max(y_max, all_unassigned_scores["subtype_score"].max())
    x_max = max(x_max * 1.1, 1.0)
    y_max = max(y_max * 1.1, 1.0)

    # Add threshold lines
    ax.axhline(y=0.8, color="#27ae60", linestyle="--", alpha=0.7, label="Subtype threshold (0.8)")
    ax.axvline(x=0.5, color="#e74c3c", linestyle="--", alpha=0.7, label="Root threshold (0.5)")

    # Shade rescue zone: low root (< 0.5) AND high subtype (> 0.8)
    ax.fill_between(
        [-0.1, 0.5],              # x: from left edge to root threshold
        [0.8, 0.8],               # bottom: subtype threshold
        [y_max, y_max],           # top: extend to y-axis max
        alpha=0.1,
        color="#27ae60",
        label="Rescue zone",
    )

    ax.set_xlabel("Root Score (parent lineage)", fontsize=11)
    ax.set_ylabel("Subtype Score", fontsize=11)
    ax.set_title(
        "Orphan Detection: Root vs Subtype Scores\n"
        "(Immune subtypes kept Unassigned - CD45 should be present)",
        fontsize=12,
    )
    ax.set_xlim(-0.1, x_max)
    ax.set_ylim(-0.1, y_max)
    ax.legend(loc="lower right", fontsize=9)

    # Add annotations (positioned in middle of each zone)
    zone_mid_y = (0.8 + y_max) / 2  # middle of "high subtype" region
    ax.annotate(
        "RESCUE ZONE\n(non-immune orphans)",
        xy=(0.2, zone_mid_y),
        fontsize=10,
        ha="center",
        color="#27ae60",
        fontweight="bold",
    )
    ax.annotate(
        "NORMAL\n(passed root gate)",
        xy=(0.75, zone_mid_y),
        fontsize=10,
        ha="center",
        color="#3498db",
    )

    plt.tight_layout()
    return save_figure(fig, output_path, dpi=dpi)
