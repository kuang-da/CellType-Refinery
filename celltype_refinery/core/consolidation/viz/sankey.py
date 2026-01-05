"""Sankey diagram visualizations for label flow.

Provides:
- Static Sankey diagram (matplotlib)
- Interactive Sankey diagram (plotly)
"""

from pathlib import Path
from typing import Dict, List, Optional, Union
import logging

import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


def plot_sankey_label_flow(
    mapping_df: pd.DataFrame,
    output_path: Union[str, Path],
    source_col: str = "original_label",
    target_col: str = "final_label",
    value_col: str = "n_cells",
    top_n: int = 15,
    dpi: int = 200,
) -> Path:
    """Plot static Sankey diagram of label flow.

    Args:
        mapping_df: DataFrame with label mapping
        output_path: Path to save figure
        source_col: Column with original labels
        target_col: Column with final labels
        value_col: Column with cell counts
        top_n: Number of top labels to show
        dpi: Figure resolution

    Returns:
        Path to saved figure
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
    from .style import set_publication_style, save_figure, get_color_palette

    set_publication_style()

    # Ensure required columns exist
    if source_col not in mapping_df.columns:
        # Try to infer from available columns
        if "assigned_label" in mapping_df.columns:
            source_col = "assigned_label"
        else:
            logger.warning(f"Source column not found: {source_col}")
            return None

    if target_col not in mapping_df.columns:
        if "final_label" in mapping_df.columns:
            target_col = "final_label"
        else:
            logger.warning(f"Target column not found: {target_col}")
            return None

    # Aggregate flows
    if value_col in mapping_df.columns:
        flows = mapping_df.groupby([source_col, target_col])[value_col].sum().reset_index()
    else:
        flows = mapping_df.groupby([source_col, target_col]).size().reset_index(name="n_cells")
        value_col = "n_cells"

    # Get top sources and targets
    source_totals = flows.groupby(source_col)[value_col].sum().nlargest(top_n)
    target_totals = flows.groupby(target_col)[value_col].sum().nlargest(top_n)

    top_sources = set(source_totals.index)
    top_targets = set(target_totals.index)

    # Filter to top labels
    flows = flows[
        flows[source_col].isin(top_sources) &
        flows[target_col].isin(top_targets)
    ]

    if len(flows) == 0:
        logger.warning("No flows to visualize after filtering")
        return None

    # Create figure using matplotlib (simplified version)
    fig, ax = plt.subplots(figsize=(14, 10))

    # Position nodes
    sources = sorted(flows[source_col].unique())
    targets = sorted(flows[target_col].unique())

    source_y = {s: i for i, s in enumerate(sources)}
    target_y = {t: i for i, t in enumerate(targets)}

    # Scale y positions
    max_y = max(len(sources), len(targets))
    source_scale = max_y / max(len(sources), 1)
    target_scale = max_y / max(len(targets), 1)

    # Get colors
    all_labels = list(set(sources) | set(targets))
    colors = get_color_palette(all_labels, palette_type="cell_type")

    # Draw flows as bezier curves
    for _, row in flows.iterrows():
        src = row[source_col]
        tgt = row[target_col]
        val = row[value_col]

        y0 = source_y[src] * source_scale
        y1 = target_y[tgt] * target_scale

        # Line width proportional to value
        width = np.sqrt(val) / 50
        width = max(0.5, min(width, 10))

        # Draw curved line
        from matplotlib.path import Path as MPath
        from matplotlib.patches import PathPatch

        verts = [
            (0.1, y0),
            (0.35, y0),
            (0.65, y1),
            (0.9, y1),
        ]
        codes = [MPath.MOVETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4]

        path = MPath(verts, codes)
        patch = PathPatch(
            path,
            facecolor="none",
            edgecolor=colors.get(src, "#bdc3c7"),
            linewidth=width,
            alpha=0.4,
        )
        ax.add_patch(patch)

    # Draw source nodes
    for src in sources:
        y = source_y[src] * source_scale
        total = flows[flows[source_col] == src][value_col].sum()

        rect = FancyBboxPatch(
            (0, y - 0.3),
            0.08,
            0.6,
            boxstyle="round,pad=0.01",
            facecolor=colors.get(src, "#bdc3c7"),
            edgecolor="white",
            linewidth=1,
        )
        ax.add_patch(rect)

        # Label
        label = src if len(src) < 20 else src[:17] + "..."
        ax.text(-0.02, y, f"{label}\n({total:,})", ha="right", va="center", fontsize=8)

    # Draw target nodes
    for tgt in targets:
        y = target_y[tgt] * target_scale
        total = flows[flows[target_col] == tgt][value_col].sum()

        rect = FancyBboxPatch(
            (0.92, y - 0.3),
            0.08,
            0.6,
            boxstyle="round,pad=0.01",
            facecolor=colors.get(tgt, "#bdc3c7"),
            edgecolor="white",
            linewidth=1,
        )
        ax.add_patch(rect)

        # Label
        label = tgt if len(tgt) < 20 else tgt[:17] + "..."
        ax.text(1.02, y, f"{label}\n({total:,})", ha="left", va="center", fontsize=8)

    ax.set_xlim(-0.3, 1.3)
    ax.set_ylim(-1, max_y + 1)
    ax.set_title("Label Flow: Before -> After Consolidation", fontsize=14, pad=20)
    ax.axis("off")

    # Add headers with column names
    ax.text(0.04, max_y + 0.5, f"BEFORE\n({source_col})", ha="center", fontsize=12, fontweight="bold")
    ax.text(0.96, max_y + 0.5, f"AFTER\n({target_col})", ha="center", fontsize=12, fontweight="bold")

    plt.tight_layout()
    return save_figure(fig, output_path, dpi=dpi)


def plot_sankey_interactive(
    mapping_df: pd.DataFrame,
    output_path: Union[str, Path],
    source_col: str = "original_label",
    target_col: str = "final_label",
    value_col: str = "n_cells",
    top_n: int = 20,
) -> Path:
    """Create interactive Sankey diagram using plotly.

    Args:
        mapping_df: DataFrame with label mapping
        output_path: Path to save HTML
        source_col: Column with original labels
        target_col: Column with final labels
        value_col: Column with cell counts
        top_n: Number of top labels to show

    Returns:
        Path to saved HTML file
    """
    try:
        import plotly.graph_objects as go
    except ImportError:
        logger.warning("plotly not available, skipping interactive Sankey")
        return None

    from .style import get_color_palette

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Ensure required columns
    if source_col not in mapping_df.columns:
        if "assigned_label" in mapping_df.columns:
            source_col = "assigned_label"
        else:
            logger.warning(f"Source column not found")
            return None

    if target_col not in mapping_df.columns:
        logger.warning(f"Target column not found")
        return None

    # Aggregate flows
    if value_col in mapping_df.columns:
        flows = mapping_df.groupby([source_col, target_col])[value_col].sum().reset_index()
    else:
        flows = mapping_df.groupby([source_col, target_col]).size().reset_index(name="n_cells")
        value_col = "n_cells"

    # Get top labels
    source_totals = flows.groupby(source_col)[value_col].sum().nlargest(top_n)
    target_totals = flows.groupby(target_col)[value_col].sum().nlargest(top_n)

    flows = flows[
        flows[source_col].isin(source_totals.index) &
        flows[target_col].isin(target_totals.index)
    ]

    # Create node list (sources + targets with unique suffix)
    sources = flows[source_col].unique().tolist()
    targets = flows[target_col].unique().tolist()

    # Add suffix to avoid collision between same labels
    all_nodes = [f"{s} (before)" for s in sources] + [f"{t} (after)" for t in targets]
    node_indices = {n: i for i, n in enumerate(all_nodes)}

    # Get colors
    all_labels = list(set(sources) | set(targets))
    colors = get_color_palette(all_labels, palette_type="cell_type")

    # Build Sankey data
    source_indices = [node_indices[f"{s} (before)"] for s in flows[source_col]]
    target_indices = [node_indices[f"{t} (after)"] for t in flows[target_col]]
    values = flows[value_col].tolist()

    # Node colors
    node_colors = (
        [colors.get(s, "#bdc3c7") for s in sources] +
        [colors.get(t, "#bdc3c7") for t in targets]
    )

    # Link colors (lighter version of source color)
    def hex_to_rgba(hex_color, alpha=0.4):
        hex_color = hex_color.lstrip("#")
        r = int(hex_color[0:2], 16)
        g = int(hex_color[2:4], 16)
        b = int(hex_color[4:6], 16)
        return f"rgba({r},{g},{b},{alpha})"

    link_colors = [hex_to_rgba(colors.get(s, "#bdc3c7"), 0.4) for s in flows[source_col]]

    # Create Sankey diagram
    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="white", width=1),
            label=all_nodes,
            color=node_colors,
            hovertemplate="%{label}<br>%{value:,} cells<extra></extra>",
        ),
        link=dict(
            source=source_indices,
            target=target_indices,
            value=values,
            color=link_colors,
            hovertemplate=(
                "%{source.label} -> %{target.label}<br>"
                "%{value:,} cells<extra></extra>"
            ),
        ),
    )])

    fig.update_layout(
        title=dict(
            text="Label Flow: Before -> After Consolidation",
            font=dict(size=16),
        ),
        font=dict(size=10),
        height=800,
        width=1200,
        # Add column name annotations
        annotations=[
            dict(
                x=0.01,
                y=1.08,
                xref="paper",
                yref="paper",
                text=f"<b>Source Column:</b> {source_col}",
                showarrow=False,
                font=dict(size=12, color="#2c3e50"),
                align="left",
            ),
            dict(
                x=0.99,
                y=1.08,
                xref="paper",
                yref="paper",
                text=f"<b>Target Column:</b> {target_col}",
                showarrow=False,
                font=dict(size=12, color="#2c3e50"),
                align="right",
            ),
        ],
    )

    fig.write_html(str(output_path))
    logger.info(f"Saved interactive Sankey to {output_path}")

    return output_path
