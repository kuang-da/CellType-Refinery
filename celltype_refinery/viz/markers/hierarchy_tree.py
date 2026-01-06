"""Tree diagram visualization for marker hierarchies.

Creates matplotlib-based hierarchical tree diagrams showing cell type
relationships and their associated markers.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Tuple

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

from .style import (
    MARKER_COLORS,
    MarkerNode,
    PathLike,
    format_marker_list,
    get_category_color,
    parse_marker_hierarchy,
    save_figure,
    set_plot_style,
)


def _compute_tree_layout(
    nodes: List[MarkerNode],
    node_lookup: Dict[str, MarkerNode],
    horizontal_spacing: float = 3.0,
    vertical_spacing: float = 2.0,
) -> Dict[str, Tuple[float, float]]:
    """
    Compute x,y positions for tree layout.

    Uses a simple recursive algorithm:
    - Y position determined by level (depth)
    - X position spreads children evenly under parent

    Returns
    -------
    Dict[str, Tuple[float, float]]
        Mapping of node path to (x, y) position.
    """
    positions: Dict[str, Tuple[float, float]] = {}

    # Find root nodes (level 0)
    roots = [n for n in nodes if n.level == 0]
    if not roots:
        return positions

    # Count leaves in each subtree for spacing
    def count_leaves(node: MarkerNode) -> int:
        if not node.children:
            return 1
        return sum(count_leaves(node_lookup[c]) for c in node.children if c in node_lookup)

    leaf_counts = {n.path: count_leaves(n) for n in nodes}

    # Assign positions
    def assign_positions(node: MarkerNode, x_start: float, x_end: float) -> None:
        y = -node.level * vertical_spacing
        x = (x_start + x_end) / 2
        positions[node.path] = (x, y)

        if node.children:
            total_leaves = sum(leaf_counts.get(c, 1) for c in node.children)
            current_x = x_start
            for child_path in node.children:
                if child_path not in node_lookup:
                    continue
                child = node_lookup[child_path]
                child_width = (leaf_counts.get(child_path, 1) / total_leaves) * (x_end - x_start)
                assign_positions(child, current_x, current_x + child_width)
                current_x += child_width

    # Lay out root nodes
    total_root_leaves = sum(leaf_counts.get(r.path, 1) for r in roots)
    total_width = total_root_leaves * horizontal_spacing
    current_x = 0.0

    for root in roots:
        root_width = (leaf_counts.get(root.path, 1) / total_root_leaves) * total_width
        assign_positions(root, current_x, current_x + root_width)
        current_x += root_width

    return positions


def plot_marker_hierarchy_tree(
    marker_map: Dict[str, Any],
    output_path: PathLike,
    *,
    title: str = "Cell Type Marker Hierarchy",
    figsize: Tuple[float, float] = (16, 12),
    show_markers: bool = True,
    show_anti_markers: bool = True,
    max_markers_display: int = 3,
    marker_color: str = "#2ecc71",
    anti_marker_color: str = "#e74c3c",
    node_color: str = "#3498db",
    font_size: int = 8,
) -> Path:
    """
    Create a hierarchical tree diagram of the marker map.

    Parameters
    ----------
    marker_map : Dict[str, Any]
        Marker hierarchy loaded from JSON file.
    output_path : PathLike
        Output file path (PNG/PDF based on extension).
    title : str
        Plot title.
    figsize : Tuple[float, float]
        Figure dimensions in inches.
    show_markers : bool
        Whether to display positive markers.
    show_anti_markers : bool
        Whether to display anti-markers.
    max_markers_display : int
        Maximum markers to show per node.
    marker_color : str
        Color for positive marker labels.
    anti_marker_color : str
        Color for anti-marker labels.
    node_color : str
        Default color for cell type nodes.
    font_size : int
        Base font size for labels.

    Returns
    -------
    Path
        Path to saved figure.
    """
    nodes, node_lookup = parse_marker_hierarchy(marker_map)
    if not nodes:
        # Create empty figure for empty marker map
        set_plot_style()
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "Empty marker map", ha="center", va="center", fontsize=14)
        ax.set_title(title)
        ax.axis("off")
        return save_figure(fig, output_path)

    positions = _compute_tree_layout(nodes, node_lookup)

    set_plot_style()
    fig, ax = plt.subplots(figsize=figsize)

    # Draw edges first (behind nodes)
    for node in nodes:
        if node.parent and node.parent in positions:
            x1, y1 = positions[node.parent]
            x2, y2 = positions[node.path]
            ax.plot([x1, x2], [y1, y2], color=MARKER_COLORS["edge_hierarchy"],
                   linewidth=1.5, zorder=1, alpha=0.7)

    # Draw nodes
    for node in nodes:
        if node.path not in positions:
            continue
        x, y = positions[node.path]

        # Get category color for this node
        color = get_category_color(node.path)

        # Draw node box
        box_width = 0.8
        box_height = 0.4
        box = FancyBboxPatch(
            (x - box_width/2, y - box_height/2),
            box_width, box_height,
            boxstyle="round,pad=0.02,rounding_size=0.1",
            facecolor=color,
            edgecolor="white",
            linewidth=2,
            alpha=0.9,
            zorder=2,
        )
        ax.add_patch(box)

        # Draw node name
        ax.text(x, y, node.name, ha="center", va="center", fontsize=font_size,
               fontweight="bold", color="white", zorder=3)

        # Draw markers below node
        marker_texts = []
        if show_markers and node.markers:
            marker_str = format_marker_list(node.markers, max_markers_display)
            marker_texts.append((marker_str, marker_color))
        if show_anti_markers and node.anti_markers:
            anti_str = format_marker_list(node.anti_markers, max_markers_display, prefix="!")
            marker_texts.append((anti_str, anti_marker_color))

        if marker_texts:
            y_offset = -box_height/2 - 0.15
            for text, color in marker_texts:
                ax.text(x, y + y_offset, text, ha="center", va="top",
                       fontsize=font_size - 1, color=color, zorder=3)
                y_offset -= 0.2

    # Set axis limits with padding
    if positions:
        xs = [p[0] for p in positions.values()]
        ys = [p[1] for p in positions.values()]
        x_margin = max(2, (max(xs) - min(xs)) * 0.1)
        y_margin = max(1, (max(ys) - min(ys)) * 0.15)
        ax.set_xlim(min(xs) - x_margin, max(xs) + x_margin)
        ax.set_ylim(min(ys) - y_margin - 1, max(ys) + y_margin)

    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.axis("off")

    # Add legend
    legend_elements = [
        mpatches.Patch(facecolor=marker_color, label="Positive markers"),
        mpatches.Patch(facecolor=anti_marker_color, label="Anti-markers"),
    ]
    ax.legend(handles=legend_elements, loc="upper right", frameon=True, fontsize=font_size)

    fig.tight_layout()
    return save_figure(fig, output_path)


__all__ = ["plot_marker_hierarchy_tree"]
