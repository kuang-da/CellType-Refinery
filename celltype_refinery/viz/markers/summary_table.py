"""Summary table visualization for marker hierarchies.

Creates a tabular view showing all cell types and their associated markers.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Tuple

import matplotlib.pyplot as plt

from .style import (
    PathLike,
    parse_marker_hierarchy,
    save_figure,
    set_plot_style,
)


def plot_marker_summary_table(
    marker_map: Dict[str, Any],
    output_path: PathLike,
    *,
    title: str = "Marker Summary by Cell Type",
    figsize: Tuple[float, float] = (14, 10),
    max_markers_display: int = 5,
) -> Path:
    """
    Create a summary table showing all cell types and their markers.

    Parameters
    ----------
    marker_map : Dict[str, Any]
        Marker hierarchy loaded from JSON file.
    output_path : PathLike
        Output file path (PNG/PDF).
    title : str
        Table title.
    figsize : Tuple[float, float]
        Figure size.
    max_markers_display : int
        Maximum markers to show per cell.

    Returns
    -------
    Path
        Path to saved figure.
    """
    nodes, _ = parse_marker_hierarchy(marker_map)

    if not nodes:
        set_plot_style()
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "Empty marker map", ha="center", va="center", fontsize=14)
        ax.set_title(title)
        ax.axis("off")
        return save_figure(fig, output_path)

    # Build table data
    rows = []
    for node in nodes:
        indent = "  " * node.level
        markers_str = ", ".join(node.markers[:max_markers_display])
        if len(node.markers) > max_markers_display:
            markers_str += f" (+{len(node.markers) - max_markers_display})"
        anti_str = ", ".join(node.anti_markers[:max_markers_display])
        if len(node.anti_markers) > max_markers_display:
            anti_str += f" (+{len(node.anti_markers) - max_markers_display})"
        rows.append([
            f"{indent}{node.name}",
            str(len(node.markers)),
            markers_str or "-",
            str(len(node.anti_markers)),
            anti_str or "-",
        ])

    columns = ["Cell Type", "#M", "Markers", "#A", "Anti-Markers"]

    set_plot_style()
    fig, ax = plt.subplots(figsize=figsize)
    ax.axis("off")

    # Create table
    table = ax.table(
        cellText=rows,
        colLabels=columns,
        cellLoc="left",
        loc="center",
        colColours=["#3498db"] * len(columns),
    )

    # Style table
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1.2, 1.5)

    # Color header
    for key, cell in table.get_celld().items():
        if key[0] == 0:  # Header row
            cell.set_text_props(color="white", fontweight="bold")
        else:
            # Alternate row colors
            if key[0] % 2 == 0:
                cell.set_facecolor("#f8f9fa")

    ax.set_title(title, fontsize=14, fontweight="bold", pad=20)
    fig.tight_layout()
    return save_figure(fig, output_path)


__all__ = ["plot_marker_summary_table"]
