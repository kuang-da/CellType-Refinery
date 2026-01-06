"""Sunburst visualization for marker hierarchies.

Creates interactive HTML sunburst diagrams using Plotly for exploring
marker hierarchy structure.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List

from .style import (
    MarkerNode,
    PathLike,
    ensure_parent,
    get_category_color,
    parse_marker_hierarchy,
)


def _to_sunburst_data(nodes: List[MarkerNode]) -> Dict[str, List]:
    """
    Convert nodes to plotly sunburst format.

    Returns
    -------
    Dict with keys: ids, labels, parents, values, colors, hover
    """
    ids: List[str] = []
    labels: List[str] = []
    parents: List[str] = []
    values: List[int] = []
    colors: List[str] = []
    hover: List[str] = []

    for node in nodes:
        ids.append(node.path)
        labels.append(node.name)
        parents.append(node.parent or "")
        # Value based on marker count (minimum 1 for visibility)
        values.append(max(1, len(node.markers) + len(node.anti_markers)))
        colors.append(get_category_color(node.path))

        # Build hover text
        hover_parts = [f"<b>{node.name}</b>"]
        if node.markers:
            hover_parts.append(f"Markers: {', '.join(node.markers[:5])}")
            if len(node.markers) > 5:
                hover_parts.append(f"  (+{len(node.markers) - 5} more)")
        if node.anti_markers:
            hover_parts.append(f"Anti: {', '.join(node.anti_markers[:5])}")
            if len(node.anti_markers) > 5:
                hover_parts.append(f"  (+{len(node.anti_markers) - 5} more)")
        if node.notes:
            hover_parts.append(f"Notes: {node.notes[:100]}...")
        hover.append("<br>".join(hover_parts))

    return {
        "ids": ids,
        "labels": labels,
        "parents": parents,
        "values": values,
        "colors": colors,
        "hover": hover,
    }


def plot_marker_hierarchy_sunburst(
    marker_map: Dict[str, Any],
    output_path: PathLike,
    *,
    title: str = "Cell Type Marker Hierarchy",
    color_by: str = "category",
) -> Path:
    """
    Create a sunburst/radial treemap of the marker hierarchy using plotly.

    Parameters
    ----------
    marker_map : Dict[str, Any]
        Marker hierarchy loaded from JSON file.
    output_path : PathLike
        Output file path (.html for interactive).
    title : str
        Plot title.
    color_by : str
        Coloring strategy: "category" (top-level), "level" (depth).

    Returns
    -------
    Path
        Path to saved HTML file.
    """
    import plotly.graph_objects as go

    nodes, _ = parse_marker_hierarchy(marker_map)
    output_path = Path(output_path)

    if not nodes:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text("<html><body><h1>Empty marker map</h1></body></html>")
        return output_path

    sunburst_data = _to_sunburst_data(nodes)

    fig = go.Figure(go.Sunburst(
        ids=sunburst_data["ids"],
        labels=sunburst_data["labels"],
        parents=sunburst_data["parents"],
        values=sunburst_data["values"],
        marker=dict(colors=sunburst_data["colors"]),
        hovertext=sunburst_data["hover"],
        hoverinfo="text",
    ))

    fig.update_layout(
        title=dict(text=title, font=dict(size=18)),
        margin=dict(t=50, l=10, r=10, b=10),
    )

    output_path = ensure_parent(output_path)
    fig.write_html(str(output_path))
    return output_path


__all__ = ["plot_marker_hierarchy_sunburst"]
