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
    StateNode,
    ensure_parent,
    get_category_color,
    parse_cell_states,
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


def _build_states_html(states: List[StateNode]) -> str:
    """
    Build HTML table for cell states sidebar.

    Parameters
    ----------
    states : List[StateNode]
        List of cell state definitions.

    Returns
    -------
    str
        HTML string for the states sidebar panel.
    """
    if not states:
        return ""

    rows = []
    for s in states:
        markers_str = ", ".join(s.markers) if s.markers else "-"
        applies = s.applies_to or "Any cell type"
        rows.append(f"""
            <tr>
                <td style="font-weight:bold; color:#16a085; padding:10px 12px; border-bottom:1px solid #ddd;">{s.name}</td>
                <td style="padding:10px 12px; border-bottom:1px solid #ddd; font-family:monospace;">{markers_str}</td>
                <td style="padding:10px 12px; border-bottom:1px solid #ddd; color:#555;">{applies}</td>
            </tr>
        """)

    return f"""
    <div style="background:#f9f9f9; border-radius:8px; padding:20px;">
        <h3 style="margin-top:0; color:#16a085; border-bottom:2px solid #16a085; padding-bottom:8px;">
            Cell States (Overlays)
        </h3>
        <p style="font-size:13px; color:#666; margin-bottom:16px;">
            States are scored separately and applied as overlay annotations after primary type assignment.
        </p>
        <table style="width:100%; border-collapse:collapse; font-size:14px;">
            <thead>
                <tr style="background:#16a085; color:white;">
                    <th style="padding:10px 12px; text-align:left; width:25%;">State</th>
                    <th style="padding:10px 12px; text-align:left; width:20%;">Markers</th>
                    <th style="padding:10px 12px; text-align:left; width:55%;">Applies To</th>
                </tr>
            </thead>
            <tbody>
                {''.join(rows)}
            </tbody>
        </table>
    </div>
    """


def plot_marker_hierarchy_sunburst(
    marker_map: Dict[str, Any],
    output_path: PathLike,
    *,
    title: str = "Cell Type Marker Hierarchy",
    color_by: str = "category",
) -> Path:
    """
    Create a sunburst/radial treemap of the marker hierarchy using plotly.

    Includes a sidebar panel showing cell states from _cell_state_markers.

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
    states = parse_cell_states(marker_map)
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
        title=dict(text=title, font=dict(size=24)),
        margin=dict(t=80, l=40, r=40, b=40),
        width=1100,
        height=900,
    )

    output_path = ensure_parent(output_path)

    # Build states sidebar HTML
    states_html = _build_states_html(states)

    # Get sunburst as embeddable HTML
    sunburst_html = fig.to_html(full_html=False, include_plotlyjs='cdn')

    # Combine into top-bottom layout (large sunburst on top, states table below)
    final_html = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>{title}</title>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            margin: 0;
            padding: 20px;
            background: #fff;
        }}
        .sunburst-panel {{
            display: flex;
            justify-content: center;
            margin-bottom: 40px;
        }}
        .states-panel {{
            max-width: 1000px;
            margin: 0 auto;
            padding: 0 20px;
        }}
    </style>
</head>
<body>
    <div class="sunburst-panel">
        {sunburst_html}
    </div>
    <div class="states-panel">
        {states_html}
    </div>
</body>
</html>"""

    output_path.write_text(final_html)
    return output_path


__all__ = ["plot_marker_hierarchy_sunburst"]
