"""Network graph visualization for marker hierarchies.

Creates force-directed network diagrams showing cell types and their
marker relationships.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np

from .style import (
    MARKER_COLORS,
    PathLike,
    get_category_color,
    parse_marker_hierarchy,
    save_figure,
    set_plot_style,
)


def _simple_spring_layout(
    edges: List[Tuple[str, str]],
    node_ids: List[str],
    n_iterations: int = 50,
    k: float = 1.0,
    seed: int = 42,
) -> Dict[str, Tuple[float, float]]:
    """
    Simple force-directed layout using Fruchterman-Reingold algorithm.

    Avoids networkx dependency for basic functionality.
    """
    rng = np.random.default_rng(seed)
    n = len(node_ids)
    if n == 0:
        return {}

    # Initialize random positions
    pos = rng.uniform(-1, 1, size=(n, 2))
    node_to_idx = {node: i for i, node in enumerate(node_ids)}

    # Build adjacency for attraction
    adj = {i: [] for i in range(n)}
    for u, v in edges:
        if u in node_to_idx and v in node_to_idx:
            i, j = node_to_idx[u], node_to_idx[v]
            adj[i].append(j)
            adj[j].append(i)

    # Optimal distance
    area = 4.0
    k_opt = k * np.sqrt(area / n)

    for iteration in range(n_iterations):
        # Repulsive forces
        disp = np.zeros((n, 2))
        for i in range(n):
            for j in range(i + 1, n):
                delta = pos[i] - pos[j]
                dist = max(np.linalg.norm(delta), 0.01)
                force = k_opt ** 2 / dist
                direction = delta / dist
                disp[i] += direction * force
                disp[j] -= direction * force

        # Attractive forces
        for i, neighbors in adj.items():
            for j in neighbors:
                delta = pos[i] - pos[j]
                dist = max(np.linalg.norm(delta), 0.01)
                force = dist ** 2 / k_opt
                direction = delta / dist
                disp[i] -= direction * force

        # Apply displacement with cooling
        temp = 1.0 - (iteration / n_iterations)
        for i in range(n):
            disp_norm = max(np.linalg.norm(disp[i]), 0.01)
            pos[i] += (disp[i] / disp_norm) * min(disp_norm, temp)

    return {node_ids[i]: (pos[i, 0], pos[i, 1]) for i in range(n)}


def plot_marker_network(
    marker_map: Dict[str, Any],
    output_path: PathLike,
    *,
    title: str = "Cell Type Marker Network",
    figsize: Tuple[float, float] = (14, 12),
    layout: str = "spring",
    show_marker_nodes: bool = True,
    marker_color: str = "#2ecc71",
    anti_marker_color: str = "#e74c3c",
    cell_type_color: str = "#3498db",
    edge_alpha: float = 0.3,
    node_size_scale: float = 1.0,
) -> Path:
    """
    Create a node-link network diagram showing cell types and markers.

    Parameters
    ----------
    marker_map : Dict[str, Any]
        Marker hierarchy.
    output_path : PathLike
        Output path.
    title : str
        Plot title.
    figsize : Tuple[float, float]
        Figure size.
    layout : str
        Layout algorithm: "spring" (force-directed).
    show_marker_nodes : bool
        If True, show markers as separate nodes connected to cell types.
        If False, show only cell type hierarchy.
    marker_color : str
        Color for positive marker nodes.
    anti_marker_color : str
        Color for anti-marker nodes.
    cell_type_color : str
        Color for cell type nodes.
    edge_alpha : float
        Edge transparency.
    node_size_scale : float
        Scaling factor for node sizes.

    Returns
    -------
    Path
        Path to saved figure.
    """
    nodes, node_lookup = parse_marker_hierarchy(marker_map)

    if not nodes:
        set_plot_style()
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "Empty marker map", ha="center", va="center", fontsize=14)
        ax.set_title(title)
        ax.axis("off")
        return save_figure(fig, output_path)

    # Build graph
    all_node_ids: List[str] = []
    node_types: Dict[str, str] = {}  # path -> "cell_type" | "marker" | "anti_marker"
    edges: List[Tuple[str, str]] = []

    # Add cell type nodes and hierarchy edges
    for node in nodes:
        all_node_ids.append(node.path)
        node_types[node.path] = "cell_type"
        if node.parent:
            edges.append((node.parent, node.path))

    # Add marker nodes if requested
    if show_marker_nodes:
        seen_markers: Dict[str, str] = {}  # marker -> first cell type that uses it
        for node in nodes:
            for m in node.markers:
                marker_id = f"M:{m}"
                if marker_id not in seen_markers:
                    all_node_ids.append(marker_id)
                    node_types[marker_id] = "marker"
                    seen_markers[marker_id] = node.path
                edges.append((node.path, marker_id))

            for m in node.anti_markers:
                marker_id = f"A:{m}"
                if marker_id not in seen_markers:
                    all_node_ids.append(marker_id)
                    node_types[marker_id] = "anti_marker"
                    seen_markers[marker_id] = node.path
                edges.append((node.path, marker_id))

    # Compute layout
    positions = _simple_spring_layout(edges, all_node_ids, n_iterations=100)

    set_plot_style()
    fig, ax = plt.subplots(figsize=figsize)

    # Draw edges
    for u, v in edges:
        if u in positions and v in positions:
            x1, y1 = positions[u]
            x2, y2 = positions[v]
            # Different edge styles for different connections
            if node_types.get(v) == "anti_marker":
                linestyle = ":"
                color = anti_marker_color
            elif node_types.get(v) == "marker":
                linestyle = "--"
                color = marker_color
            else:
                linestyle = "-"
                color = MARKER_COLORS["edge_hierarchy"]
            ax.plot([x1, x2], [y1, y2], color=color, alpha=edge_alpha,
                   linewidth=1, linestyle=linestyle, zorder=1)

    # Draw nodes
    for node_id in all_node_ids:
        if node_id not in positions:
            continue
        x, y = positions[node_id]
        ntype = node_types.get(node_id, "cell_type")

        if ntype == "cell_type":
            color = get_category_color(node_id)
            size = 300 * node_size_scale
            marker_shape = "o"
            label = node_id.split("/")[-1]
            fontsize = 7
        elif ntype == "marker":
            color = marker_color
            size = 100 * node_size_scale
            marker_shape = "s"
            label = node_id[2:]  # Remove "M:" prefix
            fontsize = 6
        else:  # anti_marker
            color = anti_marker_color
            size = 100 * node_size_scale
            marker_shape = "^"
            label = node_id[2:]  # Remove "A:" prefix
            fontsize = 6

        ax.scatter([x], [y], c=[color], s=size, marker=marker_shape,
                  edgecolors="white", linewidths=1, zorder=2, alpha=0.9)
        ax.annotate(label, (x, y), xytext=(0, -12), textcoords="offset points",
                   ha="center", va="top", fontsize=fontsize, zorder=3)

    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.axis("off")

    # Add legend
    legend_elements = [
        plt.scatter([], [], c=cell_type_color, s=100, marker="o", label="Cell Type"),
        plt.scatter([], [], c=marker_color, s=60, marker="s", label="Positive Marker"),
        plt.scatter([], [], c=anti_marker_color, s=60, marker="^", label="Anti-Marker"),
    ]
    ax.legend(handles=legend_elements, loc="upper right", frameon=True, fontsize=8)

    fig.tight_layout()
    return save_figure(fig, output_path)


__all__ = ["plot_marker_network"]
