"""Cell-type interaction network visualization."""

from __future__ import annotations

from pathlib import Path
from typing import Tuple

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd

try:
    import networkx as nx

    HAS_NETWORKX = True
except ImportError:
    HAS_NETWORKX = False


def plot_interaction_network(
    interaction_matrix: pd.DataFrame,
    output_path: Path,
    threshold: float = 0.5,
    figsize: Tuple[int, int] = (12, 12),
    dpi: int = 200,
    title: str = "Cell-Type Interaction Network",
) -> None:
    """Plot cell-type interaction network.

    Parameters
    ----------
    interaction_matrix : pd.DataFrame
        Square matrix of log2 fold enrichment scores
    output_path : Path
        Path to save figure
    threshold : float
        Minimum |score| for drawing edges
    figsize : Tuple[int, int]
        Figure size
    dpi : int
        Figure resolution
    title : str
        Plot title
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx is required for network visualization")

    if interaction_matrix.empty:
        return

    # Create directed graph
    G = nx.DiGraph()

    cell_types = interaction_matrix.index.tolist()
    G.add_nodes_from(cell_types)

    # Add edges for significant interactions
    for i, ct1 in enumerate(cell_types):
        for j, ct2 in enumerate(cell_types):
            if i != j:
                score = interaction_matrix.iloc[i, j]
                if abs(score) > threshold:
                    G.add_edge(ct1, ct2, weight=score)

    # Handle empty graph
    if len(G.edges()) == 0:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(
            0.5,
            0.5,
            f"No interactions above threshold (|score| > {threshold})",
            ha="center",
            va="center",
            fontsize=12,
        )
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        return

    fig, ax = plt.subplots(figsize=figsize)

    # Layout
    pos = nx.spring_layout(G, k=2, iterations=50, seed=42)

    # Draw nodes
    nx.draw_networkx_nodes(
        G,
        pos,
        ax=ax,
        node_size=500,
        node_color="lightblue",
        edgecolors="navy",
        linewidths=1,
    )
    nx.draw_networkx_labels(G, pos, ax=ax, font_size=8)

    # Draw edges with color based on weight
    edges = G.edges(data=True)
    edge_colors = ["red" if e[2]["weight"] > 0 else "blue" for e in edges]
    edge_widths = [min(abs(e[2]["weight"]), 3) for e in edges]

    nx.draw_networkx_edges(
        G,
        pos,
        ax=ax,
        edge_color=edge_colors,
        width=edge_widths,
        alpha=0.7,
        arrows=True,
        arrowsize=10,
        connectionstyle="arc3,rad=0.1",
    )

    # Legend
    red_patch = mpatches.Patch(color="red", label="Enriched interaction")
    blue_patch = mpatches.Patch(color="blue", label="Depleted interaction")
    ax.legend(handles=[red_patch, blue_patch], loc="upper left")

    ax.set_title(title, fontsize=14)
    ax.axis("off")

    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
