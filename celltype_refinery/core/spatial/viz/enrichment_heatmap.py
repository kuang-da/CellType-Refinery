"""Neighborhood enrichment heatmap visualization."""

from __future__ import annotations

from pathlib import Path
from typing import Tuple

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def plot_enrichment_heatmap(
    z_scores: pd.DataFrame,
    output_path: Path,
    figsize: Tuple[int, int] = (12, 10),
    vmin: float = -5.0,
    vmax: float = 5.0,
    cmap: str = "RdBu_r",
    dpi: int = 200,
    title: str = "Neighborhood Enrichment",
) -> None:
    """Plot clustered heatmap of neighborhood enrichment z-scores.

    Parameters
    ----------
    z_scores : pd.DataFrame
        Square matrix of z-scores (cell_type x cell_type)
    output_path : Path
        Path to save figure
    figsize : Tuple[int, int]
        Figure size
    vmin, vmax : float
        Color scale limits
    cmap : str
        Colormap name
    dpi : int
        Figure resolution
    title : str
        Plot title
    """
    if z_scores.empty:
        return

    # Set publication style
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
        }
    )

    # Create clustermap
    g = sns.clustermap(
        z_scores,
        cmap=cmap,
        center=0,
        vmin=vmin,
        vmax=vmax,
        figsize=figsize,
        xticklabels=True,
        yticklabels=True,
        cbar_kws={"label": "Z-score"},
        dendrogram_ratio=0.15,
        cbar_pos=(0.02, 0.8, 0.03, 0.15),
    )

    g.fig.suptitle(title, y=1.02, fontsize=14)

    # Rotate labels for readability
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha="right")
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)

    g.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close()
