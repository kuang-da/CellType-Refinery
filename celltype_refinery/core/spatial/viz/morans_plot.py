"""Moran's I comparison visualization."""

from __future__ import annotations

from pathlib import Path
from typing import Tuple

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def plot_morans_comparison(
    morans_global: pd.DataFrame,
    output_path: Path,
    figsize: Tuple[int, int] = (12, 6),
    dpi: int = 200,
    title: str = "Spatial Autocorrelation by Cell Type",
    top_n: int = 20,
) -> None:
    """Plot Moran's I comparison by cell type.

    Parameters
    ----------
    morans_global : pd.DataFrame
        Aggregated Moran's I by cell type
        Columns: cell_type, morans_i_mean, morans_i_std, n_samples
    output_path : Path
        Path to save figure
    figsize : Tuple[int, int]
        Figure size
    dpi : int
        Figure resolution
    title : str
        Plot title
    top_n : int
        Number of cell types to show (by absolute Moran's I)
    """
    if morans_global.empty:
        return

    # Sort by mean Moran's I (descending = most clustered)
    df = morans_global.sort_values("morans_i_mean", ascending=False).head(top_n)

    fig, ax = plt.subplots(figsize=figsize)

    # Bar plot with error bars
    colors = ["steelblue" if v >= 0 else "coral" for v in df["morans_i_mean"]]

    bars = ax.bar(
        range(len(df)),
        df["morans_i_mean"],
        yerr=df["morans_i_std"],
        color=colors,
        edgecolor="black",
        linewidth=0.5,
        capsize=3,
        error_kw={"elinewidth": 1, "capthick": 1},
    )

    # Reference line at 0
    ax.axhline(y=0, color="red", linestyle="--", alpha=0.7, linewidth=1)

    # Thresholds
    ax.axhline(y=0.3, color="gray", linestyle=":", alpha=0.5)
    ax.axhline(y=-0.2, color="gray", linestyle=":", alpha=0.5)

    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(df["cell_type"], rotation=45, ha="right", fontsize=9)

    ax.set_xlabel("Cell Type")
    ax.set_ylabel("Mean Moran's I")
    ax.set_title(title, fontsize=14)

    # Add annotation
    ax.text(
        0.98,
        0.95,
        "I > 0.3: clustering\nI < -0.2: dispersion",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
    )

    plt.tight_layout()
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
