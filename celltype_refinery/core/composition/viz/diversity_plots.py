"""Diversity visualization functions.

Provides:
- Diversity metrics by region (boxplots)
- Diversity metrics by sample (scatter/bar)
- Evenness distribution histogram
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Tuple
import logging

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def _sort_regions(regions, region_order: Optional[List[str]] = None) -> List[str]:
    """Sort regions by specified order or alphabetically."""
    if not region_order:
        return sorted(regions)

    order_map = {r.lower(): i for i, r in enumerate(region_order)}

    def get_order(r):
        r_lower = r.lower() if isinstance(r, str) else str(r).lower()
        return order_map.get(r_lower, len(region_order))

    return sorted(regions, key=get_order)


def _set_style():
    """Set publication-quality plot style."""
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "axes.grid": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "font.size": 10,
        "axes.labelsize": 11,
        "axes.titlesize": 12,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
    })


def _save_figure(fig, output_path: Path, dpi: int = 200) -> Path:
    """Save figure and close."""
    import matplotlib.pyplot as plt

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(output_path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    return output_path


def plot_diversity_by_region(
    diversity_by_sample: pd.DataFrame,
    output_path: Path,
    region_col: str = "region",
    region_order: Optional[List[str]] = None,
    dpi: int = 200,
    figsize: Tuple[int, int] = (12, 5),
) -> Optional[Path]:
    """Plot diversity metrics by region as boxplots.

    Creates a 2-panel figure showing Shannon entropy and Simpson index
    distributions across regions.

    Parameters
    ----------
    diversity_by_sample : pd.DataFrame
        Output from compute_diversity_by_group() with region column
    output_path : Path
        Output file path
    region_col : str
        Column with region labels
    region_order : List[str], optional
        Canonical region order
    dpi : int
        Figure resolution
    figsize : Tuple[int, int]
        Figure size

    Returns
    -------
    Optional[Path]
        Path to saved figure, or None if region column not found
    """
    import matplotlib.pyplot as plt

    _set_style()

    if region_col not in diversity_by_sample.columns:
        logger.warning(f"Region column '{region_col}' not found, skipping plot")
        return None

    fig, axes = plt.subplots(1, 2, figsize=figsize)

    # Sort regions
    sorted_regions = _sort_regions(diversity_by_sample[region_col].unique(), region_order)

    # Try seaborn for nicer boxplots
    try:
        import seaborn as sns

        # Shannon entropy
        sns.boxplot(
            data=diversity_by_sample,
            x=region_col,
            y="shannon_entropy",
            hue=region_col,
            order=sorted_regions,
            hue_order=sorted_regions,
            ax=axes[0],
            palette="Set2",
            legend=False,
        )
        sns.stripplot(
            data=diversity_by_sample,
            x=region_col,
            y="shannon_entropy",
            order=sorted_regions,
            ax=axes[0],
            color="black",
            alpha=0.5,
            size=4,
        )
        axes[0].set_title("Shannon Entropy by Region")
        axes[0].set_xlabel("Region")
        axes[0].set_ylabel("Shannon Entropy (nats)")
        plt.setp(axes[0].get_xticklabels(), rotation=45, ha="right")

        # Simpson index
        sns.boxplot(
            data=diversity_by_sample,
            x=region_col,
            y="simpson_index",
            hue=region_col,
            order=sorted_regions,
            hue_order=sorted_regions,
            ax=axes[1],
            palette="Set2",
            legend=False,
        )
        sns.stripplot(
            data=diversity_by_sample,
            x=region_col,
            y="simpson_index",
            order=sorted_regions,
            ax=axes[1],
            color="black",
            alpha=0.5,
            size=4,
        )
        axes[1].set_title("Simpson Diversity by Region")
        axes[1].set_xlabel("Region")
        axes[1].set_ylabel("Simpson Index")
        plt.setp(axes[1].get_xticklabels(), rotation=45, ha="right")

    except ImportError:
        # Fallback without seaborn
        data_shannon = [
            diversity_by_sample[diversity_by_sample[region_col] == r]["shannon_entropy"]
            for r in sorted_regions
        ]
        axes[0].boxplot(data_shannon, labels=sorted_regions)
        axes[0].set_title("Shannon Entropy by Region")
        axes[0].set_xlabel("Region")
        axes[0].set_ylabel("Shannon Entropy (nats)")
        axes[0].set_xticklabels(sorted_regions, rotation=45, ha="right")

        data_simpson = [
            diversity_by_sample[diversity_by_sample[region_col] == r]["simpson_index"]
            for r in sorted_regions
        ]
        axes[1].boxplot(data_simpson, labels=sorted_regions)
        axes[1].set_title("Simpson Diversity by Region")
        axes[1].set_xlabel("Region")
        axes[1].set_ylabel("Simpson Index")
        axes[1].set_xticklabels(sorted_regions, rotation=45, ha="right")

    plt.tight_layout()
    return _save_figure(fig, output_path, dpi)


def plot_diversity_by_sample(
    diversity_by_sample: pd.DataFrame,
    output_path: Path,
    sample_col: str = "sample_id",
    region_col: str = "region",
    dpi: int = 200,
    figsize: Tuple[int, int] = (14, 6),
) -> Path:
    """Plot diversity metrics per sample.

    Creates a 2-panel figure showing Shannon entropy and Simpson index
    per sample, colored by region if available.

    Parameters
    ----------
    diversity_by_sample : pd.DataFrame
        Output from compute_diversity_by_group()
    output_path : Path
        Output file path
    sample_col : str
        Column with sample IDs
    region_col : str
        Column with region labels (for coloring)
    dpi : int
        Figure resolution
    figsize : Tuple[int, int]
        Figure size

    Returns
    -------
    Path
        Path to saved figure
    """
    import matplotlib.pyplot as plt

    _set_style()

    # Sort by Shannon entropy
    df = diversity_by_sample.sort_values("shannon_entropy", ascending=False).reset_index(drop=True)

    fig, axes = plt.subplots(1, 2, figsize=figsize)

    # Determine colors
    if region_col in df.columns:
        regions = df[region_col].unique()
        color_values = [tuple(c) for c in plt.cm.Set2(np.linspace(0, 1, len(regions)))]
        region_colors = dict(zip(regions, color_values))
        colors = [region_colors[r] for r in df[region_col]]
    else:
        colors = "steelblue"

    # Shannon entropy
    axes[0].bar(range(len(df)), df["shannon_entropy"], color=colors)
    axes[0].set_xlabel("Sample (sorted by entropy)")
    axes[0].set_ylabel("Shannon Entropy (nats)")
    axes[0].set_title("Shannon Entropy per Sample")
    axes[0].axhline(
        df["shannon_entropy"].mean(),
        color="red",
        linestyle="--",
        label=f"Mean: {df['shannon_entropy'].mean():.2f}",
    )
    axes[0].legend()

    # Simpson index
    axes[1].bar(range(len(df)), df["simpson_index"], color=colors)
    axes[1].set_xlabel("Sample (sorted by entropy)")
    axes[1].set_ylabel("Simpson Index")
    axes[1].set_title("Simpson Diversity per Sample")
    axes[1].axhline(
        df["simpson_index"].mean(),
        color="red",
        linestyle="--",
        label=f"Mean: {df['simpson_index'].mean():.2f}",
    )
    axes[1].legend()

    # Add region legend if available
    if region_col in df.columns:
        from matplotlib.patches import Patch
        legend_patches = [
            Patch(facecolor=region_colors[r], label=r)
            for r in regions
        ]
        fig.legend(
            handles=legend_patches,
            title="Region",
            loc="upper right",
            bbox_to_anchor=(0.99, 0.99),
        )

    plt.tight_layout()
    return _save_figure(fig, output_path, dpi)


def plot_evenness_distribution(
    diversity_by_sample: pd.DataFrame,
    output_path: Path,
    dpi: int = 200,
    figsize: Tuple[int, int] = (8, 6),
) -> Path:
    """Plot evenness distribution as histogram.

    Parameters
    ----------
    diversity_by_sample : pd.DataFrame
        Output from compute_diversity_by_group()
    output_path : Path
        Output file path
    dpi : int
        Figure resolution
    figsize : Tuple[int, int]
        Figure size

    Returns
    -------
    Path
        Path to saved figure
    """
    import matplotlib.pyplot as plt

    _set_style()

    fig, ax = plt.subplots(figsize=figsize)

    ax.hist(
        diversity_by_sample["evenness"],
        bins=20,
        color="steelblue",
        edgecolor="white",
        alpha=0.7,
    )

    mean_evenness = diversity_by_sample["evenness"].mean()
    ax.axvline(mean_evenness, color="red", linestyle="--", linewidth=2)
    ax.text(
        mean_evenness + 0.02,
        ax.get_ylim()[1] * 0.9,
        f"Mean: {mean_evenness:.3f}",
        color="red",
    )

    ax.set_xlabel("Evenness (Pielou's J)")
    ax.set_ylabel("Number of Samples")
    ax.set_title("Distribution of Cell Type Evenness")
    ax.set_xlim(0, 1)

    plt.tight_layout()
    return _save_figure(fig, output_path, dpi)
