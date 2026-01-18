"""Local diversity distribution visualization."""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Tuple

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def plot_diversity_distribution(
    diversity_by_sample: pd.DataFrame,
    output_path: Path,
    figsize: Tuple[int, int] = (10, 6),
    dpi: int = 200,
    title: str = "Local Cell-Type Diversity",
    organ: Optional[str] = None,
    region_order: Optional[List[str]] = None,
) -> None:
    """Plot distribution of local diversity by region.

    Parameters
    ----------
    diversity_by_sample : pd.DataFrame
        Per-sample diversity statistics
        Columns: sample_id, local_diversity_mean, local_diversity_std, ...
        Optional: region
    output_path : Path
        Path to save figure
    figsize : Tuple[int, int]
        Figure size
    dpi : int
        Figure resolution
    title : str
        Plot title
    organ : str, optional
        Organ name for region ordering (e.g., 'fallopian_tube', 'uterus')
    region_order : List[str], optional
        Explicit region order (takes precedence over organ)
    """
    if diversity_by_sample.empty:
        return

    fig, ax = plt.subplots(figsize=figsize)

    if "region" in diversity_by_sample.columns:
        # Boxplot by region
        regions = list(diversity_by_sample["region"].unique())

        # Determine region order
        if region_order is not None:
            # Use explicit region order
            ordered_regions = [r for r in region_order if r in regions]
            other_regions = [r for r in regions if r not in region_order]
            ordered_regions = ordered_regions + sorted(other_regions)
        elif organ:
            # Use organ config for region ordering
            try:
                from celltype_refinery.config import get_organ_config
                organ_config = get_organ_config(organ)
                ordered_regions = organ_config.sort_regions(regions)
            except (ImportError, ValueError):
                # Fallback to alphabetical
                ordered_regions = sorted(regions)
        else:
            # Default: alphabetical order
            ordered_regions = sorted(regions)

        sns.boxplot(
            data=diversity_by_sample,
            x="region",
            y="local_diversity_mean",
            order=ordered_regions,
            ax=ax,
            palette="Set2",
        )

        sns.stripplot(
            data=diversity_by_sample,
            x="region",
            y="local_diversity_mean",
            order=ordered_regions,
            ax=ax,
            color="black",
            size=4,
            alpha=0.5,
        )

        ax.set_xlabel("Region")
    else:
        # Histogram if no region
        sns.histplot(
            data=diversity_by_sample,
            x="local_diversity_mean",
            ax=ax,
            bins=20,
            color="steelblue",
            edgecolor="black",
        )

        ax.set_xlabel("Mean Local Diversity")

    ax.set_ylabel("Mean Local Diversity (Shannon H)")
    ax.set_title(title, fontsize=14)

    # Add annotation with global mean
    global_mean = diversity_by_sample["local_diversity_mean"].mean()
    ax.axhline(
        y=global_mean,
        color="red",
        linestyle="--",
        alpha=0.7,
        label=f"Global mean: {global_mean:.2f}",
    )
    ax.legend(loc="upper right")

    plt.tight_layout()
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
