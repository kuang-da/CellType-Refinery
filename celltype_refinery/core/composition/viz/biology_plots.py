"""Biology metrics visualization functions.

Provides tissue-configurable biology metrics plots:
- Multi-panel biology metrics plot
- Individual metric plots by region
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


# Default metrics to plot (can be customized)
DEFAULT_BIOLOGY_METRICS = [
    ("epithelial_pct", "Epithelial (%)"),
    ("stromal_pct", "Stromal (%)"),
    ("smooth_muscle_pct", "Smooth Muscle (%)"),
    ("immune_pct", "Immune (%)"),
    ("endothelial_pct", "Endothelial (%)"),
]


def plot_biology_metrics(
    biology_by_sample: pd.DataFrame,
    output_path: Path,
    region_col: str = "region",
    metrics: Optional[List[Tuple[str, str]]] = None,
    region_order: Optional[List[str]] = None,
    dpi: int = 200,
    figsize: Optional[Tuple[int, int]] = None,
    title: str = "Biology Metrics by Region",
) -> Path:
    """Create multi-panel plot of biology metrics.

    Parameters
    ----------
    biology_by_sample : pd.DataFrame
        Output from compute_biology_metrics()
    output_path : Path
        Output file path
    region_col : str
        Column with region labels
    metrics : List[Tuple[str, str]], optional
        List of (column_name, display_label) tuples.
        If None, auto-detects available metrics.
    region_order : List[str], optional
        Canonical region order
    dpi : int
        Figure resolution
    figsize : Tuple[int, int], optional
        Figure size (auto-computed if None)
    title : str
        Overall plot title

    Returns
    -------
    Path
        Path to saved figure
    """
    import matplotlib.pyplot as plt

    _set_style()

    # Auto-detect metrics if not provided
    if metrics is None:
        # Find columns ending in _pct or _ratio
        available_metrics = []
        for col in biology_by_sample.columns:
            if col.endswith("_pct"):
                label = col.replace("_pct", "").replace("_", " ").title() + " (%)"
                available_metrics.append((col, label))
            elif col.endswith("_ratio"):
                label = col.replace("_ratio", "").replace("_", " ").title() + " Ratio"
                available_metrics.append((col, label))

        if not available_metrics:
            available_metrics = DEFAULT_BIOLOGY_METRICS

        metrics = available_metrics[:6]  # Limit to 6 for 2x3 grid

    # Filter to available metrics
    available_cols = set(biology_by_sample.columns)
    metrics = [(col, label) for col, label in metrics if col in available_cols]

    if not metrics:
        logger.warning("No biology metrics available to plot")
        # Create a simple placeholder figure
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.text(0.5, 0.5, "No biology metrics available", ha="center", va="center")
        ax.set_axis_off()
        return _save_figure(fig, output_path, dpi)

    # Determine grid size
    n_metrics = len(metrics)
    if n_metrics <= 3:
        nrows, ncols = 1, n_metrics
    elif n_metrics <= 6:
        nrows, ncols = 2, 3
    else:
        nrows, ncols = 3, 3

    # Auto-compute figure size
    if figsize is None:
        figsize = (ncols * 4.5, nrows * 4)

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
    if n_metrics == 1:
        axes = [axes]
    else:
        axes = axes.flatten()

    has_region = region_col in biology_by_sample.columns

    for i, (metric, label) in enumerate(metrics):
        ax = axes[i]

        if metric not in biology_by_sample.columns:
            ax.text(0.5, 0.5, "Data not available", ha="center", va="center")
            ax.set_title(label)
            continue

        if has_region:
            # Boxplot by region
            sorted_regions = _sort_regions(
                biology_by_sample[region_col].unique(),
                region_order
            )
            try:
                import seaborn as sns
                sns.boxplot(
                    data=biology_by_sample,
                    x=region_col,
                    y=metric,
                    hue=region_col,
                    order=sorted_regions,
                    hue_order=sorted_regions,
                    ax=ax,
                    palette="Set2",
                    legend=False,
                )
                sns.stripplot(
                    data=biology_by_sample,
                    x=region_col,
                    y=metric,
                    order=sorted_regions,
                    ax=ax,
                    color="black",
                    alpha=0.5,
                    size=4,
                )
            except ImportError:
                data = [
                    biology_by_sample[biology_by_sample[region_col] == r][metric]
                    for r in sorted_regions
                ]
                ax.boxplot(data, labels=sorted_regions)

            plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
            ax.set_xlabel("")
        else:
            # Bar plot by sample (if no region)
            ax.bar(
                range(len(biology_by_sample)),
                biology_by_sample[metric],
                color="steelblue",
            )
            ax.set_xlabel("Sample")
            ax.set_xticks([])

        ax.set_title(label)
        ax.set_ylabel("")

    # Hide unused axes
    for i in range(len(metrics), len(axes)):
        axes[i].set_visible(False)

    plt.suptitle(title, fontsize=14, y=1.02)
    plt.tight_layout()
    return _save_figure(fig, output_path, dpi)


def plot_ratio_by_region(
    biology_by_sample: pd.DataFrame,
    output_path: Path,
    ratio_col: str,
    ratio_label: str,
    region_col: str = "region",
    region_order: Optional[List[str]] = None,
    reference_lines: Optional[List[Tuple[float, str, str]]] = None,
    dpi: int = 200,
    figsize: Tuple[int, int] = (10, 6),
) -> Optional[Path]:
    """Plot a ratio metric by region.

    Parameters
    ----------
    biology_by_sample : pd.DataFrame
        Output from compute_biology_metrics()
    output_path : Path
        Output file path
    ratio_col : str
        Column name for the ratio
    ratio_label : str
        Display label for the ratio
    region_col : str
        Column with region labels
    region_order : List[str], optional
        Canonical region order
    reference_lines : List[Tuple[float, str, str]], optional
        List of (value, color, label) for reference lines
    dpi : int
        Figure resolution
    figsize : Tuple[int, int]
        Figure size

    Returns
    -------
    Optional[Path]
        Path to saved figure
    """
    import matplotlib.pyplot as plt

    _set_style()

    if region_col not in biology_by_sample.columns:
        logger.warning(f"Region column '{region_col}' not found")
        return None

    if ratio_col not in biology_by_sample.columns:
        logger.warning(f"{ratio_col} column not found")
        return None

    fig, ax = plt.subplots(figsize=figsize)

    sorted_regions = _sort_regions(
        biology_by_sample[region_col].unique(),
        region_order
    )

    try:
        import seaborn as sns
        sns.boxplot(
            data=biology_by_sample,
            x=region_col,
            y=ratio_col,
            order=sorted_regions,
            ax=ax,
            palette="Set2",
        )
        sns.stripplot(
            data=biology_by_sample,
            x=region_col,
            y=ratio_col,
            order=sorted_regions,
            ax=ax,
            color="black",
            alpha=0.5,
            size=5,
        )
    except ImportError:
        data = [
            biology_by_sample[biology_by_sample[region_col] == r][ratio_col]
            for r in sorted_regions
        ]
        ax.boxplot(data, labels=sorted_regions)

    ax.set_xlabel("Region")
    ax.set_ylabel(ratio_label)
    ax.set_title(f"{ratio_label} by Region")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")

    # Add reference lines
    if reference_lines:
        for value, color, label in reference_lines:
            ax.axhline(value, color=color, linestyle="--", alpha=0.5, label=label)
        ax.legend()

    plt.tight_layout()
    return _save_figure(fig, output_path, dpi)
