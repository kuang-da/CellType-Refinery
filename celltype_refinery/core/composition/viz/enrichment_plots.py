"""Enrichment visualization functions.

Provides:
- Enrichment heatmap (cell types x regions)
- Top enrichments bar chart
- Volcano plot
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


def plot_enrichment_heatmap(
    enrichment_df: pd.DataFrame,
    output_path: Path,
    value_col: str = "fold_change",
    region_order: Optional[List[str]] = None,
    dpi: int = 200,
    figsize: Optional[Tuple[int, int]] = None,
    top_n: int = 20,
) -> Optional[Path]:
    """Plot regional enrichment as heatmap.

    Shows fold change values with significance markers.

    Parameters
    ----------
    enrichment_df : pd.DataFrame
        Output from compute_regional_enrichment()
    output_path : Path
        Output file path
    value_col : str
        Column to use for heatmap values
    region_order : List[str], optional
        Canonical region order
    dpi : int
        Figure resolution
    figsize : Optional[Tuple[int, int]]
        Figure size (auto-computed if None)
    top_n : int
        Number of top cell types to show

    Returns
    -------
    Optional[Path]
        Path to saved figure
    """
    import matplotlib.pyplot as plt

    _set_style()

    if enrichment_df.empty:
        logger.warning("No enrichment data to plot")
        return None

    # Filter to significant or top fold changes
    sig_df = enrichment_df[enrichment_df["significant"]].copy()

    if len(sig_df) == 0:
        # Fall back to top fold changes
        sig_df = enrichment_df.nlargest(min(50, len(enrichment_df)), "fold_change")

    # Get top cell types by max fold change
    top_cell_types = (
        sig_df.groupby("cell_type")["fold_change"]
        .max()
        .nlargest(top_n)
        .index.tolist()
    )

    # Filter data
    plot_df = enrichment_df[enrichment_df["cell_type"].isin(top_cell_types)].copy()

    # Create pivot table
    pivot = plot_df.pivot_table(
        index="cell_type",
        columns="region",
        values=value_col,
        fill_value=1.0,  # No change
    )

    # Reorder columns by region order
    sorted_regions = _sort_regions(pivot.columns, region_order)
    pivot = pivot[sorted_regions]

    # Create significance mask
    sig_pivot = plot_df.pivot_table(
        index="cell_type",
        columns="region",
        values="significant",
        fill_value=False,
    )
    sig_pivot = sig_pivot[sorted_regions]

    # Auto-compute figure size
    if figsize is None:
        n_types = len(pivot)
        n_regions = len(pivot.columns)
        figsize = (max(8, n_regions * 1.5), max(6, n_types * 0.4))

    fig, ax = plt.subplots(figsize=figsize)

    # Create heatmap
    try:
        import seaborn as sns

        # Log transform fold change for better visualization
        plot_values = np.log2(pivot.values.clip(min=0.1))

        sns.heatmap(
            plot_values,
            ax=ax,
            cmap="RdBu_r",
            center=0,
            xticklabels=pivot.columns,
            yticklabels=pivot.index,
            cbar_kws={"label": "log2(Fold Change)"},
            linewidths=0.5,
        )

        # Add significance markers
        for i, cell_type in enumerate(pivot.index):
            for j, region in enumerate(pivot.columns):
                if sig_pivot.loc[cell_type, region]:
                    ax.text(
                        j + 0.5, i + 0.5, "*",
                        ha="center", va="center",
                        fontsize=12, color="black", fontweight="bold",
                    )

    except ImportError:
        # Fallback without seaborn
        plot_values = np.log2(pivot.values.clip(min=0.1))
        im = ax.imshow(plot_values, aspect="auto", cmap="RdBu_r")
        ax.set_xticks(range(len(pivot.columns)))
        ax.set_xticklabels(pivot.columns)
        ax.set_yticks(range(len(pivot.index)))
        ax.set_yticklabels(pivot.index)
        plt.colorbar(im, ax=ax, label="log2(Fold Change)")

    ax.set_xlabel("Region")
    ax.set_ylabel("Cell Type")
    ax.set_title("Regional Enrichment (* = significant)")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")

    plt.tight_layout()
    return _save_figure(fig, output_path, dpi)


def plot_top_enrichments(
    enrichment_df: pd.DataFrame,
    output_path: Path,
    top_n: int = 15,
    dpi: int = 200,
    figsize: Tuple[int, int] = (10, 8),
) -> Optional[Path]:
    """Plot top significant enrichments as horizontal bar chart.

    Parameters
    ----------
    enrichment_df : pd.DataFrame
        Output from compute_regional_enrichment()
    output_path : Path
        Output file path
    top_n : int
        Number of top enrichments to show
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

    if enrichment_df.empty:
        logger.warning("No enrichment data to plot")
        return None

    # Get significant enrichments
    sig_df = enrichment_df[enrichment_df["significant"]].copy()

    if len(sig_df) == 0:
        logger.warning("No significant enrichments to plot")
        return None

    # Get top by fold change
    top_df = sig_df.nlargest(top_n, "fold_change")

    # Create labels
    top_df["label"] = top_df["cell_type"] + " (" + top_df["region"] + ")"

    fig, ax = plt.subplots(figsize=figsize)

    # Color by region
    regions = top_df["region"].unique()
    region_colors = dict(zip(regions, plt.cm.Set2(np.linspace(0, 1, len(regions)))))
    colors = top_df["region"].map(region_colors)

    y_pos = np.arange(len(top_df))
    bars = ax.barh(y_pos, top_df["fold_change"], color=colors)

    # Add significance stars
    for i, (_, row) in enumerate(top_df.iterrows()):
        pval_adj = row["pvalue_adj"]
        if pval_adj < 0.001:
            stars = "***"
        elif pval_adj < 0.01:
            stars = "**"
        else:
            stars = "*"
        ax.text(
            row["fold_change"] + 0.1,
            i,
            stars,
            va="center",
            fontsize=10,
        )

    ax.set_yticks(y_pos)
    ax.set_yticklabels(top_df["label"])
    ax.invert_yaxis()
    ax.set_xlabel("Fold Change")
    ax.set_title(f"Top {len(top_df)} Significant Regional Enrichments")

    # Add legend for regions
    from matplotlib.patches import Patch
    legend_patches = [
        Patch(facecolor=region_colors[r], label=r)
        for r in regions
    ]
    ax.legend(handles=legend_patches, title="Region", loc="lower right")

    plt.tight_layout()
    return _save_figure(fig, output_path, dpi)


def plot_enrichment_volcano(
    enrichment_df: pd.DataFrame,
    output_path: Path,
    dpi: int = 200,
    figsize: Tuple[int, int] = (10, 8),
) -> Optional[Path]:
    """Plot enrichment results as volcano plot.

    X-axis: log2(fold change)
    Y-axis: -log10(p-value adjusted)

    Parameters
    ----------
    enrichment_df : pd.DataFrame
        Output from compute_regional_enrichment()
    output_path : Path
        Output file path
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

    if enrichment_df.empty:
        logger.warning("No enrichment data to plot")
        return None

    df = enrichment_df.copy()

    # Compute volcano coordinates
    df["log2_fc"] = np.log2(df["fold_change"].clip(min=0.1))
    df["neg_log10_pval"] = -np.log10(df["pvalue_adj"].clip(min=1e-10))

    fig, ax = plt.subplots(figsize=figsize)

    # Color by significance
    colors = np.where(df["significant"], "#e74c3c", "#bdc3c7")

    ax.scatter(
        df["log2_fc"],
        df["neg_log10_pval"],
        c=colors,
        alpha=0.6,
        s=50,
    )

    # Add threshold lines
    ax.axhline(-np.log10(0.05), color="gray", linestyle="--", alpha=0.5)
    ax.axvline(np.log2(1.5), color="gray", linestyle="--", alpha=0.5)
    ax.axvline(-np.log2(1.5), color="gray", linestyle="--", alpha=0.5)

    # Label top significant points
    top_sig = df[df["significant"]].nlargest(10, "neg_log10_pval")
    for _, row in top_sig.iterrows():
        ax.annotate(
            f"{row['cell_type']}\n({row['region']})",
            (row["log2_fc"], row["neg_log10_pval"]),
            fontsize=7,
            alpha=0.8,
        )

    ax.set_xlabel("log2(Fold Change)")
    ax.set_ylabel("-log10(Adjusted P-value)")
    ax.set_title("Regional Enrichment Volcano Plot")

    plt.tight_layout()
    return _save_figure(fig, output_path, dpi)
