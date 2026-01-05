"""Composition visualization functions.

Provides:
- Global composition bar chart
- Regional composition stacked bars
- Donor composition boxplots
- Composition heatmap with clustering
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Tuple
import logging

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def _sort_regions(regions, region_order: Optional[List[str]] = None) -> List[str]:
    """Sort regions by specified order or alphabetically.

    Parameters
    ----------
    regions : iterable
        Regions to sort
    region_order : List[str], optional
        Canonical order. If None, sorts alphabetically.

    Returns
    -------
    List[str]
        Sorted regions
    """
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


def plot_composition_global(
    adata,
    output_path: Path,
    cell_type_col: str = "cell_type_phenocycler",
    top_n: int = 15,
    dpi: int = 200,
    figsize: Tuple[int, int] = (12, 8),
) -> Path:
    """Plot global cell type composition as horizontal bar chart.

    Parameters
    ----------
    adata : AnnData
        Annotated data
    output_path : Path
        Output file path
    cell_type_col : str
        Column with cell type labels
    top_n : int
        Number of top cell types to show
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

    # Compute composition
    counts = adata.obs[cell_type_col].value_counts()
    total = counts.sum()
    proportions = (counts / total * 100).head(top_n)

    # Group remaining as "Other"
    if len(counts) > top_n:
        other_count = counts.iloc[top_n:].sum()
        other_pct = other_count / total * 100
        proportions["Other"] = other_pct

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Horizontal bar chart
    y_pos = np.arange(len(proportions))
    colors = plt.cm.tab20(np.linspace(0, 1, len(proportions)))

    bars = ax.barh(y_pos, proportions.values, color=colors)

    # Add value labels
    for bar, pct in zip(bars, proportions.values):
        ax.text(
            bar.get_width() + 0.5,
            bar.get_y() + bar.get_height() / 2,
            f"{pct:.1f}%",
            va="center",
            fontsize=9,
        )

    ax.set_yticks(y_pos)
    ax.set_yticklabels(proportions.index)
    ax.invert_yaxis()
    ax.set_xlabel("Proportion (%)")
    ax.set_title(f"Cell Type Composition (n={total:,} cells)")
    ax.set_xlim(0, max(proportions.values) * 1.15)

    return _save_figure(fig, output_path, dpi)


def plot_composition_by_region(
    adata,
    output_path: Path,
    cell_type_col: str = "cell_type_phenocycler",
    region_col: str = "region",
    top_n: int = 12,
    dpi: int = 200,
    figsize: Tuple[int, int] = (12, 7),
    region_order: Optional[List[str]] = None,
) -> Optional[Path]:
    """Plot cell type composition by region as stacked bar chart.

    Parameters
    ----------
    adata : AnnData
        Annotated data
    output_path : Path
        Output file path
    cell_type_col : str
        Column with cell type labels
    region_col : str
        Column with region labels
    top_n : int
        Number of top cell types to show
    dpi : int
        Figure resolution
    figsize : Tuple[int, int]
        Figure size
    region_order : List[str], optional
        Canonical region ordering

    Returns
    -------
    Optional[Path]
        Path to saved figure, or None if region column not found
    """
    import matplotlib.pyplot as plt

    _set_style()

    if region_col not in adata.obs.columns:
        logger.warning(f"Region column '{region_col}' not found, skipping plot")
        return None

    # Compute composition per region
    df = adata.obs[[cell_type_col, region_col]].copy()
    composition = df.groupby([region_col, cell_type_col], observed=True).size().unstack(fill_value=0)

    # Normalize to percentages
    composition_pct = composition.div(composition.sum(axis=1), axis=0) * 100

    # Get top cell types globally
    global_counts = adata.obs[cell_type_col].value_counts()
    top_types = global_counts.head(top_n).index.tolist()

    # Group others
    other_types = [c for c in composition_pct.columns if c not in top_types]
    if other_types:
        composition_pct["Other"] = composition_pct[other_types].sum(axis=1)
        composition_pct = composition_pct.drop(columns=other_types)

    # Reorder columns
    col_order = [c for c in top_types if c in composition_pct.columns]
    if "Other" in composition_pct.columns:
        col_order.append("Other")
    composition_pct = composition_pct[col_order]

    # Reorder rows by region order
    sorted_regions = _sort_regions(composition_pct.index, region_order)
    composition_pct = composition_pct.reindex(sorted_regions)

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    colors = plt.cm.tab20(np.linspace(0, 1, len(col_order)))

    composition_pct.plot(
        kind="bar",
        stacked=True,
        ax=ax,
        color=colors,
        width=0.8,
        edgecolor="white",
        linewidth=0.5,
    )

    ax.set_xlabel("Region")
    ax.set_ylabel("Proportion (%)")
    ax.set_title("Cell Type Composition by Region")
    ax.legend(
        title="Cell Type",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        fontsize=8,
    )
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")

    plt.tight_layout()
    return _save_figure(fig, output_path, dpi)


def plot_composition_by_donor(
    adata,
    output_path: Path,
    cell_type_col: str = "cell_type_phenocycler",
    donor_col: str = "donor",
    top_n: int = 10,
    dpi: int = 200,
    figsize: Tuple[int, int] = (14, 7),
) -> Optional[Path]:
    """Plot cell type composition by donor with boxplots.

    Parameters
    ----------
    adata : AnnData
        Annotated data
    output_path : Path
        Output file path
    cell_type_col : str
        Column with cell type labels
    donor_col : str
        Column with donor labels
    top_n : int
        Number of top cell types to show
    dpi : int
        Figure resolution
    figsize : Tuple[int, int]
        Figure size

    Returns
    -------
    Optional[Path]
        Path to saved figure, or None if donor column not found
    """
    import matplotlib.pyplot as plt

    _set_style()

    if donor_col not in adata.obs.columns:
        logger.warning(f"Donor column '{donor_col}' not found, skipping plot")
        return None

    # Compute composition per donor
    df = adata.obs[[cell_type_col, donor_col]].copy()
    composition = df.groupby([donor_col, cell_type_col], observed=True).size().unstack(fill_value=0)
    composition_pct = composition.div(composition.sum(axis=1), axis=0) * 100

    # Get top cell types
    global_counts = adata.obs[cell_type_col].value_counts()
    top_types = global_counts.head(top_n).index.tolist()

    # Prepare data for boxplot
    plot_data = []
    for ct in top_types:
        if ct in composition_pct.columns:
            for donor, pct in composition_pct[ct].items():
                plot_data.append({
                    "cell_type": ct,
                    "donor": donor,
                    "percentage": pct,
                })

    plot_df = pd.DataFrame(plot_data)

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Try seaborn, fallback to matplotlib
    try:
        import seaborn as sns
        colors = plt.cm.tab20(np.linspace(0, 1, len(top_types)))
        palette = {ct: colors[i] for i, ct in enumerate(top_types)}

        sns.boxplot(
            data=plot_df,
            x="cell_type",
            y="percentage",
            hue="cell_type",
            palette=palette,
            ax=ax,
            legend=False,
        )
        sns.stripplot(
            data=plot_df,
            x="cell_type",
            y="percentage",
            color="black",
            alpha=0.5,
            size=4,
            ax=ax,
        )
    except ImportError:
        # Fallback without seaborn
        positions = range(len(top_types))
        for i, ct in enumerate(top_types):
            ct_data = plot_df[plot_df["cell_type"] == ct]["percentage"]
            if len(ct_data) > 0:
                ax.boxplot(ct_data, positions=[i], widths=0.6)
        ax.set_xticks(positions)
        ax.set_xticklabels(top_types)

    ax.set_xlabel("Cell Type")
    ax.set_ylabel("Proportion (%)")
    ax.set_title(f"Cell Type Composition Across Donors (n={composition.shape[0]} donors)")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")

    plt.tight_layout()
    return _save_figure(fig, output_path, dpi)


def plot_composition_heatmap(
    adata,
    output_path: Path,
    cell_type_col: str = "cell_type_phenocycler",
    groupby_col: str = "sample_id",
    region_col: str = "region",
    top_n: int = 20,
    dpi: int = 150,
    figsize: Optional[Tuple[int, int]] = None,
) -> Optional[Path]:
    """Plot composition heatmap with hierarchical clustering.

    Parameters
    ----------
    adata : AnnData
        Annotated data
    output_path : Path
        Output file path
    cell_type_col : str
        Column with cell type labels
    groupby_col : str
        Column to group by (sample_id, donor, region)
    region_col : str
        Column with region labels for row annotation
    top_n : int
        Number of top cell types to show
    dpi : int
        Figure resolution
    figsize : Optional[Tuple[int, int]]
        Figure size (auto-computed if None)

    Returns
    -------
    Optional[Path]
        Path to saved figure
    """
    import matplotlib.pyplot as plt

    _set_style()

    # Find suitable groupby column
    if groupby_col not in adata.obs.columns:
        for fallback in ["sample_id", "donor", "region"]:
            if fallback in adata.obs.columns:
                groupby_col = fallback
                break
        else:
            logger.warning("No suitable groupby column found for heatmap")
            return None

    # Compute composition matrix
    df = adata.obs[[cell_type_col, groupby_col]].copy()
    composition = df.groupby([groupby_col, cell_type_col], observed=True).size().unstack(fill_value=0)
    composition_pct = composition.div(composition.sum(axis=1), axis=0) * 100

    # Get top cell types
    global_counts = adata.obs[cell_type_col].value_counts()
    top_types = global_counts.head(top_n).index.tolist()
    top_types = [t for t in top_types if t in composition_pct.columns]

    composition_pct = composition_pct[top_types]

    # Auto-compute figure size
    if figsize is None:
        n_samples = composition_pct.shape[0]
        n_types = composition_pct.shape[1]
        figsize = (max(10, n_types * 0.5), max(6, n_samples * 0.25))

    # Try seaborn clustermap
    try:
        import seaborn as sns
        from matplotlib.patches import Patch

        # Build row colors if region available
        row_colors = None
        region_palette = None
        unique_regions = None
        if region_col in adata.obs.columns and groupby_col == "sample_id":
            sample_region = adata.obs.groupby(groupby_col, observed=True)[region_col].first()
            unique_regions = sample_region.unique()
            colors = [tuple(c) for c in plt.cm.Set2(np.linspace(0, 1, len(unique_regions)))]
            region_palette = dict(zip(unique_regions, colors))
            row_colors = composition_pct.index.map(
                lambda x: region_palette.get(sample_region.get(x, "unknown"), (0.74, 0.76, 0.78, 1.0))
            )
            row_colors = pd.Series(row_colors, index=composition_pct.index, name="Region")

        # Adjust figure size
        figsize = (figsize[0] + 1.0, figsize[1] + 0.5)

        g = sns.clustermap(
            composition_pct,
            cmap="YlOrRd",
            figsize=figsize,
            xticklabels=True,
            yticklabels=True,
            linewidths=0.5,
            row_colors=row_colors,
            row_cluster=True,
            col_cluster=False,
            dendrogram_ratio=(0.001, 0.001),
            cbar_pos=(-0.05, 0.82, 0.015, 0.12),
            cbar_kws={"label": "Proportion (%)"},
        )

        # Hide the dendrograms
        g.ax_row_dendrogram.clear()
        g.ax_row_dendrogram.set_axis_off()
        g.ax_col_dendrogram.clear()
        g.ax_col_dendrogram.set_axis_off()

        # Rotate x-axis labels
        g.ax_heatmap.set_xticklabels(
            g.ax_heatmap.get_xticklabels(),
            rotation=45,
            ha="right",
            fontsize=9
        )

        g.ax_heatmap.set_xlabel("Cell Type", labelpad=10)
        g.ax_heatmap.set_ylabel(groupby_col.replace("_", " ").title())

        # Add region legend if available
        if region_palette is not None and unique_regions is not None:
            legend_handles = [
                Patch(facecolor=region_palette[region], edgecolor="black", linewidth=0.5, label=region)
                for region in unique_regions
            ]
            g.ax_heatmap.legend(
                handles=legend_handles,
                title="Region",
                loc="upper center",
                bbox_to_anchor=(0.5, -0.18),
                ncol=min(len(unique_regions), 6),
                fontsize=9,
                frameon=True,
            )

        g.savefig(output_path, dpi=dpi, bbox_inches="tight", pad_inches=0.1)
        plt.close(g.fig)

        return Path(output_path)

    except ImportError:
        # Fallback to simple heatmap
        fig, ax = plt.subplots(figsize=figsize)

        im = ax.imshow(composition_pct.values, aspect="auto", cmap="YlOrRd")
        ax.set_xticks(range(len(top_types)))
        ax.set_xticklabels(top_types, rotation=45, ha="right")
        ax.set_yticks(range(len(composition_pct)))
        ax.set_yticklabels(composition_pct.index)
        ax.set_xlabel("Cell Type")
        ax.set_ylabel(groupby_col.replace("_", " ").title())
        ax.set_title("Cell Type Composition Heatmap")

        plt.colorbar(im, ax=ax, label="Proportion (%)")
        plt.tight_layout()

        return _save_figure(fig, output_path, dpi)
