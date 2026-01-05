"""Composition visualization functions for consolidation results.

Provides:
- Global composition stacked bar chart
- Regional composition comparison
- Donor-level composition
- Before/after consolidation comparison
- Composition heatmap with clustering
"""

from pathlib import Path
from typing import Dict, List, Optional, Union
import logging

import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


def plot_composition_global(
    adata,
    output_path: Union[str, Path],
    cell_type_col: str = "cell_type_phenocycler",
    top_n: int = 15,
    dpi: int = 200,
) -> Path:
    """Plot global cell type composition as stacked bar chart.

    Args:
        adata: AnnData with cell type annotations
        output_path: Path to save figure
        cell_type_col: Column with cell type labels
        top_n: Number of top cell types to show (rest grouped as "Other")
        dpi: Figure resolution

    Returns:
        Path to saved figure
    """
    import matplotlib.pyplot as plt
    from .style import set_publication_style, save_figure, get_color_palette

    set_publication_style()

    # Compute composition
    counts = adata.obs[cell_type_col].value_counts()
    total = counts.sum()
    proportions = (counts / total * 100).head(top_n)

    # Group remaining as "Other"
    if len(counts) > top_n:
        other_count = counts.iloc[top_n:].sum()
        other_pct = other_count / total * 100
        proportions["Other"] = other_pct

    # Get colors
    labels = proportions.index.tolist()
    colors = get_color_palette(labels, palette_type="cell_type")

    # Create figure
    fig, ax = plt.subplots(figsize=(12, 6))

    # Horizontal bar chart
    y_pos = np.arange(len(labels))
    bars = ax.barh(y_pos, proportions.values, color=[colors[l] for l in labels])

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
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("Proportion (%)")
    ax.set_title(f"Cell Type Composition (n={total:,} cells)")
    ax.set_xlim(0, max(proportions.values) * 1.15)

    return save_figure(fig, output_path, dpi=dpi)


def plot_composition_by_region(
    adata,
    output_path: Union[str, Path],
    cell_type_col: str = "cell_type_phenocycler",
    region_col: str = "region",
    top_n: int = 12,
    dpi: int = 200,
) -> Path:
    """Plot cell type composition by anatomical region.

    Args:
        adata: AnnData with cell type and region annotations
        output_path: Path to save figure
        cell_type_col: Column with cell type labels
        region_col: Column with region labels
        top_n: Number of top cell types to show
        dpi: Figure resolution

    Returns:
        Path to saved figure
    """
    import matplotlib.pyplot as plt
    from .style import set_publication_style, save_figure, get_color_palette

    set_publication_style()

    if region_col not in adata.obs.columns:
        logger.warning(f"Region column '{region_col}' not found, skipping")
        return None

    # Compute composition per region
    df = adata.obs[[cell_type_col, region_col]].copy()
    composition = df.groupby([region_col, cell_type_col]).size().unstack(fill_value=0)

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

    # Get colors
    colors = get_color_palette(col_order, palette_type="cell_type")
    color_list = [colors[c] for c in col_order]

    # Create stacked bar chart
    fig, ax = plt.subplots(figsize=(10, 6))

    composition_pct.plot(
        kind="bar",
        stacked=True,
        ax=ax,
        color=color_list,
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
    return save_figure(fig, output_path, dpi=dpi)


def plot_composition_by_donor(
    adata,
    output_path: Union[str, Path],
    cell_type_col: str = "cell_type_phenocycler",
    donor_col: str = "donor",
    top_n: int = 10,
    dpi: int = 200,
) -> Path:
    """Plot cell type composition by donor with boxplots.

    Args:
        adata: AnnData with cell type and donor annotations
        output_path: Path to save figure
        cell_type_col: Column with cell type labels
        donor_col: Column with donor labels
        top_n: Number of top cell types to show
        dpi: Figure resolution

    Returns:
        Path to saved figure
    """
    import matplotlib.pyplot as plt
    from .style import set_publication_style, save_figure, get_color_palette

    set_publication_style()

    if donor_col not in adata.obs.columns:
        logger.warning(f"Donor column '{donor_col}' not found, skipping")
        return None

    # Compute composition per donor
    df = adata.obs[[cell_type_col, donor_col]].copy()
    composition = df.groupby([donor_col, cell_type_col]).size().unstack(fill_value=0)
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

    # Get colors
    colors = get_color_palette(top_types, palette_type="cell_type")

    # Create boxplot
    fig, ax = plt.subplots(figsize=(14, 6))

    try:
        import seaborn as sns
        sns.boxplot(
            data=plot_df,
            x="cell_type",
            y="percentage",
            palette=colors,
            ax=ax,
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
            ax.boxplot(
                ct_data,
                positions=[i],
                patch_artist=True,
                boxprops=dict(facecolor=colors.get(ct, "#bdc3c7")),
            )

    ax.set_xlabel("Cell Type")
    ax.set_ylabel("Proportion (%)")
    ax.set_title(f"Cell Type Composition Across Donors (n={composition.shape[0]} donors)")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")

    plt.tight_layout()
    return save_figure(fig, output_path, dpi=dpi)


def plot_composition_before_after(
    adata,
    mapping_df: pd.DataFrame,
    output_path: Union[str, Path],
    before_col: str = "cell_type_lvl1",
    after_col: str = "cell_type_phenocycler",
    top_n: int = 12,
    dpi: int = 200,
) -> Path:
    """Plot before/after consolidation composition comparison.

    Args:
        adata: AnnData with both before and after columns
        mapping_df: DataFrame with label mapping (for n_cells aggregation)
        output_path: Path to save figure
        before_col: Column with original labels
        after_col: Column with final labels
        top_n: Number of top cell types to show
        dpi: Figure resolution

    Returns:
        Path to saved figure
    """
    import matplotlib.pyplot as plt
    from .style import set_publication_style, save_figure, get_color_palette

    set_publication_style()

    # Compute before/after compositions
    if before_col not in adata.obs.columns:
        before_col = "cell_type_lvl1" if "cell_type_lvl1" in adata.obs.columns else after_col

    before_counts = adata.obs[before_col].value_counts()
    after_counts = adata.obs[after_col].value_counts()

    total = before_counts.sum()
    before_pct = (before_counts / total * 100)
    after_pct = (after_counts / total * 100)

    # Get union of top types
    all_types = set(before_pct.head(top_n).index) | set(after_pct.head(top_n).index)
    all_types = sorted(all_types, key=lambda x: after_pct.get(x, 0), reverse=True)[:top_n]

    # Create comparison dataframe
    comparison = pd.DataFrame({
        "Before (Stage I)": [before_pct.get(t, 0) for t in all_types],
        "After (Stage N)": [after_pct.get(t, 0) for t in all_types],
    }, index=all_types)

    # Create grouped bar chart
    fig, ax = plt.subplots(figsize=(14, 6))

    x = np.arange(len(all_types))
    width = 0.35

    bars1 = ax.bar(x - width/2, comparison["Before (Stage I)"], width,
                   label="Before (Stage I)", color="#3498db", alpha=0.8)
    bars2 = ax.bar(x + width/2, comparison["After (Stage N)"], width,
                   label="After (Stage N)", color="#27ae60", alpha=0.8)

    ax.set_xlabel("Cell Type")
    ax.set_ylabel("Proportion (%)")
    ax.set_title("Composition: Before vs After Consolidation")
    ax.set_xticks(x)
    ax.set_xticklabels(all_types, rotation=45, ha="right")
    ax.legend()

    # Add change annotations
    for i, ct in enumerate(all_types):
        before_val = comparison.loc[ct, "Before (Stage I)"]
        after_val = comparison.loc[ct, "After (Stage N)"]
        diff = after_val - before_val

        if abs(diff) > 0.5:  # Only annotate significant changes
            color = "#27ae60" if diff > 0 else "#e74c3c"
            symbol = "+" if diff > 0 else ""
            ax.annotate(
                f"{symbol}{diff:.1f}%",
                xy=(i + width/2, after_val),
                xytext=(0, 5),
                textcoords="offset points",
                ha="center",
                fontsize=7,
                color=color,
            )

    plt.tight_layout()
    return save_figure(fig, output_path, dpi=dpi)


def plot_composition_heatmap(
    adata,
    output_path: Union[str, Path],
    cell_type_col: str = "cell_type_phenocycler",
    groupby_col: str = "sample_id",
    region_col: str = "region",
    top_n: int = 20,
    dpi: int = 150,
) -> Path:
    """Plot composition heatmap with hierarchical clustering.

    Args:
        adata: AnnData with cell type annotations
        output_path: Path to save figure
        cell_type_col: Column with cell type labels
        groupby_col: Column to group by (sample_id, donor, region)
        region_col: Column with region labels for annotation
        top_n: Number of top cell types to show
        dpi: Figure resolution

    Returns:
        Path to saved figure
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from .style import set_publication_style, save_figure

    set_publication_style()

    if groupby_col not in adata.obs.columns:
        # Try fallbacks
        for fallback in ["sample_id", "donor", "region"]:
            if fallback in adata.obs.columns:
                groupby_col = fallback
                break
        else:
            logger.warning("No groupby column found for heatmap")
            return None

    # Compute composition matrix
    df = adata.obs[[cell_type_col, groupby_col]].copy()
    if region_col in adata.obs.columns:
        df[region_col] = adata.obs[region_col]
    composition = df.groupby([groupby_col, cell_type_col], observed=False).size().unstack(fill_value=0)
    composition_pct = composition.div(composition.sum(axis=1), axis=0) * 100

    # Get top cell types
    global_counts = adata.obs[cell_type_col].value_counts()
    top_types = global_counts.head(top_n).index.tolist()
    top_types = [t for t in top_types if t in composition_pct.columns]

    composition_pct = composition_pct[top_types]

    # Build region color mapping for row annotation
    row_colors = None
    region_palette = None
    if region_col in adata.obs.columns and groupby_col == "sample_id":
        # Get sample_id -> region mapping
        sample_region = adata.obs.groupby(groupby_col, observed=False)[region_col].first()

        # Define region color palette
        region_palette = {
            "isthmus": "#e74c3c",      # Red
            "ampulla": "#3498db",       # Blue
            "infundibulum": "#2ecc71",  # Green
            "fimbriae": "#9b59b6",      # Purple
        }
        # Add fallback for unknown regions
        unique_regions = sample_region.unique()
        for i, r in enumerate(unique_regions):
            if r not in region_palette:
                colors = ["#f39c12", "#1abc9c", "#e67e22", "#34495e"]
                region_palette[r] = colors[i % len(colors)]

        # Create row_colors Series aligned with composition_pct index
        row_colors = composition_pct.index.map(
            lambda x: region_palette.get(sample_region.get(x, "unknown"), "#bdc3c7")
        )
        row_colors = pd.Series(row_colors, index=composition_pct.index, name="Region")

    # Try seaborn clustermap
    try:
        import seaborn as sns

        # Compute figure size based on data
        n_samples = composition_pct.shape[0]
        n_types = composition_pct.shape[1]
        figsize = (max(12, n_types * 0.5), max(8, n_samples * 0.3))

        g = sns.clustermap(
            composition_pct,
            cmap="YlOrRd",
            figsize=figsize,
            xticklabels=True,
            yticklabels=True,
            linewidths=0.5,
            row_colors=row_colors,
            cbar_kws={"label": "Proportion (%)"},
        )
        g.ax_heatmap.set_xlabel("Cell Type")
        g.ax_heatmap.set_ylabel(groupby_col.replace("_", " ").title())
        g.fig.suptitle("Cell Type Composition Heatmap", y=1.02)

        # Add region legend if we have region colors
        if region_palette is not None and row_colors is not None:
            # Get unique regions actually in the data
            sample_region = adata.obs.groupby(groupby_col, observed=False)[region_col].first()
            present_regions = sorted(sample_region.unique())
            legend_patches = [
                Patch(facecolor=region_palette[r], label=r.title())
                for r in present_regions if r in region_palette
            ]
            # Position legend to the right of the heatmap
            g.ax_heatmap.legend(
                handles=legend_patches,
                title="Region",
                loc="upper left",
                bbox_to_anchor=(1.15, 1.0),
                frameon=True,
            )

        plt.tight_layout()
        g.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(g.fig)

        return Path(output_path)

    except ImportError:
        # Fallback to simple heatmap
        fig, ax = plt.subplots(figsize=(12, 8))

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

        return save_figure(fig, output_path, dpi=dpi)
