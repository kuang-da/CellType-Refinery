"""Spatial tissue map visualizations.

Provides:
- Tissue map colored by cell type
- Tissue map highlighting rescued orphans
- Regional comparison spatial plots
"""

from pathlib import Path
from typing import Dict, List, Optional, Union
import logging

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def plot_spatial_cell_type(
    adata,
    output_path: Union[str, Path],
    cell_type_col: str = "cell_type_phenocycler",
    x_col: str = "x",
    y_col: str = "y",
    sample_id: Optional[str] = None,
    top_n: int = 15,
    point_size: float = 0.3,
    alpha: float = 0.7,
    dpi: int = 200,
) -> Path:
    """Plot spatial tissue map colored by cell type.

    Args:
        adata: AnnData with spatial coordinates
        output_path: Path to save figure
        cell_type_col: Column with cell type labels
        x_col: Column with x coordinates
        y_col: Column with y coordinates
        sample_id: Optional sample to subset
        top_n: Number of top cell types to show
        point_size: Size of scatter points
        alpha: Transparency
        dpi: Figure resolution

    Returns:
        Path to saved figure
    """
    import matplotlib.pyplot as plt
    from .style import set_publication_style, save_figure, get_color_palette

    set_publication_style()

    # Check required columns
    if x_col not in adata.obs.columns or y_col not in adata.obs.columns:
        logger.warning(f"Spatial coordinates not found ({x_col}, {y_col})")
        return None

    # Subset to sample if specified
    if sample_id is not None and "sample_id" in adata.obs.columns:
        mask = adata.obs["sample_id"] == sample_id
        adata = adata[mask]

    x = adata.obs[x_col].values
    y = adata.obs[y_col].values
    cell_types = adata.obs[cell_type_col].values

    # Get top cell types
    counts = adata.obs[cell_type_col].value_counts()
    top_types = counts.head(top_n).index.tolist()
    colors = get_color_palette(top_types, palette_type="cell_type")

    # Create figure
    fig, ax = plt.subplots(figsize=(12, 10))

    # Plot "Other" first (background)
    other_mask = ~np.isin(cell_types, top_types)
    if other_mask.sum() > 0:
        ax.scatter(
            x[other_mask],
            y[other_mask],
            c="#d5d8dc",
            s=point_size,
            alpha=alpha * 0.3,
            rasterized=True,
        )

    # Plot each cell type
    for ct in reversed(top_types):
        mask = cell_types == ct
        if mask.sum() > 0:
            ax.scatter(
                x[mask],
                y[mask],
                c=colors.get(ct, "#bdc3c7"),
                s=point_size,
                alpha=alpha,
                label=f"{ct} ({mask.sum():,})",
                rasterized=True,
            )

    ax.set_xlabel("X (pixels)")
    ax.set_ylabel("Y (pixels)")
    title = f"Spatial: Cell Types (n={len(adata):,})"
    if sample_id:
        title += f" - {sample_id}"
    ax.set_title(title)

    # Invert y-axis (image coordinates)
    ax.invert_yaxis()

    ax.legend(
        loc="upper left",
        bbox_to_anchor=(1.02, 1),
        fontsize=8,
        markerscale=5,
    )

    ax.set_aspect("equal")
    plt.tight_layout()
    return save_figure(fig, output_path, dpi=dpi)


def plot_spatial_orphan_highlight(
    adata,
    orphan_df: pd.DataFrame,
    output_path: Union[str, Path],
    cell_type_col: str = "cell_type_phenocycler",
    cluster_col: str = "cluster_lvl1",
    x_col: str = "x",
    y_col: str = "y",
    point_size: float = 0.5,
    alpha: float = 0.8,
    dpi: int = 200,
) -> Path:
    """Plot spatial map highlighting rescued orphans.

    Args:
        adata: AnnData with spatial coordinates
        orphan_df: DataFrame with orphan candidates
        output_path: Path to save figure
        cell_type_col: Column with cell type labels
        cluster_col: Column with cluster IDs
        x_col: Column with x coordinates
        y_col: Column with y coordinates
        point_size: Size of scatter points
        alpha: Transparency
        dpi: Figure resolution

    Returns:
        Path to saved figure
    """
    import matplotlib.pyplot as plt
    from .style import set_publication_style, save_figure, get_color_palette

    set_publication_style()

    if x_col not in adata.obs.columns or y_col not in adata.obs.columns:
        logger.warning(f"Spatial coordinates not found ({x_col}, {y_col})")
        return None

    if orphan_df is None or len(orphan_df) == 0:
        logger.warning("No orphan data for spatial visualization")
        return None

    x = adata.obs[x_col].values
    y = adata.obs[y_col].values

    # Identify rescued orphan cells
    action_col = "action" if "action" in orphan_df.columns else "rescue_action"
    rescued_clusters = set()
    if action_col in orphan_df.columns:
        rescued = orphan_df[orphan_df[action_col] == "rescue"]
        if "cluster_id" in rescued.columns:
            rescued_clusters = set(rescued["cluster_id"].astype(str))

    # Create mask for rescued cells
    if cluster_col in adata.obs.columns and rescued_clusters:
        is_rescued = adata.obs[cluster_col].astype(str).isin(rescued_clusters)
    else:
        # Fallback: check for (orphan) suffix
        is_rescued = adata.obs[cell_type_col].str.contains(r"\(orphan\)", regex=True, na=False)

    # Create figure
    fig, ax = plt.subplots(figsize=(12, 10))

    # Plot non-rescued cells in gray
    non_rescued = ~is_rescued
    ax.scatter(
        x[non_rescued],
        y[non_rescued],
        c="#d5d8dc",
        s=point_size * 0.5,
        alpha=alpha * 0.3,
        label=f"Non-rescued ({non_rescued.sum():,})",
        rasterized=True,
    )

    # Get rescued cell types for coloring
    rescued_types = adata.obs.loc[is_rescued, cell_type_col].unique()
    colors = get_color_palette(list(rescued_types), palette_type="cell_type")

    # Plot rescued cells by type
    for ct in rescued_types:
        mask = is_rescued & (adata.obs[cell_type_col] == ct)
        if mask.sum() > 0:
            ax.scatter(
                x[mask],
                y[mask],
                c=colors.get(ct, "#27ae60"),
                s=point_size * 2,
                alpha=alpha,
                label=f"{ct} ({mask.sum():,})",
                rasterized=True,
            )

    ax.set_xlabel("X (pixels)")
    ax.set_ylabel("Y (pixels)")
    ax.set_title(f"Spatial: Rescued Orphans Highlighted (n={is_rescued.sum():,} rescued)")
    ax.invert_yaxis()

    ax.legend(
        loc="upper left",
        bbox_to_anchor=(1.02, 1),
        fontsize=8,
        markerscale=5,
    )

    ax.set_aspect("equal")
    plt.tight_layout()
    return save_figure(fig, output_path, dpi=dpi)


def plot_spatial_by_region(
    adata,
    output_path: Union[str, Path],
    cell_type_col: str = "cell_type_phenocycler",
    region_col: str = "region",
    x_col: str = "x",
    y_col: str = "y",
    top_n: int = 10,
    point_size: float = 0.2,
    alpha: float = 0.6,
    dpi: int = 200,
) -> Path:
    """Plot spatial maps faceted by region.

    Args:
        adata: AnnData with spatial coordinates and region labels
        output_path: Path to save figure
        cell_type_col: Column with cell type labels
        region_col: Column with region labels
        x_col: Column with x coordinates
        y_col: Column with y coordinates
        top_n: Number of top cell types per region
        point_size: Size of scatter points
        alpha: Transparency
        dpi: Figure resolution

    Returns:
        Path to saved figure
    """
    import matplotlib.pyplot as plt
    from .style import set_publication_style, save_figure, get_color_palette

    set_publication_style()

    if region_col not in adata.obs.columns:
        logger.warning(f"Region column '{region_col}' not found")
        return None

    if x_col not in adata.obs.columns or y_col not in adata.obs.columns:
        logger.warning(f"Spatial coordinates not found")
        return None

    regions = adata.obs[region_col].unique()
    n_regions = len(regions)

    if n_regions == 0:
        return None

    # Get global top cell types for consistent coloring
    global_counts = adata.obs[cell_type_col].value_counts()
    global_top = global_counts.head(top_n).index.tolist()
    colors = get_color_palette(global_top, palette_type="cell_type")

    # Create faceted figure
    n_cols = min(2, n_regions)
    n_rows = (n_regions + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(8 * n_cols, 6 * n_rows))

    if n_regions == 1:
        axes = [axes]
    else:
        axes = axes.flatten()

    for i, region in enumerate(sorted(regions)):
        ax = axes[i]
        mask = adata.obs[region_col] == region
        subset = adata[mask]

        x = subset.obs[x_col].values
        y = subset.obs[y_col].values
        cell_types = subset.obs[cell_type_col].values

        # Plot by cell type
        for ct in global_top:
            ct_mask = cell_types == ct
            if ct_mask.sum() > 0:
                ax.scatter(
                    x[ct_mask],
                    y[ct_mask],
                    c=colors.get(ct, "#bdc3c7"),
                    s=point_size,
                    alpha=alpha,
                    rasterized=True,
                )

        ax.set_title(f"{region} (n={len(subset):,})")
        ax.invert_yaxis()
        ax.set_aspect("equal")
        ax.set_xlabel("X")
        ax.set_ylabel("Y")

    # Hide unused axes
    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    # Add shared legend
    from matplotlib.patches import Patch
    legend_handles = [Patch(color=colors.get(ct, "#bdc3c7"), label=ct) for ct in global_top]
    fig.legend(
        handles=legend_handles,
        loc="center right",
        bbox_to_anchor=(1.15, 0.5),
        fontsize=8,
    )

    fig.suptitle("Spatial Distribution by Region", fontsize=14, y=1.02)
    plt.tight_layout()
    return save_figure(fig, output_path, dpi=dpi)
