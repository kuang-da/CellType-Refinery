"""IEL (Intraepithelial Immune) rescue visualization functions.

Provides:
- Rescue summary bar chart (rescued IELs by type)
- CD45 vs Veto marker score scatter plot
- Spatial IEL highlight map
"""

from pathlib import Path
from typing import Dict, List, Optional, Union
import logging

import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


# IEL type colors
IEL_TYPE_COLORS: Dict[str, str] = {
    "Intraepithelial Lymphocytes": "#9b59b6",   # Purple (lymphoid)
    "Tissue-Resident Macrophages": "#e74c3c",   # Red (myeloid)
    "Immune (intraepithelial)": "#3498db",      # Blue (generic)
}

# Veto marker colors
VETO_MARKER_COLORS: Dict[str, str] = {
    "Pan-Cytokeratin": "#f1c40f",  # Yellow (epithelial)
    "E-cadherin": "#e67e22",       # Orange (epithelial junction)
}


def plot_iel_rescue_summary(
    iel_df: pd.DataFrame,
    output_path: Union[str, Path],
    dpi: int = 200,
) -> Optional[Path]:
    """Plot horizontal bar chart of rescued IELs by type.

    Args:
        iel_df: DataFrame with IEL candidates
            Expected columns: iel_type, n_cells, final_label
        output_path: Path to save figure
        dpi: Figure resolution

    Returns:
        Path to saved figure or None if no data
    """
    import matplotlib.pyplot as plt
    from .style import set_publication_style, save_figure

    set_publication_style()

    if iel_df is None or len(iel_df) == 0:
        logger.warning("No IEL data to visualize")
        return None

    # Check for required columns
    if "iel_type" not in iel_df.columns:
        logger.warning("Column 'iel_type' not found in iel_df")
        return None

    # Aggregate by IEL type
    n_cells_col = "n_cells" if "n_cells" in iel_df.columns else None
    if n_cells_col:
        counts = iel_df.groupby("iel_type")[n_cells_col].sum().sort_values(ascending=True)
    else:
        counts = iel_df["iel_type"].value_counts().sort_values(ascending=True)

    if len(counts) == 0:
        logger.warning("No IEL types to visualize")
        return None

    # Create horizontal bar chart
    fig, ax = plt.subplots(figsize=(10, max(4, len(counts) * 0.8)))

    colors = [IEL_TYPE_COLORS.get(t, "#95a5a6") for t in counts.index]
    bars = ax.barh(range(len(counts)), counts.values, color=colors, edgecolor="white", linewidth=0.5)

    # Get cluster counts for secondary labels
    cluster_counts = iel_df.groupby("iel_type").size() if n_cells_col else None

    # Add value labels with cluster counts at end of bars
    for i, (bar, count, iel_type) in enumerate(zip(bars, counts.values, counts.index)):
        if count >= 1000:
            label = f"{count/1000:.1f}K"
        else:
            label = str(count)

        x_pos = bar.get_width() + max(counts) * 0.02
        y_center = bar.get_y() + bar.get_height() / 2

        # Cell count (main label) - offset up slightly
        ax.text(
            x_pos,
            y_center - 0.12,
            label,
            va="center",
            fontsize=10,
            fontweight="bold",
        )

        # Cluster count (secondary label) - offset down slightly
        if cluster_counts is not None:
            n_clusters = cluster_counts.get(iel_type, 0)
            ax.text(
                x_pos,
                y_center + 0.12,
                f"({n_clusters} clusters)",
                va="center",
                fontsize=8,
                color="#7f8c8d",
            )

    ax.set_yticks(range(len(counts)))
    ax.set_yticklabels(counts.index, fontsize=11)
    ax.set_xlabel("Number of Cells", fontsize=11)
    ax.set_title(f"IEL Rescue Summary\n(n={counts.sum():,} cells from {len(iel_df)} clusters)", fontsize=12)
    ax.set_xlim(0, max(counts) * 1.25)

    plt.tight_layout()
    return save_figure(fig, output_path, dpi=dpi)


def plot_iel_score_scatter(
    iel_df: pd.DataFrame,
    output_path: Union[str, Path],
    dpi: int = 200,
) -> Optional[Path]:
    """Plot scatter of CD45 pos_frac vs veto marker pos_frac for IEL candidates.

    Shows the immune vs epithelial signal balance that characterizes IELs.

    Args:
        iel_df: DataFrame with IEL candidates
            Expected columns: cd45_pos_frac, veto_marker_pos_frac, veto_marker,
                            iel_type, n_cells, cluster_id
        output_path: Path to save figure
        dpi: Figure resolution

    Returns:
        Path to saved figure or None if no data
    """
    import matplotlib.pyplot as plt
    from .style import set_publication_style, save_figure

    set_publication_style()

    if iel_df is None or len(iel_df) == 0:
        logger.warning("No IEL data to visualize")
        return None

    # Check required columns
    required = ["cd45_pos_frac", "veto_marker_pos_frac"]
    for col in required:
        if col not in iel_df.columns:
            logger.warning(f"Column '{col}' not found in iel_df")
            return None

    # Create scatter plot
    fig, ax = plt.subplots(figsize=(10, 8))

    # Plot by IEL type
    iel_types = iel_df["iel_type"].unique() if "iel_type" in iel_df.columns else ["IEL"]

    for iel_type in iel_types:
        if "iel_type" in iel_df.columns:
            mask = iel_df["iel_type"] == iel_type
            subset = iel_df[mask]
        else:
            subset = iel_df

        color = IEL_TYPE_COLORS.get(iel_type, "#95a5a6")

        # Size by n_cells if available
        if "n_cells" in subset.columns:
            sizes = np.clip(subset["n_cells"] / 20, 50, 500)
        else:
            sizes = 100

        ax.scatter(
            subset["cd45_pos_frac"],
            subset["veto_marker_pos_frac"],
            c=color,
            s=sizes,
            alpha=0.7,
            label=f"{iel_type} ({len(subset)})",
            edgecolors="white",
            linewidth=0.5,
        )

        # Add cluster labels for larger clusters
        if "cluster_id" in subset.columns and "n_cells" in subset.columns:
            for _, row in subset.iterrows():
                if row["n_cells"] > 300:
                    ax.annotate(
                        row["cluster_id"],
                        (row["cd45_pos_frac"], row["veto_marker_pos_frac"]),
                        fontsize=8,
                        alpha=0.8,
                        xytext=(4, 4),
                        textcoords="offset points",
                    )

    # Add threshold lines
    ax.axvline(x=0.30, color="#27ae60", linestyle="--", alpha=0.7, linewidth=1.5,
               label="CD45 threshold (0.30)")
    ax.axhline(y=0.20, color="#e74c3c", linestyle="--", alpha=0.7, linewidth=1.5,
               label="Veto threshold (0.20)")

    # Shade the IEL zone: CD45 >= 0.30 AND veto > 0.20
    ax.fill_between(
        [0.30, 1.05],
        [0.20, 0.20],
        [1.05, 1.05],
        alpha=0.1,
        color="#9b59b6",
        label="IEL zone",
    )

    # Add zone annotations
    ax.annotate(
        "IEL ZONE\n(CD45+ & Epithelial+)",
        xy=(0.65, 0.6),
        fontsize=11,
        ha="center",
        color="#9b59b6",
        fontweight="bold",
    )
    ax.annotate(
        "Normal Immune\n(no epithelial veto)",
        xy=(0.65, 0.10),
        fontsize=9,
        ha="center",
        color="#27ae60",
        alpha=0.7,
    )
    ax.annotate(
        "Low CD45\n(not rescued)",
        xy=(0.15, 0.6),
        fontsize=9,
        ha="center",
        color="#e74c3c",
        alpha=0.7,
    )

    ax.set_xlabel("CD45 Positive Fraction (Immune Signal)", fontsize=11)
    ax.set_ylabel("Veto Marker Positive Fraction (Epithelial Signal)", fontsize=11)
    ax.set_title(
        "IEL Detection: Immune vs Epithelial Signal\n"
        "(Cells with both CD45+ and epithelial marker+)",
        fontsize=12,
    )
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.legend(loc="lower right", fontsize=9)

    # Add veto marker breakdown if available
    if "veto_marker" in iel_df.columns:
        veto_counts = iel_df["veto_marker"].value_counts()
        veto_text = "\n".join([f"{m}: {c}" for m, c in veto_counts.items()])
        ax.text(
            0.02, 0.98,
            f"Veto Markers:\n{veto_text}",
            transform=ax.transAxes,
            fontsize=9,
            va="top",
            ha="left",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
        )

    plt.tight_layout()
    return save_figure(fig, output_path, dpi=dpi)


def plot_iel_marker_heatmap(
    iel_df: pd.DataFrame,
    output_path: Union[str, Path],
    dpi: int = 200,
    cluster_rows: bool = True,
) -> Optional[Path]:
    """Plot heatmap of marker scores for IEL candidates.

    Shows CD45, lymphoid, myeloid scores and veto marker for each cluster.

    Args:
        iel_df: DataFrame with IEL candidates
            Expected columns: cluster_id, cd45_pos_frac, lymphoid_score,
                            myeloid_score, veto_marker_pos_frac, iel_type
        output_path: Path to save figure
        dpi: Figure resolution
        cluster_rows: If True, sort rows by hierarchical clustering

    Returns:
        Path to saved figure or None if no data
    """
    import matplotlib.pyplot as plt
    from .style import set_publication_style, save_figure

    set_publication_style()

    if iel_df is None or len(iel_df) == 0:
        logger.warning("No IEL data to visualize")
        return None

    # Check required columns
    required = ["cluster_id", "cd45_pos_frac"]
    for col in required:
        if col not in iel_df.columns:
            logger.warning(f"Column '{col}' not found in iel_df")
            return None

    # Prepare data for heatmap
    df = iel_df.copy()

    # Select columns for heatmap
    heatmap_cols = []
    col_labels = []

    if "cd45_pos_frac" in df.columns:
        heatmap_cols.append("cd45_pos_frac")
        col_labels.append("CD45\npos_frac")

    if "lymphoid_score" in df.columns:
        heatmap_cols.append("lymphoid_score")
        col_labels.append("Lymphoid\nscore")

    if "myeloid_score" in df.columns:
        heatmap_cols.append("myeloid_score")
        col_labels.append("Myeloid\nscore")

    if "veto_marker_pos_frac" in df.columns:
        heatmap_cols.append("veto_marker_pos_frac")
        col_labels.append("Veto\npos_frac")

    if len(heatmap_cols) < 2:
        logger.warning("Not enough columns for heatmap")
        return None

    # Create heatmap matrix
    matrix = df[heatmap_cols].values

    # Apply hierarchical clustering to rows if requested
    if cluster_rows and len(df) > 2:
        try:
            from scipy.cluster.hierarchy import linkage, leaves_list
            from scipy.spatial.distance import pdist

            # Compute pairwise distances and hierarchical clustering
            if matrix.shape[0] > 1:
                dist_matrix = pdist(matrix, metric="euclidean")
                linkage_matrix = linkage(dist_matrix, method="ward")
                row_order = leaves_list(linkage_matrix)

                # Reorder dataframe and matrix
                df = df.iloc[row_order].reset_index(drop=True)
                matrix = matrix[row_order]
        except Exception as e:
            logger.warning(f"Hierarchical clustering failed: {e}, falling back to n_cells sort")
            df = df.sort_values("n_cells", ascending=False) if "n_cells" in df.columns else df
            matrix = df[heatmap_cols].values
    else:
        # Fall back to sorting by n_cells
        df = df.sort_values("n_cells", ascending=False) if "n_cells" in df.columns else df
        matrix = df[heatmap_cols].values

    # Create row labels
    row_labels = []
    for _, row in df.iterrows():
        label = str(row["cluster_id"])
        if "n_cells" in row:
            label += f" ({row['n_cells']})"
        row_labels.append(label)

    # Create figure with extra space for IEL type column on left and legend at bottom
    fig, (ax_type, ax_heat) = plt.subplots(
        1, 2,
        figsize=(10, max(5.5, len(df) * 0.55 + 2)),
        gridspec_kw={"width_ratios": [0.08, 1], "wspace": 0.02},
    )

    # Plot IEL type column on the left
    if "iel_type" in df.columns:
        type_colors = [IEL_TYPE_COLORS.get(t, "#95a5a6") for t in df["iel_type"]]
        type_matrix = np.arange(len(df)).reshape(-1, 1)  # Just for positioning

        for i, color in enumerate(type_colors):
            ax_type.add_patch(plt.Rectangle(
                (0, i - 0.5), 1, 1,
                facecolor=color,
                edgecolor="white",
                linewidth=0.5,
            ))

        ax_type.set_xlim(0, 1)
        ax_type.set_ylim(-0.5, len(df) - 0.5)
        ax_type.set_xticks([0.5])
        ax_type.set_xticklabels(["IEL\nType"], fontsize=9)
        ax_type.set_yticks([])
        ax_type.invert_yaxis()
        ax_type.spines["top"].set_visible(False)
        ax_type.spines["right"].set_visible(False)
        ax_type.spines["bottom"].set_visible(False)
        ax_type.spines["left"].set_visible(False)
    else:
        ax_type.axis("off")

    # Plot main heatmap using pcolormesh for seamless cells
    # pcolormesh needs coordinates for cell edges
    x_edges = np.arange(matrix.shape[1] + 1) - 0.5
    y_edges = np.arange(matrix.shape[0] + 1) - 0.5
    im = ax_heat.pcolormesh(
        x_edges, y_edges, matrix,
        cmap="RdYlGn", vmin=-2, vmax=2,
        edgecolors="face",  # No visible edges between cells
        linewidth=0,
    )

    # Set axis limits to match data
    ax_heat.set_xlim(-0.5, matrix.shape[1] - 0.5)
    ax_heat.set_ylim(matrix.shape[0] - 0.5, -0.5)  # Inverted for top-to-bottom

    # Remove ticks and spines
    ax_heat.tick_params(axis="both", which="both", length=0)
    for spine in ax_heat.spines.values():
        spine.set_visible(False)

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax_heat, shrink=0.8, pad=0.02)
    cbar.set_label("Score / Fraction", fontsize=10)

    # Set ticks
    ax_heat.set_xticks(range(len(col_labels)))
    ax_heat.set_xticklabels(col_labels, fontsize=10)
    ax_heat.set_yticks(range(len(row_labels)))
    ax_heat.set_yticklabels(row_labels, fontsize=9)

    # Add value annotations
    for i in range(len(df)):
        for j in range(len(heatmap_cols)):
            val = matrix[i, j]
            text_color = "white" if abs(val) > 1 else "black"
            ax_heat.text(j, i, f"{val:.2f}", ha="center", va="center",
                        fontsize=8, color=text_color)

    ax_heat.set_xlabel("Marker", fontsize=11)

    # Add title to figure
    fig.suptitle("IEL Candidate Marker Profile", fontsize=12, y=0.98)

    # Add legend at bottom
    if "iel_type" in df.columns:
        from matplotlib.patches import Patch
        unique_types = df["iel_type"].unique()
        # Sort to ensure consistent order
        type_order = ["Intraepithelial Lymphocytes", "Tissue-Resident Macrophages", "Immune (intraepithelial)"]
        sorted_types = [t for t in type_order if t in unique_types]

        legend_handles = [
            Patch(facecolor=IEL_TYPE_COLORS.get(t, "#95a5a6"), edgecolor="white", label=t)
            for t in sorted_types
        ]
        fig.legend(
            handles=legend_handles,
            loc="lower center",
            ncol=len(sorted_types),
            fontsize=9,
            frameon=True,
            bbox_to_anchor=(0.55, 0.01),
        )

    # Use subplots_adjust for better control over spacing
    plt.subplots_adjust(left=0.12, right=0.92, top=0.92, bottom=0.18)
    return save_figure(fig, output_path, dpi=dpi)


def plot_spatial_iel_highlight(
    adata,
    iel_df: pd.DataFrame,
    output_path: Union[str, Path],
    sample_id: Optional[str] = None,
    dpi: int = 200,
) -> Optional[Path]:
    """Plot spatial tissue map highlighting IEL cells.

    Args:
        adata: AnnData with spatial coordinates (x, y in obs)
        iel_df: DataFrame with IEL candidates
            Expected columns: cluster_id, iel_type
        output_path: Path to save figure
        sample_id: Optional sample ID to filter (plots first sample if None)
        dpi: Figure resolution

    Returns:
        Path to saved figure or None if no data
    """
    import matplotlib.pyplot as plt
    from .style import set_publication_style, save_figure

    set_publication_style()

    if adata is None or len(adata) == 0:
        logger.warning("No AnnData to visualize")
        return None

    if iel_df is None or len(iel_df) == 0:
        logger.warning("No IEL data to visualize")
        return None

    # Check for spatial coordinates
    if "x" not in adata.obs.columns or "y" not in adata.obs.columns:
        logger.warning("Spatial coordinates (x, y) not found in adata.obs")
        return None

    # Get cluster column
    cluster_col = "cluster_lvl1" if "cluster_lvl1" in adata.obs.columns else "cluster_lvl0"
    if cluster_col not in adata.obs.columns:
        logger.warning(f"Cluster column not found in adata.obs")
        return None

    # Filter to sample if specified
    if sample_id is not None and "sample_id" in adata.obs.columns:
        sample_mask = adata.obs["sample_id"] == sample_id
        if sample_mask.sum() == 0:
            logger.warning(f"Sample {sample_id} not found")
            return None
        adata_plot = adata[sample_mask]
    elif "sample_id" in adata.obs.columns:
        # Use first sample
        first_sample = adata.obs["sample_id"].iloc[0]
        sample_mask = adata.obs["sample_id"] == first_sample
        adata_plot = adata[sample_mask]
        sample_id = first_sample
    else:
        adata_plot = adata

    # Get IEL cluster IDs
    iel_cluster_ids = set(iel_df["cluster_id"].astype(str))

    # Create IEL type mapping
    iel_type_map = dict(zip(
        iel_df["cluster_id"].astype(str),
        iel_df["iel_type"] if "iel_type" in iel_df.columns else ["IEL"] * len(iel_df)
    ))

    # Assign colors
    cluster_ids = adata_plot.obs[cluster_col].astype(str)
    colors = []
    is_iel = []

    for cid in cluster_ids:
        if cid in iel_cluster_ids:
            iel_type = iel_type_map.get(cid, "IEL")
            colors.append(IEL_TYPE_COLORS.get(iel_type, "#9b59b6"))
            is_iel.append(True)
        else:
            colors.append("#e0e0e0")  # Light gray for non-IEL
            is_iel.append(False)

    is_iel = np.array(is_iel)

    # Create figure
    fig, ax = plt.subplots(figsize=(12, 10))

    # Plot non-IEL cells first (background)
    ax.scatter(
        adata_plot.obs.loc[~is_iel, "x"],
        adata_plot.obs.loc[~is_iel, "y"],
        c="#e0e0e0",
        s=1,
        alpha=0.3,
        rasterized=True,
    )

    # Plot IEL cells on top
    if is_iel.sum() > 0:
        iel_colors = [c for c, iel in zip(colors, is_iel) if iel]
        ax.scatter(
            adata_plot.obs.loc[is_iel, "x"],
            adata_plot.obs.loc[is_iel, "y"],
            c=iel_colors,
            s=8,
            alpha=0.8,
            edgecolors="white",
            linewidth=0.2,
            rasterized=True,
        )

    # Add legend
    from matplotlib.patches import Patch
    legend_handles = [
        Patch(facecolor="#e0e0e0", edgecolor="none", label="Other cells"),
    ]
    for iel_type, color in IEL_TYPE_COLORS.items():
        if iel_type in iel_type_map.values():
            n_cells = sum(1 for t in iel_type_map.values() if t == iel_type)
            legend_handles.append(
                Patch(facecolor=color, edgecolor="white", label=f"{iel_type}")
            )

    ax.legend(handles=legend_handles, loc="upper right", fontsize=9)

    # Set labels
    title = f"Spatial IEL Distribution"
    if sample_id:
        title += f"\n(Sample: {sample_id})"
    ax.set_title(title, fontsize=12)
    ax.set_xlabel("X coordinate (pixels)", fontsize=10)
    ax.set_ylabel("Y coordinate (pixels)", fontsize=10)

    # Equal aspect ratio
    ax.set_aspect("equal")

    # Invert y-axis for proper tissue orientation
    ax.invert_yaxis()

    plt.tight_layout()
    return save_figure(fig, output_path, dpi=dpi)


def summarize_iel_statistics(iel_df: pd.DataFrame) -> Dict:
    """Compute summary statistics for IEL rescue.

    Args:
        iel_df: DataFrame with IEL candidates

    Returns:
        Dict with summary statistics
    """
    if iel_df is None or len(iel_df) == 0:
        return {
            "n_clusters": 0,
            "n_cells": 0,
            "by_type": {},
            "by_veto_marker": {},
            "mean_cd45": None,
            "mean_lymphoid": None,
            "mean_myeloid": None,
        }

    stats = {
        "n_clusters": len(iel_df),
        "n_cells": int(iel_df["n_cells"].sum()) if "n_cells" in iel_df.columns else len(iel_df),
    }

    # By IEL type
    if "iel_type" in iel_df.columns and "n_cells" in iel_df.columns:
        by_type = iel_df.groupby("iel_type").agg({
            "cluster_id": "count",
            "n_cells": "sum",
        }).rename(columns={"cluster_id": "n_clusters"})
        stats["by_type"] = by_type.to_dict("index")
    else:
        stats["by_type"] = {}

    # By veto marker
    if "veto_marker" in iel_df.columns and "n_cells" in iel_df.columns:
        by_veto = iel_df.groupby("veto_marker").agg({
            "cluster_id": "count",
            "n_cells": "sum",
        }).rename(columns={"cluster_id": "n_clusters"})
        stats["by_veto_marker"] = by_veto.to_dict("index")
    else:
        stats["by_veto_marker"] = {}

    # Mean scores
    if "cd45_pos_frac" in iel_df.columns:
        stats["mean_cd45"] = float(iel_df["cd45_pos_frac"].mean())
    if "lymphoid_score" in iel_df.columns:
        stats["mean_lymphoid"] = float(iel_df["lymphoid_score"].mean())
    if "myeloid_score" in iel_df.columns:
        stats["mean_myeloid"] = float(iel_df["myeloid_score"].mean())

    return stats
