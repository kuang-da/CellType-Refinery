"""Style utilities and color schemes for consolidation visualizations.

This module provides consistent styling across all visualizations:
- Color palettes for cell types, categories, confidence bands
- Matplotlib style configuration
- Figure saving utilities
"""

from pathlib import Path
from typing import Dict, List, Optional, Union
import logging

logger = logging.getLogger(__name__)

# =============================================================================
# Color Palettes
# =============================================================================

# Category colors (root lineages)
CATEGORY_COLORS: Dict[str, str] = {
    "Epithelium": "#f1c40f",           # Yellow
    "Immune Cells": "#9b59b6",          # Purple
    "Endothelium": "#e74c3c",           # Red
    "Mesenchymal Cells": "#3498db",     # Blue
    "Misc": "#95a5a6",                  # Gray
    "Unassigned": "#7f8c8d",            # Dark gray
    "Hybrid": "#1abc9c",                # Teal
}

# Confidence band colors (traffic light scheme)
CONFIDENCE_COLORS: Dict[str, str] = {
    "high": "#27ae60",          # Green
    "medium": "#f39c12",        # Orange
    "low": "#e67e22",           # Dark orange
    "very_low": "#c0392b",      # Red
}

# Rescue status colors
RESCUE_STATUS_COLORS: Dict[str, str] = {
    "original": "#3498db",      # Blue - unchanged
    "rescued": "#27ae60",       # Green - rescued orphan
    "overridden": "#9b59b6",    # Purple - manual override
    "relabeled": "#f39c12",     # Orange - relabeled
    "unassigned": "#7f8c8d",    # Gray - still unassigned
}

# Label category colors
LABEL_CATEGORY_COLORS: Dict[str, str] = {
    "subtype": "#27ae60",       # Green - specific subtype
    "root": "#3498db",          # Blue - root level
    "hybrid": "#1abc9c",        # Teal - hybrid X~Y
    "unassigned": "#7f8c8d",    # Gray - unassigned
    "rescued": "#f39c12",       # Orange - rescued orphan
}

# Extended cell type colors (for specific subtypes)
CELL_TYPE_COLORS: Dict[str, str] = {
    # Epithelium subtypes
    "Ciliated Epithelium": "#f4d03f",
    "Glandular Epithelium": "#f5b041",
    "Secretory Epithelium": "#eb984e",
    "Epithelium": "#f1c40f",

    # Endothelium subtypes
    "Vascular Endothelium": "#e74c3c",
    "Lymphatic Endothelium": "#ec7063",
    "Endothelium": "#cb4335",

    # Mesenchymal subtypes
    "Smooth Muscle Cells": "#5dade2",
    "Fibroblasts": "#85c1e9",
    "Pericytes": "#aed6f1",
    "Mesenchymal Cells": "#3498db",

    # Immune subtypes
    "T Cells": "#a569bd",
    "B Cells": "#bb8fce",
    "Macrophages": "#d7bde2",
    "Myeloids": "#8e44ad",
    "Immune Cells": "#9b59b6",

    # Misc
    "Proliferating Cells": "#58d68d",
    "Misc": "#95a5a6",

    # Special
    "Unassigned": "#7f8c8d",
}

# Orphan action colors
ORPHAN_ACTION_COLORS: Dict[str, str] = {
    "rescue": "#27ae60",           # Green
    "keep_unassigned": "#e74c3c",  # Red
    "expert_review": "#f39c12",    # Orange
}

# Orphan plausibility colors
ORPHAN_PLAUSIBILITY_COLORS: Dict[str, str] = {
    "high": "#27ae60",         # Green
    "medium": "#f39c12",       # Orange
    "suspicious": "#e74c3c",   # Red
    "uncertain": "#f39c12",    # Orange
}

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


# =============================================================================
# Style Configuration
# =============================================================================

def set_publication_style():
    """Set matplotlib style for publication-quality figures."""
    import matplotlib.pyplot as plt

    # Try seaborn style first
    try:
        plt.style.use("seaborn-v0_8-whitegrid")
    except OSError:
        try:
            plt.style.use("seaborn-whitegrid")
        except OSError:
            pass

    # Custom rcParams for publication
    plt.rcParams.update({
        "font.size": 10,
        "axes.titlesize": 12,
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "figure.titlesize": 14,
        "figure.dpi": 100,
        "savefig.dpi": 200,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.1,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": True,
        "legend.framealpha": 0.9,
    })


def get_color_palette(
    labels: List[str],
    palette_type: str = "cell_type",
    default_color: str = "#bdc3c7",
) -> Dict[str, str]:
    """Get color mapping for a list of labels.

    Args:
        labels: List of label strings
        palette_type: One of "cell_type", "category", "confidence", "rescue"
        default_color: Color for labels not in predefined palette

    Returns:
        Dict mapping labels to hex colors
    """
    # Select base palette
    if palette_type == "cell_type":
        base_palette = CELL_TYPE_COLORS
    elif palette_type == "category":
        base_palette = CATEGORY_COLORS
    elif palette_type == "confidence":
        base_palette = CONFIDENCE_COLORS
    elif palette_type == "rescue":
        base_palette = RESCUE_STATUS_COLORS
    elif palette_type == "label_category":
        base_palette = LABEL_CATEGORY_COLORS
    elif palette_type == "orphan_action":
        base_palette = ORPHAN_ACTION_COLORS
    else:
        base_palette = {}

    # Build color mapping
    colors = {}
    fallback_colors = _get_fallback_palette(len(labels))
    fallback_idx = 0

    for label in labels:
        if label in base_palette:
            colors[label] = base_palette[label]
        else:
            # Try to match category
            matched = False
            for cat_name, cat_color in CATEGORY_COLORS.items():
                if cat_name.lower() in label.lower():
                    colors[label] = cat_color
                    matched = True
                    break

            if not matched:
                # Use fallback color
                if fallback_idx < len(fallback_colors):
                    colors[label] = fallback_colors[fallback_idx]
                    fallback_idx += 1
                else:
                    colors[label] = default_color

    return colors


def _get_fallback_palette(n: int) -> List[str]:
    """Get a fallback color palette for n items."""
    import matplotlib.pyplot as plt

    if n <= 10:
        cmap = plt.cm.get_cmap("tab10")
    elif n <= 20:
        cmap = plt.cm.get_cmap("tab20")
    else:
        cmap = plt.cm.get_cmap("viridis")

    colors = []
    for i in range(n):
        rgba = cmap(i / max(n - 1, 1))
        hex_color = "#{:02x}{:02x}{:02x}".format(
            int(rgba[0] * 255),
            int(rgba[1] * 255),
            int(rgba[2] * 255),
        )
        colors.append(hex_color)

    return colors


# =============================================================================
# Figure Utilities
# =============================================================================

def save_figure(
    fig,
    output_path: Union[str, Path],
    dpi: int = 200,
    close: bool = True,
) -> Path:
    """Save matplotlib figure with consistent settings.

    Args:
        fig: Matplotlib figure object
        output_path: Path to save figure
        dpi: Resolution
        close: Whether to close figure after saving

    Returns:
        Path to saved figure
    """
    import matplotlib.pyplot as plt

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(
        output_path,
        dpi=dpi,
        bbox_inches="tight",
        facecolor="white",
        edgecolor="none",
    )

    if close:
        plt.close(fig)

    logger.debug(f"Saved figure to {output_path}")
    return output_path


def create_figure(
    nrows: int = 1,
    ncols: int = 1,
    figsize: Optional[tuple] = None,
    **kwargs,
):
    """Create matplotlib figure with consistent styling.

    Args:
        nrows: Number of subplot rows
        ncols: Number of subplot columns
        figsize: Figure size (width, height) in inches
        **kwargs: Additional arguments to plt.subplots

    Returns:
        Tuple of (figure, axes)
    """
    import matplotlib.pyplot as plt

    set_publication_style()

    if figsize is None:
        # Auto-size based on subplots
        figsize = (4 * ncols + 1, 3 * nrows + 0.5)

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, **kwargs)
    return fig, axes


def add_cell_count_legend(
    ax,
    labels: List[str],
    counts: Dict[str, int],
    colors: Dict[str, str],
    loc: str = "upper right",
    ncol: int = 1,
    fontsize: int = 8,
):
    """Add legend with cell counts to axis.

    Args:
        ax: Matplotlib axis
        labels: List of labels in display order
        counts: Dict mapping label to cell count
        colors: Dict mapping label to color
        loc: Legend location
        ncol: Number of columns
        fontsize: Font size for legend
    """
    from matplotlib.patches import Patch

    handles = []
    legend_labels = []

    for label in labels:
        count = counts.get(label, 0)
        color = colors.get(label, "#bdc3c7")

        handles.append(Patch(facecolor=color, edgecolor="none"))
        legend_labels.append(f"{label} ({count:,})")

    ax.legend(
        handles,
        legend_labels,
        loc=loc,
        ncol=ncol,
        fontsize=fontsize,
        frameon=True,
        framealpha=0.9,
    )


def truncate_labels(labels: List[str], max_length: int = 25) -> List[str]:
    """Truncate long labels for display.

    Args:
        labels: List of label strings
        max_length: Maximum length before truncation

    Returns:
        List of truncated labels
    """
    return [
        label[:max_length-3] + "..." if len(label) > max_length else label
        for label in labels
    ]
