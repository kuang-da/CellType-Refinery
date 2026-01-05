"""UMAP visualization functions for consolidation results.

Provides:
- UMAP computation (if missing)
- UMAP colored by cell type
- UMAP colored by confidence band
- UMAP colored by rescue status
- UMAP colored by category
"""

from pathlib import Path
from typing import Dict, Optional, Union
import logging

import numpy as np

logger = logging.getLogger(__name__)


def compute_umap_if_missing(
    adata,
    use_rep: str = "X_pca",
    n_neighbors: int = 15,
    min_dist: float = 0.5,
    random_state: int = 42,
    use_gpu: bool = True,
) -> bool:
    """Compute UMAP embedding if not present.

    Uses GPU acceleration via rapids-singlecell when available,
    falls back to CPU scanpy otherwise.

    Args:
        adata: AnnData object
        use_rep: Representation to use (X_pca recommended)
        n_neighbors: Number of neighbors for UMAP
        min_dist: Minimum distance for UMAP
        random_state: Random seed for reproducibility
        use_gpu: Try to use GPU acceleration (default: True)

    Returns:
        True if UMAP was computed, False if already present
    """
    if "X_umap" in adata.obsm:
        logger.info("UMAP already present in adata.obsm")
        return False

    n_cells = len(adata)
    logger.info(f"Computing UMAP for {n_cells:,} cells from {use_rep}...")

    # Try GPU first if requested
    if use_gpu:
        try:
            import rapids_singlecell as rsc
            import cupy as cp

            logger.info("Using GPU-accelerated UMAP (rapids-singlecell)...")

            # Ensure neighbors are computed (GPU)
            if "neighbors" not in adata.uns:
                logger.info("Computing neighbors graph on GPU...")
                rsc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=n_neighbors)

            # Compute UMAP on GPU
            rsc.tl.umap(adata, min_dist=min_dist, random_state=random_state)

            # Ensure result is on CPU (numpy array)
            if hasattr(adata.obsm["X_umap"], "get"):
                adata.obsm["X_umap"] = adata.obsm["X_umap"].get()

            logger.info(f"GPU UMAP computed: shape {adata.obsm['X_umap'].shape}")
            return True

        except ImportError:
            logger.warning("rapids-singlecell not available, falling back to CPU")
        except Exception as e:
            logger.warning(f"GPU UMAP failed: {e}, falling back to CPU")
            # Clear any partial GPU state
            try:
                import gc
                gc.collect()
            except:
                pass

    # Fallback to CPU scanpy
    try:
        import scanpy as sc

        logger.info("Using CPU UMAP (scanpy)...")
        logger.warning(f"CPU UMAP for {n_cells:,} cells may take 10-30 minutes")

        # Ensure neighbors are computed
        if "neighbors" not in adata.uns:
            logger.info("Computing neighbors graph on CPU...")
            sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=n_neighbors)

        # Compute UMAP
        sc.tl.umap(adata, min_dist=min_dist, random_state=random_state)

        logger.info(f"CPU UMAP computed: shape {adata.obsm['X_umap'].shape}")
        return True

    except ImportError:
        logger.error("scanpy not available for UMAP computation")
        return False
    except Exception as e:
        logger.error(f"UMAP computation failed: {e}")
        return False


def plot_umap_by_cell_type(
    adata,
    output_path: Union[str, Path],
    cell_type_col: str = "cell_type_phenocycler",
    top_n: int = 20,
    point_size: float = 0.5,
    alpha: float = 0.5,
    dpi: int = 200,
) -> Path:
    """Plot UMAP colored by cell type.

    Args:
        adata: AnnData with X_umap and cell type column
        output_path: Path to save figure
        cell_type_col: Column with cell type labels
        top_n: Number of top cell types to show individually
        point_size: Size of scatter points
        alpha: Transparency
        dpi: Figure resolution

    Returns:
        Path to saved figure
    """
    import matplotlib.pyplot as plt
    from .style import set_publication_style, save_figure, get_color_palette

    set_publication_style()

    if "X_umap" not in adata.obsm:
        logger.warning("X_umap not found, skipping UMAP plot")
        return None

    umap_coords = adata.obsm["X_umap"]
    cell_types = adata.obs[cell_type_col].values

    # Get top cell types
    counts = adata.obs[cell_type_col].value_counts()
    top_types = counts.head(top_n).index.tolist()

    # Get colors
    colors = get_color_palette(top_types, palette_type="cell_type")

    # Create figure
    fig, ax = plt.subplots(figsize=(12, 10))

    # Plot "Other" first (background)
    other_mask = ~np.isin(cell_types, top_types)
    if other_mask.sum() > 0:
        ax.scatter(
            umap_coords[other_mask, 0],
            umap_coords[other_mask, 1],
            c="#d5d8dc",
            s=point_size,
            alpha=alpha * 0.5,
            label=f"Other ({other_mask.sum():,})",
            rasterized=True,
        )

    # Plot each cell type
    for ct in reversed(top_types):  # Reverse so smallest are on top
        mask = cell_types == ct
        count = mask.sum()
        if count > 0:
            ax.scatter(
                umap_coords[mask, 0],
                umap_coords[mask, 1],
                c=colors.get(ct, "#bdc3c7"),
                s=point_size,
                alpha=alpha,
                label=f"{ct} ({count:,})",
                rasterized=True,
            )

    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title(f"UMAP: Cell Types (n={len(adata):,})")
    ax.legend(
        loc="upper left",
        bbox_to_anchor=(1.02, 1),
        fontsize=8,
        markerscale=5,
    )

    plt.tight_layout()
    return save_figure(fig, output_path, dpi=dpi)


def plot_umap_by_confidence(
    adata,
    output_path: Union[str, Path],
    confidence_col: str = "confidence_band",
    point_size: float = 0.5,
    alpha: float = 0.5,
    dpi: int = 200,
) -> Path:
    """Plot UMAP colored by confidence band.

    Args:
        adata: AnnData with X_umap and confidence column
        output_path: Path to save figure
        confidence_col: Column with confidence bands
        point_size: Size of scatter points
        alpha: Transparency
        dpi: Figure resolution

    Returns:
        Path to saved figure
    """
    import matplotlib.pyplot as plt
    from .style import set_publication_style, save_figure, CONFIDENCE_COLORS

    set_publication_style()

    if "X_umap" not in adata.obsm:
        logger.warning("X_umap not found, skipping UMAP plot")
        return None

    if confidence_col not in adata.obs.columns:
        logger.warning(f"Confidence column '{confidence_col}' not found")
        return None

    umap_coords = adata.obsm["X_umap"]
    confidence = adata.obs[confidence_col].values

    # Create figure
    fig, ax = plt.subplots(figsize=(12, 10))

    # Order: very_low first (background), high last (foreground)
    order = ["very_low", "low", "medium", "high"]

    for band in order:
        mask = confidence == band
        count = mask.sum()
        if count > 0:
            ax.scatter(
                umap_coords[mask, 0],
                umap_coords[mask, 1],
                c=CONFIDENCE_COLORS.get(band, "#bdc3c7"),
                s=point_size,
                alpha=alpha,
                label=f"{band} ({count:,})",
                rasterized=True,
            )

    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title(f"UMAP: Confidence Bands (n={len(adata):,})")
    ax.legend(
        loc="upper left",
        bbox_to_anchor=(1.02, 1),
        fontsize=10,
        markerscale=5,
    )

    plt.tight_layout()
    return save_figure(fig, output_path, dpi=dpi)


def plot_umap_by_rescue_status(
    adata,
    output_path: Union[str, Path],
    rescue_col: str = "rescue_status",
    point_size: float = 0.5,
    alpha: float = 0.5,
    dpi: int = 200,
) -> Path:
    """Plot UMAP colored by rescue status.

    Args:
        adata: AnnData with X_umap and rescue status column
        output_path: Path to save figure
        rescue_col: Column with rescue status
        point_size: Size of scatter points
        alpha: Transparency
        dpi: Figure resolution

    Returns:
        Path to saved figure
    """
    import matplotlib.pyplot as plt
    from .style import set_publication_style, save_figure, RESCUE_STATUS_COLORS

    set_publication_style()

    if "X_umap" not in adata.obsm:
        logger.warning("X_umap not found, skipping UMAP plot")
        return None

    # Try to infer rescue status if column doesn't exist
    if rescue_col not in adata.obs.columns:
        # Check if we can infer from cell_type_phenocycler
        if "cell_type_phenocycler" in adata.obs.columns:
            # Cells with "(orphan)" suffix are rescued
            adata.obs[rescue_col] = "original"
            orphan_mask = adata.obs["cell_type_phenocycler"].str.contains(r"\(orphan\)", regex=True, na=False)
            adata.obs.loc[orphan_mask, rescue_col] = "rescued"

            # Check for unassigned
            unassigned_mask = adata.obs["cell_type_phenocycler"] == "Unassigned"
            adata.obs.loc[unassigned_mask, rescue_col] = "unassigned"
        else:
            logger.warning(f"Rescue status column '{rescue_col}' not found")
            return None

    umap_coords = adata.obsm["X_umap"]
    status = adata.obs[rescue_col].values

    # Create figure
    fig, ax = plt.subplots(figsize=(12, 10))

    # Order: original first, rescued highlighted on top
    order = ["original", "unassigned", "overridden", "relabeled", "rescued"]

    for stat in order:
        mask = status == stat
        count = mask.sum()
        if count > 0:
            ax.scatter(
                umap_coords[mask, 0],
                umap_coords[mask, 1],
                c=RESCUE_STATUS_COLORS.get(stat, "#bdc3c7"),
                s=point_size if stat == "original" else point_size * 2,
                alpha=alpha if stat == "original" else alpha * 1.5,
                label=f"{stat.title()} ({count:,})",
                rasterized=True,
            )

    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title(f"UMAP: Rescue Status (n={len(adata):,})")
    ax.legend(
        loc="upper left",
        bbox_to_anchor=(1.02, 1),
        fontsize=10,
        markerscale=5,
    )

    plt.tight_layout()
    return save_figure(fig, output_path, dpi=dpi)


def plot_umap_by_category(
    adata,
    output_path: Union[str, Path],
    category_col: str = "label_category",
    cell_type_col: str = "cell_type_phenocycler",
    point_size: float = 0.5,
    alpha: float = 0.5,
    dpi: int = 200,
) -> Path:
    """Plot UMAP colored by label category (subtype, root, hybrid, etc).

    Args:
        adata: AnnData with X_umap
        output_path: Path to save figure
        category_col: Column with label categories
        cell_type_col: Column with cell type labels (used to infer category)
        point_size: Size of scatter points
        alpha: Transparency
        dpi: Figure resolution

    Returns:
        Path to saved figure
    """
    import matplotlib.pyplot as plt
    from .style import set_publication_style, save_figure, LABEL_CATEGORY_COLORS

    set_publication_style()

    if "X_umap" not in adata.obsm:
        logger.warning("X_umap not found, skipping UMAP plot")
        return None

    # Infer category if not present
    if category_col not in adata.obs.columns:
        if cell_type_col in adata.obs.columns:
            try:
                from celltype_refinery.core.consolidation.rules import classify_label, LabelCategory

                categories = []
                for label in adata.obs[cell_type_col]:
                    cat = classify_label(label)
                    # Check for rescued orphans
                    if "(orphan)" in str(label):
                        categories.append("rescued")
                    else:
                        categories.append(cat.value)
                adata.obs[category_col] = categories
            except ImportError:
                # Fallback to simple heuristic
                categories = []
                for label in adata.obs[cell_type_col]:
                    label_str = str(label)
                    if "(orphan)" in label_str:
                        categories.append("rescued")
                    elif label_str == "Unassigned":
                        categories.append("unassigned")
                    elif "~" in label_str:
                        categories.append("hybrid")
                    else:
                        categories.append("subtype")
                adata.obs[category_col] = categories
        else:
            logger.warning(f"Cannot infer category without {cell_type_col}")
            return None

    umap_coords = adata.obsm["X_umap"]
    category = adata.obs[category_col].values

    # Create figure
    fig, ax = plt.subplots(figsize=(12, 10))

    # Order for layering
    order = ["subtype", "root", "hybrid", "unassigned", "rescued"]

    for cat in order:
        mask = category == cat
        count = mask.sum()
        if count > 0:
            ax.scatter(
                umap_coords[mask, 0],
                umap_coords[mask, 1],
                c=LABEL_CATEGORY_COLORS.get(cat, "#bdc3c7"),
                s=point_size,
                alpha=alpha,
                label=f"{cat.title()} ({count:,})",
                rasterized=True,
            )

    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title(f"UMAP: Label Categories (n={len(adata):,})")
    ax.legend(
        loc="upper left",
        bbox_to_anchor=(1.02, 1),
        fontsize=10,
        markerscale=5,
    )

    plt.tight_layout()
    return save_figure(fig, output_path, dpi=dpi)
