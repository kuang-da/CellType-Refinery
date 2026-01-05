"""Visualization module for Stage N consolidation results.

This module provides comprehensive visualizations for cell-type consolidation:

Static Figures (PNG):
- Composition charts (global, by region, by donor, before/after)
- Orphan rescue visualizations (bar chart, pie chart, scatter)
- UMAP plots (by cell type, confidence, rescue status)
- Spatial tissue maps
- Sankey label flow diagrams

Interactive (HTML):
- Consolidation dashboard
- Interactive UMAP explorer
- Sankey with hover/filter
- Spatial tissue explorer

Usage:
    from celltype_refinery.core.consolidation.viz import generate_all_figures

    generate_all_figures(
        adata=consolidated_adata,
        mapping_df=mapping_table,
        orphan_df=orphan_candidates,
        output_dir=Path("out/stage_n/figures"),
        compute_umap=True,
    )
"""

from .style import (
    CELL_TYPE_COLORS,
    CATEGORY_COLORS,
    CONFIDENCE_COLORS,
    RESCUE_STATUS_COLORS,
    IEL_TYPE_COLORS,
    VETO_MARKER_COLORS,
    set_publication_style,
    get_color_palette,
    save_figure,
)

from .composition import (
    plot_composition_global,
    plot_composition_by_region,
    plot_composition_by_donor,
    plot_composition_before_after,
    plot_composition_heatmap,
)

from .orphan_viz import (
    plot_orphan_rescue_summary,
    plot_orphan_decision_pie,
    plot_orphan_score_scatter,
)

from .iel_viz import (
    plot_iel_rescue_summary,
    plot_iel_score_scatter,
    plot_iel_marker_heatmap,
    plot_spatial_iel_highlight,
    summarize_iel_statistics,
)

from .umap_plots import (
    compute_umap_if_missing,
    plot_umap_by_cell_type,
    plot_umap_by_confidence,
    plot_umap_by_rescue_status,
    plot_umap_by_category,
)

from .spatial import (
    plot_spatial_cell_type,
    plot_spatial_orphan_highlight,
    plot_spatial_by_region,
)

from .sankey import (
    plot_sankey_label_flow,
    plot_sankey_interactive,
)

# NOTE: dashboard module removed - better visualization in composition analysis

__all__ = [
    # Style utilities
    "CELL_TYPE_COLORS",
    "CATEGORY_COLORS",
    "CONFIDENCE_COLORS",
    "RESCUE_STATUS_COLORS",
    "IEL_TYPE_COLORS",
    "VETO_MARKER_COLORS",
    "set_publication_style",
    "get_color_palette",
    "save_figure",
    # Composition
    "plot_composition_global",
    "plot_composition_by_region",
    "plot_composition_by_donor",
    "plot_composition_before_after",
    "plot_composition_heatmap",
    # Orphan
    "plot_orphan_rescue_summary",
    "plot_orphan_decision_pie",
    "plot_orphan_score_scatter",
    # IEL
    "plot_iel_rescue_summary",
    "plot_iel_score_scatter",
    "plot_iel_marker_heatmap",
    "plot_spatial_iel_highlight",
    "summarize_iel_statistics",
    # UMAP
    "compute_umap_if_missing",
    "plot_umap_by_cell_type",
    "plot_umap_by_confidence",
    "plot_umap_by_rescue_status",
    "plot_umap_by_category",
    # Spatial
    "plot_spatial_cell_type",
    "plot_spatial_orphan_highlight",
    "plot_spatial_by_region",
    # Sankey
    "plot_sankey_label_flow",
    "plot_sankey_interactive",
]


def generate_all_figures(
    adata,
    mapping_df,
    orphan_df=None,
    iel_df=None,
    marker_scores=None,
    diagnostic_df=None,
    output_dir=None,
    compute_umap: bool = True,
    use_gpu: bool = True,
    static_only: bool = False,
    interactive_only: bool = False,
    dpi: int = 200,
):
    """Generate all visualization outputs.

    Args:
        adata: AnnData with cell_type_phenocycler column
        mapping_df: DataFrame with cluster_id, original_label, final_label, n_cells
        orphan_df: DataFrame with orphan candidates (optional)
        iel_df: DataFrame with IEL candidates (optional)
            Expected columns: cluster_id, n_cells, cd45_pos_frac, veto_marker,
                            veto_marker_pos_frac, lymphoid_score, myeloid_score,
                            iel_type, final_label
        marker_scores: DataFrame with marker scores for all clusters (optional)
            Used to compute scores for non-orphan Unassigned clusters in scatter plot.
        diagnostic_df: DataFrame with diagnostic report (optional)
            Used to identify Unassigned clusters for scatter plot.
        output_dir: Output directory (Path or str)
        compute_umap: Whether to compute UMAP if missing
        use_gpu: Use GPU acceleration for UMAP (default: True)
        static_only: Only generate static PNG figures
        interactive_only: Only generate interactive HTML
        dpi: DPI for static figures

    Returns:
        Dict mapping figure names to output paths
    """
    from pathlib import Path
    import logging

    logger = logging.getLogger(__name__)

    if output_dir is None:
        raise ValueError("output_dir is required")

    output_dir = Path(output_dir)
    figures_dir = output_dir / "figures"
    interactive_dir = output_dir / "interactive"

    figures_dir.mkdir(parents=True, exist_ok=True)
    interactive_dir.mkdir(parents=True, exist_ok=True)

    outputs = {}

    # Compute UMAP if needed
    if compute_umap:
        logger.info("Computing UMAP if missing (use_gpu=%s)...", use_gpu)
        compute_umap_if_missing(adata, use_gpu=use_gpu)

    # Static figures
    if not interactive_only:
        logger.info("Generating static figures...")

        # Composition - only before/after comparison (others moved to Composition module)
        # Note: composition_global, composition_by_region, composition_by_donor,
        # and composition_heatmap are now generated by composition module
        # to avoid redundancy. Run composition analysis after Stage N for these plots.
        outputs["composition_before_after"] = plot_composition_before_after(
            adata, mapping_df, figures_dir / "composition_before_after.png", dpi=dpi
        )

        # Orphan visualizations
        if orphan_df is not None and len(orphan_df) > 0:
            outputs["orphan_rescue_summary"] = plot_orphan_rescue_summary(
                orphan_df, figures_dir / "orphan_rescue_summary.png", dpi=dpi
            )
            outputs["orphan_decision_pie"] = plot_orphan_decision_pie(
                orphan_df, figures_dir / "orphan_decision_pie.png", dpi=dpi
            )

            # Compute all unassigned cluster scores for complete scatter plot
            all_unassigned_scores = None
            if marker_scores is not None and diagnostic_df is not None:
                try:
                    from celltype_refinery.core.consolidation.orphan_detection import compute_all_unassigned_scores

                    # Get all Unassigned cluster IDs
                    unassigned_mask = diagnostic_df["assigned_label"].str.contains("Unassigned", na=False)
                    unassigned_ids = diagnostic_df.loc[unassigned_mask, "cluster_id"].astype(str).tolist()

                    if len(unassigned_ids) > 0:
                        all_unassigned_scores = compute_all_unassigned_scores(
                            marker_scores, unassigned_ids
                        )
                        logger.info(
                            f"Computed scores for {len(all_unassigned_scores)} Unassigned clusters "
                            f"({len(orphan_df)} are orphan candidates)"
                        )
                except ImportError:
                    logger.warning("Could not import compute_all_unassigned_scores, scatter plot will be limited")

            outputs["orphan_score_scatter"] = plot_orphan_score_scatter(
                orphan_df,
                figures_dir / "orphan_score_scatter.png",
                all_unassigned_scores=all_unassigned_scores,
                dpi=dpi,
            )

        # IEL visualizations
        if iel_df is not None and len(iel_df) > 0:
            outputs["iel_rescue_summary"] = plot_iel_rescue_summary(
                iel_df, figures_dir / "iel_rescue_summary.png", dpi=dpi
            )
            outputs["iel_score_scatter"] = plot_iel_score_scatter(
                iel_df, figures_dir / "iel_score_scatter.png", dpi=dpi
            )
            outputs["iel_marker_heatmap"] = plot_iel_marker_heatmap(
                iel_df, figures_dir / "iel_marker_heatmap.png", dpi=dpi
            )

        # UMAP plots
        if "X_umap" in adata.obsm:
            outputs["umap_cell_type"] = plot_umap_by_cell_type(
                adata, figures_dir / "umap_cell_type.png", dpi=dpi
            )
            outputs["umap_confidence"] = plot_umap_by_confidence(
                adata, figures_dir / "umap_confidence.png", dpi=dpi
            )
            outputs["umap_rescue_status"] = plot_umap_by_rescue_status(
                adata, figures_dir / "umap_rescue_status.png", dpi=dpi
            )
            outputs["umap_category"] = plot_umap_by_category(
                adata, figures_dir / "umap_category.png", dpi=dpi
            )

        # Spatial plots
        if "x" in adata.obs.columns and "y" in adata.obs.columns:
            outputs["spatial_cell_type"] = plot_spatial_cell_type(
                adata, figures_dir / "spatial_cell_type.png", dpi=dpi
            )
            if orphan_df is not None:
                outputs["spatial_orphan_highlight"] = plot_spatial_orphan_highlight(
                    adata, orphan_df, figures_dir / "spatial_orphan_highlight.png", dpi=dpi
                )
            if iel_df is not None and len(iel_df) > 0:
                outputs["spatial_iel_highlight"] = plot_spatial_iel_highlight(
                    adata, iel_df, figures_dir / "spatial_iel_highlight.png", dpi=dpi
                )

        # Sankey
        outputs["sankey_label_flow"] = plot_sankey_label_flow(
            mapping_df, figures_dir / "sankey_label_flow.png", dpi=dpi
        )

    # Interactive figures
    if not static_only:
        logger.info("Generating interactive figures...")

        outputs["sankey_interactive"] = plot_sankey_interactive(
            mapping_df, interactive_dir / "sankey_interactive.html"
        )

        # NOTE: consolidation_dashboard removed - better visualization in composition analysis

    logger.info(f"Generated {len(outputs)} visualization outputs")
    return outputs
