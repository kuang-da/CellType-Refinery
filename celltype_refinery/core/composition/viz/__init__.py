"""Visualization module for composition analysis.

Provides functions to generate:
- Composition plots (global, by region, by donor, heatmap)
- Diversity plots (by region, by sample)
- Biology plots (multi-panel, tissue-configurable)
- Enrichment plots (heatmap)
- Interactive dashboard (requires plotly)
"""

from pathlib import Path
from typing import Dict, Optional
import logging

from ..config import CompositionConfig
from ..engine import CompositionResult, MultiColumnResult

logger = logging.getLogger(__name__)


def generate_all_figures(
    adata,
    result: CompositionResult,
    output_dir: Path,
    cell_type_col: str = "cell_type_phenocycler",
    config: Optional[CompositionConfig] = None,
) -> Dict[str, Optional[Path]]:
    """Generate all composition visualization figures.

    Parameters
    ----------
    adata : AnnData
        Annotated data
    result : CompositionResult
        Composition analysis result
    output_dir : Path
        Output directory for figures
    cell_type_col : str
        Column with cell type labels
    config : CompositionConfig, optional
        Configuration for visualization settings

    Returns
    -------
    Dict[str, Optional[Path]]
        Dictionary mapping figure names to file paths
    """
    from .composition_plots import (
        plot_composition_global,
        plot_composition_by_region,
        plot_composition_by_donor,
        plot_composition_heatmap,
    )
    from .diversity_plots import (
        plot_diversity_by_region,
        plot_diversity_by_sample,
    )
    from .biology_plots import (
        plot_biology_metrics,
    )
    from .enrichment_plots import (
        plot_enrichment_heatmap,
    )

    if config is None:
        config = CompositionConfig.default()

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    dpi = config.visualization.dpi
    top_n = config.visualization.top_n_cell_types

    outputs: Dict[str, Optional[Path]] = {}

    # Composition plots
    try:
        outputs["composition_global"] = plot_composition_global(
            adata,
            output_dir / "composition_global.png",
            cell_type_col=cell_type_col,
            top_n=top_n,
            dpi=dpi,
        )
    except Exception as e:
        logger.warning(f"Failed to generate composition_global: {e}")
        outputs["composition_global"] = None

    try:
        outputs["composition_by_region"] = plot_composition_by_region(
            adata,
            output_dir / "composition_by_region.png",
            cell_type_col=cell_type_col,
            top_n=top_n,
            dpi=dpi,
        )
    except Exception as e:
        logger.warning(f"Failed to generate composition_by_region: {e}")
        outputs["composition_by_region"] = None

    try:
        outputs["composition_by_donor"] = plot_composition_by_donor(
            adata,
            output_dir / "composition_by_donor.png",
            cell_type_col=cell_type_col,
            top_n=top_n,
            dpi=dpi,
        )
    except Exception as e:
        logger.warning(f"Failed to generate composition_by_donor: {e}")
        outputs["composition_by_donor"] = None

    try:
        outputs["composition_heatmap"] = plot_composition_heatmap(
            adata,
            output_dir / "composition_heatmap.png",
            cell_type_col=cell_type_col,
            top_n=top_n,
            dpi=dpi,
        )
    except Exception as e:
        logger.warning(f"Failed to generate composition_heatmap: {e}")
        outputs["composition_heatmap"] = None

    # Diversity plots
    if result.diversity_by_sample is not None and not result.diversity_by_sample.empty:
        try:
            outputs["diversity_by_region"] = plot_diversity_by_region(
                result.diversity_by_sample,
                output_dir / "diversity_by_region.png",
                dpi=dpi,
            )
        except Exception as e:
            logger.warning(f"Failed to generate diversity_by_region: {e}")
            outputs["diversity_by_region"] = None

        try:
            outputs["diversity_by_sample"] = plot_diversity_by_sample(
                result.diversity_by_sample,
                output_dir / "diversity_by_sample.png",
                dpi=dpi,
            )
        except Exception as e:
            logger.warning(f"Failed to generate diversity_by_sample: {e}")
            outputs["diversity_by_sample"] = None

    # Biology plots
    if result.biology_by_sample is not None and not result.biology_by_sample.empty:
        try:
            outputs["biology_metrics"] = plot_biology_metrics(
                result.biology_by_sample,
                output_dir / "biology_metrics.png",
                dpi=dpi,
            )
        except Exception as e:
            logger.warning(f"Failed to generate biology_metrics: {e}")
            outputs["biology_metrics"] = None

    # Enrichment plots
    if result.enrichment is not None and not result.enrichment.empty:
        try:
            outputs["enrichment_heatmap"] = plot_enrichment_heatmap(
                result.enrichment,
                output_dir / "enrichment_heatmap.png",
                dpi=dpi,
            )
        except Exception as e:
            logger.warning(f"Failed to generate enrichment_heatmap: {e}")
            outputs["enrichment_heatmap"] = None

    return outputs


def generate_dashboard(
    output_dir: Path,
    output_path: Optional[Path] = None,
    title: str = "Cell-Type Composition Analysis",
) -> Optional[Path]:
    """Generate interactive HTML dashboard from composition CSVs.

    This function reads the CSV outputs from composition analysis and
    generates a comprehensive interactive HTML dashboard using Plotly.

    Parameters
    ----------
    output_dir : Path
        Directory containing composition output CSVs
    output_path : Path, optional
        Path to save HTML. If None, saves to output_dir/dashboard.html
    title : str
        Dashboard title

    Returns
    -------
    Optional[Path]
        Path to saved HTML file, or None if plotly not available
    """
    from .dashboard import generate_composition_dashboard

    return generate_composition_dashboard(
        output_dir=output_dir,
        output_path=output_path,
        title=title,
    )


def generate_all_figures_multi(
    adata,
    multi_result: MultiColumnResult,
    output_dir: Path,
    config: Optional[CompositionConfig] = None,
) -> Dict[str, Dict[str, Optional[Path]]]:
    """Generate figures for all processed columns.

    Creates figures in subdirectories:
        output_dir/cell_type_phenocycler/figures/
        output_dir/cell_type_broad/figures/

    Parameters
    ----------
    adata : AnnData
        Annotated data
    multi_result : MultiColumnResult
        Multi-column analysis results
    output_dir : Path
        Output directory
    config : CompositionConfig, optional
        Visualization configuration

    Returns
    -------
    Dict[str, Dict[str, Optional[Path]]]
        Nested dict: {column: {figure_name: path}}
    """
    if config is None:
        config = CompositionConfig.default()

    output_dir = Path(output_dir)
    all_outputs: Dict[str, Dict[str, Optional[Path]]] = {}

    for col in multi_result.columns_analyzed:
        result = multi_result.results[col]
        fig_dir = output_dir / col / "figures"

        logger.info(f"Generating figures for column: {col} -> {fig_dir}")

        try:
            outputs = generate_all_figures(
                adata=adata,
                result=result,
                output_dir=fig_dir,
                cell_type_col=col,
                config=config,
            )
            all_outputs[col] = outputs
        except Exception as e:
            logger.warning(f"Figure generation failed for {col}: {e}")
            all_outputs[col] = {}

    logger.info(
        f"Multi-column figure generation complete: "
        f"{len(multi_result.columns_analyzed)} columns"
    )
    return all_outputs


__all__ = [
    "generate_all_figures",
    "generate_all_figures_multi",
    "generate_dashboard",
]
