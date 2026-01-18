"""Visualization submodule for spatial analysis.

Provides plotting functions for:
- Neighborhood enrichment heatmaps
- Cell-type interaction networks
- Moran's I comparisons
- Local diversity distributions
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, TYPE_CHECKING
import logging

if TYPE_CHECKING:
    from ..config import SpatialConfig
    from ..engine import SpatialResult

logger = logging.getLogger(__name__)


def generate_all_visualizations(
    result: "SpatialResult",
    output_dir: Path,
    config: "SpatialConfig",
    organ: Optional[str] = None,
) -> None:
    """Generate all visualizations for spatial analysis results.

    Parameters
    ----------
    result : SpatialResult
        Analysis results
    output_dir : Path
        Output directory (figures will be saved in output_dir/figures/)
    config : SpatialConfig
        Configuration with visualization settings
    organ : str, optional
        Organ name for region ordering (e.g., 'fallopian_tube', 'uterus')
    """
    output_dir = Path(output_dir)
    figures_dir = output_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    viz_config = config.visualization

    # Neighborhood enrichment heatmap
    if result.neighborhood_enrichment is not None:
        try:
            from .enrichment_heatmap import plot_enrichment_heatmap

            plot_enrichment_heatmap(
                z_scores=result.neighborhood_enrichment.z_scores,
                output_path=figures_dir / "neighborhood_enrichment_heatmap.png",
                figsize=viz_config.figsize_heatmap,
                vmin=viz_config.z_score_vmin,
                vmax=viz_config.z_score_vmax,
                cmap=viz_config.cmap_enrichment,
                dpi=viz_config.dpi,
            )
            logger.info("Generated enrichment heatmap")
        except Exception as e:
            logger.warning(f"Failed to generate enrichment heatmap: {e}")

    # Interaction network
    if result.interaction_scores is not None:
        try:
            from .interaction_network import plot_interaction_network

            plot_interaction_network(
                interaction_matrix=result.interaction_scores.interaction_matrix,
                output_path=figures_dir / "interaction_network.png",
                threshold=viz_config.interaction_threshold,
                figsize=viz_config.figsize_network,
                dpi=viz_config.dpi,
            )
            logger.info("Generated interaction network")
        except Exception as e:
            logger.warning(f"Failed to generate interaction network: {e}")

    # Moran's I comparison
    if result.morans is not None:
        try:
            from .morans_plot import plot_morans_comparison

            plot_morans_comparison(
                morans_global=result.morans.morans_global,
                output_path=figures_dir / "morans_i_comparison.png",
                figsize=viz_config.figsize_bar,
                dpi=viz_config.dpi,
            )
            logger.info("Generated Moran's I comparison")
        except Exception as e:
            logger.warning(f"Failed to generate Moran's I plot: {e}")

    # Local diversity distribution
    if result.local_diversity is not None:
        try:
            from .diversity_plot import plot_diversity_distribution

            plot_diversity_distribution(
                diversity_by_sample=result.local_diversity.diversity_by_sample,
                output_path=figures_dir / "local_diversity_distribution.png",
                figsize=viz_config.figsize_bar,
                dpi=viz_config.dpi,
                organ=organ,
            )
            logger.info("Generated diversity distribution")
        except Exception as e:
            logger.warning(f"Failed to generate diversity plot: {e}")


__all__ = ["generate_all_visualizations"]
