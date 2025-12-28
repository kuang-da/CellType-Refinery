"""Export functions for spatial analysis results.

Provides CSV/JSON export with provenance tracking.
"""

from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Union
import json
import logging

import pandas as pd

from .config import SpatialConfig
from .engine import SpatialResult, MultiColumnSpatialResult

logger = logging.getLogger(__name__)


def export_single_result(
    result: SpatialResult,
    output_dir: Path,
    config: SpatialConfig,
) -> None:
    """Export results from single-column analysis.

    Parameters
    ----------
    result : SpatialResult
        Analysis results
    output_dir : Path
        Output directory
    config : SpatialConfig
        Configuration used
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    enrichment_dir = output_dir / "neighborhood_enrichment"
    interaction_dir = output_dir / "interaction_scores"
    heterogeneity_dir = output_dir / "heterogeneity"
    figures_dir = output_dir / "figures"

    for d in [enrichment_dir, interaction_dir, heterogeneity_dir, figures_dir]:
        d.mkdir(exist_ok=True)

    if result.neighborhood_enrichment:
        ne = result.neighborhood_enrichment
        ne.z_scores.to_csv(enrichment_dir / "enrichment_matrix.csv")
        ne.p_values_raw.to_csv(enrichment_dir / "enrichment_pvalues_raw.csv")
        ne.p_values_adj.to_csv(enrichment_dir / "enrichment_pvalues_adj.csv")
        ne.enrichment_pairs.to_csv(enrichment_dir / "enrichment_pairs.csv", index=False)
        ne.sample_stats.to_csv(enrichment_dir / "sample_stats.csv", index=False)
        logger.info(f"Exported neighborhood enrichment to {enrichment_dir}")

    if result.interaction_scores:
        inter = result.interaction_scores
        inter.interaction_matrix.to_csv(interaction_dir / "interaction_matrix.csv")
        inter.significant_interactions.to_csv(
            interaction_dir / "significant_interactions.csv", index=False
        )
        logger.info(f"Exported interaction scores to {interaction_dir}")

    if result.morans:
        morans = result.morans
        morans.morans_by_sample.to_csv(
            heterogeneity_dir / "morans_i_by_sample.csv", index=False
        )
        morans.morans_global.to_csv(
            heterogeneity_dir / "morans_i_global.csv", index=False
        )
        logger.info(f"Exported Moran's I to {heterogeneity_dir}")

    if result.local_diversity:
        div = result.local_diversity
        div.diversity_by_sample.to_csv(
            heterogeneity_dir / "local_diversity.csv", index=False
        )
        if not div.diversity_by_region.empty:
            div.diversity_by_region.to_csv(
                heterogeneity_dir / "local_diversity_by_region.csv", index=False
            )
        logger.info(f"Exported local diversity to {heterogeneity_dir}")

    summary = result.summary_dict()
    with open(output_dir / "spatial_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    logger.info(f"Exported summary to {output_dir / 'spatial_summary.json'}")


def export_multi_result(
    result: MultiColumnSpatialResult,
    output_dir: Path,
    config: SpatialConfig,
) -> None:
    """Export results from multi-column analysis.

    Parameters
    ----------
    result : MultiColumnSpatialResult
        Analysis results
    output_dir : Path
        Output directory
    config : SpatialConfig
        Configuration used
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for col, col_result in result.results.items():
        col_dir_name = col.replace("/", "_").replace(" ", "_")
        col_dir = output_dir / col_dir_name

        export_single_result(col_result, col_dir, config)

    summary = result.summary_dict()
    with open(output_dir / "spatial_summary_multi.json", "w") as f:
        json.dump(summary, f, indent=2)

    logger.info(f"Exported multi-column summary to {output_dir}")


def create_provenance(
    config: SpatialConfig,
    input_path: Path,
    graphs_dir: Path,
    output_dir: Path,
    result: Union[SpatialResult, MultiColumnSpatialResult],
    cell_type_cols: List[str],
) -> Dict[str, Any]:
    """Create provenance record for audit trail.

    Parameters
    ----------
    config : SpatialConfig
        Configuration used
    input_path : Path
        Path to input AnnData
    graphs_dir : Path
        Path to graphs directory
    output_dir : Path
        Path to output directory
    result : SpatialResult or MultiColumnSpatialResult
        Analysis result
    cell_type_cols : List[str]
        Cell type columns analyzed

    Returns
    -------
    Dict[str, Any]
        Provenance record
    """
    is_multi = isinstance(result, MultiColumnSpatialResult)

    provenance = {
        "timestamp": datetime.now().isoformat(),
        "module": "celltype_refinery.core.spatial",
        "version": "1.0.0",
        "inputs": {
            "adata_path": str(input_path),
            "graphs_dir": str(graphs_dir),
        },
        "outputs": {
            "output_dir": str(output_dir),
        },
        "config": config.to_dict(),
        "cell_type_cols": cell_type_cols,
        "execution": {
            "multi_column_mode": is_multi,
            "success": result.success,
        },
    }

    if is_multi:
        provenance["execution"]["columns_processed"] = result.columns_processed
        provenance["execution"]["columns_skipped"] = result.columns_skipped
        provenance["execution"]["total_time_seconds"] = round(
            result.total_execution_time_seconds, 2
        )
    else:
        provenance["execution"]["n_samples"] = result.n_samples
        provenance["execution"]["n_cells"] = result.n_cells_total
        provenance["execution"]["time_seconds"] = round(
            result.execution_time_seconds, 2
        )

    return provenance


def export_provenance(
    provenance: Dict[str, Any],
    output_dir: Path,
) -> None:
    """Export provenance to JSON file.

    Parameters
    ----------
    provenance : Dict
        Provenance record
    output_dir : Path
        Output directory
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    with open(output_dir / "provenance.json", "w") as f:
        json.dump(provenance, f, indent=2)

    logger.info(f"Exported provenance to {output_dir / 'provenance.json'}")


def export_all(
    result: Union[SpatialResult, MultiColumnSpatialResult],
    output_dir: Path,
    config: SpatialConfig,
    input_path: Path,
    graphs_dir: Path,
    cell_type_cols: List[str],
) -> None:
    """Export all results and provenance.

    Parameters
    ----------
    result : SpatialResult or MultiColumnSpatialResult
        Analysis result
    output_dir : Path
        Output directory
    config : SpatialConfig
        Configuration used
    input_path : Path
        Input AnnData path
    graphs_dir : Path
        Graphs directory path
    cell_type_cols : List[str]
        Cell type columns analyzed
    """
    output_dir = Path(output_dir)

    if isinstance(result, MultiColumnSpatialResult):
        export_multi_result(result, output_dir, config)
    else:
        export_single_result(result, output_dir, config)

    provenance = create_provenance(
        config=config,
        input_path=input_path,
        graphs_dir=graphs_dir,
        output_dir=output_dir,
        result=result,
        cell_type_cols=cell_type_cols,
    )
    export_provenance(provenance, output_dir)
