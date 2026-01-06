"""SpatialEngine - main orchestrator for spatial analysis.

Coordinates all spatial analysis components:
- Neighborhood enrichment (permutation-based)
- Cell-type interaction scores
- Moran's I spatial autocorrelation
- Local diversity indices
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional
import logging
import time

import numpy as np
import pandas as pd
from scipy import sparse

try:
    import scanpy as sc
except ImportError:
    sc = None

from .config import SpatialConfig
from .graph_io import load_spatial_graph, validate_graphs_dir
from .neighborhood import (
    NeighborhoodEnrichmentResult,
    compute_neighborhood_enrichment_sample,
    aggregate_enrichment_results,
)
from .interaction import (
    InteractionResult,
    compute_interactions_for_sample,
    aggregate_interaction_results,
)
from .morans import (
    MoransResult,
    compute_morans_for_sample,
    aggregate_morans_results,
)
from .diversity import (
    LocalDiversityResult,
    compute_diversity_for_sample,
    aggregate_diversity_results,
)
from .parallel import SampleInput

logger = logging.getLogger(__name__)


@dataclass
class SpatialResult:
    """Result from spatial analysis.

    Attributes
    ----------
    success : bool
        Whether analysis completed successfully
    n_cells_total : int
        Total cells analyzed
    n_samples : int
        Number of samples with valid graphs
    n_cell_types : int
        Number of distinct cell types
    cell_type_col : str
        Column used for cell type labels
    neighborhood_enrichment : NeighborhoodEnrichmentResult
        Neighborhood enrichment results
    interaction_scores : InteractionResult
        Cell-type interaction results
    morans : MoransResult
        Moran's I results
    local_diversity : LocalDiversityResult
        Local diversity results
    execution_time_seconds : float
        Total execution time
    warnings : List[str]
        Warnings encountered
    errors : List[str]
        Errors encountered
    """

    success: bool = True
    n_cells_total: int = 0
    n_samples: int = 0
    n_samples_skipped: int = 0
    n_cell_types: int = 0
    cell_type_col: str = ""

    neighborhood_enrichment: Optional[NeighborhoodEnrichmentResult] = None
    interaction_scores: Optional[InteractionResult] = None
    morans: Optional[MoransResult] = None
    local_diversity: Optional[LocalDiversityResult] = None

    execution_time_seconds: float = 0.0
    warnings: List[str] = field(default_factory=list)
    errors: List[str] = field(default_factory=list)

    def summary_dict(self) -> Dict[str, Any]:
        """Return summary dictionary for JSON export."""
        result = {
            "success": self.success,
            "n_cells_total": self.n_cells_total,
            "n_samples": self.n_samples,
            "n_samples_skipped": self.n_samples_skipped,
            "n_cell_types": self.n_cell_types,
            "cell_type_col": self.cell_type_col,
            "execution_time_seconds": round(self.execution_time_seconds, 2),
            "n_warnings": len(self.warnings),
            "n_errors": len(self.errors),
        }

        if self.neighborhood_enrichment:
            result["neighborhood_enrichment"] = self.neighborhood_enrichment.to_dict()
        if self.interaction_scores:
            result["interaction_scores"] = self.interaction_scores.to_dict()
        if self.morans:
            result["morans"] = self.morans.to_dict()
        if self.local_diversity:
            result["local_diversity"] = self.local_diversity.to_dict()

        return result


@dataclass
class MultiColumnSpatialResult:
    """Result from multi-column spatial analysis."""

    results: Dict[str, SpatialResult] = field(default_factory=dict)
    columns_processed: List[str] = field(default_factory=list)
    columns_skipped: List[str] = field(default_factory=list)
    total_execution_time_seconds: float = 0.0
    warnings: List[str] = field(default_factory=list)
    errors: List[str] = field(default_factory=list)

    @property
    def success(self) -> bool:
        return len(self.columns_processed) > 0

    def get_result(self, col: str) -> Optional[SpatialResult]:
        return self.results.get(col)

    def summary_dict(self) -> Dict[str, Any]:
        return {
            "columns_processed": self.columns_processed,
            "columns_skipped": self.columns_skipped,
            "per_column_summary": {
                col: result.summary_dict() for col, result in self.results.items()
            },
            "total_execution_time_seconds": round(self.total_execution_time_seconds, 2),
            "n_warnings": len(self.warnings),
            "n_errors": len(self.errors),
        }


class SpatialEngine:
    """Engine for comprehensive spatial analysis.

    Parameters
    ----------
    config : SpatialConfig, optional
        Configuration for analysis

    Example
    -------
    >>> from celltype_refinery.core.spatial import SpatialEngine, SpatialConfig
    >>> config = SpatialConfig.default()
    >>> engine = SpatialEngine(config)
    >>> result = engine.execute(adata, graphs_dir)
    """

    def __init__(self, config: Optional[SpatialConfig] = None):
        self.config = config or SpatialConfig.default()

    def execute(
        self,
        adata,
        graphs_dir: Path,
        cell_type_col: Optional[str] = None,
        skip_enrichment: bool = False,
        skip_interaction: bool = False,
        skip_morans: bool = False,
        skip_diversity: bool = False,
    ) -> SpatialResult:
        """Execute spatial analysis for a single cell type column.

        Parameters
        ----------
        adata : AnnData
            Annotated data with cell type labels
        graphs_dir : Path
            Directory containing spatial graphs (NPZ files)
        cell_type_col : str, optional
            Column name for cell types
        skip_enrichment : bool
            Skip neighborhood enrichment analysis
        skip_interaction : bool
            Skip interaction score analysis
        skip_morans : bool
            Skip Moran's I analysis
        skip_diversity : bool
            Skip local diversity analysis

        Returns
        -------
        SpatialResult
            Analysis results
        """
        start_time = time.time()
        result = SpatialResult()
        graphs_dir = Path(graphs_dir)

        ct_col = cell_type_col or self.config.cell_type_col
        result.cell_type_col = ct_col

        if ct_col not in adata.obs.columns:
            result.success = False
            result.errors.append(f"Cell type column '{ct_col}' not found")
            return result

        sample_col = self.config.sample_col
        if sample_col not in adata.obs.columns:
            result.success = False
            result.errors.append(f"Sample column '{sample_col}' not found")
            return result

        unique_samples = sorted(adata.obs[sample_col].unique().tolist())
        unique_types = sorted(adata.obs[ct_col].dropna().unique().tolist())

        type_counts = adata.obs[ct_col].value_counts()
        valid_types = type_counts[
            type_counts >= self.config.min_cells_per_type
        ].index.tolist()

        logger.info(f"Cell type column: {ct_col}")
        logger.info(f"Samples: {len(unique_samples)}")
        logger.info(f"Cell types: {len(unique_types)}")
        logger.info(f"Valid cell types (>={self.config.min_cells_per_type}): {len(valid_types)}")

        result.n_cell_types = len(valid_types)

        available, missing = validate_graphs_dir(graphs_dir, unique_samples)
        if missing:
            result.warnings.append(
                f"{len(missing)} samples missing graphs: {missing[:5]}..."
            )

        if not available:
            result.success = False
            result.errors.append("No spatial graphs found")
            return result

        samples = []
        region_col = self.config.region_col

        for sample_id in available:
            mask = adata.obs[sample_col] == sample_id
            n_cells = mask.sum()

            if n_cells < self.config.min_cells:
                result.n_samples_skipped += 1
                continue

            graph = load_spatial_graph(graphs_dir, sample_id)
            if graph is None:
                result.n_samples_skipped += 1
                continue

            cell_types = adata.obs.loc[mask, ct_col].values
            region = None
            if region_col in adata.obs.columns:
                region = adata.obs.loc[mask, region_col].iloc[0]

            samples.append(SampleInput(
                sample_id=sample_id,
                graph=graph,
                cell_types=cell_types,
                region=region,
            ))

        result.n_samples = len(samples)
        result.n_cells_total = sum(s.n_cells for s in samples)

        logger.info(f"Processing {len(samples)} samples ({result.n_cells_total:,} cells)")

        # Build analysis plan description
        analyses = []
        if not skip_enrichment:
            analyses.append(f"enrichment ({self.config.permutation.n_permutations} perms)")
        if not skip_interaction:
            analyses.append("interaction")
        if not skip_morans:
            analyses.append("morans")
        if not skip_diversity:
            analyses.append("diversity")
        logger.info(f"Analyses: {', '.join(analyses)}")

        enrichment_results = []
        interaction_results = []
        morans_results = []
        diversity_results = []
        sample_ids = []

        n_samples = len(samples)
        sample_start_time = time.time()

        for idx, sample in enumerate(samples, 1):
            sample_ids.append(sample.sample_id)
            iter_start = time.time()

            # Progress logging
            elapsed = time.time() - sample_start_time
            if idx > 1:
                avg_per_sample = elapsed / (idx - 1)
                remaining = avg_per_sample * (n_samples - idx + 1)
                logger.info(
                    f"[{idx}/{n_samples}] {sample.sample_id} ({sample.n_cells:,} cells) | "
                    f"elapsed: {elapsed:.0f}s, ETA: {remaining:.0f}s"
                )
            else:
                logger.info(f"[{idx}/{n_samples}] {sample.sample_id} ({sample.n_cells:,} cells)")

            if not skip_enrichment:
                try:
                    enrich_start = time.time()
                    z, p, n_c, n_e = compute_neighborhood_enrichment_sample(
                        graph=sample.graph,
                        cell_types=sample.cell_types,
                        unique_types=valid_types,
                        n_permutations=self.config.permutation.n_permutations,
                        n_threads=self.config.permutation.n_threads,
                        seed=self.config.permutation.seed,
                    )
                    enrichment_results.append((z, p, n_c, n_e))
                    logger.debug(f"  └─ enrichment: {time.time() - enrich_start:.1f}s")
                except Exception as e:
                    logger.warning(f"Enrichment failed for {sample.sample_id}: {e}")

            if not skip_interaction:
                try:
                    inter_start = time.time()
                    interactions, _, _ = compute_interactions_for_sample(
                        sample_id=sample.sample_id,
                        graph=sample.graph,
                        cell_types=sample.cell_types,
                        valid_types=valid_types,
                    )
                    interaction_results.append(interactions)
                    logger.debug(f"  └─ interaction: {time.time() - inter_start:.1f}s")
                except Exception as e:
                    logger.warning(f"Interaction failed for {sample.sample_id}: {e}")

            if not skip_morans:
                try:
                    morans_start = time.time()
                    morans = compute_morans_for_sample(
                        sample_id=sample.sample_id,
                        graph=sample.graph,
                        cell_types=sample.cell_types,
                        valid_types=valid_types,
                        min_cells_per_type=self.config.min_cells_per_type,
                    )
                    morans_results.append(morans)
                    logger.debug(f"  └─ morans: {time.time() - morans_start:.1f}s")
                except Exception as e:
                    logger.warning(f"Moran's I failed for {sample.sample_id}: {e}")

            if not skip_diversity:
                try:
                    div_start = time.time()
                    diversity = compute_diversity_for_sample(
                        sample_id=sample.sample_id,
                        graph=sample.graph,
                        cell_types=sample.cell_types,
                        unique_types=valid_types,
                        region=sample.region,
                    )
                    diversity_results.append(diversity)
                    logger.debug(f"  └─ diversity: {time.time() - div_start:.1f}s")
                except Exception as e:
                    logger.warning(f"Diversity failed for {sample.sample_id}: {e}")

        # Log total sample processing time
        total_sample_time = time.time() - sample_start_time
        logger.info(f"Sample processing complete: {total_sample_time:.1f}s ({total_sample_time/n_samples:.1f}s/sample avg)")

        # Aggregation phase
        logger.info("Aggregating results...")

        if enrichment_results:
            result.neighborhood_enrichment = aggregate_enrichment_results(
                sample_results=enrichment_results,
                sample_ids=sample_ids[: len(enrichment_results)],
                valid_types=valid_types,
                correction_method=self.config.correction.method,
                alpha=self.config.correction.alpha,
                n_permutations=self.config.permutation.n_permutations,
            )

        if interaction_results:
            result.interaction_scores = aggregate_interaction_results(
                sample_interactions=interaction_results,
                valid_types=valid_types,
                threshold=self.config.visualization.interaction_threshold,
            )

        if morans_results:
            result.morans = aggregate_morans_results(morans_results)

        if diversity_results:
            result.local_diversity = aggregate_diversity_results(diversity_results)

        result.execution_time_seconds = time.time() - start_time
        logger.info(f"Completed in {result.execution_time_seconds:.1f}s")

        return result

    def execute_multi(
        self,
        adata,
        graphs_dir: Path,
        cell_type_cols: Optional[List[str]] = None,
        **kwargs,
    ) -> MultiColumnSpatialResult:
        """Execute spatial analysis for multiple cell type columns.

        Parameters
        ----------
        adata : AnnData
            Annotated data
        graphs_dir : Path
            Directory with spatial graphs
        cell_type_cols : List[str], optional
            Columns to analyze
        **kwargs
            Passed to execute()

        Returns
        -------
        MultiColumnSpatialResult
            Results for all columns
        """
        start_time = time.time()
        result = MultiColumnSpatialResult()

        cols = cell_type_cols or self.config.cell_type_cols

        for col in cols:
            if col not in adata.obs.columns:
                result.columns_skipped.append(col)
                result.warnings.append(f"Column '{col}' not found in adata.obs")
                continue

            logger.info(f"\n{'='*60}\nProcessing column: {col}\n{'='*60}")

            col_result = self.execute(
                adata=adata,
                graphs_dir=graphs_dir,
                cell_type_col=col,
                **kwargs,
            )

            result.results[col] = col_result
            result.columns_processed.append(col)

        result.total_execution_time_seconds = time.time() - start_time

        return result
