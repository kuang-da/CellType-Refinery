"""Composition analysis engine.

This module provides the CompositionEngine class that orchestrates
cell-type composition analysis.
"""

from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd

try:
    import scanpy as sc
except ImportError:
    sc = None

from .aggregation import (
    compute_composition_by_sample,
    aggregate_by_region,
    aggregate_by_donor,
    compute_global_summary,
    create_composition_wide,
)
from .biology import (
    TissueBiologyMetrics,
    GenericBiologyMetrics,
    compute_biology_metrics,
)
from .config import CompositionConfig
from .diversity import (
    compute_diversity_by_group,
    compute_diversity_summary,
)
from .enrichment import (
    compute_regional_enrichment,
    summarize_enrichments_by_region,
    summarize_enrichments_by_cell_type,
)


@dataclass
class CompositionResult:
    """Result of composition analysis.

    Attributes
    ----------
    composition_by_sample : pd.DataFrame
        Per-sample composition
    composition_by_region : pd.DataFrame
        Regional aggregation
    composition_by_donor : pd.DataFrame
        Per-donor aggregation
    composition_global : pd.DataFrame
        Global summary
    composition_wide : pd.DataFrame
        Wide-format matrix (samples x cell types)
    diversity_by_sample : pd.DataFrame
        Diversity metrics per sample
    diversity_summary : pd.DataFrame
        Diversity summary statistics
    biology_by_sample : pd.DataFrame, optional
        Biology metrics per sample
    biology_by_region : pd.DataFrame, optional
        Biology metrics per region
    enrichment : pd.DataFrame, optional
        Regional enrichment results
    enrichment_summary_region : pd.DataFrame, optional
        Enrichment summary by region
    enrichment_summary_cell_type : pd.DataFrame, optional
        Enrichment summary by cell type
    config : CompositionConfig
        Configuration used
    cell_type_column : str
        Cell type column analyzed
    provenance : Dict[str, Any]
        Execution provenance
    """

    composition_by_sample: pd.DataFrame
    composition_by_region: pd.DataFrame
    composition_by_donor: pd.DataFrame
    composition_global: pd.DataFrame
    composition_wide: pd.DataFrame
    diversity_by_sample: pd.DataFrame
    diversity_summary: pd.DataFrame
    biology_by_sample: Optional[pd.DataFrame] = None
    biology_by_region: Optional[pd.DataFrame] = None
    enrichment: Optional[pd.DataFrame] = None
    enrichment_summary_region: Optional[pd.DataFrame] = None
    enrichment_summary_cell_type: Optional[pd.DataFrame] = None
    config: CompositionConfig = field(default_factory=CompositionConfig)
    cell_type_column: str = "cell_type"
    provenance: Dict[str, Any] = field(default_factory=dict)


@dataclass
class MultiColumnResult:
    """Result of multi-column composition analysis.

    Attributes
    ----------
    results : Dict[str, CompositionResult]
        Results for each cell type column
    columns_analyzed : List[str]
        List of columns analyzed
    provenance : Dict[str, Any]
        Execution provenance
    """

    results: Dict[str, CompositionResult]
    columns_analyzed: List[str]
    provenance: Dict[str, Any] = field(default_factory=dict)


class CompositionEngine:
    """Engine for cell-type composition analysis.

    This engine computes composition statistics, diversity metrics,
    biology metrics, and regional enrichment tests.

    Parameters
    ----------
    config : CompositionConfig, optional
        Configuration. Uses defaults if not provided.
    biology_metrics : TissueBiologyMetrics, optional
        Custom tissue-specific biology metrics. Uses GenericBiologyMetrics if None.

    Example
    -------
    >>> from celltype_refinery.core.composition import CompositionEngine
    >>> engine = CompositionEngine()
    >>> result = engine.execute(adata, cell_type_col="cell_type")
    >>> print(f"Processed {len(result.composition_global)} cell types")
    """

    def __init__(
        self,
        config: Optional[CompositionConfig] = None,
        biology_metrics: Optional[TissueBiologyMetrics] = None,
    ):
        self.config = config or CompositionConfig.default()
        self.biology_metrics = biology_metrics or GenericBiologyMetrics()

    def execute(
        self,
        adata: "sc.AnnData",
        cell_type_col: str = "cell_type",
    ) -> CompositionResult:
        """Execute composition analysis.

        Parameters
        ----------
        adata : AnnData
            AnnData object with cell annotations
        cell_type_col : str
            Column with cell type labels

        Returns
        -------
        CompositionResult
            Composition analysis results
        """
        start_time = datetime.now()

        # Get configuration
        config = self.config
        sample_col = config.sample_column
        region_col = config.region_column
        donor_col = config.donor_column

        # Compute composition by sample
        composition_by_sample = compute_composition_by_sample(
            adata,
            cell_type_col=cell_type_col,
            sample_col=sample_col,
            region_col=region_col,
            donor_col=donor_col,
        )

        # Regional aggregation
        composition_by_region = aggregate_by_region(
            composition_by_sample,
            region_col=region_col,
            region_order=config.region_order if config.region_order else None,
        ) if region_col in adata.obs.columns else pd.DataFrame()

        # Donor aggregation
        composition_by_donor = aggregate_by_donor(
            composition_by_sample,
            donor_col=donor_col,
        ) if donor_col in adata.obs.columns else pd.DataFrame()

        # Global summary
        composition_global = compute_global_summary(composition_by_sample)

        # Wide format
        composition_wide = create_composition_wide(
            composition_by_sample,
            index_col=sample_col,
        )

        # Diversity metrics
        df = adata.obs[[sample_col, cell_type_col]].copy()
        diversity_by_sample = compute_diversity_by_group(
            df,
            group_col=sample_col,
            cell_type_col=cell_type_col,
            min_cells=config.diversity.min_cells_per_sample,
        )

        # Add region and donor columns to diversity_by_sample
        # Create sample-level metadata lookup
        sample_metadata_cols = [sample_col]
        if region_col in adata.obs.columns:
            sample_metadata_cols.append(region_col)
        if donor_col in adata.obs.columns:
            sample_metadata_cols.append(donor_col)

        if len(sample_metadata_cols) > 1:
            sample_metadata = (
                adata.obs[sample_metadata_cols]
                .drop_duplicates()
                .set_index(sample_col)
            )
            diversity_by_sample = diversity_by_sample.merge(
                sample_metadata,
                left_on=sample_col,
                right_index=True,
                how="left",
            )

        diversity_summary = compute_diversity_summary(diversity_by_sample)

        # Biology metrics (optional)
        biology_by_sample = None
        biology_by_region = None

        if not config.skip_biology:
            df_with_region = adata.obs[[sample_col, cell_type_col]].copy()
            if region_col in adata.obs.columns:
                df_with_region[region_col] = adata.obs[region_col]

            biology_by_sample = compute_biology_metrics(
                df_with_region,
                group_col=sample_col,
                cell_type_col=cell_type_col,
                patterns=config.patterns,
                metrics_impl=self.biology_metrics,
            )

            if region_col in adata.obs.columns:
                biology_by_region = compute_biology_metrics(
                    df_with_region,
                    group_col=region_col,
                    cell_type_col=cell_type_col,
                    patterns=config.patterns,
                    metrics_impl=self.biology_metrics,
                )

        # Regional enrichment (optional)
        enrichment = None
        enrichment_summary_region = None
        enrichment_summary_cell_type = None

        if not config.skip_enrichment and region_col in adata.obs.columns:
            enrichment = compute_regional_enrichment(
                composition_by_sample,
                region_col=region_col,
                cell_type_col="cell_type",
                value_col="proportion",
                sample_col=sample_col,
                min_samples_per_region=config.enrichment.min_samples_per_region,
                correction_method=config.enrichment.correction_method,
                alpha=config.enrichment.alpha,
            )

            if len(enrichment) > 0:
                enrichment_summary_region = summarize_enrichments_by_region(enrichment)
                enrichment_summary_cell_type = summarize_enrichments_by_cell_type(enrichment)

        end_time = datetime.now()

        # Build provenance
        provenance = {
            "timestamp": start_time.isoformat(),
            "duration_seconds": (end_time - start_time).total_seconds(),
            "n_cells": adata.n_obs,
            "n_samples": adata.obs[sample_col].nunique(),
            "n_cell_types": composition_global["cell_type"].nunique() if len(composition_global) > 0 else 0,
            "cell_type_column": cell_type_col,
            "config": config.to_dict(),
        }

        return CompositionResult(
            composition_by_sample=composition_by_sample,
            composition_by_region=composition_by_region,
            composition_by_donor=composition_by_donor,
            composition_global=composition_global,
            composition_wide=composition_wide,
            diversity_by_sample=diversity_by_sample,
            diversity_summary=diversity_summary,
            biology_by_sample=biology_by_sample,
            biology_by_region=biology_by_region,
            enrichment=enrichment,
            enrichment_summary_region=enrichment_summary_region,
            enrichment_summary_cell_type=enrichment_summary_cell_type,
            config=config,
            cell_type_column=cell_type_col,
            provenance=provenance,
        )

    def execute_multi(
        self,
        adata: "sc.AnnData",
        cell_type_columns: Optional[List[str]] = None,
    ) -> MultiColumnResult:
        """Execute composition analysis for multiple cell type columns.

        Parameters
        ----------
        adata : AnnData
            AnnData object with cell annotations
        cell_type_columns : List[str], optional
            Cell type columns to analyze. Uses config if not provided.

        Returns
        -------
        MultiColumnResult
            Results for each cell type column
        """
        if cell_type_columns is None:
            cell_type_columns = self.config.cell_type_columns

        results = {}
        for col in cell_type_columns:
            if col not in adata.obs.columns:
                continue
            results[col] = self.execute(adata, cell_type_col=col)

        return MultiColumnResult(
            results=results,
            columns_analyzed=list(results.keys()),
            provenance={
                "timestamp": datetime.now().isoformat(),
                "n_columns": len(results),
                "columns": list(results.keys()),
            },
        )
