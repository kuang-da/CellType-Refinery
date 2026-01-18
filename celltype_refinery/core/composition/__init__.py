"""Composition analysis module for cell-type statistics.

This module provides tools for analyzing cell-type composition across
samples, regions, and donors. It computes diversity metrics, biology
metrics, and regional enrichment tests.

Key Features:
- Composition aggregation: By sample, region, donor, global
- Diversity metrics: Shannon entropy, Simpson index, Pielou's evenness
- Biology metrics: Tissue-configurable metrics via abstract base class
- Regional enrichment: Mann-Whitney U tests with multiple testing correction
- Multi-column support: Analyze multiple cell type columns

Example Usage
-------------
Basic composition analysis:

    >>> from celltype_refinery.core.composition import CompositionEngine
    >>> engine = CompositionEngine()
    >>> result = engine.execute(adata, cell_type_col="cell_type")
    >>> print(f"Processed {len(result.composition_global)} cell types")

With configuration:

    >>> from celltype_refinery.core.composition import (
    ...     CompositionEngine,
    ...     CompositionConfig,
    ... )
    >>> config = CompositionConfig.from_yaml("composition.yaml")
    >>> engine = CompositionEngine(config=config)
    >>> result = engine.execute(adata)

With custom tissue biology metrics:

    >>> from celltype_refinery.core.composition import (
    ...     CompositionEngine,
    ...     TissueBiologyMetrics,
    ... )
    >>> class MyTissueMetrics(TissueBiologyMetrics):
    ...     # Implement tissue-specific metrics
    ...     pass
    >>> engine = CompositionEngine(biology_metrics=MyTissueMetrics())
    >>> result = engine.execute(adata)
"""

# Configuration
from .config import (
    CompositionConfig,
    PatternConfig,
    EnrichmentConfig,
    DiversityConfig,
    VisualizationConfig,
)

# Diversity metrics
from .diversity import (
    compute_shannon_entropy,
    compute_simpson_index,
    compute_evenness,
    compute_diversity_by_group,
    compute_diversity_summary,
    compute_diversity_from_composition,
)

# Aggregation
from .aggregation import (
    compute_composition_by_group,
    compute_composition_by_sample,
    aggregate_by_region,
    aggregate_by_donor,
    compute_global_summary,
    create_composition_wide,
    compute_rare_cell_types,
    compute_cell_type_co_occurrence,
    compute_composition_differences,
)

# Biology metrics
from .biology import (
    BiologyMetric,
    BiologyResult,
    TissueBiologyMetrics,
    GenericBiologyMetrics,
    classify_cell_type,
    count_by_pattern,
    compute_ratio,
    compute_percentage,
    compute_biology_metrics,
    get_default_biology_metrics,
    # Organ registry
    get_organ_metrics,
    register_organ_metrics,
    list_available_organs,
)

# Organ-specific biology metrics (auto-registers on import)
from .biology_ft import FallopianTubeMetrics, validate_ft_expectations
from .biology_uterus import UterusMetrics, validate_uterus_expectations

# Regional enrichment
from .enrichment import (
    compute_regional_enrichment,
    get_significant_enrichments,
    get_significant_depletions,
    summarize_enrichments_by_region,
    summarize_enrichments_by_cell_type,
    compute_pairwise_enrichment,
)

# Engine
from .engine import (
    CompositionEngine,
    CompositionResult,
    MultiColumnResult,
)

# Export functions
from .export import (
    export_composition_by_sample,
    export_composition_by_region,
    export_composition_by_donor,
    export_composition_global,
    export_composition_wide,
    export_diversity_by_sample,
    export_diversity_summary,
    export_biology_by_sample,
    export_biology_by_region,
    export_enrichment,
    export_enrichment_summary_region,
    export_enrichment_summary_cell_type,
    export_provenance,
    export_summary_json,
    export_all,
    export_all_multi,
)

# Visualization
from .viz import (
    generate_all_figures,
    generate_all_figures_multi,
    generate_dashboard,
)

__all__ = [
    # Configuration
    "CompositionConfig",
    "PatternConfig",
    "EnrichmentConfig",
    "DiversityConfig",
    "VisualizationConfig",
    # Diversity
    "compute_shannon_entropy",
    "compute_simpson_index",
    "compute_evenness",
    "compute_diversity_by_group",
    "compute_diversity_summary",
    "compute_diversity_from_composition",
    # Aggregation
    "compute_composition_by_group",
    "compute_composition_by_sample",
    "aggregate_by_region",
    "aggregate_by_donor",
    "compute_global_summary",
    "create_composition_wide",
    "compute_rare_cell_types",
    "compute_cell_type_co_occurrence",
    "compute_composition_differences",
    # Biology
    "BiologyMetric",
    "BiologyResult",
    "TissueBiologyMetrics",
    "GenericBiologyMetrics",
    "classify_cell_type",
    "count_by_pattern",
    "compute_ratio",
    "compute_percentage",
    "compute_biology_metrics",
    "get_default_biology_metrics",
    # Organ registry
    "get_organ_metrics",
    "register_organ_metrics",
    "list_available_organs",
    # Organ-specific biology
    "FallopianTubeMetrics",
    "validate_ft_expectations",
    "UterusMetrics",
    "validate_uterus_expectations",
    # Enrichment
    "compute_regional_enrichment",
    "get_significant_enrichments",
    "get_significant_depletions",
    "summarize_enrichments_by_region",
    "summarize_enrichments_by_cell_type",
    "compute_pairwise_enrichment",
    # Engine
    "CompositionEngine",
    "CompositionResult",
    "MultiColumnResult",
    # Export
    "export_composition_by_sample",
    "export_composition_by_region",
    "export_composition_by_donor",
    "export_composition_global",
    "export_composition_wide",
    "export_diversity_by_sample",
    "export_diversity_summary",
    "export_biology_by_sample",
    "export_biology_by_region",
    "export_enrichment",
    "export_enrichment_summary_region",
    "export_enrichment_summary_cell_type",
    "export_provenance",
    "export_summary_json",
    "export_all",
    "export_all_multi",
    # Visualization
    "generate_all_figures",
    "generate_all_figures_multi",
    "generate_dashboard",
]
