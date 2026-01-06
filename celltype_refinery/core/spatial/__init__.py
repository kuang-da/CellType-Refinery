"""Spatial analysis module for cell-type neighborhood analysis.

This module provides comprehensive spatial analysis including:
- Neighborhood enrichment analysis (permutation-based)
- Cell-type interaction scores (log2 fold enrichment)
- Moran's I spatial autocorrelation
- Local diversity indices (Shannon entropy)

Example Usage
-------------
Basic spatial analysis:

    >>> from celltype_refinery.core.spatial import SpatialEngine, SpatialConfig
    >>> config = SpatialConfig.default()
    >>> engine = SpatialEngine(config)
    >>> result = engine.execute(adata, graphs_dir)

With custom configuration:

    >>> from celltype_refinery.core.spatial import SpatialConfig
    >>> config = SpatialConfig.from_yaml("spatial.yaml")
    >>> engine = SpatialEngine(config)
    >>> result = engine.execute(adata, graphs_dir)
"""

__version__ = "1.0.0"

# Configuration
from .config import (
    SpatialConfig,
    PermutationConfig,
    CorrectionConfig,
    VisualizationConfig,
)

# Graph I/O
from .graph_io import (
    load_spatial_graph,
    get_graph_stats,
    list_available_graphs,
    validate_graphs_dir,
    load_graphs_for_samples,
)

# Parallel execution
from .parallel import (
    SampleInput,
    compute_single_permutation,
    batch_permutations_parallel,
    run_permutations_sequential,
    process_samples_adaptive,
    get_optimal_thread_allocation,
)

# Engine
from .engine import (
    SpatialEngine,
    SpatialResult,
    MultiColumnSpatialResult,
)

# Neighborhood Enrichment
from .neighborhood import (
    NeighborhoodEnrichmentResult,
    compute_enrichment_zscores,
    compute_neighborhood_enrichment_sample,
    apply_fdr_correction,
    aggregate_enrichment_results,
)

# Interaction
from .interaction import (
    InteractionResult,
    compute_neighbor_counts,
    compute_neighbor_counts_vectorized,
    compute_interaction_scores,
    compute_interactions_for_sample,
    aggregate_interaction_results,
)

# Moran's I
from .morans import (
    MoransResult,
    compute_morans_i,
    compute_morans_i_for_type,
    compute_morans_for_sample,
    aggregate_morans_results,
)

# Local Diversity
from .diversity import (
    LocalDiversityResult,
    compute_local_diversity,
    compute_diversity_for_sample,
    aggregate_diversity_results,
)

# Export
from .export import (
    export_single_result,
    export_multi_result,
    create_provenance,
    export_provenance,
    export_all,
)

# Visualization
from .viz import generate_all_visualizations

__all__ = [
    # Version
    "__version__",
    # Config
    "SpatialConfig",
    "PermutationConfig",
    "CorrectionConfig",
    "VisualizationConfig",
    # Graph I/O
    "load_spatial_graph",
    "get_graph_stats",
    "list_available_graphs",
    "validate_graphs_dir",
    "load_graphs_for_samples",
    # Parallel
    "SampleInput",
    "compute_single_permutation",
    "batch_permutations_parallel",
    "run_permutations_sequential",
    "process_samples_adaptive",
    "get_optimal_thread_allocation",
    # Engine
    "SpatialEngine",
    "SpatialResult",
    "MultiColumnSpatialResult",
    # Neighborhood Enrichment
    "NeighborhoodEnrichmentResult",
    "compute_enrichment_zscores",
    "compute_neighborhood_enrichment_sample",
    "apply_fdr_correction",
    "aggregate_enrichment_results",
    # Interaction
    "InteractionResult",
    "compute_neighbor_counts",
    "compute_neighbor_counts_vectorized",
    "compute_interaction_scores",
    "compute_interactions_for_sample",
    "aggregate_interaction_results",
    # Moran's I
    "MoransResult",
    "compute_morans_i",
    "compute_morans_i_for_type",
    "compute_morans_for_sample",
    "aggregate_morans_results",
    # Local Diversity
    "LocalDiversityResult",
    "compute_local_diversity",
    "compute_diversity_for_sample",
    "aggregate_diversity_results",
    # Export
    "export_single_result",
    "export_multi_result",
    "create_provenance",
    "export_provenance",
    "export_all",
    # Visualization
    "generate_all_visualizations",
]
