"""Annotation module for marker-based cell-type assignment.

Provides hierarchical gating and marker scoring for cell-type annotation.

Key Features:
- Hierarchical marker scoring with multi-level gating
- Enhanced annotation exports (44 columns) with lineage tracking
- Stage H format exports (24 columns) for backward compatibility
- Cluster label mapping with score normalization
- Groups structure export for organized review by lineage
"""

from .assignment import annotate_obs
from .engine import AnnotationEngine, AnnotationParams, AnnotationResult
from .export import (
    # Column detection
    detect_cell_type_column,
    detect_cluster_column,
    # Core exports
    export_cluster_annotations_stage_h_format,
    export_cluster_label_mapping,
    export_enhanced_annotations,
    export_marker_scores,
    get_hierarchical_runner_up,
    run_annotation_exports,
    # Composition exports
    export_composition_stats,
    export_marker_scores_from_adata,
    export_cluster_annotations_simple,
    # Review exports
    export_review_summary,
    export_workflow_state,
    run_review_exports,
    # Group structure exports (NEW)
    GroupConfig,
    ClusterGroup,
    DEFAULT_GROUP_ORDER,
    build_groups_from_annotations,
    export_groups_structure,
    export_groups_derived_yaml,
    summarize_groups,
    get_root_group,
)
from .gating import (
    BASE_GATING_PARAMS,
    DEFAULT_GATING_PARAMS,
    assign_labels_hierarchical,
    merge_gating_params,
    extract_gating_params_from_marker_map,
    log_gating_params,
)
from .marker_loading import (
    MarkerSet,
    canonicalize_marker,
    compute_marker_document_frequency,
    load_marker_sets,
)
from .scoring import (
    ScoringContext,
    build_de_rank_lookup,
    compute_marker_idf,
    compute_marker_scores,
    compute_marker_scores_parallel,
    extract_existing_de_results,
)

__all__ = [
    # Engine
    "AnnotationEngine",
    "AnnotationParams",
    "AnnotationResult",
    # Marker loading
    "MarkerSet",
    "canonicalize_marker",
    "load_marker_sets",
    "compute_marker_document_frequency",
    # Scoring
    "ScoringContext",
    "compute_marker_scores",
    "compute_marker_scores_parallel",
    "compute_marker_idf",
    "extract_existing_de_results",
    "build_de_rank_lookup",
    # Gating
    "BASE_GATING_PARAMS",
    "DEFAULT_GATING_PARAMS",
    "assign_labels_hierarchical",
    "merge_gating_params",
    "extract_gating_params_from_marker_map",
    "log_gating_params",
    # Assignment
    "annotate_obs",
    # Column detection
    "detect_cell_type_column",
    "detect_cluster_column",
    # Core exports
    "export_enhanced_annotations",
    "export_cluster_annotations_stage_h_format",
    "export_marker_scores",
    "export_cluster_label_mapping",
    "get_hierarchical_runner_up",
    "run_annotation_exports",
    # Composition exports
    "export_composition_stats",
    "export_marker_scores_from_adata",
    "export_cluster_annotations_simple",
    # Review exports
    "export_review_summary",
    "export_workflow_state",
    "run_review_exports",
    # Group structure exports (NEW)
    "GroupConfig",
    "ClusterGroup",
    "DEFAULT_GROUP_ORDER",
    "build_groups_from_annotations",
    "export_groups_structure",
    "export_groups_derived_yaml",
    "summarize_groups",
    "get_root_group",
]
