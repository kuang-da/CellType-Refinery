"""Annotation module for marker-based cell-type assignment.

Provides hierarchical gating and marker scoring for cell-type annotation.
"""

from .assignment import annotate_obs
from .engine import AnnotationEngine, AnnotationParams, AnnotationResult
from .gating import (
    DEFAULT_GATING_PARAMS,
    assign_labels_hierarchical,
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
    "DEFAULT_GATING_PARAMS",
    "assign_labels_hierarchical",
    # Assignment
    "annotate_obs",
]
