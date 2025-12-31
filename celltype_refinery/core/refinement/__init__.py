"""
Cell-type annotation refinement module.

This module provides tools for iterative refinement of cell-type annotations:
- RefinementEngine: Executes refinement plans on AnnData objects
- RefinePlan/RefineOp: Dataclasses for refinement operations
- AutoPolicy/ManualPolicy: Generate plans from criteria or YAML config
- InputAdapter: Normalize column names across different input formats
- Parallel subclustering for performance
- Cluster metrics for heterogeneity-based refinement decisions

Example workflow:
    >>> from celltype_refinery.core.refinement import (
    ...     RefinementEngine, AutoPolicy, InputAdapter, RefinePlan
    ... )
    >>>
    >>> # Normalize input data
    >>> adapter = InputAdapter()
    >>> annotations = adapter.adapt_cluster_annotations(raw_annotations)
    >>> scores = adapter.adapt_marker_scores(raw_scores)
    >>>
    >>> # Generate plan using automatic criteria
    >>> policy = AutoPolicy(score_threshold=1.0)
    >>> plan = policy.generate_plan(annotations, scores)
    >>>
    >>> # Execute plan
    >>> engine = RefinementEngine(plan, adata)
    >>> result = engine.execute()
    >>> print(f"Modified {result.n_cells_modified} cells")
"""

# Plan and operation dataclasses
from .plan import (
    RefineOp,
    RefinePlan,
    SubclusterOp,
    MergeOp,
    OverrideOp,
    RelabelOp,
    RescoreOp,
)

# Engine for executing plans
from .engine import (
    RefinementEngine,
    EngineResult,
    DERunner,
    default_scanpy_de_runner,
)

# Schema and validation
from .schema import (
    CanonicalSchema,
    InputAdapter,
    StageHAdapter,  # Backward compat alias
    normalize_inputs,
)
from .validation import (
    ValidationError,
    ValidationResult,
    MissingColumnError,
    InvalidClusterIdError,
    MissingProvenanceError,
    PathSeparatorError,
    TypeMismatchError,
)

# Operations
from .operations import (
    run_subcluster,
    create_hierarchical_labels,
    validate_cluster_exists,
    get_cluster_cells,
    get_current_iteration,
    initialize_iteration_level,
)

# Parallel processing
from .parallel import (
    SubclusterWorkItem,
    SubclusterResult,
    extract_work_items,
    worker_subcluster,
    run_subclusters_parallel,
)

# Provenance tracking
from .provenance import (
    RefinementProvenance,
    create_provenance,
    detect_parent_stage,
    store_provenance,
    get_provenance,
    get_provenance_chain,
)

# Policies
from .policies import (
    AutoPolicy,
    AutoPolicyConfig,
    ManualPolicy,
    ManualPolicyConfig,
)

# Cluster metrics for heterogeneity detection
from .cluster_metrics import (
    compute_cluster_marker_heterogeneity,
    merge_heterogeneity_with_annotations,
)

__all__ = [
    # Plan classes
    "RefineOp",
    "RefinePlan",
    "SubclusterOp",
    "MergeOp",
    "OverrideOp",
    "RelabelOp",
    "RescoreOp",
    # Engine
    "RefinementEngine",
    "EngineResult",
    "DERunner",
    "default_scanpy_de_runner",
    # Schema
    "CanonicalSchema",
    "InputAdapter",
    "StageHAdapter",
    "normalize_inputs",
    # Validation
    "ValidationError",
    "ValidationResult",
    "MissingColumnError",
    "InvalidClusterIdError",
    "MissingProvenanceError",
    "PathSeparatorError",
    "TypeMismatchError",
    # Operations
    "run_subcluster",
    "create_hierarchical_labels",
    "validate_cluster_exists",
    "get_cluster_cells",
    "get_current_iteration",
    "initialize_iteration_level",
    # Parallel
    "SubclusterWorkItem",
    "SubclusterResult",
    "extract_work_items",
    "worker_subcluster",
    "run_subclusters_parallel",
    # Provenance
    "RefinementProvenance",
    "create_provenance",
    "detect_parent_stage",
    "store_provenance",
    "get_provenance",
    "get_provenance_chain",
    # Policies
    "AutoPolicy",
    "AutoPolicyConfig",
    "ManualPolicy",
    "ManualPolicyConfig",
    # Cluster metrics
    "compute_cluster_marker_heterogeneity",
    "merge_heterogeneity_with_annotations",
]
