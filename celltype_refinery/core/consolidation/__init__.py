"""Consolidation module for final cell-type label assignment.

This module provides the ConsolidationEngine for producing final, production-ready
cell-type labels from refinement output.

Key Features:
- Orphan rescue: Detect and rescue biologically plausible orphaned subtypes
- IEL rescue: Rescue intraepithelial immune cells (optional)
- Manual overrides: Apply expert-curated overrides via YAML config
- Relabel rules: Apply global label transformations
- Two-level harmonization: Map to fine (10-class) and broad (6-class) vocabularies
- Full audit trail: Provenance tracking for reproducibility

Example Usage
-------------
Basic consolidation:

    >>> from celltype_refinery.core.consolidation import ConsolidationEngine
    >>> engine = ConsolidationEngine()
    >>> result = engine.execute(adata, diagnostic_report, marker_scores)
    >>> print(f"Processed {result.n_cells_total:,} cells")

With configuration:

    >>> from celltype_refinery.core.consolidation import (
    ...     ConsolidationEngine,
    ...     ConsolidationConfig,
    ... )
    >>> config = ConsolidationConfig.from_yaml("consolidation.yaml")
    >>> engine = ConsolidationEngine(config=config, enable_orphan_rescue=True)
    >>> result = engine.execute(adata, diagnostic_report, marker_scores)

With harmonization:

    >>> from celltype_refinery.core.consolidation import (
    ...     ConsolidationEngine,
    ...     HarmonizeConfig,
    ... )
    >>> harmonize_config = HarmonizeConfig.default()
    >>> engine = ConsolidationEngine(harmonize_config=harmonize_config)
    >>> result = engine.execute(adata, diagnostic_report, marker_scores)
    >>> # adata.obs now has cell_type_fine and cell_type_broad columns
"""

# Configuration
from .config import (
    ConsolidationConfig,
    OverrideRule,
    RelabelRule,
    OrphanRescueConfig,
    OrphanRule,
    IELRescueConfig,
)

# Rules and classification
from .rules import (
    LabelCategory,
    ROOT_TYPES,
    INTERMEDIATE_TYPES,
    classify_label,
    select_final_label,
    apply_override,
    apply_relabel_rules,
    build_reason_chain,
    normalize_hybrid_label,
    simplify_label,
    is_valid_cell_type,
)

# Orphan detection and rescue
from .orphan_detection import (
    OrphanCandidate,
    OrphanPlausibility,
    OrphanAction,
    DEFAULT_RESCUE_RULES,
    detect_orphaned_subtypes,
    apply_orphan_rescue,
    get_orphan_rescue_map,
    summarize_orphans,
    compute_all_unassigned_scores,
)

# IEL rescue
from .iel_rescue import (
    IELCandidate,
    IELType,
    IELRescueConfig as IELConfig,
    detect_iel_candidates,
    apply_iel_rescue,
    get_iel_rescue_map,
    summarize_iel_candidates,
    summarize_iel_by_type,
)

# Harmonization
from .harmonize import (
    HarmonizeConfig,
    LabelInfo,
    FINE_VOCAB,
    BROAD_VOCAB,
    parse_label,
    map_to_fine,
    map_to_broad,
    harmonize_single_label,
    harmonize_adata,
)

# Engine
from .engine import (
    ConsolidationEngine,
    ConsolidationResult,
)

# Export functions
from .export import (
    export_consolidation_summary,
    export_consolidation_summary_detailed,
    export_confidence_breakdown,
    export_override_log,
    export_orphan_report,
    export_orphan_summary,
    export_iel_report,
    export_iel_summary,
    export_mapping_table,
    export_provenance_json,
    export_harmonized_summary,
    export_all,
)

__all__ = [
    # Configuration
    "ConsolidationConfig",
    "OverrideRule",
    "RelabelRule",
    "OrphanRescueConfig",
    "OrphanRule",
    "IELRescueConfig",
    # Rules
    "LabelCategory",
    "ROOT_TYPES",
    "INTERMEDIATE_TYPES",
    "classify_label",
    "select_final_label",
    "apply_override",
    "apply_relabel_rules",
    "build_reason_chain",
    "normalize_hybrid_label",
    "simplify_label",
    "is_valid_cell_type",
    # Orphan detection
    "OrphanCandidate",
    "OrphanPlausibility",
    "OrphanAction",
    "DEFAULT_RESCUE_RULES",
    "detect_orphaned_subtypes",
    "apply_orphan_rescue",
    "get_orphan_rescue_map",
    "summarize_orphans",
    "compute_all_unassigned_scores",
    # IEL rescue
    "IELCandidate",
    "IELType",
    "IELConfig",
    "detect_iel_candidates",
    "apply_iel_rescue",
    "get_iel_rescue_map",
    "summarize_iel_candidates",
    "summarize_iel_by_type",
    # Harmonization
    "HarmonizeConfig",
    "LabelInfo",
    "FINE_VOCAB",
    "BROAD_VOCAB",
    "parse_label",
    "map_to_fine",
    "map_to_broad",
    "harmonize_single_label",
    "harmonize_adata",
    # Engine
    "ConsolidationEngine",
    "ConsolidationResult",
    # Export functions
    "export_consolidation_summary",
    "export_consolidation_summary_detailed",
    "export_confidence_breakdown",
    "export_override_log",
    "export_orphan_report",
    "export_orphan_summary",
    "export_iel_report",
    "export_iel_summary",
    "export_mapping_table",
    "export_provenance_json",
    "export_harmonized_summary",
    "export_all",
]
