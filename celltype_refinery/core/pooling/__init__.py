"""
Pooling module for Stage J lineage consolidation.

This module provides:
- Relabeling rules for pooling clusters by lineage
- PoolingEngine for executing pooling operations
- Export utilities for Stage I-compatible outputs
"""

from .rules import (
    load_root_types_from_marker_map,
    canonicalize_ambiguous_label,
    classify_cluster_for_pooling,
    build_pool_mapping,
    get_pool_summary,
)
from .engine import PoolingEngine, PoolingResult

__all__ = [
    "load_root_types_from_marker_map",
    "canonicalize_ambiguous_label",
    "classify_cluster_for_pooling",
    "build_pool_mapping",
    "get_pool_summary",
    "PoolingEngine",
    "PoolingResult",
]
