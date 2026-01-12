"""Scope compiler for marker map state rules."""

from .compiler import ScopeCompiler, compile_marker_map
from .config import ScopeCompilerConfig
from .normalizer import to_title_case
from .resolver import ScopeResolution, resolve_scope
from .tree import MarkerMapNode, MarkerMapTree
from .validator import ValidationError, validate_scope

__all__ = [
    "ScopeCompilerConfig",
    "MarkerMapTree",
    "MarkerMapNode",
    "resolve_scope",
    "ScopeResolution",
    "ScopeCompiler",
    "compile_marker_map",
    "to_title_case",
    "validate_scope",
    "ValidationError",
]
