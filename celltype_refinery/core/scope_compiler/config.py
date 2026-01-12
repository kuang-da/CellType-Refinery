"""Configuration for scope compilation."""

from dataclasses import dataclass


@dataclass
class ScopeCompilerConfig:
    """Configuration for scope compilation."""

    normalize_labels: bool = True  # Apply Title Case normalization
    include_roots: bool = True  # Add scope_resolved_roots field
    validate_paths: bool = True  # Validate all paths exist
    fail_on_empty: bool = True  # Fail if resolved set is empty
    output_version: str = "v5"  # Version string for metadata
