"""Main scope compiler class."""

import copy
import json
from datetime import date
from pathlib import Path
from typing import Any

from .config import ScopeCompilerConfig
from .resolver import resolve_scope
from .tree import MarkerMapTree
from .validator import ValidationError, validate_marker_map


class ScopeCompiler:
    """Compiles scope paths in marker map state rules."""

    def __init__(self, config: ScopeCompilerConfig | None = None):
        self.config = config or ScopeCompilerConfig()

    def compile(self, marker_map: dict[str, Any]) -> dict[str, Any]:
        """
        Compile all state rule scopes in a marker map.

        Returns new marker map dict with compiled scope fields added.
        """
        # Build tree from TYPE hierarchy
        tree = MarkerMapTree(marker_map)

        # Validate if enabled
        if self.config.validate_paths:
            errors = validate_marker_map(marker_map, tree)
            if errors:
                raise ValidationError("Validation failed:\n" + "\n".join(errors))

        # Deep copy to avoid modifying original
        result = copy.deepcopy(marker_map)

        # Compile scopes in each rule
        cell_state_markers = result.get("_cell_state_markers", {})
        states = cell_state_markers.get("states", {})

        for state_name, state_def in states.items():
            rules = state_def.get("rules", [])
            for rule in rules:
                scope = rule.get("scope", {})
                resolution = resolve_scope(
                    scope, tree, normalize=self.config.normalize_labels
                )

                # Validate non-empty (unless wildcard)
                if (
                    self.config.fail_on_empty
                    and not resolution.is_wildcard
                    and resolution.label_count == 0
                ):
                    rule_id = rule.get("rule_id", "unknown")
                    raise ValidationError(
                        f"State '{state_name}', rule '{rule_id}': Empty resolved label set"
                    )

                # Add compiled fields
                rule["scope_resolved_labels"] = resolution.labels
                rule["scope_resolved_label_count"] = resolution.label_count
                if self.config.include_roots:
                    rule["scope_resolved_roots"] = resolution.roots

        # Update metadata
        self._update_metadata(result)

        return result

    def _update_metadata(self, marker_map: dict[str, Any]) -> None:
        """Update marker map metadata for new version."""
        metadata = marker_map.setdefault("_marker_map_metadata", {})

        # Store old version for changelog
        old_version = metadata.get("version", "unknown")

        # Update version and date
        metadata["version"] = self.config.output_version
        metadata["created_at"] = date.today().isoformat()

        # Add changelog
        changelog_key = f"key_changes_from_{old_version}"
        metadata[changelog_key] = [
            "Compiled/resolved each state rule's scope paths into explicit `scope_resolved_labels` lists",
            "Added `scope_resolved_label_count` (integer count, -1 for '*' wildcard)",
            "Added `scope_resolved_roots` (root lineages for debugging/filtering)",
            "Labels are Title Case normalized for consistent matching",
        ]

    def compile_file(self, input_path: Path, output_path: Path) -> None:
        """Load, compile, and save marker map."""
        with open(input_path, "r") as f:
            marker_map = json.load(f)

        compiled = self.compile(marker_map)

        with open(output_path, "w") as f:
            json.dump(compiled, f, indent=2)

        print(f"Compiled {input_path} -> {output_path}")


def compile_marker_map(marker_map: dict[str, Any], **config_kwargs) -> dict[str, Any]:
    """Convenience function to compile a marker map."""
    config = ScopeCompilerConfig(**config_kwargs)
    compiler = ScopeCompiler(config)
    return compiler.compile(marker_map)
