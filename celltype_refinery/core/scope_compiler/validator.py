"""Validation logic for scope definitions."""

from typing import Any

from .tree import MarkerMapTree


class ValidationError(Exception):
    """Raised when scope validation fails."""

    pass


def validate_scope(scope: dict[str, Any], tree: MarkerMapTree) -> list[str]:
    """
    Validate a scope definition.

    Returns list of error messages (empty if valid).
    """
    errors = []
    include_paths = scope.get("include", [])
    exclude_paths = scope.get("exclude", [])

    # Skip validation for wildcard
    if "*" in include_paths or scope.get("match") == "all":
        return errors

    # Validate include paths exist
    for path in include_paths:
        try:
            tree.find_node(path)
        except KeyError:
            errors.append(f"Include path not found: {path}")

    # Validate exclude paths exist
    for path in exclude_paths:
        try:
            tree.find_node(path)
        except KeyError:
            errors.append(f"Exclude path not found: {path}")

    return errors


def validate_marker_map(marker_map: dict[str, Any], tree: MarkerMapTree) -> list[str]:
    """
    Validate all scopes in a marker map.

    Returns list of all error messages.
    """
    all_errors = []

    cell_state_markers = marker_map.get("_cell_state_markers", {})
    states = cell_state_markers.get("states", {})

    for state_name, state_def in states.items():
        rules = state_def.get("rules", [])
        for i, rule in enumerate(rules):
            rule_id = rule.get("rule_id", f"rule_{i}")
            scope = rule.get("scope", {})
            errors = validate_scope(scope, tree)
            for error in errors:
                all_errors.append(f"State '{state_name}', rule '{rule_id}': {error}")

    return all_errors
