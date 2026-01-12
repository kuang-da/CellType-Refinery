"""JSON export utilities."""

import json
from pathlib import Path
from typing import Any


def export_json(
    marker_map: dict[str, Any], output_path: Path, indent: int = 2
) -> None:
    """Export marker map to JSON with pretty formatting."""
    with open(output_path, "w") as f:
        json.dump(marker_map, f, indent=indent, ensure_ascii=False)


def format_scope_summary(marker_map: dict[str, Any]) -> str:
    """Format a human-readable summary of compiled scopes."""
    lines = ["Compiled Scope Summary", "=" * 50, ""]

    cell_state_markers = marker_map.get("_cell_state_markers", {})
    states = cell_state_markers.get("states", {})

    for state_name, state_def in states.items():
        lines.append(f"State: {state_name}")
        lines.append("-" * 30)

        rules = state_def.get("rules", [])
        for rule in rules:
            rule_id = rule.get("rule_id", "unknown")
            labels = rule.get("scope_resolved_labels", [])
            count = rule.get("scope_resolved_label_count", 0)
            roots = rule.get("scope_resolved_roots", [])

            lines.append(f"  Rule: {rule_id}")
            lines.append(f"    Labels ({count}): {labels}")
            lines.append(f"    Roots: {roots}")
        lines.append("")

    return "\n".join(lines)
