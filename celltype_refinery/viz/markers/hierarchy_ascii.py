"""ASCII tree visualization for marker hierarchies.

Generates console-friendly tree representation of marker map structure.
"""

from __future__ import annotations

from typing import Any, Dict, List, Tuple

from .style import (
    MarkerNode,
    format_marker_list,
    skip_metadata_keys,
)


def print_marker_hierarchy_ascii(
    marker_map: Dict[str, Any],
    *,
    show_markers: bool = True,
    show_anti_markers: bool = True,
    max_markers_display: int = 3,
    indent: str = "  ",
    branch_chars: Tuple[str, str, str, str] = ("├── ", "└── ", "│   ", "    "),
) -> str:
    """
    Generate ASCII tree representation of marker hierarchy for console output.

    Parameters
    ----------
    marker_map : Dict[str, Any]
        Marker hierarchy loaded from JSON file.
    show_markers : bool
        Whether to show positive markers.
    show_anti_markers : bool
        Whether to show anti-markers (prefixed with "!").
    max_markers_display : int
        Maximum markers per node before truncating.
    indent : str
        Base indentation string.
    branch_chars : Tuple[str, str, str, str]
        Characters for tree branches: (branch, last_branch, pipe, space).

    Returns
    -------
    str
        ASCII tree representation.

    Example Output
    --------------
    Cell Type Marker Hierarchy
    ├── Epithelium [E-cadherin, Pan-CK] [!CD45, !CD31]
    │   ├── Ciliated Epithelium [FOXJ1, DNALI1]
    │   └── Glandular Epithelium [PAX8, EPCAM]
    │       └── Peg Cells [CD44]
    └── Immune Cells [CD45]
        └── T Cells [CD3e]
    """
    if not marker_map:
        return ""

    branch, last_branch, pipe, space = branch_chars
    lines: List[str] = ["Cell Type Marker Hierarchy"]

    def format_node_line(node: MarkerNode) -> str:
        """Format a single node's display text."""
        parts = [node.name]
        if show_markers and node.markers:
            parts.append(format_marker_list(node.markers, max_markers_display))
        if show_anti_markers and node.anti_markers:
            parts.append(format_marker_list(node.anti_markers, max_markers_display, prefix="!"))
        return " ".join(parts)

    def render_subtree(
        data: Dict[str, Any],
        prefix: str = "",
        parent_path: str = "",
        level: int = 0,
    ) -> None:
        """Recursively render subtree."""
        cleaned = skip_metadata_keys(data)
        items = list(cleaned.items())

        for i, (name, node_data) in enumerate(items):
            if not isinstance(node_data, dict):
                continue

            is_last = i == len(items) - 1
            connector = last_branch if is_last else branch
            child_prefix = space if is_last else pipe

            path = f"{parent_path}/{name}" if parent_path else name
            node = MarkerNode(
                name=name,
                level=level,
                path=path,
                markers=node_data.get("markers", []) if isinstance(node_data.get("markers"), list) else [],
                anti_markers=node_data.get("Anti_markers", []) if isinstance(node_data.get("Anti_markers"), list) else [],
            )

            lines.append(f"{prefix}{connector}{format_node_line(node)}")

            # Recurse into subtypes
            if "subtypes" in node_data and isinstance(node_data["subtypes"], dict):
                render_subtree(
                    node_data["subtypes"],
                    prefix=prefix + child_prefix,
                    parent_path=path,
                    level=level + 1,
                )

    render_subtree(marker_map)
    return "\n".join(lines)


__all__ = ["print_marker_hierarchy_ascii"]
