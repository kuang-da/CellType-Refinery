"""Marker set loading and resolution for cell-type annotation.

This module handles loading marker maps from JSON files and resolving
marker names against available panel markers.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np


@dataclass(frozen=True)
class MarkerSet:
    """A set of markers defining a cell type.

    Attributes:
        label: Cell type name (e.g., "Ciliated Epithelium")
        path: Tuple of labels from root to this node (e.g., ("Epithelium", "Ciliated Epithelium"))
        level: Hierarchy depth (0=root, 1=first level, etc.)
        markers: Original marker names from JSON
        resolved_markers: Markers that exist in the panel
        missing_markers: Markers defined but not in panel
        anti_markers: Markers that should be absent
        resolved_anti_markers: Anti-markers that exist in panel
        missing_anti_markers: Anti-markers not in panel
        gating_overrides: Optional per-cell-type gating threshold overrides
    """
    label: str
    path: Tuple[str, ...]
    level: int
    markers: Tuple[str, ...]
    resolved_markers: Tuple[str, ...]
    missing_markers: Tuple[str, ...]
    anti_markers: Tuple[str, ...] = ()
    resolved_anti_markers: Tuple[str, ...] = ()
    missing_anti_markers: Tuple[str, ...] = ()
    gating_overrides: Optional[Dict[str, float]] = None


def canonicalize_marker(marker: str) -> str:
    """Normalize marker name for consistent lookup.

    Lowercases and strips whitespace for case-insensitive matching.
    """
    return marker.strip().lower()


def _collect_child_nodes(
    node: dict,
    current_path: Tuple[str, ...],
    level: int,
    var_names: Sequence[str],
    result: List[MarkerSet],
    logger: logging.Logger,
) -> None:
    """Recursively collect child nodes from the marker hierarchy.

    Args:
        node: Current node in the hierarchy
        current_path: Path from root to current node
        level: Current depth in hierarchy
        var_names: Available marker names in the panel
        result: List to append MarkerSet objects to
        logger: Logger instance
    """
    # Collect children from multiple sources:
    # 1. Explicit "subtypes" key
    # 2. Explicit "children" key
    # 3. Any dict values that aren't reserved keys
    children: Dict[str, dict] = {}

    if isinstance(node.get("subtypes"), dict):
        children.update(node["subtypes"])
    if isinstance(node.get("children"), dict):
        children.update(node["children"])

    # Also check for inline child definitions (dict values not in reserved keys)
    reserved_keys = {"markers", "subtypes", "children", "Anti_markers", "anti_markers", "gating_overrides"}
    for key, value in node.items():
        if key in reserved_keys:
            continue
        if isinstance(value, dict):
            children[key] = value

    if not children:
        return

    for child_name, child_node in children.items():
        child_path = current_path + (child_name,)
        markers_raw = child_node.get("markers", [])
        anti_raw = child_node.get("anti_markers", [])
        gating_overrides = child_node.get("gating_overrides", None)

        # Resolve markers
        resolved = []
        missing = []
        for marker in markers_raw:
            canon = canonicalize_marker(marker)
            matched = next(
                (v for v in var_names if canonicalize_marker(v) == canon),
                None
            )
            if matched:
                resolved.append(matched)
            else:
                missing.append(marker)

        # Resolve anti-markers
        anti_resolved = []
        anti_missing = []
        for anti in anti_raw:
            canon = canonicalize_marker(anti)
            matched = next(
                (v for v in var_names if canonicalize_marker(v) == canon),
                None
            )
            if matched:
                anti_resolved.append(matched)
            else:
                anti_missing.append(anti)

        if missing:
            logger.debug(
                "Marker set '%s': missing %d/%d markers: %s",
                child_name,
                len(missing),
                len(markers_raw),
                missing,
            )

        result.append(
            MarkerSet(
                label=child_name,
                path=child_path,
                level=level + 1,
                markers=tuple(markers_raw),
                resolved_markers=tuple(resolved),
                missing_markers=tuple(missing),
                anti_markers=tuple(anti_raw),
                resolved_anti_markers=tuple(anti_resolved),
                missing_anti_markers=tuple(anti_missing),
                gating_overrides=gating_overrides,
            )
        )

        # Recurse into children
        _collect_child_nodes(
            child_node, child_path, level + 1, var_names, result, logger
        )


def load_marker_sets(
    marker_map: Union[Dict, Path, str],
    var_names: Sequence[str],
    logger: Optional[logging.Logger] = None,
) -> List[MarkerSet]:
    """Load and resolve marker sets from a marker map.

    The marker map is a hierarchical JSON structure defining cell types
    and their characteristic markers. This function:
    1. Parses the JSON (if path provided)
    2. Traverses the hierarchy
    3. Resolves marker names against the panel
    4. Returns a flat list of MarkerSet objects

    Args:
        marker_map: Either a dict (parsed JSON), Path, or string path to JSON file
        var_names: List of marker names available in the panel
        logger: Optional logger instance

    Returns:
        List of MarkerSet objects, one per cell type in the hierarchy

    Raises:
        FileNotFoundError: If marker_map path doesn't exist
        json.JSONDecodeError: If marker_map file is not valid JSON
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    # Load JSON if path provided
    if isinstance(marker_map, (str, Path)):
        path = Path(marker_map)
        if not path.exists():
            raise FileNotFoundError(f"Marker map not found: {path}")
        with open(path, "r") as f:
            marker_map = json.load(f)

    result: List[MarkerSet] = []

    # Process each root category
    for root_name, root_node in marker_map.items():
        # Skip metadata keys
        if root_name.startswith("_"):
            continue

        markers_raw = root_node.get("markers", [])
        anti_raw = root_node.get("anti_markers", [])
        gating_overrides = root_node.get("gating_overrides", None)

        # Resolve markers
        resolved = []
        missing = []
        for marker in markers_raw:
            canon = canonicalize_marker(marker)
            matched = next(
                (v for v in var_names if canonicalize_marker(v) == canon),
                None
            )
            if matched:
                resolved.append(matched)
            else:
                missing.append(marker)

        # Resolve anti-markers
        anti_resolved = []
        anti_missing = []
        for anti in anti_raw:
            canon = canonicalize_marker(anti)
            matched = next(
                (v for v in var_names if canonicalize_marker(v) == canon),
                None
            )
            if matched:
                anti_resolved.append(matched)
            else:
                anti_missing.append(anti)

        if missing:
            logger.debug(
                "Root '%s': missing %d/%d markers: %s",
                root_name,
                len(missing),
                len(markers_raw),
                missing,
            )

        result.append(
            MarkerSet(
                label=root_name,
                path=(root_name,),
                level=0,
                markers=tuple(markers_raw),
                resolved_markers=tuple(resolved),
                missing_markers=tuple(missing),
                anti_markers=tuple(anti_raw),
                resolved_anti_markers=tuple(anti_resolved),
                missing_anti_markers=tuple(anti_missing),
                gating_overrides=gating_overrides,
            )
        )

        # Process children
        _collect_child_nodes(
            root_node, (root_name,), 0, var_names, result, logger
        )

    # Summary logging
    n_roots = sum(1 for ms in result if ms.level == 0)
    n_total = len(result)
    n_with_markers = sum(1 for ms in result if ms.resolved_markers)

    logger.info(
        "Loaded %d marker sets (%d roots, %d with resolved markers)",
        n_total,
        n_roots,
        n_with_markers,
    )

    return result


def compute_marker_document_frequency(
    marker_sets: List[MarkerSet],
) -> Dict[str, int]:
    """Compute document frequency for each marker.

    Document frequency = number of cell types that use a marker.
    Used for IDF weighting and commonness penalty.

    Args:
        marker_sets: List of MarkerSet objects

    Returns:
        Dict mapping marker name → count of cell types using it
    """
    doc_freq: Dict[str, int] = {}
    for mset in marker_sets:
        for marker in mset.resolved_markers:
            canon = canonicalize_marker(marker)
            doc_freq[canon] = doc_freq.get(canon, 0) + 1
            # Also store original name
            doc_freq[marker] = doc_freq.get(marker, 0) + 1
    return doc_freq
