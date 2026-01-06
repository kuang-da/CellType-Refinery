"""Shared styling constants and utilities for marker visualization.

This module provides:
- Color palettes for marker and cell type visualization
- Data structures for marker hierarchy nodes
- Helper functions for parsing marker maps
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import matplotlib.pyplot as plt

PathLike = str | Path


# ============================================================================
# Color Constants
# ============================================================================

MARKER_COLORS = {
    "positive": "#2ecc71",       # Green - positive markers
    "anti": "#e74c3c",           # Red - anti-markers (exclusion)
    "cell_type": "#3498db",      # Blue - cell type nodes
    "edge_hierarchy": "#95a5a6", # Gray - parent-child edges
    "edge_marker": "#bdc3c7",    # Light gray - cell-type to marker edges
}

CATEGORY_COLORS = {
    "Epithelium": "#f1c40f",          # Yellow
    "Immune Cells": "#9b59b6",        # Purple
    "Endothelium": "#e74c3c",         # Red
    "Mesenchymal Cells": "#3498db",   # Blue
    "Misc": "#95a5a6",                # Gray
}


# ============================================================================
# Data Structures
# ============================================================================

@dataclass
class MarkerNode:
    """Internal representation of a marker hierarchy node."""
    name: str
    level: int
    path: str
    markers: List[str] = field(default_factory=list)
    anti_markers: List[str] = field(default_factory=list)
    notes: Optional[str] = None
    gating_overrides: Optional[Dict[str, Any]] = None
    parent: Optional[str] = None
    children: List[str] = field(default_factory=list)


# ============================================================================
# Matplotlib Style Setup
# ============================================================================

_DEFAULT_STYLE = {
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "legend.frameon": False,
}


def set_plot_style() -> None:
    """Apply a lightweight matplotlib style suitable for reports."""
    plt.style.use("seaborn-v0_8" if "seaborn-v0_8" in plt.style.available else "default")
    plt.rcParams.update(_DEFAULT_STYLE)


def ensure_parent(path: PathLike) -> Path:
    """Ensure parent directory exists and return Path object."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def save_figure(fig: plt.Figure, path: PathLike, *, dpi: int = 200) -> Path:
    """Save figure to disk and close it."""
    path = ensure_parent(path)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return path


# ============================================================================
# Marker Map Parsing Utilities
# ============================================================================

def skip_metadata_keys(marker_map: Dict[str, Any]) -> Dict[str, Any]:
    """Filter out metadata keys (starting with '_') from marker map."""
    return {k: v for k, v in marker_map.items() if not k.startswith("_")}


def parse_marker_hierarchy(
    marker_map: Dict[str, Any],
    parent_path: str = "",
    level: int = 0,
) -> Tuple[List[MarkerNode], Dict[str, MarkerNode]]:
    """
    Recursively parse marker JSON into flat list and lookup dict of nodes.

    Parameters
    ----------
    marker_map : Dict[str, Any]
        Marker hierarchy dictionary (or subtypes dict).
    parent_path : str
        Path of parent node (empty for root).
    level : int
        Current depth in hierarchy.

    Returns
    -------
    Tuple[List[MarkerNode], Dict[str, MarkerNode]]
        - nodes: Flat list of all MarkerNode objects
        - node_lookup: Dict mapping path -> MarkerNode
    """
    nodes: List[MarkerNode] = []
    node_lookup: Dict[str, MarkerNode] = {}
    cleaned = skip_metadata_keys(marker_map)

    for name, data in cleaned.items():
        if not isinstance(data, dict):
            continue

        path = f"{parent_path}/{name}" if parent_path else name

        node = MarkerNode(
            name=name,
            level=level,
            path=path,
            markers=data.get("markers", []) if isinstance(data.get("markers"), list) else [],
            anti_markers=data.get("Anti_markers", []) if isinstance(data.get("Anti_markers"), list) else [],
            notes=data.get("notes"),
            gating_overrides=data.get("gating_overrides"),
            parent=parent_path or None,
            children=[],
        )
        nodes.append(node)
        node_lookup[path] = node

        # Recurse into subtypes
        if "subtypes" in data and isinstance(data["subtypes"], dict):
            child_nodes, child_lookup = parse_marker_hierarchy(
                data["subtypes"],
                parent_path=path,
                level=level + 1
            )
            # Update children list for this node
            node.children = [c.path for c in child_nodes if c.level == level + 1]
            nodes.extend(child_nodes)
            node_lookup.update(child_lookup)

    return nodes, node_lookup


def format_marker_list(
    markers: List[str],
    max_display: int = 5,
    prefix: str = "",
) -> str:
    """
    Format marker list with optional truncation.

    Parameters
    ----------
    markers : List[str]
        List of marker names.
    max_display : int
        Maximum markers to show before truncating.
    prefix : str
        Prefix to add to each marker (e.g., "!" for anti-markers).

    Returns
    -------
    str
        Formatted marker string like "[M1, M2, +3 more]"
    """
    if not markers:
        return ""
    display = markers[:max_display]
    formatted = ", ".join(f"{prefix}{m}" for m in display)
    if len(markers) > max_display:
        formatted += f", +{len(markers) - max_display} more"
    return f"[{formatted}]"


def get_category_color(path: str) -> str:
    """Get color for a node based on its top-level category."""
    root = path.split("/")[0]
    return CATEGORY_COLORS.get(root, "#95a5a6")


__all__ = [
    # Constants
    "MARKER_COLORS",
    "CATEGORY_COLORS",
    # Data structures
    "MarkerNode",
    # Style functions
    "set_plot_style",
    "ensure_parent",
    "save_figure",
    # Parsing utilities
    "skip_metadata_keys",
    "parse_marker_hierarchy",
    "format_marker_list",
    "get_category_color",
    # Types
    "PathLike",
]
