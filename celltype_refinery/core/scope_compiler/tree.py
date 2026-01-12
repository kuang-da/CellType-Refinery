"""Marker map tree representation."""

from dataclasses import dataclass, field
from typing import Any


@dataclass
class MarkerMapNode:
    """A node in the marker map hierarchy."""

    name: str  # e.g., "T Cells"
    path: str  # e.g., "Immune Cells/Lymphoids/T Cells"
    root: str  # e.g., "Immune Cells"
    level: int  # 0 = root, 1 = first level, etc.
    children: list["MarkerMapNode"] = field(default_factory=list)

    def get_all_descendants(self, include_self: bool = True) -> list[str]:
        """Return all descendant names (including internal nodes)."""
        result = [self.name] if include_self else []
        for child in self.children:
            result.extend(child.get_all_descendants(include_self=True))
        return result


class MarkerMapTree:
    """Hierarchical tree built from marker map JSON."""

    # Keys that are metadata, not cell types
    METADATA_KEYS = {"_marker_map_metadata", "_cell_state_markers", "_gating_params"}

    def __init__(self, marker_map: dict[str, Any]):
        self.roots: list[MarkerMapNode] = []
        self.path_to_node: dict[str, MarkerMapNode] = {}
        self._build_tree(marker_map)

    def _build_tree(self, marker_map: dict[str, Any]) -> None:
        """Build tree from marker map JSON."""
        for key, value in marker_map.items():
            if key in self.METADATA_KEYS or key.startswith("_"):
                continue
            if not isinstance(value, dict):
                continue
            # Skip deprecated/empty entries
            if value.get("subtypes") == {} and "markers" not in value:
                continue
            root_node = self._build_node(key, key, key, 0, value)
            self.roots.append(root_node)

    def _build_node(
        self, name: str, path: str, root: str, level: int, data: dict
    ) -> MarkerMapNode:
        """Recursively build a node and its children."""
        node = MarkerMapNode(name=name, path=path, root=root, level=level)
        self.path_to_node[path] = node

        subtypes = data.get("subtypes", {})
        for child_name, child_data in subtypes.items():
            if isinstance(child_data, dict):
                child_path = f"{path}/{child_name}"
                child_node = self._build_node(
                    child_name, child_path, root, level + 1, child_data
                )
                node.children.append(child_node)

        return node

    def find_node(self, path: str) -> MarkerMapNode:
        """Find node by path. Raises KeyError if not found."""
        if path not in self.path_to_node:
            raise KeyError(f"Path not found in marker map: {path}")
        return self.path_to_node[path]

    def get_all_labels(self) -> list[str]:
        """Return all cell type labels in the tree."""
        labels = []
        for root in self.roots:
            labels.extend(root.get_all_descendants(include_self=True))
        return labels
