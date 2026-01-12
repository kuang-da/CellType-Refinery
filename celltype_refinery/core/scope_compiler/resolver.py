"""Scope resolution logic."""

from dataclasses import dataclass
from typing import Any

from .normalizer import to_title_case
from .tree import MarkerMapTree


@dataclass
class ScopeResolution:
    """Result of resolving a scope definition."""

    labels: list[str]  # Sorted, Title Case normalized
    label_count: int  # -1 for wildcard
    roots: list[str]  # Root lineages
    is_wildcard: bool  # True if scope is "*"


def resolve_scope(
    scope: dict[str, Any], tree: MarkerMapTree, normalize: bool = True
) -> ScopeResolution:
    """
    Resolve a scope definition to explicit labels.

    Args:
        scope: The scope dict from a rule (contains include, match, exclude, exclude_match)
        tree: The MarkerMapTree built from the marker map
        normalize: Whether to apply Title Case normalization

    Returns:
        ScopeResolution with resolved labels, count, and roots
    """
    include_paths = scope.get("include", [])
    match_mode = scope.get("match", "exact")
    exclude_paths = scope.get("exclude", [])
    exclude_match = scope.get("exclude_match", "exact")

    # Case 1: Wildcard ("*" or match="all")
    if "*" in include_paths or match_mode == "all":
        return ScopeResolution(
            labels=["*"], label_count=-1, roots=["*"], is_wildcard=True
        )

    # Case 2: Specific paths
    include_labels: set[str] = set()
    include_roots: set[str] = set()

    for path in include_paths:
        node = tree.find_node(path)  # Raises KeyError if not found
        if match_mode == "descendants":
            labels = node.get_all_descendants(include_self=True)
        else:  # "exact"
            labels = [node.name]
        include_labels.update(labels)
        include_roots.add(node.root)

    # Apply exclusions
    exclude_labels: set[str] = set()
    for path in exclude_paths:
        node = tree.find_node(path)
        if exclude_match == "descendants":
            exclude_labels.update(node.get_all_descendants(include_self=True))
        else:
            exclude_labels.add(node.name)

    final_labels = include_labels - exclude_labels

    # Normalize and sort
    if normalize:
        final_labels = {to_title_case(label) for label in final_labels}
        include_roots = {to_title_case(root) for root in include_roots}

    return ScopeResolution(
        labels=sorted(final_labels),
        label_count=len(final_labels),
        roots=sorted(include_roots),
        is_wildcard=False,
    )
