"""Decision rules for label classification and consolidation.

This module provides functions to:
- Classify labels as subtype/root/hybrid/unassigned
- Apply decision logic for final label selection
- Apply manual overrides and relabel rules
"""

from enum import Enum
from typing import Dict, Optional, Set, Tuple

from .config import OverrideRule, RelabelRule


class LabelCategory(Enum):
    """Categories for cell type labels."""

    SUBTYPE = "subtype"  # Specific cell type (e.g., Ciliated Epithelium)
    ROOT = "root"  # Root-level type (e.g., Epithelium)
    HYBRID = "hybrid"  # Ambiguous between lineages (e.g., Epithelium~Endothelium)
    UNASSIGNED = "unassigned"  # No assignment
    MIXED = "mixed"  # Mixed population within a lineage (e.g., Myeloids)


# Default root cell types
ROOT_TYPES: Set[str] = {
    "Epithelium",
    "Endothelium",
    "Immune Cells",
    "Mesenchymal Cells",
    "Misc",
}

# Intermediate types (not leaf, not root)
INTERMEDIATE_TYPES: Set[str] = {
    "Myeloids",
    "Lymphoids",
    "Granulocytes",
    "Monocytes",
    "Glandular Epithelium",  # Has Peg Cells subtype
    "T Cells",  # Has Activated T Cells subtype
}


def classify_label(label: str, root_types: Optional[Set[str]] = None) -> LabelCategory:
    """Classify a cell type label into categories.

    Parameters
    ----------
    label : str
        The cell type label to classify
    root_types : Set[str], optional
        Set of root type names. If None, uses default ROOT_TYPES.

    Returns
    -------
    LabelCategory
        The category of the label
    """
    if root_types is None:
        root_types = ROOT_TYPES

    # Handle edge cases
    if not label or label.lower() in ("", "nan", "none"):
        return LabelCategory.UNASSIGNED

    label = str(label).strip()

    # Check for Unassigned
    if label.lower() == "unassigned" or label.startswith("No root gate"):
        return LabelCategory.UNASSIGNED

    # Check for hybrid (contains ~)
    if "~" in label:
        return LabelCategory.HYBRID

    # Check for root types
    if label in root_types:
        return LabelCategory.ROOT

    # Check for intermediate/mixed types
    if label in INTERMEDIATE_TYPES:
        return LabelCategory.MIXED

    # Check if label contains "Mixed population"
    if "mixed population" in label.lower():
        return LabelCategory.MIXED

    # Otherwise it's a subtype
    return LabelCategory.SUBTYPE


def select_final_label(
    cluster_id: str,
    assigned_label: str,
    confidence_band: str,
    recommendation: str,
    root_types: Optional[Set[str]] = None,
    root_fail_reasons: Optional[str] = None,
) -> Tuple[str, str]:
    """Apply decision logic to select final label.

    Decision rules:
    - Clusters with any confidence -> use assigned_label as-is
    - Root labels -> keep as root (honest uncertainty)
    - Hybrid labels (X~Y) -> keep hybrid notation
    - Unassigned -> keep as "Unassigned"

    Parameters
    ----------
    cluster_id : str
        The cluster ID
    assigned_label : str
        The label from diagnostic report
    confidence_band : str
        Confidence level (high, medium, low, very_low)
    recommendation : str
        Recommendation from diagnostic (SKIP, SUBCLUSTER, RELABEL)
    root_types : Set[str], optional
        Set of root type names
    root_fail_reasons : str, optional
        Semicolon-separated list of root:reason pairs

    Returns
    -------
    Tuple[str, str]
        (final_label, decision_reason)
    """
    category = classify_label(assigned_label, root_types)

    if category == LabelCategory.UNASSIGNED:
        if root_fail_reasons:
            return "Unassigned", f"No root gate passed: {root_fail_reasons}"
        return "Unassigned", "No root gate passed"

    if category == LabelCategory.HYBRID:
        return assigned_label, "Hybrid population (ambiguous root)"

    if category == LabelCategory.ROOT:
        return assigned_label, "Root-level assignment (subtype unclear)"

    if category == LabelCategory.MIXED:
        return assigned_label, "Mixed population within lineage"

    # Subtype - use as-is
    return assigned_label, f"Subtype assignment ({confidence_band} confidence)"


def apply_override(
    cluster_id: str,
    current_label: str,
    override_map: Dict[str, OverrideRule],
) -> Tuple[str, Optional[str]]:
    """Apply manual override if exists.

    Parameters
    ----------
    cluster_id : str
        The cluster ID to check
    current_label : str
        The current label before override
    override_map : Dict[str, OverrideRule]
        Mapping from cluster_id to override rule

    Returns
    -------
    Tuple[str, Optional[str]]
        (final_label, override_reason or None if no override)
    """
    if cluster_id in override_map:
        override = override_map[cluster_id]
        reason = override.reason or "Manual override"
        return override.final_label, reason

    return current_label, None


def apply_relabel_rules(
    label: str,
    relabel_map: Dict[str, RelabelRule],
) -> Tuple[str, Optional[str]]:
    """Apply global relabel rules.

    Parameters
    ----------
    label : str
        The current label to check
    relabel_map : Dict[str, RelabelRule]
        Mapping from from_label to relabel rule

    Returns
    -------
    Tuple[str, Optional[str]]
        (final_label, relabel_reason or None if no relabel)
    """
    if label in relabel_map:
        rule = relabel_map[label]
        reason = rule.reason or f"Relabeled from '{label}'"
        return rule.to_label, reason

    return label, None


def build_reason_chain(
    base_reason: str,
    override_reason: Optional[str] = None,
    relabel_reason: Optional[str] = None,
    orphan_reason: Optional[str] = None,
) -> str:
    """Build a reason chain from multiple decision steps.

    Parameters
    ----------
    base_reason : str
        The base decision reason
    override_reason : str, optional
        Reason from manual override
    relabel_reason : str, optional
        Reason from relabel rule
    orphan_reason : str, optional
        Reason from orphan rescue

    Returns
    -------
    str
        Combined reason string
    """
    reasons = [base_reason]

    if orphan_reason:
        reasons.append(f"Orphan rescue: {orphan_reason}")

    if override_reason:
        reasons.append(f"Override: {override_reason}")

    if relabel_reason:
        reasons.append(f"Relabel: {relabel_reason}")

    return " | ".join(reasons)


def normalize_hybrid_label(label: str) -> str:
    """Normalize hybrid label to consistent format.

    Sorts components alphabetically: "B~A" -> "A~B"

    Parameters
    ----------
    label : str
        The hybrid label

    Returns
    -------
    str
        Normalized label with components sorted
    """
    if "~" not in label:
        return label

    parts = label.split("~")
    sorted_parts = sorted(parts)
    return "~".join(sorted_parts)


def simplify_label(label: str) -> str:
    """Simplify label for atlas output.

    1. Strip "(orphan)" suffix
    2. Sort hybrid components alphabetically

    Parameters
    ----------
    label : str
        The label to simplify

    Returns
    -------
    str
        Simplified label ready for atlas output
    """
    # Strip orphan suffix
    if label.endswith(" (orphan)"):
        label = label[:-9].strip()

    # Normalize hybrid (sort components)
    return normalize_hybrid_label(label)


def is_valid_cell_type(label: str) -> bool:
    """Check if a label is a valid cell type (not empty/nan/error).

    Parameters
    ----------
    label : str
        The label to check

    Returns
    -------
    bool
        True if valid
    """
    if not label:
        return False

    label_lower = str(label).lower().strip()
    invalid = {"", "nan", "none", "null", "na", "n/a"}

    return label_lower not in invalid
