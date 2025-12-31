"""Export functions for annotation results.

Provides export utilities for cluster annotations and enhanced annotations
with multi-level lineage tracking.

Key Functions:
- export_enhanced_annotations: Create 44-column enhanced annotations
- export_cluster_annotations_stage_h_format: Create Stage H format (24 columns)
- get_hierarchical_runner_up: Get runner-up label from decision tree
- export_marker_scores: Export marker scores DataFrame to CSV
- export_composition_stats: Export composition statistics
- export_review_summary: Export review summary JSON
- export_workflow_state: Export workflow state YAML
- generate_enhanced_annotations_44col: Create 44-column format from 24-col input
- run_review_exports: High-level orchestrator for review exports

Group Structure Export (NEW):
- GroupConfig: Configuration for lineage grouping (customizable group order)
- ClusterGroup: Data class representing a group of clusters by lineage
- build_groups_from_annotations: Derive groups from root_label column
- export_groups_structure: Create groups/ directory with per-group mappings
- export_groups_derived_yaml: Export groups_derived.yaml configuration
- summarize_groups: Format group summary for logging
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple, TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    import scanpy as sc


# =============================================================================
# Lineage Grouping Constants
# =============================================================================

LINEAGE_GROUPS = {
    'Epithelium': [
        'Epithelium', 'Ciliated Epithelium', 'Glandular Epithelium',
        'Secretory Epithelium', 'Peg Cells'
    ],
    'Endothelium': ['Endothelium', 'Vascular Endothelium', 'Lymphatic Endothelium'],
    'Mesenchymal Cells': [
        'Mesenchymal Cells', 'Smooth Muscle Cells', 'Fibroblasts', 'Pericytes'
    ],
    'Immune Cells': [
        'Immune Cells', 'T Cells', 'B Cells', 'Macrophages', 'Dendritic Cells',
        'Myeloids', 'Lymphoids', 'NK Cells', 'Natural-Killer (NK) Cells',
        'Neutrophils', 'Monocytes', 'Granulocytes', 'Activated T Cells'
    ],
    'Unassigned': ['Unassigned'],
}

# Default group ordering for review consistency (common tissue lineages)
DEFAULT_GROUP_ORDER: List[str] = [
    "Epithelium",
    "Endothelium",
    "Mesenchymal Cells",
    "Immune Cells",
    "Hybrids",
    "Unassigned",
]


# =============================================================================
# Group Configuration and Data Classes
# =============================================================================


@dataclass
class GroupConfig:
    """Configuration for lineage grouping.

    Allows customization of group ordering and naming conventions for
    different tissue types. Default configuration matches common spatial
    biology tissue lineages.

    Attributes
    ----------
    group_order : List[str]
        Ordered list of group names for consistent output ordering.
        Groups not in this list will be appended at the end.
    hybrid_group_name : str
        Name for the group containing ambiguous/hybrid clusters.
    unassigned_group_name : str
        Name for the group containing unassigned clusters.
    normalize_names : bool
        Whether to normalize group names (lowercase, replace spaces with underscores).

    Examples
    --------
    >>> config = GroupConfig()
    >>> config.group_order
    ['Epithelium', 'Endothelium', 'Mesenchymal Cells', 'Immune Cells', 'Hybrids', 'Unassigned']

    >>> # Custom ordering for a specific tissue
    >>> config = GroupConfig(group_order=['Neurons', 'Glia', 'Immune Cells', 'Unassigned'])
    """

    group_order: List[str] = field(default_factory=lambda: DEFAULT_GROUP_ORDER.copy())
    hybrid_group_name: str = "Hybrids"
    unassigned_group_name: str = "Unassigned"
    normalize_names: bool = True

    def get_dir_name(self, index: int, group_name: str) -> str:
        """Generate directory name for a group.

        Parameters
        ----------
        index : int
            1-based index of the group
        group_name : str
            Name of the group

        Returns
        -------
        str
            Directory name like 'g1_epithelium'
        """
        if self.normalize_names:
            safe_name = group_name.lower().replace(' ', '_').replace('-', '_')
        else:
            safe_name = group_name.replace(' ', '_')
        return f"g{index}_{safe_name}"


@dataclass
class ClusterGroup:
    """A group of clusters organized by lineage.

    Represents clusters that share a common root lineage (e.g., all Epithelium
    clusters) for organized review and processing.

    Attributes
    ----------
    name : str
        Group name (e.g., "Epithelium", "Hybrids")
    cluster_ids : Set[str]
        All cluster IDs in this group
    focus_label : str
        Label to use for focused processing (same as name for non-hybrids)
    subcluster_ids : Set[str]
        Only clusters with SUBCLUSTER recommendation (subset of cluster_ids)
    n_cells : int
        Total cells across all clusters in this group

    Examples
    --------
    >>> group = ClusterGroup(name="Epithelium", cluster_ids={"0", "1", "5"})
    >>> group.n_clusters
    3
    >>> group.has_work
    False  # No subclusters defined yet
    """

    name: str
    cluster_ids: Set[str] = field(default_factory=set)
    focus_label: str = ""
    subcluster_ids: Set[str] = field(default_factory=set)
    n_cells: int = 0

    def __post_init__(self):
        if not self.focus_label:
            self.focus_label = self.name

    @property
    def has_work(self) -> bool:
        """Return True if there are clusters to subcluster."""
        return len(self.subcluster_ids) > 0

    @property
    def n_clusters(self) -> int:
        """Return number of clusters in this group."""
        return len(self.cluster_ids)

    @property
    def n_subcluster(self) -> int:
        """Return number of clusters to subcluster."""
        return len(self.subcluster_ids)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for YAML/JSON export."""
        return {
            "name": self.name,
            "cluster_ids": sorted(self.cluster_ids),
            "focus_label": self.focus_label,
            "subcluster_ids": sorted(self.subcluster_ids),
            "n_clusters": self.n_clusters,
            "n_subcluster": self.n_subcluster,
            "n_cells": self.n_cells,
        }


# =============================================================================
# Group Building Functions
# =============================================================================


def get_root_group(
    row: pd.Series,
    config: Optional[GroupConfig] = None,
) -> str:
    """Map a cluster annotation row to its group name.

    Uses the is_ambiguous_root column to identify hybrids, otherwise
    uses the root_label column directly.

    Parameters
    ----------
    row : pd.Series
        A row from cluster_annotations DataFrame.
        Expected columns: root_label, is_ambiguous_root (optional)
    config : GroupConfig, optional
        Group configuration. Uses defaults if not provided.

    Returns
    -------
    str
        Group name (one of the configured group names)

    Examples
    --------
    >>> row = pd.Series({'is_ambiguous_root': True, 'root_label': 'Endo~Mesen'})
    >>> get_root_group(row)
    'Hybrids'

    >>> row = pd.Series({'is_ambiguous_root': False, 'root_label': 'Epithelium'})
    >>> get_root_group(row)
    'Epithelium'
    """
    if config is None:
        config = GroupConfig()

    # Check for hybrid using is_ambiguous_root column
    is_ambiguous = row.get("is_ambiguous_root", False)
    if is_ambiguous is True or str(is_ambiguous).lower() == "true":
        return config.hybrid_group_name

    # Check assigned_label for hybrid marker (~)
    assigned_label = str(row.get("assigned_label", ""))
    if "~" in assigned_label:
        return config.hybrid_group_name

    # Get root_label
    root_label = str(row.get("root_label", ""))

    # Handle Unassigned
    if root_label.lower() in ("unassigned", "unknown", "", "nan"):
        return config.unassigned_group_name

    # Map to known groups using keyword matching
    root_label_lower = root_label.lower()

    if "epithelium" in root_label_lower or "epithelial" in root_label_lower:
        return "Epithelium"
    if "endothelium" in root_label_lower or "endothelial" in root_label_lower:
        return "Endothelium"
    if "mesenchym" in root_label_lower or "stromal" in root_label_lower:
        return "Mesenchymal Cells"
    if "immune" in root_label_lower:
        return "Immune Cells"

    # Return as-is if it matches a known group in config
    if root_label in config.group_order:
        return root_label

    # Default to Unassigned for unknown types
    return config.unassigned_group_name


def build_groups_from_annotations(
    cluster_annotations: pd.DataFrame,
    diagnostic_report: Optional[pd.DataFrame] = None,
    config: Optional[GroupConfig] = None,
    logger: Optional[logging.Logger] = None,
) -> List[ClusterGroup]:
    """Build cluster groups from cluster annotations.

    Derives groups dynamically from the root_label column of cluster_annotations,
    using is_ambiguous_root to identify hybrid clusters.

    Parameters
    ----------
    cluster_annotations : pd.DataFrame
        DataFrame from cluster_annotations.csv.
        Required columns: cluster_id, root_label
        Optional columns: is_ambiguous_root, n_cells, assigned_label
    diagnostic_report : pd.DataFrame, optional
        DataFrame from diagnostic_report.csv.
        If provided, filters subcluster_ids to only SUBCLUSTER recommendations.
        Required columns: cluster_id, recommendation
    config : GroupConfig, optional
        Group configuration. Uses defaults if not provided.
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    List[ClusterGroup]
        List of ClusterGroup objects, ordered by config.group_order

    Examples
    --------
    >>> annotations = pd.read_csv("cluster_annotations.csv")
    >>> groups = build_groups_from_annotations(annotations)
    >>> for g in groups:
    ...     print(f"{g.name}: {g.n_clusters} clusters, {g.n_cells} cells")
    Epithelium: 17 clusters, 450000 cells
    Endothelium: 5 clusters, 120000 cells
    ...
    """
    logger = logger or logging.getLogger(__name__)

    if config is None:
        config = GroupConfig()

    # Ensure cluster_id is string
    annotations = cluster_annotations.copy()
    annotations["cluster_id"] = annotations["cluster_id"].astype(str)

    # Build set of SUBCLUSTER clusters from diagnostic report
    subcluster_clusters: Set[str] = set()
    if diagnostic_report is not None and not diagnostic_report.empty:
        diag = diagnostic_report.copy()
        diag["cluster_id"] = diag["cluster_id"].astype(str)
        if "recommendation" in diag.columns:
            subcluster_mask = diag["recommendation"] == "SUBCLUSTER"
            subcluster_clusters = set(diag.loc[subcluster_mask, "cluster_id"].tolist())

    # Group clusters by root group
    groups_dict: Dict[str, ClusterGroup] = {}

    for _, row in annotations.iterrows():
        cluster_id = str(row["cluster_id"])
        group_name = get_root_group(row, config)
        n_cells = int(row.get("n_cells", 0))

        if group_name not in groups_dict:
            groups_dict[group_name] = ClusterGroup(
                name=group_name,
                cluster_ids=set(),
                focus_label=group_name,
                subcluster_ids=set(),
                n_cells=0,
            )

        group = groups_dict[group_name]
        group.cluster_ids.add(cluster_id)
        group.n_cells += n_cells

        # Add to subcluster_ids if:
        # - No diagnostic report provided (subcluster all), OR
        # - Cluster is in SUBCLUSTER set
        if diagnostic_report is None or cluster_id in subcluster_clusters:
            group.subcluster_ids.add(cluster_id)

    # Order groups according to config.group_order
    ordered_groups: List[ClusterGroup] = []
    for group_name in config.group_order:
        if group_name in groups_dict:
            ordered_groups.append(groups_dict[group_name])

    # Add any unexpected groups at the end (sorted alphabetically)
    unexpected_groups = sorted([
        name for name in groups_dict.keys()
        if name not in config.group_order
    ])
    for group_name in unexpected_groups:
        ordered_groups.append(groups_dict[group_name])

    logger.debug(
        "Built %d groups from %d clusters",
        len(ordered_groups),
        len(annotations)
    )

    return ordered_groups


def summarize_groups(groups: List[ClusterGroup]) -> str:
    """Format a summary of groups for logging.

    Parameters
    ----------
    groups : List[ClusterGroup]
        List of ClusterGroup objects

    Returns
    -------
    str
        Formatted summary string
    """
    lines = ["Group Summary:", "-" * 60]
    total_clusters = 0
    total_cells = 0
    total_subcluster = 0

    for i, group in enumerate(groups, 1):
        prefix = f"g{i}"
        safe_name = group.name.lower().replace(' ', '_')
        status = "WORK" if group.has_work else "skip"
        lines.append(
            f"  {prefix}_{safe_name:20s}: "
            f"{group.n_clusters:4d} clusters, "
            f"{group.n_cells:8,d} cells, "
            f"{group.n_subcluster:3d} SUBCLUSTER [{status}]"
        )
        total_clusters += group.n_clusters
        total_cells += group.n_cells
        total_subcluster += group.n_subcluster

    lines.append("-" * 60)
    lines.append(
        f"  Total: {total_clusters} clusters, "
        f"{total_cells:,} cells, "
        f"{total_subcluster} to SUBCLUSTER"
    )

    return "\n".join(lines)


# =============================================================================
# Group Export Functions
# =============================================================================


def export_groups_derived_yaml(
    groups: List[ClusterGroup],
    output_path: Path,
    marker_map_path: Optional[str] = None,
    iteration: int = 1,
    logger: Optional[logging.Logger] = None,
) -> Path:
    """Export groups_derived.yaml configuration file.

    Parameters
    ----------
    groups : List[ClusterGroup]
        List of ClusterGroup objects from build_groups_from_annotations()
    output_path : Path
        Output file path for groups_derived.yaml
    marker_map_path : str, optional
        Path to marker map (included in metadata)
    iteration : int
        Current iteration number
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    Path
        Path to the exported YAML file
    """
    import yaml

    logger = logger or logging.getLogger(__name__)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    groups_data = {
        "version": "2.0",
        "iteration": iteration,
        "generated_at": datetime.now().isoformat(),
    }

    if marker_map_path:
        groups_data["marker_map"] = str(marker_map_path)

    # Add groups list
    groups_data["groups"] = [g.to_dict() for g in groups]

    with open(output_path, 'w') as f:
        yaml.safe_dump(groups_data, f, default_flow_style=False, sort_keys=False)

    logger.info("Exported groups_derived.yaml: %d groups → %s", len(groups), output_path)
    return output_path


def export_groups_structure(
    cluster_annotations: pd.DataFrame,
    cluster_label_mapping: pd.DataFrame,
    output_dir: Path,
    diagnostic_report: Optional[pd.DataFrame] = None,
    config: Optional[GroupConfig] = None,
    marker_map_path: Optional[str] = None,
    iteration: int = 1,
    logger: Optional[logging.Logger] = None,
) -> Dict[str, Any]:
    """Export groups/ directory structure from annotations.

    Creates per-group subdirectories with cluster_label_mapping.csv files
    and generates groups_derived.yaml for workflow tracking.

    This matches the reference implementation from ft/src/workflow/stage_i_unified.py,
    where the same global mapping is copied to each group directory for
    organizational review purposes.

    Parameters
    ----------
    cluster_annotations : pd.DataFrame
        DataFrame from cluster_annotations.csv.
        Required columns: cluster_id, root_label, n_cells
        Optional columns: is_ambiguous_root, assigned_label
    cluster_label_mapping : pd.DataFrame
        DataFrame from cluster_label_mapping.csv.
        Required columns: cluster_id, cell_type (or assigned_label), n_cells
    output_dir : Path
        Output directory (groups/ will be created as subdirectory)
    diagnostic_report : pd.DataFrame, optional
        DataFrame from diagnostic_report.csv for SUBCLUSTER filtering
    config : GroupConfig, optional
        Group configuration. Uses defaults if not provided.
    marker_map_path : str, optional
        Path to marker map (included in groups_derived.yaml)
    iteration : int
        Current iteration number
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    Dict[str, Any]
        Result dictionary with:
        - groups_dir: Path to groups/ directory
        - groups_derived_path: Path to groups_derived.yaml
        - group_dirs: Dict mapping group name to directory path
        - n_groups: Number of groups created
        - groups: List of ClusterGroup objects

    Examples
    --------
    >>> annotations = pd.read_csv("cluster_annotations.csv")
    >>> mapping = pd.read_csv("cluster_label_mapping.csv")
    >>> result = export_groups_structure(annotations, mapping, Path("output"))
    >>> print(result['n_groups'])
    6
    >>> print(result['group_dirs'].keys())
    dict_keys(['Epithelium', 'Endothelium', 'Mesenchymal Cells', ...])
    """
    logger = logger or logging.getLogger(__name__)

    if config is None:
        config = GroupConfig()

    output_dir = Path(output_dir)
    groups_dir = output_dir / "groups"
    groups_dir.mkdir(parents=True, exist_ok=True)

    # Build groups from annotations
    groups = build_groups_from_annotations(
        cluster_annotations=cluster_annotations,
        diagnostic_report=diagnostic_report,
        config=config,
        logger=logger,
    )

    if not groups:
        logger.warning("No groups derived from annotations")
        return {
            "groups_dir": groups_dir,
            "groups_derived_path": None,
            "group_dirs": {},
            "n_groups": 0,
            "groups": [],
        }

    # Log group summary
    logger.info(summarize_groups(groups))

    # Create per-group directories and export mapping files
    group_dirs: Dict[str, Path] = {}

    for i, group in enumerate(groups, 1):
        dir_name = config.get_dir_name(i, group.name)
        group_dir = groups_dir / dir_name
        group_dir.mkdir(parents=True, exist_ok=True)

        # Copy the global mapping to each group directory
        # (matches reference behavior - same file in each group for organization)
        mapping_path = group_dir / "cluster_label_mapping.csv"
        cluster_label_mapping.to_csv(mapping_path, index=False)

        group_dirs[group.name] = group_dir

        logger.info(
            "  Created %s/cluster_label_mapping.csv (%d clusters)",
            dir_name,
            group.n_clusters,
        )

    # Export groups_derived.yaml
    groups_derived_path = output_dir / "groups_derived.yaml"
    export_groups_derived_yaml(
        groups=groups,
        output_path=groups_derived_path,
        marker_map_path=marker_map_path,
        iteration=iteration,
        logger=logger,
    )

    logger.info(
        "Groups structure complete: %d groups → %s",
        len(groups),
        groups_dir,
    )

    return {
        "groups_dir": groups_dir,
        "groups_derived_path": groups_derived_path,
        "group_dirs": group_dirs,
        "n_groups": len(groups),
        "groups": groups,
    }


# =============================================================================
# Column Detection Utilities
# =============================================================================


def detect_cell_type_column(adata: "sc.AnnData") -> str:
    """Detect best cell type column with fallback chain.

    Priority: cell_type_curated > cell_type_lvl1 > cell_type_lvl0 >
              cell_type_auto > assigned_label > cell_type

    Parameters
    ----------
    adata : sc.AnnData
        AnnData object

    Returns
    -------
    str
        Name of detected cell type column

    Raises
    ------
    ValueError
        If no cell type column found
    """
    priority = [
        'cell_type_curated',
        'cell_type_lvl1',
        'cell_type_lvl0',
        'cell_type_auto',
        'assigned_label',
        'cell_type',
    ]
    for col in priority:
        if col in adata.obs.columns:
            return col
    raise ValueError("No cell type column found in adata.obs")


def detect_cluster_column(adata: "sc.AnnData") -> str:
    """Detect best cluster column with fallback chain.

    Priority: cluster_id > cluster_lvl1 > cluster_lvl0 > leiden

    Parameters
    ----------
    adata : sc.AnnData
        AnnData object

    Returns
    -------
    str
        Name of detected cluster column
    """
    priority = ['cluster_id', 'cluster_lvl1', 'cluster_lvl0', 'leiden']
    for col in priority:
        if col in adata.obs.columns:
            return col
    return 'cluster_id'  # Default


# =============================================================================
# Helper Functions
# =============================================================================


def _sort_cluster_key(cluster_id: str) -> Tuple:
    """Sort key for cluster IDs that handles subclusters (e.g., '1:0', '18:2').

    Examples:
        "1" -> (1, -1)
        "1:0" -> (1, 0)
        "18:2" -> (18, 2)
    """
    if ":" in cluster_id:
        parent, sub = cluster_id.split(":", 1)
        try:
            return (int(parent), int(sub))
        except ValueError:
            return (float("inf"), cluster_id)
    else:
        try:
            return (int(cluster_id), -1)
        except ValueError:
            return (float("inf"), cluster_id)


def _classify_confidence_band(score: float) -> str:
    """Classify score into confidence bands.

    Parameters
    ----------
    score : float
        Marker score

    Returns
    -------
    str
        Confidence level: "high", "medium", "low", or "very_low"

    Note
    ----
    Thresholds aligned with reference ft pipeline:
    - score >= 2.0: "high"
    - score >= 1.0: "medium"
    - score >= 0.5: "low"
    - else: "very_low"
    """
    if score >= 2.0:
        return "high"
    elif score >= 1.0:
        return "medium"
    elif score >= 0.5:
        return "low"
    else:
        return "very_low"


# =============================================================================
# Runner-Up Detection
# =============================================================================


def get_hierarchical_runner_up(
    cluster_id: str,
    decision_steps: pd.DataFrame,
    cluster_annotations: pd.DataFrame,
    marker_scores: Optional[pd.DataFrame] = None,
) -> Tuple[str, float, float]:
    """Get hierarchical runner-up: the losing sibling at the tightest decision point.

    Instead of using flat score ranking (which can incorrectly identify the assigned
    label as runner-up), this function finds the actual competing sibling at the
    decision step where the margin was smallest.

    For `ambiguous_siblings` stop reason, extracts runner-up from the composition
    JSON stored in annotations (which contains `_runner_up_by_score` at the correct
    stopping depth).

    Parameters
    ----------
    cluster_id : str
        The cluster ID to look up
    decision_steps : pd.DataFrame
        DataFrame with columns: cluster_id, step_idx, parent_label, child_label,
        child_passed_gate, child_score, selected, margin_to_runner_up
    cluster_annotations : pd.DataFrame
        DataFrame with columns: cluster_id, min_margin_along_path, assigned_label,
        stop_reason, composition (JSON), etc.
    marker_scores : pd.DataFrame, optional
        DataFrame with marker scores to look up runner-up score if not in composition

    Returns
    -------
    Tuple[str, float, float]
        (runner_up_label, runner_up_score, min_margin)
        - runner_up_label: The losing sibling at the tightest decision point
        - runner_up_score: Score of the runner-up
        - min_margin: The minimum margin along the hierarchical path (gap)

    Notes
    -----
    Returns ("", 0.0, NaN) when no meaningful competition exists (e.g., single
    candidate at each level). This matches reference behavior for byte-for-byte
    compatibility.
    """
    # Handle empty inputs
    if cluster_annotations is None or cluster_annotations.empty:
        return "", 0.0, float("nan")

    # Get annotation row for this cluster
    ann = cluster_annotations[cluster_annotations["cluster_id"].astype(str) == str(cluster_id)]
    if ann.empty:
        return "", 0.0, float("nan")

    ann_row = ann.iloc[0]
    min_margin = ann_row.get("min_margin_along_path", float("nan"))
    if pd.isna(min_margin):
        min_margin = float("nan")

    stop_reason = str(ann_row.get("stop_reason", ""))

    # For ambiguous_siblings, extract runner-up from composition JSON
    # The composition JSON stores _runner_up_by_score at the correct stopping depth
    if stop_reason == "ambiguous_siblings":
        composition_str = ann_row.get("composition", "")
        if composition_str and isinstance(composition_str, str) and composition_str.strip():
            try:
                composition = json.loads(composition_str)
                runner_up_label = composition.get("_runner_up_by_score", "")
                gap = composition.get("_gap", min_margin)

                # Get runner-up score from composition or marker_scores
                runner_up_score = 0.0
                if runner_up_label:
                    # Try to get score from marker_scores if available
                    if marker_scores is not None and not marker_scores.empty:
                        score_row = marker_scores[
                            (marker_scores["cluster_id"].astype(str) == str(cluster_id)) &
                            (marker_scores["label"] == runner_up_label)
                        ]
                        if not score_row.empty:
                            runner_up_score = float(score_row.iloc[0].get("score", 0.0))
                    # Alternative: infer from gap and winner score
                    if runner_up_score == 0.0 and not pd.isna(gap):
                        top_label = composition.get("_top_by_score", "")
                        if top_label and marker_scores is not None and not marker_scores.empty:
                            top_row = marker_scores[
                                (marker_scores["cluster_id"].astype(str) == str(cluster_id)) &
                                (marker_scores["label"] == top_label)
                            ]
                            if not top_row.empty:
                                top_score = float(top_row.iloc[0].get("score", 0.0))
                                runner_up_score = top_score - float(gap)

                if runner_up_label:
                    return (
                        str(runner_up_label),
                        float(runner_up_score),
                        float(gap) if not pd.isna(gap) else float(min_margin),
                    )
            except (json.JSONDecodeError, TypeError, KeyError):
                pass  # Fall through to decision_steps logic

    # Standard logic: use decision_steps
    if decision_steps is None or decision_steps.empty:
        return "", 0.0, min_margin

    # Find decision steps for this cluster
    cluster_steps = decision_steps[decision_steps["cluster_id"].astype(str) == str(cluster_id)]
    if cluster_steps.empty:
        return "", 0.0, min_margin

    if "margin_to_runner_up" not in cluster_steps.columns:
        return "", 0.0, min_margin

    # For ambiguous_siblings, the winner step is marked selected=False but has margin_to_runner_up
    # We need to consider ALL steps with valid margins, not just selected ones
    steps_with_margins = cluster_steps[
        cluster_steps["margin_to_runner_up"].notna() &
        (cluster_steps["margin_to_runner_up"] != np.inf)
    ]
    if steps_with_margins.empty:
        return "", 0.0, min_margin

    # Find the step with minimum margin (tightest decision)
    min_idx = steps_with_margins["margin_to_runner_up"].idxmin()
    min_step = steps_with_margins.loc[min_idx]
    step_margin = float(min_step["margin_to_runner_up"])

    # Get all siblings at that step (same parent, passed gate)
    siblings = cluster_steps[
        (cluster_steps["step_idx"] == min_step["step_idx"]) &
        (cluster_steps["parent_label"] == min_step["parent_label"]) &
        (cluster_steps["child_passed_gate"] == True)
    ]

    if len(siblings) < 2:
        return "", 0.0, min_margin

    # Sort by score descending; the winner is the one with highest score
    siblings = siblings.sort_values("child_score", ascending=False)
    runner = siblings.iloc[1]  # Second-best is the runner-up

    return (
        str(runner.get("child_label", "")),
        float(runner.get("child_score", 0.0)),
        float(step_margin),  # Use the actual gap at this step, not min_margin_along_path
    )


# =============================================================================
# Enhanced Annotations Export
# =============================================================================


def export_enhanced_annotations(
    adata: "sc.AnnData",
    output_path: Path,
    marker_scores: Optional[pd.DataFrame] = None,
    decision_steps: Optional[pd.DataFrame] = None,
    cluster_annotations: Optional[pd.DataFrame] = None,
    max_iteration: int = 5,
    include_runner_up: bool = True,
    include_top_markers: bool = True,
    top_n_markers: int = 3,
    logger: Optional[logging.Logger] = None,
) -> pd.DataFrame:
    """Export enhanced cluster annotations with multi-level lineage tracking.

    Creates a comprehensive CSV with 44 columns tracking:
    - Cluster provenance (origin, iteration created)
    - Cell type assignments at each iteration level
    - Marker scores and confidence at each level
    - Runner-up labels and margins
    - Top expressing markers
    - Regional distribution

    Parameters
    ----------
    adata : sc.AnnData
        AnnData with cluster annotations
    output_path : Path
        Output CSV path
    marker_scores : pd.DataFrame, optional
        Marker scores DataFrame. Falls back to adata.uns if not provided.
    decision_steps : pd.DataFrame, optional
        Hierarchical decision steps for runner-up detection
    cluster_annotations : pd.DataFrame, optional
        Pre-computed cluster annotations. Falls back to adata.uns.
    max_iteration : int
        Maximum iteration level to include (default: 5)
    include_runner_up : bool
        Whether to include runner-up labels (default: True)
    include_top_markers : bool
        Whether to include top marker columns (default: True)
    top_n_markers : int
        Number of top markers to include (default: 3)
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    pd.DataFrame
        Enhanced annotations DataFrame with 44 columns
    """
    logger = logger or logging.getLogger(__name__)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Get marker scores from adata.uns if not provided
    if marker_scores is None:
        if "marker_scores_refined" in adata.uns:
            marker_scores = adata.uns["marker_scores_refined"]
            if isinstance(marker_scores, dict):
                marker_scores = pd.DataFrame(marker_scores)
        else:
            marker_scores = pd.DataFrame()

    # Determine cluster column (use highest available level)
    cluster_col = None
    for level in range(max_iteration, -1, -1):
        col = f"cluster_lvl{level}"
        if col in adata.obs:
            cluster_col = col
            break
    if cluster_col is None:
        cluster_col = "cluster_lvl0"

    # Get unique cluster IDs
    all_cluster_ids = sorted(
        adata.obs[cluster_col].astype(str).unique(),
        key=_sort_cluster_key
    )

    logger.info("Generating enhanced annotations for %d clusters", len(all_cluster_ids))

    # Build score lookup
    score_lookup = {}
    if not marker_scores.empty and "cluster_id" in marker_scores.columns:
        for cid in marker_scores["cluster_id"].unique():
            cid_str = str(cid)
            cluster_scores = marker_scores[marker_scores["cluster_id"].astype(str) == cid_str]
            if not cluster_scores.empty:
                # Get all scores for this cluster
                score_lookup[cid_str] = cluster_scores

    # Build records
    records = []
    total_cells = len(adata)

    for cluster_id in all_cluster_ids:
        cluster_id_str = str(cluster_id)
        mask = adata.obs[cluster_col].astype(str) == cluster_id_str
        n_cells = mask.sum()

        if n_cells == 0:
            continue

        # Determine origin cluster and iteration created
        if "cluster_lvl0" in adata.obs:
            origin_cluster = adata.obs.loc[mask, "cluster_lvl0"].astype(str).iloc[0]
        else:
            origin_cluster = cluster_id_str.split(":")[0] if ":" in cluster_id_str else cluster_id_str

        iteration_created = cluster_id_str.count(":")

        # Build record with base columns
        record = {
            "cluster_id": cluster_id_str,
            "origin_cluster": origin_cluster,
            "iteration_created": iteration_created,
            "n_cells": n_cells,
            "proportion": round(n_cells / total_cells, 6),
        }

        # Add cell type and score columns for each iteration level
        for level in range(max_iteration + 1):
            level_suffix = f"_lvl{level}"
            ct_col = f"cell_type_lvl{level}" if level > 0 else "cell_type_lvl0"
            cluster_lvl_col = f"cluster_lvl{level}"

            # Get cell type for this level
            cell_type = ""
            if ct_col in adata.obs:
                ct_values = adata.obs.loc[mask, ct_col].value_counts()
                if len(ct_values) > 0:
                    cell_type = str(ct_values.index[0])

            record[f"cell_type{level_suffix}"] = cell_type

            # Get score for this level (from marker_scores or adata.uns)
            score = 0.0
            assigned_path = ""
            stop_reason = ""
            confidence = ""
            coverage = 0.0
            is_ambiguous_root = False
            root_label = ""
            root_fail_reasons = ""
            min_margin = None
            decision_trace = ""

            # Look for hierarchical annotation data in adata.uns
            hier_key = f"hierarchical_annotations_lvl{level}"
            if hier_key in adata.uns:
                hier_df = adata.uns[hier_key]
                if isinstance(hier_df, pd.DataFrame) and not hier_df.empty:
                    row = hier_df[hier_df["cluster_id"].astype(str) == cluster_id_str]
                    if not row.empty:
                        row = row.iloc[0]
                        score = float(row.get("score", 0))
                        assigned_path = str(row.get("assigned_path", ""))
                        stop_reason = str(row.get("stop_reason", ""))
                        confidence = str(row.get("confidence_level", ""))
                        coverage = float(row.get("coverage", 0))
                        is_ambiguous_root = bool(row.get("is_ambiguous_root", False))
                        root_label = str(row.get("root_label", ""))
                        root_fail_reasons = str(row.get("root_fail_reasons", ""))
                        min_margin = row.get("min_margin_along_path")
                        decision_trace = str(row.get("decision_trace", ""))

            # Fallback to marker_scores if no hierarchical data
            if score == 0.0 and cluster_id_str in score_lookup:
                cluster_scores = score_lookup[cluster_id_str]
                if cell_type and not cluster_scores.empty:
                    label_scores = cluster_scores[cluster_scores["label"] == cell_type]
                    if not label_scores.empty:
                        score = float(label_scores.iloc[0].get("score", 0))
                        coverage = float(label_scores.iloc[0].get("coverage", 0))

            record[f"score{level_suffix}"] = round(score, 3)
            record[f"assigned_path{level_suffix}"] = assigned_path
            record[f"stop_reason{level_suffix}"] = stop_reason
            record[f"confidence{level_suffix}"] = confidence or _classify_confidence_band(round(float(score), 3))
            record[f"coverage{level_suffix}"] = round(coverage, 3)
            record[f"is_ambiguous_root{level_suffix}"] = is_ambiguous_root
            record[f"root_label{level_suffix}"] = root_label
            record[f"root_fail_reasons{level_suffix}"] = root_fail_reasons
            if min_margin is not None and not pd.isna(min_margin):
                record[f"min_margin{level_suffix}"] = round(float(min_margin), 3)
            else:
                record[f"min_margin{level_suffix}"] = None
            record[f"decision_trace{level_suffix}"] = decision_trace

        # Add runner-up information (based on highest level)
        if include_runner_up and cluster_annotations is not None:
            runner_up, runner_up_score, margin = get_hierarchical_runner_up(
                cluster_id_str, decision_steps, cluster_annotations, marker_scores
            )
            record["runner_up_label"] = runner_up
            record["runner_up_score"] = round(runner_up_score, 3)
            record["margin_to_runner_up"] = round(margin, 3) if margin != float("inf") else None

        # Add top markers
        if include_top_markers and cluster_id_str in score_lookup:
            cluster_scores = score_lookup[cluster_id_str]
            if not cluster_scores.empty and "label" in cluster_scores.columns:
                top_scores = cluster_scores.nlargest(top_n_markers, "score")
                for i, (_, row) in enumerate(top_scores.iterrows(), 1):
                    record[f"top_marker_{i}"] = str(row.get("label", ""))
                    record[f"top_marker_{i}_score"] = round(float(row.get("score", 0)), 3)

        # Add regional distribution (if region column exists)
        if "region" in adata.obs:
            region_counts = adata.obs.loc[mask, "region"].value_counts()
            for region, count in region_counts.items():
                record[f"region_{region}_count"] = count
                record[f"region_{region}_pct"] = round(100 * count / n_cells, 1)

        # Add mean enrichment and positive fraction
        if cluster_id_str in score_lookup:
            cluster_scores = score_lookup[cluster_id_str]
            # Find the assigned label for this cluster
            assigned_label = record.get("cell_type_lvl1", record.get("cell_type_lvl0", ""))
            if assigned_label and not cluster_scores.empty:
                label_scores = cluster_scores[cluster_scores["label"] == assigned_label]
                if not label_scores.empty:
                    record["mean_enrichment"] = round(
                        float(label_scores.iloc[0].get("mean_enrichment", 0)), 3
                    )
                    record["mean_positive_fraction"] = round(
                        float(label_scores.iloc[0].get("mean_positive_fraction", 0)), 3
                    )

        records.append(record)

    df = pd.DataFrame.from_records(records)
    df.to_csv(output_path, index=False)

    logger.info("Exported enhanced annotations: %d clusters, %d columns → %s",
                len(df), len(df.columns), output_path)
    return df


# =============================================================================
# Stage H Format Export
# =============================================================================


def export_cluster_annotations_stage_h_format(
    adata: "sc.AnnData",
    enhanced_annotations: pd.DataFrame,
    output_path: Path,
    iteration: int = 1,
    logger: Optional[logging.Logger] = None,
) -> pd.DataFrame:
    """Export cluster annotations in Stage H format with 24 columns.

    Provides backward-compatible output matching the original Stage H format,
    derived from enhanced annotations at a specific iteration level.

    Parameters
    ----------
    adata : sc.AnnData
        AnnData with cluster annotations
    enhanced_annotations : pd.DataFrame
        Enhanced annotations DataFrame from export_enhanced_annotations()
    output_path : Path
        Output CSV path
    iteration : int
        Iteration level to export (default: 1)
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    pd.DataFrame
        Stage H format annotations with 24 columns:
        - cluster_id, origin_cluster, iteration_created, proportion
        - mean_enrichment, mean_positive_fraction, n_cells
        - assigned_label, assigned_path, assigned_level, assigned_score
        - root_label, confidence, min_margin_along_path, margin_is_infinite
        - stop_reason, stopped_before_leaf, composition, decision_trace
        - coverage, resolved_markers, is_ambiguous_root
        - ambiguous_root_candidates, root_fail_reasons, confidence_level
    """
    logger = logger or logging.getLogger(__name__)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if enhanced_annotations.empty:
        logger.warning("Empty enhanced annotations; creating empty Stage H format file")
        df = pd.DataFrame()
        df.to_csv(output_path, index=False)
        return df

    records = []
    for _, row in enhanced_annotations.iterrows():
        cluster_id = str(row.get("cluster_id", ""))

        # Get values from the specified iteration level
        cell_type = row.get(f"cell_type_lvl{iteration}", row.get("cell_type_lvl1", ""))
        score = row.get(f"score_lvl{iteration}", row.get("score_lvl1", 0))
        path = row.get(f"assigned_path_lvl{iteration}", row.get("assigned_path_lvl1", ""))
        level = row.get(f"assigned_level_lvl{iteration}", row.get("assigned_level_lvl1", -1))

        # Derive is_ambiguous_root
        is_ambiguous = row.get(f"is_ambiguous_root_lvl{iteration}", False)
        if not is_ambiguous and "~" in str(cell_type):
            is_ambiguous = True

        # Derive root_label
        root_label = row.get(f"root_label_lvl{iteration}", "")
        if not root_label:
            if "~" in str(cell_type):
                root_label = cell_type
            elif "/" in str(path):
                root_label = str(path).split("/")[0].strip()
            else:
                root_label = cell_type

        records.append({
            # Provenance columns
            "cluster_id": cluster_id,
            "origin_cluster": row.get("origin_cluster", str(cluster_id).split(":")[0]),
            "iteration_created": row.get("iteration_created", str(cluster_id).count(":")),
            "proportion": row.get("proportion", 0),
            "mean_enrichment": row.get("mean_enrichment", 0),
            "mean_positive_fraction": row.get("mean_positive_fraction", 0),

            # Stage H format columns
            "n_cells": row.get("n_cells", 0),
            "assigned_label": cell_type,
            "assigned_path": path if path else cell_type,
            "assigned_level": int(level) if level is not None and level != "" else -1,
            "assigned_score": round(float(score) if score else 0.0, 3),
            "root_label": root_label,
            "confidence": row.get(f"confidence_lvl{iteration}", score),
            "min_margin_along_path": row.get(f"min_margin_lvl{iteration}"),
            "margin_is_infinite": row.get(f"margin_is_infinite_lvl{iteration}", False),
            "stop_reason": row.get(f"stop_reason_lvl{iteration}", ""),
            "stopped_before_leaf": row.get(f"stopped_before_leaf_lvl{iteration}", True),
            "composition": row.get(f"composition_lvl{iteration}", ""),
            "decision_trace": row.get(f"decision_trace_lvl{iteration}", ""),
            "coverage": row.get(f"coverage_lvl{iteration}", 0),
            "resolved_markers": row.get(f"resolved_markers_lvl{iteration}", ""),
            "is_ambiguous_root": is_ambiguous,
            "ambiguous_root_candidates": row.get(f"ambiguous_root_candidates_lvl{iteration}", ""),
            "root_fail_reasons": row.get(f"root_fail_reasons_lvl{iteration}", ""),
            "confidence_level": row.get("confidence_level", _classify_confidence_band(round(float(score), 3) if score else 0)),
        })

    df = pd.DataFrame.from_records(records)
    df.to_csv(output_path, index=False)

    logger.info("Exported Stage H format annotations: %d clusters → %s", len(df), output_path)
    return df


# =============================================================================
# Marker Scores Export
# =============================================================================


def export_marker_scores(
    adata: "sc.AnnData",
    output_path: Path,
    logger: Optional[logging.Logger] = None,
) -> Path:
    """Export marker_scores_refined from adata.uns to CSV.

    Parameters
    ----------
    adata : sc.AnnData
        AnnData with marker scores in uns["marker_scores_refined"]
    output_path : Path
        Output CSV path
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    Path
        Path to the exported CSV file
    """
    logger = logger or logging.getLogger(__name__)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if "marker_scores_refined" in adata.uns:
        scores = adata.uns["marker_scores_refined"]
        if isinstance(scores, pd.DataFrame):
            scores.to_csv(output_path, index=False)
            logger.info("Exported marker scores: %d rows → %s", len(scores), output_path)
        elif isinstance(scores, (dict, list)):
            df = pd.DataFrame(scores)
            df.to_csv(output_path, index=False)
            logger.info("Exported marker scores: %d rows → %s", len(df), output_path)
        else:
            logger.warning("marker_scores_refined has unexpected type: %s", type(scores))
            pd.DataFrame().to_csv(output_path, index=False)
    else:
        logger.warning("No marker_scores_refined in adata.uns; created empty file")
        pd.DataFrame().to_csv(output_path, index=False)

    return output_path


# =============================================================================
# Cluster Label Mapping Export
# =============================================================================


def export_cluster_label_mapping(
    adata: "sc.AnnData",
    output_path: Path,
    marker_scores: Optional[pd.DataFrame] = None,
    cluster_col: str = "cluster_lvl1",
    label_col: str = "cell_type_lvl1",
    logger: Optional[logging.Logger] = None,
) -> pd.DataFrame:
    """Export cluster to label mapping with scores.

    Creates a simple mapping table with cluster_id, assigned_label, assigned_score,
    and n_cells. Sets assigned_score=0.0 for Unassigned clusters.

    Parameters
    ----------
    adata : sc.AnnData
        AnnData with cluster annotations
    output_path : Path
        Output CSV path
    marker_scores : pd.DataFrame, optional
        Marker scores for score lookup
    cluster_col : str
        Cluster column name (default: "cluster_lvl1")
    label_col : str
        Label column name (default: "cell_type_lvl1")
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    pd.DataFrame
        Cluster label mapping with columns: cluster_id, assigned_label, assigned_score, n_cells
    """
    logger = logger or logging.getLogger(__name__)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Get marker scores from adata.uns if not provided
    if marker_scores is None:
        if "marker_scores_refined" in adata.uns:
            marker_scores = adata.uns["marker_scores_refined"]
            if isinstance(marker_scores, dict):
                marker_scores = pd.DataFrame(marker_scores)
        else:
            marker_scores = pd.DataFrame()

    # Build score lookup (max score per cluster)
    score_lookup = {}
    if not marker_scores.empty and "cluster_id" in marker_scores.columns and "score" in marker_scores.columns:
        best_scores = marker_scores.loc[marker_scores.groupby("cluster_id")["score"].idxmax()]
        for _, row in best_scores.iterrows():
            score_lookup[str(row["cluster_id"])] = float(row["score"])

    # Get cluster-label mapping from adata.obs
    if cluster_col not in adata.obs or label_col not in adata.obs:
        logger.warning("Required columns not found: %s or %s", cluster_col, label_col)
        df = pd.DataFrame()
        df.to_csv(output_path, index=False)
        return df

    cluster_label_df = (
        adata.obs.groupby(cluster_col)
        .agg({label_col: "first"})
        .reset_index()
        .rename(columns={cluster_col: "cluster_id", label_col: "assigned_label"})
    )

    # Add n_cells
    cell_counts = adata.obs[cluster_col].value_counts()
    cluster_label_df["n_cells"] = cluster_label_df["cluster_id"].map(cell_counts)

    # Add assigned_score
    cluster_label_df["assigned_score"] = cluster_label_df["cluster_id"].astype(str).map(
        lambda x: score_lookup.get(x, 0.0)
    )

    # Set score=0.0 for Unassigned clusters (no root gate passed = no confidence)
    # This matches the reference pipeline semantics: assigned_score indicates
    # confidence in the assigned label, not best possible match potential
    unassigned_mask = cluster_label_df["assigned_label"] == "Unassigned"
    cluster_label_df.loc[unassigned_mask, "assigned_score"] = 0.0

    # Round scores
    cluster_label_df["assigned_score"] = cluster_label_df["assigned_score"].round(3)

    # Sort by cluster_id
    cluster_label_df = cluster_label_df.sort_values(
        "cluster_id",
        key=lambda x: x.map(lambda c: _sort_cluster_key(str(c)))
    ).reset_index(drop=True)

    # Reorder columns
    cluster_label_df = cluster_label_df[["cluster_id", "assigned_label", "assigned_score", "n_cells"]]

    cluster_label_df.to_csv(output_path, index=False)
    logger.info("Exported cluster label mapping: %d clusters → %s", len(cluster_label_df), output_path)

    return cluster_label_df


# =============================================================================
# High-Level Export Function
# =============================================================================


def run_annotation_exports(
    adata: "sc.AnnData",
    output_dir: Path,
    marker_scores: Optional[pd.DataFrame] = None,
    decision_steps: Optional[pd.DataFrame] = None,
    cluster_annotations: Optional[pd.DataFrame] = None,
    iteration: int = 1,
    cluster_col: str = "cluster_lvl1",
    label_col: str = "cell_type_lvl1",
    logger: Optional[logging.Logger] = None,
) -> Dict[str, Path]:
    """Run all annotation exports.

    Exports enhanced annotations, Stage H format, marker scores, and cluster mapping.

    Parameters
    ----------
    adata : sc.AnnData
        AnnData with annotations
    output_dir : Path
        Output directory
    marker_scores : pd.DataFrame, optional
        Marker scores DataFrame
    decision_steps : pd.DataFrame, optional
        Decision steps for runner-up detection
    cluster_annotations : pd.DataFrame, optional
        Pre-computed cluster annotations
    iteration : int
        Iteration level for Stage H format (default: 1)
    cluster_col : str
        Cluster column name
    label_col : str
        Label column name
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    Dict[str, Path]
        Mapping from export type to output path
    """
    logger = logger or logging.getLogger(__name__)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    output_paths = {}

    # Export enhanced annotations
    enhanced_path = output_dir / "cluster_annotations_enhanced.csv"
    enhanced_df = export_enhanced_annotations(
        adata=adata,
        output_path=enhanced_path,
        marker_scores=marker_scores,
        decision_steps=decision_steps,
        cluster_annotations=cluster_annotations,
        max_iteration=5,
        logger=logger,
    )
    output_paths["enhanced"] = enhanced_path

    # Export Stage H format
    stage_h_path = output_dir / "cluster_annotations.csv"
    export_cluster_annotations_stage_h_format(
        adata=adata,
        enhanced_annotations=enhanced_df,
        output_path=stage_h_path,
        iteration=iteration,
        logger=logger,
    )
    output_paths["stage_h"] = stage_h_path

    # Export marker scores
    marker_scores_path = output_dir / "marker_scores.csv"
    export_marker_scores(adata, marker_scores_path, logger)
    output_paths["marker_scores"] = marker_scores_path

    # Export cluster label mapping
    mapping_path = output_dir / "cluster_label_mapping.csv"
    export_cluster_label_mapping(
        adata=adata,
        output_path=mapping_path,
        marker_scores=marker_scores,
        cluster_col=cluster_col,
        label_col=label_col,
        logger=logger,
    )
    output_paths["mapping"] = mapping_path

    logger.info("Completed annotation exports: %d files → %s", len(output_paths), output_dir)
    return output_paths


# =============================================================================
# Composition Statistics Export
# =============================================================================


def export_composition_stats(
    adata: "sc.AnnData",
    output_dir: Path,
    cell_type_col: str,
    logger: Optional[logging.Logger] = None,
) -> Dict[str, Path]:
    """Export composition statistics (global, by sample, by region, by donor).

    Parameters
    ----------
    adata : sc.AnnData
        AnnData with cell type annotations
    output_dir : Path
        Output directory
    cell_type_col : str
        Column name containing cell type labels
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    Dict[str, Path]
        Mapping from export type to output path
    """
    logger = logger or logging.getLogger(__name__)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_paths = {}

    sample_col = 'sample_id' if 'sample_id' in adata.obs.columns else None
    region_col = 'region' if 'region' in adata.obs.columns else None
    donor_col = 'donor' if 'donor' in adata.obs.columns else None

    # Global composition (matching reference format: proportion as %, total_cells column)
    logger.info('Generating global composition...')
    total_cells = len(adata)
    comp_global = adata.obs[cell_type_col].value_counts().reset_index()
    comp_global.columns = ['cell_type', 'n_cells']
    comp_global['proportion'] = round(comp_global['n_cells'] / total_cells * 100, 2)
    comp_global['total_cells'] = total_cells
    comp_global = comp_global.sort_values('n_cells', ascending=False)
    global_path = output_dir / 'composition_global_refined.csv'
    comp_global.to_csv(global_path, index=False)
    output_paths['global'] = global_path
    logger.info(f'  Wrote composition_global_refined.csv ({len(comp_global)} cell types)')

    # Composition by sample (matching reference format exactly)
    # Uses value_counts() approach to match reference order (descending by count,
    # then first-occurrence order for ties)
    if sample_col:
        logger.info('Generating composition by sample...')
        records = []
        # Iterate samples in sorted order (matches reference after sort)
        for sample_val in sorted(adata.obs[sample_col].unique()):
            mask = adata.obs[sample_col] == sample_val
            subset = adata.obs.loc[mask]
            sample_total = len(subset)
            # value_counts() returns descending by count, ties in first-occurrence order
            counts = subset[cell_type_col].value_counts()
            for cell_type, count in counts.items():
                records.append({
                    'sample_id': sample_val,
                    'cell_type': cell_type,
                    'n_cells': int(count),
                    'proportion': round(100 * count / sample_total, 2),
                    'sample_id_total': sample_total,
                })
        comp_sample_long = pd.DataFrame.from_records(records)
        sample_path = output_dir / 'composition_by_sample_refined.csv'
        comp_sample_long.to_csv(sample_path, index=False)
        output_paths['by_sample'] = sample_path
        logger.info(f'  Wrote composition_by_sample_refined.csv ({len(comp_sample_long)} rows)')

    # Composition by region (matching reference format exactly)
    if region_col:
        logger.info('Generating composition by region...')
        # Long format for internal calculations
        comp_region_long = adata.obs.groupby([region_col, cell_type_col]).size().reset_index(name='n_cells')
        comp_region_long.columns = ['region', 'cell_type', 'n_cells']
        # Add region totals
        region_totals = comp_region_long.groupby('region')['n_cells'].transform('sum')
        comp_region_long['proportion'] = round(
            comp_region_long['n_cells'] / region_totals * 100, 2
        )
        comp_region_long['region_total'] = region_totals
        # Sort by region and n_cells descending
        comp_region_long = comp_region_long.sort_values(
            ['region', 'n_cells'], ascending=[True, False]
        )
        region_path = output_dir / 'composition_by_region_refined.csv'
        comp_region_long.to_csv(region_path, index=False)
        output_paths['by_region'] = region_path
        logger.info(f'  Wrote composition_by_region_refined.csv ({len(comp_region_long)} rows)')

        # === REFERENCE FORMAT PIVOT TABLES ===
        # Use reference region order (anatomical: proximal to distal)
        region_order = ["isthmus", "ampulla", "infundibulum", "fimbriae"]
        available_regions = comp_region_long["region"].unique().tolist()
        region_order = [r for r in region_order if r in available_regions]
        # Add any regions not in the default order
        for r in available_regions:
            if r not in region_order:
                region_order.append(r)

        # === TABLE 1: Raw Counts (cell types as rows, regions as columns) ===
        counts_pivot = comp_region_long.pivot_table(
            index="cell_type",
            columns="region",
            values="n_cells",
            aggfunc="sum"
        ).fillna(0).astype(int)
        # Reorder columns
        counts_pivot = counts_pivot[[r for r in region_order if r in counts_pivot.columns]]
        # Add Total column
        counts_pivot["Total"] = counts_pivot.sum(axis=1)
        # Sort by total descending
        counts_pivot = counts_pivot.sort_values("Total", ascending=False)
        # Add TOTAL row
        col_totals = counts_pivot.sum(axis=0)
        col_totals.name = "TOTAL"
        counts_pivot = pd.concat([counts_pivot, col_totals.to_frame().T])

        counts_path = output_dir / 'composition_counts_by_region.csv'
        counts_pivot.to_csv(counts_path)
        output_paths['counts_by_region'] = counts_path
        logger.info('  Wrote composition_counts_by_region.csv (reference format)')

        # === TABLE 2: % Across Regions (row sums to 100%) ===
        across_pivot = comp_region_long.pivot_table(
            index="cell_type",
            columns="region",
            values="n_cells",
            aggfunc="sum"
        ).fillna(0)
        across_pivot = across_pivot[[r for r in region_order if r in across_pivot.columns]]
        row_totals = across_pivot.sum(axis=1)
        across_pct = across_pivot.div(row_totals, axis=0) * 100
        # Sort by total cells descending
        across_pct["_total"] = row_totals
        across_pct = across_pct.sort_values("_total", ascending=False)
        across_pct = across_pct.drop("_total", axis=1)
        across_pct = across_pct.round(1)

        across_path = output_dir / 'composition_pct_across_regions.csv'
        across_pct.to_csv(across_path)
        output_paths['pct_across_regions'] = across_path
        logger.info('  Wrote composition_pct_across_regions.csv (rows sum to 100%)')

        # === TABLE 3: % Within Regions (column sums to 100%) ===
        within_pivot = comp_region_long.pivot_table(
            index="cell_type",
            columns="region",
            values="n_cells",
            aggfunc="sum"
        ).fillna(0)
        within_pivot = within_pivot[[r for r in region_order if r in within_pivot.columns]]
        col_totals_within = within_pivot.sum(axis=0)
        within_pct = within_pivot.div(col_totals_within, axis=1) * 100
        # Sort by average proportion descending
        within_pct["_avg"] = within_pct.mean(axis=1)
        within_pct = within_pct.sort_values("_avg", ascending=False)
        within_pct = within_pct.drop("_avg", axis=1)
        within_pct = within_pct.round(2)

        within_path = output_dir / 'composition_pct_within_regions.csv'
        within_pct.to_csv(within_path)
        output_paths['pct_within_regions'] = within_path
        logger.info('  Wrote composition_pct_within_regions.csv (columns sum to 100%)')

    # Composition by donor (matching reference format exactly)
    if donor_col:
        logger.info('Generating composition by donor...')
        comp_donor = adata.obs.groupby([donor_col, cell_type_col]).size().unstack(fill_value=0)
        comp_donor_long = comp_donor.stack().reset_index()
        comp_donor_long.columns = ['donor', 'cell_type', 'n_cells']
        # Add donor totals
        donor_totals = comp_donor_long.groupby('donor')['n_cells'].transform('sum')
        comp_donor_long['proportion'] = round(
            comp_donor_long['n_cells'] / donor_totals * 100, 2
        )
        comp_donor_long['donor_total'] = donor_totals
        # Sort by donor and n_cells descending
        comp_donor_long = comp_donor_long.sort_values(
            ['donor', 'n_cells'], ascending=[True, False]
        )
        donor_path = output_dir / 'composition_by_donor_refined.csv'
        comp_donor_long.to_csv(donor_path, index=False)
        output_paths['by_donor'] = donor_path
        logger.info(f'  Wrote composition_by_donor_refined.csv ({len(comp_donor_long)} rows)')

    return output_paths


# =============================================================================
# Marker Scores Export from AnnData
# =============================================================================


def export_marker_scores_from_adata(
    adata: "sc.AnnData",
    output_dir: Path,
    logger: Optional[logging.Logger] = None,
) -> Optional[pd.DataFrame]:
    """Export marker scores from adata.uns.

    Parameters
    ----------
    adata : sc.AnnData
        AnnData with marker scores in uns
    output_dir : Path
        Output directory
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    pd.DataFrame or None
        Marker scores DataFrame if found, else None
    """
    logger = logger or logging.getLogger(__name__)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info('Extracting marker scores...')

    marker_scores_key = None
    for key in ['marker_scores_refined', 'marker_scores', 'scores']:
        if key in adata.uns:
            marker_scores_key = key
            break

    if not marker_scores_key:
        logger.warning('  No marker scores found in adata.uns')
        return None

    scores_data = adata.uns[marker_scores_key]
    if isinstance(scores_data, pd.DataFrame):
        scores_df = scores_data
    elif isinstance(scores_data, list):
        scores_df = pd.DataFrame(scores_data)
    elif isinstance(scores_data, dict):
        # Convert dict format to DataFrame
        records = []
        for cluster_id, labels in scores_data.items():
            if isinstance(labels, dict):
                for label, score in labels.items():
                    records.append({
                        'cluster_id': cluster_id,
                        'label': label,
                        'score': score,
                    })
        scores_df = pd.DataFrame(records)
    else:
        logger.warning(f'  Unexpected marker scores format: {type(scores_data)}')
        return None

    if len(scores_df) > 0:
        scores_path = output_dir / 'marker_scores.csv'
        scores_df.to_csv(scores_path, index=False)
        logger.info(f'  Wrote marker_scores.csv ({len(scores_df)} records)')
        return scores_df

    return None


# =============================================================================
# Simple Cluster Annotations Export
# =============================================================================


def export_cluster_annotations_simple(
    adata: "sc.AnnData",
    output_dir: Path,
    cell_type_col: str,
    cluster_col: str,
    scores_df: Optional[pd.DataFrame],
    logger: Optional[logging.Logger] = None,
) -> None:
    """Export simple cluster annotations with label mapping.

    Parameters
    ----------
    adata : sc.AnnData
        AnnData with cluster annotations
    output_dir : Path
        Output directory
    cell_type_col : str
        Column name containing cell type labels
    cluster_col : str
        Column name containing cluster IDs
    scores_df : pd.DataFrame, optional
        Marker scores DataFrame
    logger : logging.Logger, optional
        Logger instance
    """
    logger = logger or logging.getLogger(__name__)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info('Generating cluster annotations...')

    if cluster_col not in adata.obs.columns:
        logger.warning(f'  Cluster column {cluster_col} not found')
        return

    # Get unique cluster -> cell type mapping
    cluster_mapping = adata.obs.groupby(cluster_col)[cell_type_col].agg(
        lambda x: x.value_counts().index[0] if len(x) > 0 else 'Unknown'
    ).reset_index()
    cluster_mapping.columns = ['cluster_id', 'assigned_label']

    # Add cell counts
    cluster_counts = adata.obs[cluster_col].value_counts().reset_index()
    cluster_counts.columns = ['cluster_id', 'n_cells']
    cluster_mapping = cluster_mapping.merge(cluster_counts, on='cluster_id')

    # Add scores if available
    # IMPORTANT: Use the ASSIGNED LABEL's score, not the highest-scoring label
    # (Fixed 2025-12-30: was incorrectly using idxmax which gave highest score)
    if (scores_df is not None and
            'cluster_id' in scores_df.columns and
            'score' in scores_df.columns):
        scores_df_copy = scores_df.copy()
        scores_df_copy['cluster_id'] = scores_df_copy['cluster_id'].astype(str)
        scores_df_copy['label'] = scores_df_copy['label'].astype(str)
        cluster_mapping['cluster_id'] = cluster_mapping['cluster_id'].astype(str)

        # Build lookup: (cluster_id, label) -> score
        score_lookup = {}
        for _, row in scores_df_copy.iterrows():
            key = (row['cluster_id'], row['label'])
            score_lookup[key] = row['score']

        # Build fallback: cluster_id -> highest score
        best_idx = scores_df_copy.groupby('cluster_id')['score'].idxmax()
        fallback_lookup = dict(zip(
            scores_df_copy.loc[best_idx, 'cluster_id'],
            scores_df_copy.loc[best_idx, 'score']
        ))

        # Get score for assigned_label, fall back to highest score
        def get_score(row):
            key = (str(row['cluster_id']), str(row['assigned_label']))
            if key in score_lookup:
                return score_lookup[key]
            return fallback_lookup.get(str(row['cluster_id']), 0.0)

        cluster_mapping['assigned_score'] = cluster_mapping.apply(get_score, axis=1)

        # Set score=0.0 for Unassigned clusters
        unassigned_mask = cluster_mapping['assigned_label'] == 'Unassigned'
        cluster_mapping.loc[unassigned_mask, 'assigned_score'] = 0.0

        # Round assigned_score to 3 decimal places (ISSUE-003f fix)
        cluster_mapping['assigned_score'] = cluster_mapping['assigned_score'].apply(
            lambda x: round(float(x), 3) if pd.notna(x) else x
        )

    # Sort by cluster_id using reference sort key for consistent ordering
    cluster_mapping = cluster_mapping.sort_values(
        'cluster_id',
        key=lambda x: x.map(lambda c: _sort_cluster_key(str(c)))
    ).reset_index(drop=True)

    # Save cluster annotations (simple format)
    mapping_path = output_dir / 'cluster_annotations.csv'
    cluster_mapping.to_csv(mapping_path, index=False)
    logger.info(f'  Wrote cluster_annotations.csv ({len(cluster_mapping)} clusters)')

    # Enhanced annotations with regional distribution
    region_col = 'region' if 'region' in adata.obs.columns else None
    enhanced = cluster_mapping.copy()

    if region_col:
        # Add regional distribution
        region_dist = adata.obs.groupby([cluster_col, region_col]).size().unstack(fill_value=0)
        region_dist_pct = region_dist.div(region_dist.sum(axis=1), axis=0) * 100
        region_cols = [f'pct_{r}' for r in region_dist.columns]
        region_dist_pct.columns = region_cols
        region_dist_pct = region_dist_pct.reset_index()
        region_dist_pct[cluster_col] = region_dist_pct[cluster_col].astype(str)
        enhanced = enhanced.merge(
            region_dist_pct,
            left_on='cluster_id',
            right_on=cluster_col,
            how='left',
        )
        if cluster_col != 'cluster_id':
            enhanced = enhanced.drop(columns=[cluster_col])

    # Ensure consistent sort order (merge may change order)
    enhanced = enhanced.sort_values(
        'cluster_id',
        key=lambda x: x.map(lambda c: _sort_cluster_key(str(c)))
    ).reset_index(drop=True)

    enhanced_path = output_dir / 'cluster_annotations_enhanced.csv'
    enhanced.to_csv(enhanced_path, index=False)
    logger.info('  Wrote cluster_annotations_enhanced.csv')


# =============================================================================
# Review Summary Export
# =============================================================================


def export_review_summary(
    adata: "sc.AnnData",
    output_dir: Path,
    cell_type_col: str,
    cluster_col: str,
    comp_global: pd.DataFrame,
    workflow_state: Optional[Dict[str, Any]] = None,
    logger: Optional[logging.Logger] = None,
) -> Path:
    """Export review summary JSON with lineage groupings.

    Matches reference schema from ft/src/workflow/stage_i_unified.py:
    - generated_at, total_cells, composition_by_group, workflow_state

    Parameters
    ----------
    adata : sc.AnnData
        AnnData with cell type annotations
    output_dir : Path
        Output directory
    cell_type_col : str
        Column name containing cell type labels
    cluster_col : str
        Column name containing cluster IDs
    comp_global : pd.DataFrame
        Global composition DataFrame
    workflow_state : Dict[str, Any], optional
        Workflow state dictionary to embed (matches reference schema)
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    Path
        Path to the exported JSON file
    """
    logger = logger or logging.getLogger(__name__)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info('Generating review summary...')

    # Group composition by lineage (matching reference grouping logic)
    # Reference uses simple keyword matching with "Other" catch-all
    composition_by_group: Dict[str, Dict[str, int]] = {}

    for cell_type, count in adata.obs[cell_type_col].value_counts().items():
        cell_type_str = str(cell_type)
        cell_type_lower = cell_type_str.lower()
        assigned = False

        # Determine group from cell type name (matching reference logic)
        group = None
        if "epithelium" in cell_type_lower or "epithelial" in cell_type_lower:
            group = "Epithelium"
        elif "endothelium" in cell_type_lower or "endothelial" in cell_type_lower:
            group = "Endothelium"
        elif ("mesenchym" in cell_type_lower or "stromal" in cell_type_lower or
              "smooth muscle" in cell_type_lower):
            group = "Mesenchymal Cells"
        elif ("immune" in cell_type_lower or "macrophage" in cell_type_lower or
              "t cell" in cell_type_lower):
            group = "Immune Cells"
        elif "unassigned" in cell_type_lower:
            group = "Unassigned"
        elif "~" in cell_type_str:
            group = "Hybrids"

        if group:
            if group not in composition_by_group:
                composition_by_group[group] = {}
            composition_by_group[group][cell_type_str] = int(count)
            assigned = True

        if not assigned:
            # "Other" catch-all for unmatched cell types (matches reference)
            if "Other" not in composition_by_group:
                composition_by_group["Other"] = {}
            composition_by_group["Other"][cell_type_str] = int(count)

    # Build review summary matching reference schema exactly
    # Reference: ft/src/workflow/stage_i_unified.py lines 1103-1111
    review_summary: Dict[str, Any] = {
        'generated_at': datetime.now().isoformat(),
        'total_cells': int(adata.n_obs),
        'composition_by_group': composition_by_group,
    }

    # Embed workflow_state if provided (matches reference schema)
    if workflow_state is not None:
        review_summary['workflow_state'] = workflow_state

    summary_path = output_dir / 'review_summary.json'
    with open(summary_path, 'w') as f:
        json.dump(review_summary, f, indent=2)
    logger.info('  Wrote review_summary.json')

    return summary_path


# =============================================================================
# Workflow State Export
# =============================================================================


def export_workflow_state(
    output_dir: Path,
    marker_map_path: Optional[str] = None,
    input_path: Optional[str] = None,
    stage_h_dir: Optional[str] = None,
    # Annotate step
    annotate_complete: bool = False,
    annotate_timestamp: Optional[str] = None,
    annotate_output: Optional[str] = None,
    # Diagnose step
    diagnose_complete: bool = False,
    diagnose_timestamp: Optional[str] = None,
    diagnostic_summary: Optional[Dict[str, int]] = None,
    diagnostic_report_path: Optional[str] = None,
    # Execute step
    execute_complete: bool = False,
    execute_timestamp: Optional[str] = None,
    groups_executed: Optional[List[str]] = None,
    groups_skipped: Optional[List[str]] = None,
    final_output_path: Optional[str] = None,
    # Review step
    review_complete: bool = True,
    review_timestamp: Optional[str] = None,
    # Iteration tracking
    iterations: Optional[List[Dict[str, Any]]] = None,
    current_iteration: int = 0,
    stopping_reason: Optional[str] = None,
    # Created/updated timestamps
    created_at: Optional[str] = None,
    updated_at: Optional[str] = None,
    logger: Optional[logging.Logger] = None,
) -> Path:
    """Export workflow state YAML for tracking.

    Matches reference schema from ft/src/workflow/state.py WorkflowState class.
    Contains all 18 fields for full compatibility.

    Parameters
    ----------
    output_dir : Path
        Output directory
    marker_map_path : str, optional
        Path to marker map
    input_path : str, optional
        Path to input h5ad
    stage_h_dir : str, optional
        Path to Stage H directory
    annotate_complete : bool
        Whether annotate step is complete
    annotate_timestamp : str, optional
        When annotate completed
    annotate_output : str, optional
        Path to annotate output
    diagnose_complete : bool
        Whether diagnose step is complete
    diagnose_timestamp : str, optional
        When diagnose completed
    diagnostic_summary : Dict[str, int], optional
        Summary counts {SUBCLUSTER: n, RELABEL: n, SKIP: n}
    diagnostic_report_path : str, optional
        Path to diagnostic_report.csv
    execute_complete : bool
        Whether execute step is complete
    execute_timestamp : str, optional
        When execute completed
    groups_executed : List[str], optional
        List of group names that were executed
    groups_skipped : List[str], optional
        List of group names that were skipped
    final_output_path : str, optional
        Path to refined_final.h5ad
    review_complete : bool
        Whether review step is complete
    review_timestamp : str, optional
        When review completed
    iterations : List[Dict], optional
        List of iteration state dicts
    current_iteration : int
        Current iteration number
    stopping_reason : str, optional
        Reason for stopping iterations
    created_at : str, optional
        Timestamp when state was created
    updated_at : str, optional
        Timestamp when state was last updated
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    Path
        Path to the exported YAML file
    """
    import yaml

    logger = logger or logging.getLogger(__name__)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    now = datetime.now().isoformat()

    # Build workflow state matching reference schema exactly
    # Reference: ft/src/workflow/state.py WorkflowState dataclass
    workflow_state = {
        # Config
        'version': '1.0',
        'marker_map_path': marker_map_path or '',
        'input_path': input_path or '',
        'stage_h_dir': stage_h_dir or '',
        'out_dir': str(output_dir),
        'created_at': created_at or now,
        'updated_at': updated_at or now,
        # Annotate step
        'annotate_complete': annotate_complete,
        'annotate_timestamp': annotate_timestamp,
        'annotate_output': annotate_output,
        # Diagnose step
        'diagnose_complete': diagnose_complete,
        'diagnose_timestamp': diagnose_timestamp or now,
        'diagnostic_summary': diagnostic_summary or {},
        'diagnostic_report_path': diagnostic_report_path,
        # Execute step
        'execute_complete': execute_complete,
        'execute_timestamp': execute_timestamp,
        'groups_executed': groups_executed or [],
        'groups_skipped': groups_skipped or [],
        'final_output_path': final_output_path,
        # Review step
        'review_complete': review_complete,
        'review_timestamp': review_timestamp or now,
        # Iteration tracking
        'iterations': iterations or [],
        'current_iteration': current_iteration,
        'stopping_reason': stopping_reason,
    }

    state_path = output_dir / 'workflow_state.yaml'
    with open(state_path, 'w') as f:
        yaml.safe_dump(workflow_state, f, default_flow_style=False, sort_keys=False)
    logger.info('  Wrote workflow_state.yaml')

    return state_path


# =============================================================================
# 44-Column Enhanced Annotations Format
# =============================================================================


def _map_col_name_44(col: str) -> str:
    """Map enhanced column name to cluster_ann column name."""
    mapping = {
        'cell_type': 'assigned_label',
        'score': 'assigned_score',
        'assigned_path': 'assigned_path',
        'assigned_level': 'assigned_level',
        'root_label': 'root_label',
        'confidence': 'confidence',
        'min_margin_along_path': 'min_margin_along_path',
        'margin_is_infinite': 'margin_is_infinite',
        'stop_reason': 'stop_reason',
        'stopped_before_leaf': 'stopped_before_leaf',
        'composition': 'composition',
        'decision_trace': 'decision_trace',
        'coverage': 'coverage',
        'resolved_markers': 'resolved_markers',
        'is_ambiguous_root': 'is_ambiguous_root',
        'ambiguous_root_candidates': 'ambiguous_root_candidates',
    }
    return mapping.get(col, col)


def _build_reason_string(row: pd.Series) -> str:
    """Build reason string from row data."""
    label = row.get('assigned_label', '')
    score = row.get('assigned_score', 0)
    margin = row.get('min_margin_along_path', 0)

    # Handle None/NaN values
    if score is None or (isinstance(score, float) and pd.isna(score)):
        score = 0
    if margin is None or (isinstance(margin, float) and pd.isna(margin)):
        margin = 0

    if not label or label == 'Unassigned':
        return ''

    # Check if it's a path (contains ' / ')
    path = row.get('assigned_path', '')
    if ' / ' in str(path):
        parts = str(path).split(' / ')
        if len(parts) >= 2:
            return f"{parts[0]} -> {label} (score={score:.2f}, margin={margin:.2f})"

    return f"{label} (score={score:.2f})"


def generate_enhanced_annotations_44col(
    cluster_ann: pd.DataFrame,
    scores_df: Optional[pd.DataFrame],
    stage_h_dir: Optional[str] = None,
    iteration: int = 1,
    logger: Optional[logging.Logger] = None,
    decision_steps: Optional[pd.DataFrame] = None,
    cluster_annotations: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    """Generate 44-column enhanced annotations matching reference format.

    Creates hierarchical annotations with:
    - 6 basic columns: cluster_id, origin_cluster, iteration_created, n_cells, proportion, reason
    - 16 lvl0 columns: Stage H annotations (or empty if unavailable)
    - 16 lvl1 columns: Current iteration annotations
    - 3 quality metrics: mean_enrichment, mean_positive_fraction, confidence_level
    - 3 runner-up columns: runner_up_label, runner_up_score, gap

    Parameters
    ----------
    cluster_ann : pd.DataFrame
        Cluster annotations in 24-column format (from cluster_annotations_subcluster)
    scores_df : pd.DataFrame, optional
        Marker scores for runner-up calculation
    stage_h_dir : str, optional
        Path to Stage H directory for lvl0 data and decision_steps.csv
    iteration : int
        Current iteration number (default: 1)
    logger : logging.Logger, optional
        Logger instance
    decision_steps : pd.DataFrame, optional
        Combined decision steps (Stage H + current iteration) for hierarchical runner-up.
        If not provided, will try to load from stage_h_dir/decision_steps.csv
    cluster_annotations : pd.DataFrame, optional
        Combined cluster annotations for hierarchical runner-up lookup.
        If not provided, uses cluster_ann parameter.

    Returns
    -------
    pd.DataFrame
        44-column enhanced annotations DataFrame
    """
    logger = logger or logging.getLogger(__name__)

    # Use cluster_ann as cluster_annotations if not provided
    if cluster_annotations is None:
        cluster_annotations = cluster_ann

    # Hierarchical columns that exist per level
    HIER_COLS = [
        'cell_type', 'score', 'assigned_path', 'assigned_level',
        'root_label', 'confidence', 'min_margin_along_path', 'margin_is_infinite',
        'stop_reason', 'stopped_before_leaf', 'composition', 'decision_trace',
        'coverage', 'resolved_markers', 'is_ambiguous_root', 'ambiguous_root_candidates',
    ]

    # Load Stage H annotations for lvl0 data if available
    stage_h_lookup = {}
    stage_h_annotations_df = pd.DataFrame()
    stage_h_decision_steps_df = pd.DataFrame()
    if stage_h_dir:
        stage_h_path = Path(stage_h_dir) / 'cluster_annotations.csv'
        if stage_h_path.exists():
            try:
                stage_h_annotations_df = pd.read_csv(stage_h_path)
                for _, row in stage_h_annotations_df.iterrows():
                    cid = str(row.get('cluster_id', ''))
                    stage_h_lookup[cid] = row.to_dict()
                logger.info(f'  Loaded {len(stage_h_lookup)} Stage H annotations for lvl0')
            except Exception as e:
                logger.warning(f'  Failed to load Stage H annotations: {e}')

        # Load Stage H decision_steps for hierarchical runner-up calculation
        stage_h_steps_path = Path(stage_h_dir) / 'decision_steps.csv'
        if stage_h_steps_path.exists():
            try:
                stage_h_decision_steps_df = pd.read_csv(stage_h_steps_path)
                logger.info(f'  Loaded {len(stage_h_decision_steps_df)} Stage H decision steps')
            except Exception as e:
                logger.warning(f'  Failed to load Stage H decision steps: {e}')

    # Combine Stage H and current iteration data for hierarchical runner-up
    combined_annotations = pd.concat(
        [df for df in [stage_h_annotations_df, cluster_annotations] if not df.empty],
        ignore_index=True
    ) if not stage_h_annotations_df.empty or not cluster_annotations.empty else cluster_annotations

    # Combine decision_steps if provided
    if decision_steps is not None and not decision_steps.empty:
        combined_decision_steps = pd.concat(
            [df for df in [stage_h_decision_steps_df, decision_steps] if not df.empty],
            ignore_index=True
        ) if not stage_h_decision_steps_df.empty else decision_steps
    else:
        combined_decision_steps = stage_h_decision_steps_df

    records = []
    for _, row in cluster_ann.iterrows():
        cluster_id = str(row.get('cluster_id', ''))
        origin_cluster = str(row.get(
            'origin_cluster',
            cluster_id.split(':')[0] if ':' in cluster_id else cluster_id
        ))
        iteration_created = int(row.get('iteration_created', cluster_id.count(':')))

        # Build base record
        record = {
            'cluster_id': cluster_id,
            'origin_cluster': origin_cluster,
            'iteration_created': iteration_created,
            'n_cells': row.get('n_cells', 0),
            'proportion': row.get('proportion', 0),
            'reason': _build_reason_string(row),
        }

        # Add lvl0 columns (from Stage H or explicit defaults)
        # ISSUE-003l fix: use explicit defaults matching reference (-1, False)
        h_data = stage_h_lookup.get(origin_cluster, {})
        for col in HIER_COLS:
            src_col = _map_col_name_44(col)
            val = h_data.get(src_col, '')
            # Handle special cases
            if col == 'cell_type':
                val = h_data.get('assigned_label', '')
            elif col == 'score':
                val = h_data.get('assigned_score', 0.0)

            # Apply explicit defaults for missing values (matches reference)
            if val == '' or (isinstance(val, float) and pd.isna(val)):
                if col == 'assigned_level':
                    val = -1
                elif col in ('margin_is_infinite', 'stopped_before_leaf', 'is_ambiguous_root'):
                    val = False
                elif 'score' in col or 'confidence' in col or 'coverage' in col:
                    val = 0.0
                else:
                    val = ''
            record[f'{col}_lvl0'] = val

        # Add lvl1 columns (from current cluster_ann)
        for col in HIER_COLS:
            src_col = _map_col_name_44(col)
            val = row.get(src_col, '')
            # Handle special cases
            if col == 'cell_type':
                val = row.get('assigned_label', '')
            elif col == 'score':
                val = row.get('assigned_score', 0.0)
            record[f'{col}_lvl{iteration}'] = val if val != '' else (
                0.0 if 'score' in col or 'confidence' in col or 'coverage' in col else ''
            )

        # Add quality metrics
        record['mean_enrichment'] = round(
            float(row.get('mean_enrichment', 0) or 0), 3
        )
        record['mean_positive_fraction'] = round(
            float(row.get('mean_positive_fraction', 0) or 0), 3
        )
        record['confidence_level'] = row.get('confidence_level', '')

        # Add runner-up columns using hierarchical sibling-based algorithm
        # ISSUE-003j fix: use get_hierarchical_runner_up instead of global top-2
        runner_up_label, runner_up_score, gap = get_hierarchical_runner_up(
            cluster_id=cluster_id,
            decision_steps=combined_decision_steps,
            cluster_annotations=combined_annotations,
            marker_scores=scores_df,
        )
        record['runner_up_label'] = runner_up_label
        record['runner_up_score'] = round(runner_up_score, 3) if not pd.isna(runner_up_score) else 0.0
        record['gap'] = round(gap, 3) if not pd.isna(gap) else float('nan')

        records.append(record)

    df = pd.DataFrame.from_records(records)

    # Ensure column order matches reference (44 columns)
    ref_col_order = [
        'cluster_id', 'origin_cluster', 'iteration_created', 'n_cells', 'proportion', 'reason',
    ]
    # Add lvl0 columns
    for col in HIER_COLS:
        ref_col_order.append(f'{col}_lvl0')
    # Add lvl1 columns
    for col in HIER_COLS:
        ref_col_order.append(f'{col}_lvl{iteration}')
    # Add quality and runner-up columns
    ref_col_order.extend([
        'mean_enrichment', 'mean_positive_fraction', 'confidence_level',
        'runner_up_label', 'runner_up_score', 'gap',
    ])

    # Reorder columns, only keeping those that exist
    final_cols = [c for c in ref_col_order if c in df.columns]
    df = df[final_cols]

    return df


# =============================================================================
# High-Level Review Exports Orchestrator
# =============================================================================


def run_review_exports(
    adata: "sc.AnnData",
    output_dir: Path,
    cell_type_col: Optional[str] = None,
    cluster_col: Optional[str] = None,
    stage_h_dir: Optional[str] = None,
    marker_map_path: Optional[str] = None,
    input_path: Optional[str] = None,
    # Workflow state parameters (for reference schema compatibility)
    diagnose_complete: bool = True,
    diagnostic_summary: Optional[Dict[str, int]] = None,
    diagnostic_report_path: Optional[str] = None,
    execute_complete: bool = True,
    groups_executed: Optional[List[str]] = None,
    groups_skipped: Optional[List[str]] = None,
    final_output_path: Optional[str] = None,
    logger: Optional[logging.Logger] = None,
) -> Dict[str, Path]:
    """Run all review exports (composition, annotations, summary, workflow state).

    This is the high-level orchestrator that combines all review export functions.
    Exports match reference schema from ft/src/workflow/stage_i_unified.py.

    Parameters
    ----------
    adata : sc.AnnData
        AnnData with annotations
    output_dir : Path
        Output directory
    cell_type_col : str, optional
        Cell type column name (auto-detected if not specified)
    cluster_col : str, optional
        Cluster column name (auto-detected if not specified)
    stage_h_dir : str, optional
        Path to Stage H directory for workflow state
    marker_map_path : str, optional
        Path to marker map for workflow state
    input_path : str, optional
        Path to input h5ad for workflow state
    diagnose_complete : bool
        Whether diagnose step is complete (default: True)
    diagnostic_summary : Dict[str, int], optional
        Summary counts {SUBCLUSTER: n, RELABEL: n, SKIP: n}
    diagnostic_report_path : str, optional
        Path to diagnostic_report.csv
    execute_complete : bool
        Whether execute step is complete (default: True)
    groups_executed : List[str], optional
        List of group names that were executed
    groups_skipped : List[str], optional
        List of group names that were skipped
    final_output_path : str, optional
        Path to refined_final.h5ad
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    Dict[str, Path]
        Mapping from export type to output path
    """
    logger = logger or logging.getLogger(__name__)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    output_paths = {}

    # Detect columns if not specified
    if cell_type_col is None:
        cell_type_col = detect_cell_type_column(adata)
    if cluster_col is None:
        cluster_col = detect_cluster_column(adata)

    logger.info(f'Using cell type column: {cell_type_col}')
    logger.info(f'Using cluster column: {cluster_col}')

    # Export composition statistics
    comp_paths = export_composition_stats(adata, output_dir, cell_type_col, logger)
    output_paths.update(comp_paths)

    # Export marker scores
    scores_df = export_marker_scores_from_adata(adata, output_dir, logger)

    # Export cluster annotations from adata.uns if available (fast path)
    cluster_ann_key = 'cluster_annotations_subcluster'
    decision_steps_key = 'decision_steps_subcluster'

    # Extract decision_steps for hierarchical runner-up calculation (ISSUE-003j fix)
    decision_steps_df = None
    if decision_steps_key in adata.uns:
        decision_steps_df = adata.uns[decision_steps_key]
        if isinstance(decision_steps_df, dict):
            decision_steps_df = pd.DataFrame(decision_steps_df)
        logger.info(f'  Loaded {len(decision_steps_df)} decision steps from adata.uns')

    if cluster_ann_key in adata.uns:
        logger.info('Generating cluster annotations from stored data...')
        cluster_ann = adata.uns[cluster_ann_key].copy()
        total_cells = len(adata)

        # Add missing columns to match reference 24-column format
        if 'origin_cluster' not in cluster_ann.columns:
            cluster_ann['origin_cluster'] = cluster_ann['cluster_id'].astype(str).apply(
                lambda x: x.split(':')[0] if ':' in x else x
            )
        if 'iteration_created' not in cluster_ann.columns:
            cluster_ann['iteration_created'] = cluster_ann['cluster_id'].astype(str).apply(
                lambda x: x.count(':')
            )
        if 'proportion' not in cluster_ann.columns:
            cluster_ann['proportion'] = round(
                cluster_ann['n_cells'] / total_cells * 100, 2
            )

        # Add mean_enrichment and mean_positive_fraction from marker_scores
        # IMPORTANT: Use the ASSIGNED LABEL's metrics, not the highest-scoring label
        # (Fixed 2025-12-30: was incorrectly using idxmax which gave highest score's metrics)
        if scores_df is not None and not scores_df.empty:
            if 'cluster_id' in scores_df.columns and 'score' in scores_df.columns:
                # Build lookup: (cluster_id, label) -> {mean_enrichment, mean_positive_fraction}
                scores_df_copy = scores_df.copy()
                scores_df_copy['cluster_id'] = scores_df_copy['cluster_id'].astype(str)
                scores_df_copy['label'] = scores_df_copy['label'].astype(str)

                # Create lookup dictionary keyed by (cluster_id, label)
                # Note: Convert to Python float to ensure consistent round() behavior
                # (np.float64 and Python float have different rounding at boundaries like 0.8875)
                metrics_lookup = {}
                for _, row in scores_df_copy.iterrows():
                    key = (row['cluster_id'], row['label'])
                    metrics_lookup[key] = {
                        'mean_enrichment': float(row.get('mean_enrichment', 0) or 0),
                        'mean_positive_fraction': float(row.get('mean_positive_fraction', 0) or 0),
                    }

                # Also build fallback: cluster_id -> highest score's metrics
                best_idx = scores_df_copy.groupby('cluster_id')['score'].idxmax()
                fallback_lookup = {}
                for idx in best_idx:
                    row = scores_df_copy.loc[idx]
                    fallback_lookup[row['cluster_id']] = {
                        'mean_enrichment': float(row.get('mean_enrichment', 0) or 0),
                        'mean_positive_fraction': float(row.get('mean_positive_fraction', 0) or 0),
                    }

                # Apply: use assigned_label's metrics, fall back to highest score
                if 'mean_enrichment' not in cluster_ann.columns or 'mean_positive_fraction' not in cluster_ann.columns:
                    me_values = []
                    mpf_values = []
                    for _, row in cluster_ann.iterrows():
                        cid = str(row['cluster_id'])
                        label = str(row.get('assigned_label', ''))
                        key = (cid, label)

                        if key in metrics_lookup:
                            # Use assigned label's metrics (correct behavior)
                            me_values.append(metrics_lookup[key]['mean_enrichment'])
                            mpf_values.append(metrics_lookup[key]['mean_positive_fraction'])
                        elif cid in fallback_lookup:
                            # Fallback to highest score (for Unassigned or missing labels)
                            me_values.append(fallback_lookup[cid]['mean_enrichment'])
                            mpf_values.append(fallback_lookup[cid]['mean_positive_fraction'])
                        else:
                            me_values.append(0)
                            mpf_values.append(0)

                    if 'mean_enrichment' not in cluster_ann.columns:
                        cluster_ann['mean_enrichment'] = [round(v, 3) if pd.notna(v) else 0 for v in me_values]
                    if 'mean_positive_fraction' not in cluster_ann.columns:
                        cluster_ann['mean_positive_fraction'] = [round(v, 3) if pd.notna(v) else 0 for v in mpf_values]

        # Add confidence_level if not present
        # ISSUE-003i fix: round score to 3 decimals before classifying to match reference behavior
        if 'confidence_level' not in cluster_ann.columns:
            cluster_ann['confidence_level'] = cluster_ann['assigned_score'].apply(
                lambda x: _classify_confidence_band(round(float(x), 3) if pd.notna(x) else 0)
            )

        # Reorder columns to match reference format
        ref_cols = [
            'cluster_id', 'origin_cluster', 'iteration_created', 'proportion',
            'mean_enrichment', 'mean_positive_fraction', 'n_cells', 'assigned_label',
            'assigned_path', 'assigned_level', 'assigned_score', 'root_label',
            'confidence', 'min_margin_along_path', 'margin_is_infinite', 'stop_reason',
            'stopped_before_leaf', 'composition', 'decision_trace', 'coverage',
            'resolved_markers', 'is_ambiguous_root', 'ambiguous_root_candidates',
            'confidence_level',
        ]
        final_cols = [c for c in ref_cols if c in cluster_ann.columns]
        cluster_ann = cluster_ann[final_cols]

        # Round numeric columns to 3 decimal places to match reference precision
        # (ISSUE-003f fix: reference rounds these columns for consistent output)
        numeric_cols_to_round = [
            'assigned_score', 'confidence', 'min_margin_along_path', 'coverage',
            'mean_enrichment', 'mean_positive_fraction',
        ]
        for col in numeric_cols_to_round:
            if col in cluster_ann.columns:
                cluster_ann[col] = cluster_ann[col].apply(
                    lambda x: round(float(x), 3) if pd.notna(x) else x
                )

        # Fix NaN representation to match reference (ISSUE-003h)
        # 1. For Unassigned: assigned_path should be "Unassigned" not NaN
        if 'assigned_path' in cluster_ann.columns and 'assigned_label' in cluster_ann.columns:
            unassigned_mask = cluster_ann['assigned_label'] == 'Unassigned'
            cluster_ann.loc[unassigned_mask, 'assigned_path'] = 'Unassigned'

        # 2. For Unassigned: min_margin_along_path should be NaN not 0.0
        if 'min_margin_along_path' in cluster_ann.columns and 'assigned_label' in cluster_ann.columns:
            unassigned_mask = cluster_ann['assigned_label'] == 'Unassigned'
            cluster_ann.loc[unassigned_mask, 'min_margin_along_path'] = float('nan')

        # 3. For root-level (assigned_label == root_label): assigned_level should be -1 not 0
        if 'assigned_level' in cluster_ann.columns and 'assigned_label' in cluster_ann.columns and 'root_label' in cluster_ann.columns:
            root_level_mask = (cluster_ann['assigned_label'] == cluster_ann['root_label']) & (cluster_ann['assigned_level'] == 0)
            cluster_ann.loc[root_level_mask, 'assigned_level'] = -1

        # Sort by cluster_id using reference sort key for consistent ordering
        # (Handles subclusters: '1:0' -> (1, 0), deeper subclusters at end)
        cluster_ann = cluster_ann.sort_values(
            'cluster_id',
            key=lambda x: x.map(lambda c: _sort_cluster_key(str(c)))
        ).reset_index(drop=True)

        # Save cluster_annotations.csv (24 columns)
        stage_h_path = output_dir / 'cluster_annotations.csv'
        cluster_ann.to_csv(stage_h_path, index=False)
        output_paths['cluster_annotations'] = stage_h_path
        logger.info(
            f'  Wrote cluster_annotations.csv ({len(cluster_ann)} clusters, '
            f'{len(final_cols)} columns)'
        )

        # Generate enhanced version with 44-column format
        # Pass decision_steps for hierarchical runner-up calculation (ISSUE-003j fix)
        enhanced_df = generate_enhanced_annotations_44col(
            cluster_ann=cluster_ann,
            scores_df=scores_df,
            stage_h_dir=stage_h_dir,
            iteration=1,
            logger=logger,
            decision_steps=decision_steps_df,
            cluster_annotations=cluster_ann,
        )
        # Ensure consistent sort order
        if 'cluster_id' in enhanced_df.columns:
            enhanced_df = enhanced_df.sort_values(
                'cluster_id',
                key=lambda x: x.map(lambda c: _sort_cluster_key(str(c)))
            ).reset_index(drop=True)
        enhanced_path = output_dir / 'cluster_annotations_enhanced.csv'
        enhanced_df.to_csv(enhanced_path, index=False)
        output_paths['cluster_annotations_enhanced'] = enhanced_path
        logger.info(
            f'  Wrote cluster_annotations_enhanced.csv ({len(enhanced_df.columns)} columns)'
        )
    else:
        # Fallback: use simple export
        logger.info('Generating cluster annotations (simple format)...')
        export_cluster_annotations_simple(
            adata, output_dir, cell_type_col, cluster_col, scores_df, logger
        )
        output_paths['cluster_annotations'] = output_dir / 'cluster_annotations.csv'
        output_paths['cluster_annotations_enhanced'] = (
            output_dir / 'cluster_annotations_enhanced.csv'
        )

    # Build workflow state dict (matches reference schema)
    now = datetime.now().isoformat()
    workflow_state_dict = {
        'version': '1.0',
        'marker_map_path': marker_map_path or '',
        'input_path': input_path or '',
        'stage_h_dir': stage_h_dir or '',
        'out_dir': str(output_dir),
        'created_at': now,
        'updated_at': now,
        # Annotate step (not tracked in review exports, set to False)
        'annotate_complete': False,
        'annotate_timestamp': None,
        'annotate_output': None,
        # Diagnose step
        'diagnose_complete': diagnose_complete,
        'diagnose_timestamp': now,
        'diagnostic_summary': diagnostic_summary or {},
        'diagnostic_report_path': diagnostic_report_path,
        # Execute step
        'execute_complete': execute_complete,
        'execute_timestamp': now,
        'groups_executed': groups_executed or [],
        'groups_skipped': groups_skipped or [],
        'final_output_path': final_output_path,
        # Review step
        'review_complete': True,
        'review_timestamp': now,
        # Iteration tracking
        'iterations': [],
        'current_iteration': 0,
        'stopping_reason': None,
    }

    # Export review summary (with embedded workflow_state)
    comp_global = adata.obs[cell_type_col].value_counts()
    summary_path = export_review_summary(
        adata, output_dir, cell_type_col, cluster_col, comp_global,
        workflow_state=workflow_state_dict,
        logger=logger,
    )
    output_paths['review_summary'] = summary_path

    # Export workflow state (with all fields)
    state_path = export_workflow_state(
        output_dir,
        marker_map_path=marker_map_path,
        input_path=input_path,
        stage_h_dir=stage_h_dir,
        diagnose_complete=diagnose_complete,
        diagnose_timestamp=now,
        diagnostic_summary=diagnostic_summary,
        diagnostic_report_path=diagnostic_report_path,
        execute_complete=execute_complete,
        execute_timestamp=now,
        groups_executed=groups_executed,
        groups_skipped=groups_skipped,
        final_output_path=final_output_path,
        review_complete=True,
        review_timestamp=now,
        logger=logger,
    )
    output_paths['workflow_state'] = state_path

    logger.info('')
    logger.info('Review export complete!')
    logger.info(f'  Total cells: {adata.n_obs:,}')
    logger.info(f'  Cell types: {len(comp_global)}')
    logger.info(f'  Output files: {len(output_paths)}')

    return output_paths
