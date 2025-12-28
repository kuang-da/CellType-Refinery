"""
RefinePlan and RefineOp dataclasses for cell-type annotation refinement.

This module defines the data structures for representing refinement operations
and complete refinement plans that can be serialized to YAML/JSON.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field, asdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Union

import yaml


# =============================================================================
# RefineOp Base and Subclasses
# =============================================================================

@dataclass
class RefineOp:
    """Base class for refinement operations.

    Attributes:
        op_type: Type of operation (subcluster, merge, override, relabel, rescore)
        cluster_id: Target cluster ID(s) - single string or list for merge ops
        reason: Human-readable reason for this operation
        source: Origin policy ("auto" or "manual")
        priority: Ordering priority (higher = later execution)
    """
    op_type: Literal["subcluster", "merge", "override", "relabel", "rescore"]
    cluster_id: Union[str, List[str]]
    reason: str
    source: Literal["auto", "manual"] = "manual"
    priority: int = 0

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return asdict(self)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "RefineOp":
        """Create RefineOp from dictionary.

        Routes to appropriate subclass based on op_type.
        """
        op_type = data.get("op_type")
        if op_type == "subcluster":
            return SubclusterOp.from_dict(data)
        elif op_type == "merge":
            return MergeOp.from_dict(data)
        elif op_type == "override":
            return OverrideOp.from_dict(data)
        elif op_type == "relabel":
            return RelabelOp.from_dict(data)
        elif op_type == "rescore":
            return RescoreOp.from_dict(data)
        else:
            raise ValueError(f"Unknown op_type: {op_type}")


@dataclass
class SubclusterOp(RefineOp):
    """Subclustering operation to split a cluster into finer subclusters.

    Creates hierarchical labels like "3:0", "3:1", "3:2" from parent cluster "3".

    Attributes:
        resolution: Leiden clustering resolution (lower = fewer clusters)
        n_pcs: Number of principal components for PCA
        neighbors_k: Number of neighbors for k-NN graph
        focus_markers: Optional list of markers to focus on (subset features)
        min_cells: Minimum cells required to proceed with subclustering
    """
    op_type: Literal["subcluster"] = field(default="subcluster", init=False)
    resolution: float = 0.4
    n_pcs: int = 30
    neighbors_k: int = 15
    focus_markers: Optional[List[str]] = None
    min_cells: int = 100

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "SubclusterOp":
        """Create SubclusterOp from dictionary."""
        return cls(
            cluster_id=data["cluster_id"],
            reason=data.get("reason", ""),
            source=data.get("source", "manual"),
            priority=data.get("priority", 0),
            resolution=data.get("resolution", 0.4),
            n_pcs=data.get("n_pcs", 30),
            neighbors_k=data.get("neighbors_k", 15),
            focus_markers=data.get("focus_markers"),
            min_cells=data.get("min_cells", 100),
        )


@dataclass
class MergeOp(RefineOp):
    """Merge operation to combine multiple clusters into one.

    Attributes:
        source_clusters: List of cluster IDs to merge
        target_label: New label for the merged cluster
    """
    op_type: Literal["merge"] = field(default="merge", init=False)
    source_clusters: List[str] = field(default_factory=list)
    target_label: str = ""

    def __post_init__(self):
        # For merge ops, cluster_id should be the list of source clusters
        if not self.source_clusters and isinstance(self.cluster_id, list):
            self.source_clusters = self.cluster_id
        elif self.source_clusters:
            self.cluster_id = self.source_clusters

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "MergeOp":
        """Create MergeOp from dictionary."""
        source_clusters = data.get("source_clusters", [])
        if not source_clusters and isinstance(data.get("cluster_id"), list):
            source_clusters = data["cluster_id"]
        return cls(
            cluster_id=source_clusters,
            reason=data.get("reason", ""),
            source=data.get("source", "manual"),
            priority=data.get("priority", 0),
            source_clusters=source_clusters,
            target_label=data.get("target_label", ""),
        )


@dataclass
class OverrideOp(RefineOp):
    """Override operation to directly assign a cell type label to a cluster.

    Attributes:
        new_label: New cell type label to assign
    """
    op_type: Literal["override"] = field(default="override", init=False)
    new_label: str = ""

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "OverrideOp":
        """Create OverrideOp from dictionary."""
        return cls(
            cluster_id=data["cluster_id"],
            reason=data.get("reason", ""),
            source=data.get("source", "manual"),
            priority=data.get("priority", 0),
            new_label=data.get("new_label", data.get("cell_type", "")),
        )


@dataclass
class RelabelOp(RefineOp):
    """Relabel operation for instant relabeling without clustering.

    Used when a homogeneous parent cluster can be directly relabeled to
    a child type without subclustering (from automatic logic).

    Attributes:
        new_label: New cell type label
        old_label: Previous label being replaced
        confidence_score: Score supporting this relabeling decision
    """
    op_type: Literal["relabel"] = field(default="relabel", init=False)
    new_label: str = ""
    old_label: str = ""
    confidence_score: float = 0.0

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "RelabelOp":
        """Create RelabelOp from dictionary."""
        return cls(
            cluster_id=data["cluster_id"],
            reason=data.get("reason", ""),
            source=data.get("source", "auto"),
            priority=data.get("priority", 0),
            new_label=data.get("new_label", ""),
            old_label=data.get("old_label", ""),
            confidence_score=data.get("confidence_score", 0.0),
        )


@dataclass
class RescoreOp(RefineOp):
    """Rescore operation to recompute marker scores.

    Attributes:
        mode: Rescoring mode - "smart" (reuse unchanged) or "full" (recompute all)
        recompute_de: Whether to recompute differential expression tests
        target_clusters: Optional list of specific clusters to rescore (None = all modified)
    """
    op_type: Literal["rescore"] = field(default="rescore", init=False)
    mode: Literal["smart", "full"] = "smart"
    recompute_de: bool = False
    target_clusters: Optional[List[str]] = None

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "RescoreOp":
        """Create RescoreOp from dictionary."""
        return cls(
            cluster_id=data.get("cluster_id", "all"),
            reason=data.get("reason", "Rescore after modifications"),
            source=data.get("source", "auto"),
            priority=data.get("priority", 100),  # Rescore happens last
            mode=data.get("mode", "smart"),
            recompute_de=data.get("recompute_de", False),
            target_clusters=data.get("target_clusters"),
        )


# =============================================================================
# RefinePlan
# =============================================================================

@dataclass
class RefinePlan:
    """Complete refinement plan with ordered operations.

    A RefinePlan contains a list of RefineOp operations to be executed
    on an AnnData object. Plans can be serialized to YAML/JSON for
    inspection, auditing, and replay.

    Attributes:
        version: Plan schema version for compatibility checking
        operations: List of RefineOp operations to execute
        metadata: Additional metadata (timestamp, source files, etc.)

    Example:
        >>> plan = RefinePlan()
        >>> plan.add_override("3", "Ciliated_Epithelium", "High FOXJ1 expression")
        >>> plan.add_subcluster("5", resolution=0.3, reason="M1/M2 heterogeneity")
        >>> plan.to_yaml(Path("plan.yaml"))
    """
    version: str = "1.0"
    operations: List[RefineOp] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        """Initialize metadata with timestamp if not provided."""
        if "created_at" not in self.metadata:
            self.metadata["created_at"] = datetime.now().isoformat()

    # -------------------------------------------------------------------------
    # Convenience methods for adding operations
    # -------------------------------------------------------------------------

    def add_override(
        self,
        cluster_id: str,
        new_label: str,
        reason: str = "",
        source: Literal["auto", "manual"] = "manual",
    ) -> "RefinePlan":
        """Add an override operation."""
        self.operations.append(OverrideOp(
            cluster_id=cluster_id,
            new_label=new_label,
            reason=reason,
            source=source,
        ))
        return self

    def add_subcluster(
        self,
        cluster_id: str,
        resolution: float = 0.4,
        n_pcs: int = 30,
        neighbors_k: int = 15,
        focus_markers: Optional[List[str]] = None,
        min_cells: int = 100,
        reason: str = "",
        source: Literal["auto", "manual"] = "manual",
    ) -> "RefinePlan":
        """Add a subclustering operation."""
        self.operations.append(SubclusterOp(
            cluster_id=cluster_id,
            resolution=resolution,
            n_pcs=n_pcs,
            neighbors_k=neighbors_k,
            focus_markers=focus_markers,
            min_cells=min_cells,
            reason=reason,
            source=source,
        ))
        return self

    def add_merge(
        self,
        source_clusters: List[str],
        target_label: str,
        reason: str = "",
        source: Literal["auto", "manual"] = "manual",
    ) -> "RefinePlan":
        """Add a merge operation."""
        self.operations.append(MergeOp(
            cluster_id=source_clusters,
            source_clusters=source_clusters,
            target_label=target_label,
            reason=reason,
            source=source,
        ))
        return self

    def add_relabel(
        self,
        cluster_id: str,
        new_label: str,
        old_label: str,
        confidence_score: float,
        reason: str = "",
        source: Literal["auto", "manual"] = "auto",
    ) -> "RefinePlan":
        """Add a relabel operation (instant, no clustering)."""
        self.operations.append(RelabelOp(
            cluster_id=cluster_id,
            new_label=new_label,
            old_label=old_label,
            confidence_score=confidence_score,
            reason=reason,
            source=source,
        ))
        return self

    def add_rescore(
        self,
        mode: Literal["smart", "full"] = "smart",
        recompute_de: bool = False,
        target_clusters: Optional[List[str]] = None,
        reason: str = "Rescore after modifications",
        source: Literal["auto", "manual"] = "auto",
    ) -> "RefinePlan":
        """Add a rescore operation."""
        self.operations.append(RescoreOp(
            cluster_id="all",
            mode=mode,
            recompute_de=recompute_de,
            target_clusters=target_clusters,
            reason=reason,
            source=source,
        ))
        return self

    # -------------------------------------------------------------------------
    # Query methods
    # -------------------------------------------------------------------------

    def get_operations_by_type(self, op_type: str) -> List[RefineOp]:
        """Get all operations of a specific type."""
        return [op for op in self.operations if op.op_type == op_type]

    def get_affected_clusters(self) -> List[str]:
        """Get list of all cluster IDs affected by operations."""
        clusters = set()
        for op in self.operations:
            if isinstance(op.cluster_id, list):
                clusters.update(op.cluster_id)
            elif op.cluster_id and op.cluster_id != "all":
                clusters.add(op.cluster_id)
        return sorted(clusters)

    def get_operations_for_cluster(self, cluster_id: str) -> List[RefineOp]:
        """Get all operations targeting a specific cluster."""
        result = []
        for op in self.operations:
            if isinstance(op.cluster_id, list):
                if cluster_id in op.cluster_id:
                    result.append(op)
            elif op.cluster_id == cluster_id:
                result.append(op)
        return result

    @property
    def is_empty(self) -> bool:
        """Check if plan has no operations."""
        return len(self.operations) == 0

    def summary(self) -> Dict[str, int]:
        """Get summary count of operations by type."""
        summary = {}
        for op in self.operations:
            summary[op.op_type] = summary.get(op.op_type, 0) + 1
        return summary

    # -------------------------------------------------------------------------
    # Serialization
    # -------------------------------------------------------------------------

    def to_dict(self) -> Dict[str, Any]:
        """Convert plan to dictionary for serialization."""
        return {
            "version": self.version,
            "metadata": self.metadata,
            "operations": [op.to_dict() for op in self.operations],
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "RefinePlan":
        """Create RefinePlan from dictionary."""
        operations = [RefineOp.from_dict(op) for op in data.get("operations", [])]
        return cls(
            version=data.get("version", "1.0"),
            operations=operations,
            metadata=data.get("metadata", {}),
        )

    def to_yaml(self, path: Path) -> None:
        """Export plan to YAML file."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as f:
            yaml.dump(self.to_dict(), f, default_flow_style=False, sort_keys=False)

    @classmethod
    def from_yaml(cls, path: Path) -> "RefinePlan":
        """Load RefinePlan from YAML file."""
        with open(path, "r") as f:
            data = yaml.safe_load(f)
        return cls.from_dict(data or {})

    def to_json(self, path: Path, indent: int = 2) -> None:
        """Export plan to JSON file."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=indent)

    @classmethod
    def from_json(cls, path: Path) -> "RefinePlan":
        """Load RefinePlan from JSON file."""
        with open(path, "r") as f:
            data = json.load(f)
        return cls.from_dict(data)

    def content_checksum(self) -> str:
        """Compute MD5 checksum of normalized plan content.

        This hashes a deterministic JSON representation of the plan object,
        NOT the source file bytes. Two semantically identical plans will
        produce the same checksum regardless of YAML formatting or whitespace.
        """
        # Sort operations for deterministic checksum
        plan_str = json.dumps(self.to_dict(), sort_keys=True)
        return hashlib.md5(plan_str.encode()).hexdigest()

    # -------------------------------------------------------------------------
    # Sorting and ordering
    # -------------------------------------------------------------------------

    def sort_operations(self) -> "RefinePlan":
        """Sort operations in deterministic execution order.

        Order: overrides → merges → subclusters → relabels → rescore
        Within each type, sort by priority (lower first), then cluster_id.
        """
        type_order = {
            "override": 0,
            "merge": 1,
            "subcluster": 2,
            "relabel": 3,
            "rescore": 4,
        }

        def sort_key(op: RefineOp):
            type_idx = type_order.get(op.op_type, 99)
            cluster_key = str(op.cluster_id) if not isinstance(op.cluster_id, list) else ",".join(sorted(op.cluster_id))
            return (type_idx, op.priority, cluster_key)

        self.operations.sort(key=sort_key)
        return self

    def __repr__(self) -> str:
        summary = self.summary()
        ops_str = ", ".join(f"{k}={v}" for k, v in summary.items())
        return f"RefinePlan(version={self.version!r}, operations=[{ops_str}])"
