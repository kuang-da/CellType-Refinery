"""
Provenance tracking for cell-type annotation refinement.

This module provides:
- RefinementProvenance: Dataclass for tracking refinement metadata
- Functions to store and retrieve provenance from AnnData objects
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field, asdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, TYPE_CHECKING

try:
    import scanpy as sc
except ImportError:
    sc = None

if TYPE_CHECKING:
    from .plan import RefinePlan
    from .engine import EngineResult


@dataclass
class RefinementProvenance:
    """Unified provenance for refinement operations.

    This dataclass tracks all metadata needed for reproducibility:
    - Source stage and file information
    - Plan details and content checksum
    - Execution results
    - Timestamps and configuration
    """
    # Source information
    parent_stage: str  # "clustering", "refinement", or custom stage name
    parent_file: str
    iteration: int

    # Plan information
    plan_file: Optional[str] = None
    plan_content_checksum: str = ""
    plan_summary: Dict[str, int] = field(default_factory=dict)

    # Execution results
    label_key_out: str = "cell_type_curated"
    operations_summary: Dict[str, int] = field(default_factory=dict)
    clusters_modified: List[str] = field(default_factory=list)
    subclusters_created: List[str] = field(default_factory=list)
    n_cells_modified: int = 0

    # Configuration
    auto_policy_config: Optional[Dict[str, Any]] = None
    manual_config_path: Optional[str] = None

    # Metadata
    timestamp: str = ""
    gpu_used: bool = False
    execution_time_seconds: float = 0.0

    def __post_init__(self):
        if not self.timestamp:
            self.timestamp = datetime.now().isoformat()

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for storage in adata.uns."""
        return asdict(self)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "RefinementProvenance":
        """Create from dictionary."""
        return cls(**data)

    def to_uns_key(self) -> str:
        """Generate the adata.uns key for this provenance."""
        return f"refinement_iteration_{self.iteration}"


def create_provenance(
    adata: "sc.AnnData",
    plan: "RefinePlan",
    result: "EngineResult",
    label_key_out: str = "cell_type_curated",
    auto_policy_config: Optional[Dict] = None,
    manual_config_path: Optional[Path] = None,
    gpu_used: bool = False,
    input_path: Optional[str] = None,
) -> RefinementProvenance:
    """Create provenance from execution context.

    Parameters
    ----------
    adata : sc.AnnData
        AnnData object with source provenance
    plan : RefinePlan
        Executed plan
    result : EngineResult
        Execution result
    label_key_out : str
        Output label column name
    auto_policy_config : Dict, optional
        AutoPolicy configuration if used
    manual_config_path : Path, optional
        Path to manual config YAML if used
    gpu_used : bool
        Whether GPU was used
    input_path : str, optional
        Explicit input file path (preferred over adata.uns lookup)

    Returns
    -------
    RefinementProvenance
        Provenance object
    """
    # Detect parent stage
    parent_stage, iteration = detect_parent_stage(adata)

    # Get parent file - prefer explicit input_path, fall back to adata.uns lookup
    parent_file = input_path or ""
    if not parent_file:
        # Try various provenance keys
        for key in ["clustering_provenance", "stage_h_provenance", "refinement_provenance"]:
            if key in adata.uns:
                parent_file = adata.uns[key].get("output_path", "")
                if parent_file:
                    break

    return RefinementProvenance(
        parent_stage=parent_stage,
        parent_file=parent_file,
        iteration=iteration + 1,  # This is the new iteration
        plan_content_checksum=plan.content_checksum(),
        plan_summary=plan.summary(),
        label_key_out=label_key_out,
        operations_summary=result.operations_executed,
        clusters_modified=plan.get_affected_clusters(),
        subclusters_created=result.subclusters_created,
        n_cells_modified=result.n_cells_modified,
        auto_policy_config=auto_policy_config,
        manual_config_path=str(manual_config_path) if manual_config_path else None,
        gpu_used=gpu_used,
        execution_time_seconds=result.execution_time_seconds,
    )


def detect_parent_stage(adata: "sc.AnnData") -> tuple[str, int]:
    """Detect the parent stage and current iteration from provenance.

    Parameters
    ----------
    adata : sc.AnnData
        AnnData object with provenance in .uns

    Returns
    -------
    tuple[str, int]
        (parent_stage, current_iteration)
        parent_stage is "clustering", "refinement", or "unknown"
        current_iteration is 0 for clustering, or the iteration number
    """
    # Check for refinement iterations (generic key)
    refine_keys = sorted([k for k in adata.uns.keys() if k.startswith("refinement_iteration_")])
    if refine_keys:
        # Extract highest iteration number
        n = max(int(k.split("_")[-1]) for k in refine_keys)
        return ("refinement", n)

    # Check for Stage J iterations (backward compat)
    j_keys = sorted([k for k in adata.uns.keys() if k.startswith("stage_j_iteration_")])
    if j_keys:
        n = len(j_keys)
        return ("refinement", n)

    # Check for Stage I provenance (backward compat)
    if "stage_i_provenance" in adata.uns:
        iteration = adata.uns["stage_i_provenance"].get("iteration", 1)
        return ("refinement", iteration)

    # Check for Stage refine iterations (backward compat)
    stage_refine_keys = sorted([k for k in adata.uns.keys() if k.startswith("stage_refine_iteration_")])
    if stage_refine_keys:
        n = max(int(k.split("_")[-1]) for k in stage_refine_keys)
        return ("refinement", n)

    # Check for clustering provenance (generic key)
    if "clustering_provenance" in adata.uns:
        return ("clustering", 0)

    # Check for Stage H provenance (backward compat)
    if "stage_h_provenance" in adata.uns:
        return ("clustering", 0)

    # No provenance found
    return ("unknown", 0)


def store_provenance(adata: "sc.AnnData", provenance: RefinementProvenance) -> None:
    """Store provenance in AnnData.uns.

    Parameters
    ----------
    adata : sc.AnnData
        AnnData object to modify
    provenance : RefinementProvenance
        Provenance to store
    """
    key = provenance.to_uns_key()
    adata.uns[key] = provenance.to_dict()


def get_provenance(adata: "sc.AnnData", iteration: Optional[int] = None) -> Optional[RefinementProvenance]:
    """Retrieve provenance from AnnData.uns.

    Parameters
    ----------
    adata : sc.AnnData
        AnnData object
    iteration : int, optional
        Specific iteration to retrieve. If None, returns latest.

    Returns
    -------
    RefinementProvenance or None
        Provenance object if found
    """
    if iteration is not None:
        key = f"refinement_iteration_{iteration}"
        if key in adata.uns:
            return RefinementProvenance.from_dict(adata.uns[key])
        return None

    # Find latest iteration
    refine_keys = sorted([k for k in adata.uns.keys() if k.startswith("refinement_iteration_")])
    if not refine_keys:
        return None

    latest_key = refine_keys[-1]
    return RefinementProvenance.from_dict(adata.uns[latest_key])


def get_provenance_chain(adata: "sc.AnnData") -> List[Dict[str, Any]]:
    """Get full provenance chain from AnnData.

    Returns list of all provenance entries in chronological order.
    """
    chain = []

    # Check for clustering provenance (generic)
    if "clustering_provenance" in adata.uns:
        chain.append({"stage": "clustering", "data": adata.uns["clustering_provenance"]})
    # Backward compat: Stage H
    elif "stage_h_provenance" in adata.uns:
        chain.append({"stage": "clustering", "data": adata.uns["stage_h_provenance"]})

    # Refinement iterations (generic)
    for key in sorted(adata.uns.keys()):
        if key.startswith("refinement_iteration_"):
            chain.append({"stage": key, "data": adata.uns[key]})

    # Backward compat: Stage refine iterations
    for key in sorted(adata.uns.keys()):
        if key.startswith("stage_refine_iteration_"):
            chain.append({"stage": key, "data": adata.uns[key]})

    # Backward compat: Stage J iterations
    for key in sorted(adata.uns.keys()):
        if key.startswith("stage_j_iteration_"):
            chain.append({"stage": key, "data": adata.uns[key]})

    return chain
