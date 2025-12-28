"""Harmonization module for two-level cell type reporting.

This module provides functions to project cell type labels
to harmonized vocabularies for multi-omics alignment:

1. Fine vocabulary (10 classes + NA): Specific cell subtypes
2. Broad vocabulary (6 classes): Major lineage categories

Key design principle: Don't force cells into subtypes. Preserve uncertainty
explicitly via NA in fine vocab and honest broad mappings.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Set
import logging

import pandas as pd
import yaml

logger = logging.getLogger(__name__)


# Default vocabulary constants
# Note: These can be customized via HarmonizeConfig
FINE_VOCAB: Set[str] = {
    "Ciliated epithelium",
    "Glandular epithelium",
    "Lymphatic endothelium",
    "Blood endothelium",
    "Stromal fibroblast",
    "Smooth muscle",
    "Pericyte",
    "T/NK cell",
    "Macrophage",
}

BROAD_VOCAB: Set[str] = {
    "Epithelium",
    "Mesenchymal",
    "Smooth Muscle",
    "Endothelium",
    "Immune",
    "Other/Unassigned",
}


@dataclass
class LabelInfo:
    """Parsed label information.

    Attributes
    ----------
    raw_label : str
        Original label as-is
    base_label : str
        Label without (orphan) suffix
    is_orphan : bool
        True if label ends with " (orphan)"
    is_hybrid : bool
        True if label contains "~"
    hybrid_pair : str, optional
        Canonicalized pair for hybrids, e.g., "Endothelium|Immune Cells"
    root_label : str, optional
        Inferred root lineage
    annotation_level : str
        One of: leaf_subtype, intermediate, root_only, hybrid, unassigned
    """

    raw_label: str
    base_label: str
    is_orphan: bool
    is_hybrid: bool
    hybrid_pair: Optional[str]
    root_label: Optional[str]
    annotation_level: str


@dataclass
class HarmonizeConfig:
    """Configuration for label harmonization.

    Attributes
    ----------
    version : str
        Version identifier for mapping
    fine_mapping : Dict[str, str]
        Map from detailed label -> fine vocab label
    immune_to_na : List[str]
        Immune types that map to NA (strict policy)
    misc_to_na : List[str]
        Misc types that map to NA (no multi-omics equivalent)
    broad_from_fine : Dict[str, str]
        Map from fine vocab -> broad vocab
    broad_from_root : Dict[str, str]
        Map from root label -> broad vocab (for root-only/NA cases)
    root_types : Set[str]
        Set of root type names
    intermediate_types : Set[str]
        Set of intermediate type names (not leaf, not root)
    """

    version: str = "default_v1"
    fine_mapping: Dict[str, str] = field(default_factory=dict)
    immune_to_na: List[str] = field(default_factory=list)
    misc_to_na: List[str] = field(default_factory=list)
    broad_from_fine: Dict[str, str] = field(default_factory=dict)
    broad_from_root: Dict[str, str] = field(default_factory=dict)
    root_types: Set[str] = field(default_factory=lambda: {
        "Epithelium",
        "Endothelium",
        "Immune Cells",
        "Mesenchymal Cells",
        "Misc",
    })
    intermediate_types: Set[str] = field(default_factory=lambda: {
        "Myeloids",
        "Lymphoids",
        "Granulocytes",
        "Monocytes",
        "Glandular Epithelium",
        "T Cells",
    })

    @classmethod
    def from_yaml(cls, path: Path) -> "HarmonizeConfig":
        """Load configuration from YAML file.

        Parameters
        ----------
        path : Path
            Path to YAML config file

        Returns
        -------
        HarmonizeConfig
            Loaded configuration
        """
        path = Path(path)
        with open(path) as f:
            data = yaml.safe_load(f)

        return cls(
            version=data.get("version", "default_v1"),
            fine_mapping=data.get("fine_mapping", {}),
            immune_to_na=data.get("immune_to_na", []),
            misc_to_na=data.get("misc_to_na", []),
            broad_from_fine=data.get("broad_from_fine", {}),
            broad_from_root=data.get("broad_from_root", {}),
            root_types=set(data.get("root_types", [])),
            intermediate_types=set(data.get("intermediate_types", [])),
        )

    @classmethod
    def default(cls) -> "HarmonizeConfig":
        """Create default configuration with built-in mappings.

        Returns
        -------
        HarmonizeConfig
            Default configuration
        """
        return cls(
            version="default_v1",
            fine_mapping={
                # Epithelium subtypes
                "Ciliated Epithelium": "Ciliated epithelium",
                "Glandular Epithelium": "Glandular epithelium",
                "Peg Cells": "Glandular epithelium",
                "Secretory Epithelium": "Glandular epithelium",
                # Endothelium subtypes
                "Lymphatic Endothelium": "Lymphatic endothelium",
                "Vascular Endothelium": "Blood endothelium",
                # Mesenchymal subtypes
                "Fibroblasts": "Stromal fibroblast",
                # Mural/contractile cells
                "Smooth Muscle Cells": "Smooth muscle",
                "Pericytes": "Pericyte",
                # Immune subtypes
                "T Cells": "T/NK cell",
                "Activated T Cells": "T/NK cell",
                "Natural-Killer (NK) Cells": "T/NK cell",
                "NK Cells": "T/NK cell",
                "Macrophages": "Macrophage",
            },
            immune_to_na=[
                "B Cells",
                "Dendritic Cells",
                "Myeloids",
                "Monocytes",
                "Granulocytes",
                "Neutrophils",
                "Lymphoids",
            ],
            misc_to_na=[
                "Proliferating Cells",
                "Proliferating",
            ],
            broad_from_fine={
                "Ciliated epithelium": "Epithelium",
                "Glandular epithelium": "Epithelium",
                "Lymphatic endothelium": "Endothelium",
                "Blood endothelium": "Endothelium",
                "Stromal fibroblast": "Mesenchymal",
                "Smooth muscle": "Smooth Muscle",
                "Pericyte": "Smooth Muscle",
                "T/NK cell": "Immune",
                "Macrophage": "Immune",
            },
            broad_from_root={
                "Epithelium": "Epithelium",
                "Endothelium": "Endothelium",
                "Immune Cells": "Immune",
                "Mesenchymal Cells": "Mesenchymal",
                "Misc": "Other/Unassigned",
            },
        )


def parse_label(label: str, config: HarmonizeConfig) -> LabelInfo:
    """Parse a cell type label into structured information.

    Parameters
    ----------
    label : str
        The cell type label to parse
    config : HarmonizeConfig
        Configuration with root/intermediate type sets

    Returns
    -------
    LabelInfo
        Parsed label information
    """
    if not label or pd.isna(label):
        return LabelInfo(
            raw_label="",
            base_label="",
            is_orphan=False,
            is_hybrid=False,
            hybrid_pair=None,
            root_label=None,
            annotation_level="unassigned",
        )

    label = str(label).strip()
    raw_label = label

    # Check for Unassigned
    if label.lower() == "unassigned" or label.startswith("No root gate"):
        return LabelInfo(
            raw_label=raw_label,
            base_label="Unassigned",
            is_orphan=False,
            is_hybrid=False,
            hybrid_pair=None,
            root_label=None,
            annotation_level="unassigned",
        )

    # Check for orphan suffix
    is_orphan = label.endswith(" (orphan)")
    base_label = label[:-9].strip() if is_orphan else label

    # Check for hybrid (contains ~)
    is_hybrid = "~" in base_label
    hybrid_pair = None
    if is_hybrid:
        parts = sorted(base_label.split("~"))
        hybrid_pair = "|".join(parts)
        return LabelInfo(
            raw_label=raw_label,
            base_label=base_label,
            is_orphan=is_orphan,
            is_hybrid=True,
            hybrid_pair=hybrid_pair,
            root_label=None,  # Hybrids have ambiguous root
            annotation_level="hybrid",
        )

    # Determine annotation level and root label
    if base_label in config.root_types:
        return LabelInfo(
            raw_label=raw_label,
            base_label=base_label,
            is_orphan=is_orphan,
            is_hybrid=False,
            hybrid_pair=None,
            root_label=base_label,
            annotation_level="root_only",
        )

    if base_label in config.intermediate_types:
        # Infer root from intermediate type
        root_label = _infer_root_from_intermediate(base_label, config)
        return LabelInfo(
            raw_label=raw_label,
            base_label=base_label,
            is_orphan=is_orphan,
            is_hybrid=False,
            hybrid_pair=None,
            root_label=root_label,
            annotation_level="intermediate",
        )

    # Leaf subtype - infer root from fine mapping or default
    root_label = _infer_root_from_subtype(base_label, config)
    return LabelInfo(
        raw_label=raw_label,
        base_label=base_label,
        is_orphan=is_orphan,
        is_hybrid=False,
        hybrid_pair=None,
        root_label=root_label,
        annotation_level="leaf_subtype",
    )


def _infer_root_from_intermediate(label: str, config: HarmonizeConfig) -> Optional[str]:
    """Infer root lineage from intermediate type.

    Parameters
    ----------
    label : str
        Intermediate type label
    config : HarmonizeConfig
        Configuration

    Returns
    -------
    str, optional
        Inferred root label
    """
    # Known intermediate -> root mappings
    intermediate_to_root = {
        "Myeloids": "Immune Cells",
        "Lymphoids": "Immune Cells",
        "Granulocytes": "Immune Cells",
        "Monocytes": "Immune Cells",
        "Glandular Epithelium": "Epithelium",
        "T Cells": "Immune Cells",
    }
    return intermediate_to_root.get(label)


def _infer_root_from_subtype(label: str, config: HarmonizeConfig) -> Optional[str]:
    """Infer root lineage from subtype label.

    Parameters
    ----------
    label : str
        Subtype label
    config : HarmonizeConfig
        Configuration

    Returns
    -------
    str, optional
        Inferred root label
    """
    # Check if label maps to a fine vocab term, then infer root from broad mapping
    fine = config.fine_mapping.get(label)
    if fine:
        broad = config.broad_from_fine.get(fine)
        if broad:
            # Reverse lookup: broad -> root
            for root, b in config.broad_from_root.items():
                if b == broad:
                    return root

    # Check immune_to_na list
    if label in config.immune_to_na:
        return "Immune Cells"

    # Fallback: try pattern matching on label name
    label_lower = label.lower()
    if "epithelium" in label_lower or "ciliated" in label_lower:
        return "Epithelium"
    if "endothelium" in label_lower or "lymphatic" in label_lower or "vascular" in label_lower:
        return "Endothelium"
    if any(kw in label_lower for kw in ["fibroblast", "muscle", "pericyte", "stromal"]):
        return "Mesenchymal Cells"
    if any(kw in label_lower for kw in ["t cell", "b cell", "macrophage", "dendritic", "nk cell", "immune"]):
        return "Immune Cells"

    return None


def map_to_fine(
    base_label: str,
    annotation_level: str,
    config: HarmonizeConfig,
) -> Optional[str]:
    """Map a label to the fine vocabulary.

    Parameters
    ----------
    base_label : str
        Base label (without orphan suffix)
    annotation_level : str
        The annotation level (leaf_subtype, intermediate, root_only, hybrid, unassigned)
    config : HarmonizeConfig
        Configuration with mapping tables

    Returns
    -------
    str, optional
        Fine vocab label or None (for NA)
    """
    # Root-only, hybrid, and unassigned always map to NA
    if annotation_level in ("root_only", "hybrid", "unassigned"):
        return None

    # Check misc_to_na list first (no multi-omics equivalent)
    if base_label in config.misc_to_na:
        return None

    # Check direct mapping
    if base_label in config.fine_mapping:
        fine = config.fine_mapping[base_label]
        if fine in FINE_VOCAB:
            return fine
        logger.warning(f"Fine mapping '{fine}' for '{base_label}' not in FINE_VOCAB")
        return None

    # Check immune_to_na list (strict policy)
    if base_label in config.immune_to_na:
        return None

    # Intermediate types that aren't in fine_mapping map to NA
    if annotation_level == "intermediate":
        return None

    # Unknown subtype - log warning and return NA
    logger.debug(f"No fine mapping for '{base_label}', returning NA")
    return None


def map_to_broad(
    fine: Optional[str],
    root_label: Optional[str],
    annotation_level: str,
    config: HarmonizeConfig,
) -> str:
    """Map to the broad vocabulary.

    Parameters
    ----------
    fine : str, optional
        Fine vocab label (or None for NA)
    root_label : str, optional
        Root lineage label
    annotation_level : str
        The annotation level
    config : HarmonizeConfig
        Configuration with mapping tables

    Returns
    -------
    str
        Broad vocab label (always returns a valid value)
    """
    # If we have a fine label, use fine -> broad mapping
    if fine is not None:
        broad = config.broad_from_fine.get(fine)
        if broad and broad in BROAD_VOCAB:
            return broad
        logger.warning(f"Fine label '{fine}' has no broad mapping, falling back to root")

    # Hybrid always maps to Other/Unassigned
    if annotation_level == "hybrid":
        return "Other/Unassigned"

    # Unassigned always maps to Other/Unassigned
    if annotation_level == "unassigned":
        return "Other/Unassigned"

    # Use root -> broad mapping if root is known
    if root_label is not None:
        broad = config.broad_from_root.get(root_label)
        if broad and broad in BROAD_VOCAB:
            return broad
        logger.warning(f"Root label '{root_label}' has no broad mapping")

    # Fallback
    return "Other/Unassigned"


def harmonize_single_label(
    label: str,
    config: HarmonizeConfig,
    external_root_label: Optional[str] = None,
) -> Dict[str, Any]:
    """Harmonize a single label to all output fields.

    Parameters
    ----------
    label : str
        The cell type label
    config : HarmonizeConfig
        Configuration
    external_root_label : str, optional
        Root label from external source (e.g., cluster_annotations)
        If provided, overrides inferred root_label

    Returns
    -------
    Dict[str, Any]
        Dictionary with all harmonized fields
    """
    # Parse the label
    info = parse_label(label, config)

    # Use external root_label if provided and not hybrid
    root_label = info.root_label
    if external_root_label and not info.is_hybrid:
        root_label = external_root_label

    # Map to fine and broad
    fine = map_to_fine(info.base_label, info.annotation_level, config)
    broad = map_to_broad(fine, root_label, info.annotation_level, config)

    # Build mapping notes
    mapping_notes = None
    if fine is None:
        if info.annotation_level == "hybrid":
            mapping_notes = "hybrid_ambiguous"
        elif info.annotation_level == "root_only":
            mapping_notes = "root_only_no_subtype"
        elif info.annotation_level == "unassigned":
            mapping_notes = "unassigned"
        elif info.base_label in config.misc_to_na:
            mapping_notes = "misc_no_equivalent"
        elif info.base_label in config.immune_to_na:
            mapping_notes = "immune_strict_na"
        elif info.annotation_level == "intermediate":
            mapping_notes = f"intermediate_{info.base_label.lower().replace(' ', '_')}"
        else:
            mapping_notes = "unmapped_subtype"

    return {
        "cell_type_fine": fine,
        "cell_type_broad": broad,
        "root_label": root_label,
        "annotation_level": info.annotation_level,
        "is_orphan": info.is_orphan,
        "hybrid_pair": info.hybrid_pair,
        "mapping_version": config.version,
        "mapping_notes": mapping_notes,
    }


def harmonize_adata(
    adata: "sc.AnnData",
    diagnostic_report: pd.DataFrame,
    config: HarmonizeConfig,
    label_col: str = "cell_type_final",
    cluster_col: str = "cluster_lvl1",
) -> None:
    """Add harmonized columns to AnnData.obs in-place.

    Parameters
    ----------
    adata : sc.AnnData
        AnnData with cell type column
    diagnostic_report : pd.DataFrame
        Diagnostic report with cluster_id and root_label columns
    config : HarmonizeConfig
        Harmonization configuration
    label_col : str
        Column containing cell type labels
    cluster_col : str
        Column containing cluster IDs
    """
    logger.info(f"Harmonizing labels from '{label_col}' column...")

    # Build cluster -> root_label lookup from diagnostic_report
    cluster_to_root: Dict[str, str] = {}
    if "root_label" in diagnostic_report.columns:
        for _, row in diagnostic_report.iterrows():
            cid = str(row.get("cluster_id", ""))
            root = row.get("root_label")
            if cid and root and not pd.isna(root):
                cluster_to_root[cid] = str(root)
        logger.info(f"Loaded root_label for {len(cluster_to_root)} clusters from diagnostic_report")

    # Get cluster IDs if available
    cluster_ids = None
    if cluster_col in adata.obs.columns:
        cluster_ids = adata.obs[cluster_col].astype(str)

    # Initialize output columns
    n_cells = len(adata)
    results = {
        "cell_type_fine": [None] * n_cells,
        "cell_type_broad": [""] * n_cells,
        "root_label": [None] * n_cells,
        "annotation_level": [""] * n_cells,
        "is_orphan": [False] * n_cells,
        "hybrid_pair": [None] * n_cells,
        "mapping_version": [config.version] * n_cells,
        "mapping_notes": [None] * n_cells,
    }

    # Harmonize each cell
    labels = adata.obs[label_col].values
    for i, label in enumerate(labels):
        # Get external root_label if available
        ext_root = None
        if cluster_ids is not None:
            cid = cluster_ids.iloc[i]
            ext_root = cluster_to_root.get(cid)

        # Harmonize
        harmonized = harmonize_single_label(str(label) if not pd.isna(label) else "", config, ext_root)

        # Store results
        for key, value in harmonized.items():
            results[key][i] = value

    # Add columns to adata.obs
    adata.obs["cell_type_fine"] = pd.Categorical(results["cell_type_fine"])
    adata.obs["cell_type_broad"] = pd.Categorical(results["cell_type_broad"])
    adata.obs["root_label"] = results["root_label"]
    adata.obs["annotation_level"] = pd.Categorical(results["annotation_level"])
    adata.obs["is_orphan"] = results["is_orphan"]
    adata.obs["hybrid_pair"] = results["hybrid_pair"]
    adata.obs["mapping_version"] = config.version
    adata.obs["mapping_notes"] = results["mapping_notes"]

    # Log summary
    fine_counts = adata.obs["cell_type_fine"].value_counts(dropna=False)
    broad_counts = adata.obs["cell_type_broad"].value_counts()
    n_orphan = sum(results["is_orphan"])
    n_hybrid = sum(1 for h in results["hybrid_pair"] if h is not None)

    logger.info("Harmonization complete:")
    logger.info(f"  - {len(fine_counts) - 1} fine classes + NA")
    logger.info(f"  - {len(broad_counts)} broad classes")
    logger.info(f"  - {n_orphan:,} orphan cells")
    logger.info(f"  - {n_hybrid:,} hybrid cells")


# Import scanpy optionally (for type hints)
try:
    import scanpy as sc
except ImportError:
    sc = None
