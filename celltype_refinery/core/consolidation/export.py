"""Export functions for consolidation outputs.

This module provides functions to export:
- Consolidation summary (final label distribution)
- Confidence breakdown (cells by confidence band)
- Override log (audit trail of manual overrides)
- Orphan report (detected and rescued orphans)
- Provenance JSON (full audit trail)
"""

from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional
import json
import logging

import pandas as pd

try:
    import scanpy as sc
except ImportError:
    sc = None

from .config import ConsolidationConfig
from .engine import ConsolidationResult
from .orphan_detection import OrphanCandidate
from .iel_rescue import IELCandidate, summarize_iel_candidates, summarize_iel_by_type
from .rules import classify_label, simplify_label

logger = logging.getLogger(__name__)


def export_consolidation_summary(
    adata: "sc.AnnData",
    result: ConsolidationResult,
    output_path: Path,
    output_col: str = "cell_type_final",
) -> pd.DataFrame:
    """Export final label distribution summary with harmonized columns.

    Parameters
    ----------
    adata : sc.AnnData
        AnnData with final labels and harmonized columns
    result : ConsolidationResult
        Consolidation result
    output_path : Path
        Output path for CSV
    output_col : str
        Column containing final labels

    Returns
    -------
    pd.DataFrame
        Summary DataFrame with harmonized mappings
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Build mapping: simplified_label -> list of detailed labels that merged into it
    detailed_to_simplified = {}
    if result.label_distribution_detailed:
        for detailed_label in result.label_distribution_detailed.keys():
            simplified = simplify_label(detailed_label)
            if simplified not in detailed_to_simplified:
                detailed_to_simplified[simplified] = []
            if detailed_label != simplified:
                detailed_to_simplified[simplified].append(detailed_label)

    # Build lookup for harmonized values per label
    # All cells with the same label have the same harmonized values
    harmonized_cols = [
        "cell_type_fine",
        "cell_type_broad",
        "root_label",
        "annotation_level",
        "is_orphan",
        "hybrid_pair",
        "mapping_version",
        "mapping_notes",
    ]
    harmonized_lookup = {}
    for label in result.label_distribution.keys():
        mask = adata.obs[output_col] == label
        if mask.any():
            first_idx = mask.idxmax()
            harmonized_lookup[label] = {
                col: adata.obs.loc[first_idx, col] if col in adata.obs.columns else None
                for col in harmonized_cols
            }

    # Build summary rows
    rows = []
    for label, count in sorted(result.label_distribution.items(), key=lambda x: -x[1]):
        category = classify_label(label)
        proportion = count / result.n_cells_total if result.n_cells_total > 0 else 0

        row = {
            "cell_type_final": label,
            "n_cells": count,
            "proportion": round(proportion, 4),
            "percentage": round(proportion * 100, 2),
            "category": category.value,
        }

        # Add detailed variants that merged into this simplified label
        merged_variants = detailed_to_simplified.get(label, [])
        row["merged_from_detailed"] = "; ".join(sorted(merged_variants)) if merged_variants else ""

        # Add harmonized columns
        if label in harmonized_lookup:
            h = harmonized_lookup[label]
            fine = h.get("cell_type_fine")
            row["cell_type_fine"] = "NA" if pd.isna(fine) else str(fine)
            row["cell_type_broad"] = str(h.get("cell_type_broad", ""))
            row["root_label"] = str(h.get("root_label", "")) if h.get("root_label") else ""
            row["annotation_level"] = str(h.get("annotation_level", ""))
            row["is_orphan"] = bool(h.get("is_orphan", False))
            row["hybrid_pair"] = str(h.get("hybrid_pair", "")) if h.get("hybrid_pair") else ""
            row["mapping_version"] = str(h.get("mapping_version", ""))
            row["mapping_notes"] = str(h.get("mapping_notes", "")) if h.get("mapping_notes") else ""

        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(output_path, index=False)
    logger.info(f"Consolidation summary saved to {output_path}")

    return df


def export_consolidation_summary_detailed(
    result: ConsolidationResult,
    output_path: Path,
) -> pd.DataFrame:
    """Export detailed label distribution (with orphan suffix and unsorted hybrids).

    Parameters
    ----------
    result : ConsolidationResult
        Consolidation result
    output_path : Path
        Output path for CSV

    Returns
    -------
    pd.DataFrame
        Summary DataFrame with detailed labels
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if not result.label_distribution_detailed:
        logger.warning("No detailed label distribution available")
        return pd.DataFrame()

    rows = []
    for label, count in sorted(result.label_distribution_detailed.items(), key=lambda x: -x[1]):
        category = classify_label(label)
        proportion = count / result.n_cells_total if result.n_cells_total > 0 else 0
        simplified = simplify_label(label)

        row = {
            "cell_type_final_detailed": label,
            "cell_type_final": simplified,
            "n_cells": count,
            "proportion": round(proportion, 4),
            "percentage": round(proportion * 100, 2),
            "category": category.value,
            "is_orphan_label": label.endswith(" (orphan)"),
            "was_merged": label != simplified,
        }
        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(output_path, index=False)
    logger.info(f"Detailed consolidation summary saved to {output_path}")

    return df


def export_confidence_breakdown(
    result: ConsolidationResult,
    output_path: Path,
) -> pd.DataFrame:
    """Export confidence band breakdown.

    Parameters
    ----------
    result : ConsolidationResult
        Consolidation result
    output_path : Path
        Output path for CSV

    Returns
    -------
    pd.DataFrame
        Breakdown DataFrame
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    total = sum(result.confidence_breakdown.values()) if result.confidence_breakdown else 0

    for band in ["high", "medium", "low", "very_low", "unknown"]:
        count = result.confidence_breakdown.get(band, 0)
        proportion = count / total if total > 0 else 0
        rows.append({
            "confidence_band": band,
            "n_cells": count,
            "proportion": round(proportion, 4),
            "percentage": round(proportion * 100, 2),
        })

    df = pd.DataFrame(rows)
    df.to_csv(output_path, index=False)
    logger.info(f"Confidence breakdown saved to {output_path}")

    return df


def export_override_log(
    config: ConsolidationConfig,
    result: ConsolidationResult,
    output_path: Path,
) -> pd.DataFrame:
    """Export log of all manual overrides configured.

    Parameters
    ----------
    config : ConsolidationConfig
        Configuration with overrides
    result : ConsolidationResult
        Consolidation result
    output_path : Path
        Output path for CSV

    Returns
    -------
    pd.DataFrame
        Override log DataFrame
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    for override in config.overrides:
        rows.append({
            "cluster_id": override.cluster_id,
            "final_label": override.final_label,
            "reason": override.reason,
            "type": "manual_override",
        })

    for rule in config.relabel_rules:
        rows.append({
            "cluster_id": f"*{rule.from_label}*",
            "final_label": rule.to_label,
            "reason": rule.reason,
            "type": "relabel_rule",
        })

    df = pd.DataFrame(rows)
    df.to_csv(output_path, index=False)
    logger.info(f"Override log saved to {output_path}")

    return df


def export_orphan_report(
    result: ConsolidationResult,
    output_path: Path,
) -> pd.DataFrame:
    """Export orphan detection and rescue report.

    Parameters
    ----------
    result : ConsolidationResult
        Consolidation result with orphan candidates
    output_path : Path
        Output path for CSV

    Returns
    -------
    pd.DataFrame
        Orphan report DataFrame
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if not result.orphan_candidates:
        # Create empty DataFrame with expected columns
        df = pd.DataFrame(columns=[
            "cluster_id",
            "n_cells",
            "subtype",
            "subtype_score",
            "parent_root",
            "root_score",
            "plausibility",
            "action",
            "final_label",
            "flag",
        ])
        df.to_csv(output_path, index=False)
        logger.info(f"Empty orphan report saved to {output_path}")
        return df

    rows = []
    for c in result.orphan_candidates:
        rows.append({
            "cluster_id": c.cluster_id,
            "n_cells": c.n_cells,
            "subtype": c.subtype,
            "subtype_score": round(c.subtype_score, 3),
            "parent_root": c.parent_root,
            "root_score": round(c.root_score, 3),
            "plausibility": c.plausibility.value,
            "action": c.action.value,
            "final_label": c.final_label,
            "flag": c.flag,
        })

    df = pd.DataFrame(rows)
    df = df.sort_values("n_cells", ascending=False)
    df.to_csv(output_path, index=False)
    logger.info(f"Orphan report saved to {output_path} ({len(df)} candidates)")

    return df


def export_orphan_summary(
    result: ConsolidationResult,
    output_path: Path,
) -> pd.DataFrame:
    """Export orphan summary by subtype.

    Parameters
    ----------
    result : ConsolidationResult
        Consolidation result with orphan candidates
    output_path : Path
        Output path for CSV

    Returns
    -------
    pd.DataFrame
        Summary DataFrame
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if not result.orphan_candidates:
        df = pd.DataFrame(columns=[
            "subtype",
            "n_clusters",
            "n_cells",
            "mean_subtype_score",
            "mean_root_score",
            "plausibility",
            "action",
        ])
        df.to_csv(output_path, index=False)
        return df

    # Group by subtype
    from collections import defaultdict
    subtype_data: Dict[str, Dict[str, Any]] = defaultdict(
        lambda: {
            "n_clusters": 0,
            "n_cells": 0,
            "subtype_scores": [],
            "root_scores": [],
            "plausibility": None,
            "action": None,
        }
    )

    for c in result.orphan_candidates:
        d = subtype_data[c.subtype]
        d["n_clusters"] += 1
        d["n_cells"] += c.n_cells
        d["subtype_scores"].append(c.subtype_score)
        d["root_scores"].append(c.root_score)
        d["plausibility"] = c.plausibility.value
        d["action"] = c.action.value

    rows = []
    for subtype, d in subtype_data.items():
        rows.append({
            "subtype": subtype,
            "n_clusters": d["n_clusters"],
            "n_cells": d["n_cells"],
            "mean_subtype_score": round(sum(d["subtype_scores"]) / len(d["subtype_scores"]), 3),
            "mean_root_score": round(sum(d["root_scores"]) / len(d["root_scores"]), 3),
            "plausibility": d["plausibility"],
            "action": d["action"],
        })

    df = pd.DataFrame(rows)
    df = df.sort_values("n_cells", ascending=False)
    df.to_csv(output_path, index=False)
    logger.info(f"Orphan summary saved to {output_path}")

    return df


def export_iel_report(
    result: ConsolidationResult,
    output_path: Path,
) -> pd.DataFrame:
    """Export IEL (intraepithelial immune) detection and rescue report.

    Parameters
    ----------
    result : ConsolidationResult
        Consolidation result with IEL candidates
    output_path : Path
        Output path for CSV

    Returns
    -------
    pd.DataFrame
        IEL report DataFrame
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if not result.iel_candidates:
        # Create empty DataFrame with expected columns
        df = pd.DataFrame(columns=[
            "cluster_id",
            "n_cells",
            "cd45_pos_frac",
            "veto_marker",
            "veto_marker_pos_frac",
            "lymphoid_score",
            "myeloid_score",
            "iel_type",
            "final_label",
            "reason",
        ])
        df.to_csv(output_path, index=False)
        logger.info(f"Empty IEL report saved to {output_path}")
        return df

    df = summarize_iel_candidates(result.iel_candidates)
    df.to_csv(output_path, index=False)
    logger.info(f"IEL report saved to {output_path} ({len(df)} candidates)")

    return df


def export_iel_summary(
    result: ConsolidationResult,
    output_path: Path,
) -> pd.DataFrame:
    """Export IEL summary by type.

    Parameters
    ----------
    result : ConsolidationResult
        Consolidation result with IEL candidates
    output_path : Path
        Output path for CSV

    Returns
    -------
    pd.DataFrame
        Summary DataFrame by IEL type
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if not result.iel_candidates:
        df = pd.DataFrame(columns=[
            "iel_type",
            "n_clusters",
            "n_cells",
            "mean_cd45",
            "mean_lymphoid",
            "mean_myeloid",
        ])
        df.to_csv(output_path, index=False)
        return df

    df = summarize_iel_by_type(result.iel_candidates)
    df.to_csv(output_path, index=False)
    logger.info(f"IEL summary saved to {output_path}")

    return df


def export_mapping_table(
    mapping_df: pd.DataFrame,
    output_path: Path,
) -> None:
    """Export cluster to final label mapping table.

    Parameters
    ----------
    mapping_df : pd.DataFrame
        Mapping table from engine.get_mapping_table()
    output_path : Path
        Output path for CSV
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    mapping_df.to_csv(output_path, index=False)
    logger.info(f"Mapping table saved to {output_path}")


def export_provenance_json(
    result: ConsolidationResult,
    config: Optional[ConsolidationConfig],
    input_path: Path,
    output_path: Path,
) -> Dict[str, Any]:
    """Export provenance information as JSON for audit trail.

    Parameters
    ----------
    result : ConsolidationResult
        Consolidation result
    config : ConsolidationConfig, optional
        Configuration used
    input_path : Path
        Path to input AnnData
    output_path : Path
        Output path for JSON

    Returns
    -------
    Dict[str, Any]
        Provenance dictionary
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    provenance = {
        "module": "celltype_refinery.core.consolidation",
        "version": "1.0.0",
        "timestamp": datetime.now().isoformat(),
        "input": {
            "adata_path": str(input_path),
        },
        "config": {
            "version": config.version if config else None,
            "description": config.description if config else None,
            "n_overrides": len(config.overrides) if config else 0,
            "n_relabel_rules": len(config.relabel_rules) if config else 0,
            "orphan_rescue_enabled": config.orphan_rescue.enabled if config else False,
        },
        "result": {
            "success": result.success,
            "n_cells_total": result.n_cells_total,
            "n_clusters_total": result.n_clusters_total,
            "n_overrides_applied": result.n_overrides_applied,
            "n_relabels_applied": result.n_relabels_applied,
            "n_orphans_rescued": result.n_orphans_rescued,
            "n_orphans_flagged": result.n_orphans_flagged,
            "n_iel_rescued": result.n_iel_rescued,
            "execution_time_seconds": round(result.execution_time_seconds, 2),
        },
        "label_distribution": result.label_distribution,
        "category_breakdown": result.category_breakdown,
        "errors": result.errors,
        "warnings": result.warnings,
    }

    with open(output_path, "w") as f:
        json.dump(provenance, f, indent=2)

    logger.info(f"Provenance saved to {output_path}")

    return provenance


def export_harmonized_summary(
    adata: "sc.AnnData",
    output_path: Path,
) -> pd.DataFrame:
    """Export harmonized label distribution summary.

    Parameters
    ----------
    adata : sc.AnnData
        AnnData with harmonized columns
    output_path : Path
        Output path for CSV

    Returns
    -------
    pd.DataFrame
        Summary DataFrame with fine and broad distributions
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Check if harmonized columns exist
    if "cell_type_fine" not in adata.obs.columns:
        logger.warning("Harmonized columns not found in adata.obs, skipping export")
        return pd.DataFrame()

    n_cells = len(adata)

    # Fine distribution
    fine_counts = adata.obs["cell_type_fine"].value_counts(dropna=False)
    fine_rows = []
    for label, count in fine_counts.items():
        label_str = "NA" if pd.isna(label) else str(label)
        proportion = count / n_cells if n_cells > 0 else 0
        fine_rows.append({
            "level": "fine",
            "label": label_str,
            "n_cells": count,
            "proportion": round(proportion, 4),
            "percentage": round(proportion * 100, 2),
        })

    # Broad distribution
    broad_counts = adata.obs["cell_type_broad"].value_counts(dropna=False)
    broad_rows = []
    for label, count in broad_counts.items():
        label_str = str(label) if not pd.isna(label) else "Unknown"
        proportion = count / n_cells if n_cells > 0 else 0
        broad_rows.append({
            "level": "broad",
            "label": label_str,
            "n_cells": count,
            "proportion": round(proportion, 4),
            "percentage": round(proportion * 100, 2),
        })

    # Combine and sort
    df = pd.DataFrame(fine_rows + broad_rows)
    df = df.sort_values(["level", "n_cells"], ascending=[True, False])
    df.to_csv(output_path, index=False)

    logger.info(f"Harmonized summary saved to {output_path}")

    return df


def export_all(
    adata: "sc.AnnData",
    result: ConsolidationResult,
    config: Optional[ConsolidationConfig],
    mapping_df: pd.DataFrame,
    input_path: Path,
    output_dir: Path,
    output_col: str = "cell_type_final",
) -> Dict[str, Path]:
    """Export all consolidation outputs.

    Parameters
    ----------
    adata : sc.AnnData
        AnnData with final labels
    result : ConsolidationResult
        Consolidation result
    config : ConsolidationConfig, optional
        Configuration used
    mapping_df : pd.DataFrame
        Mapping table
    input_path : Path
        Path to input AnnData
    output_dir : Path
        Output directory
    output_col : str
        Column containing final labels

    Returns
    -------
    Dict[str, Path]
        Dictionary of output names to paths
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    outputs: Dict[str, Path] = {}

    # Main outputs
    outputs["summary"] = output_dir / "consolidation_summary.csv"
    export_consolidation_summary(adata, result, outputs["summary"], output_col)

    # Detailed summary (with orphan suffix and unsorted hybrids)
    if result.label_distribution_detailed:
        outputs["summary_detailed"] = output_dir / "consolidation_summary_detailed.csv"
        export_consolidation_summary_detailed(result, outputs["summary_detailed"])

    outputs["mapping"] = output_dir / "cluster_label_mapping.csv"
    export_mapping_table(mapping_df, outputs["mapping"])

    outputs["provenance"] = output_dir / "provenance.json"
    export_provenance_json(result, config, input_path, outputs["provenance"])

    # Confidence breakdown (if available)
    if result.confidence_breakdown:
        outputs["confidence"] = output_dir / "confidence_breakdown.csv"
        export_confidence_breakdown(result, outputs["confidence"])

    # Override log (if config has overrides)
    if config and (config.overrides or config.relabel_rules):
        outputs["overrides"] = output_dir / "override_log.csv"
        export_override_log(config, result, outputs["overrides"])

    # Orphan reports
    if result.orphan_candidates:
        outputs["orphan_report"] = output_dir / "orphan_candidates.csv"
        export_orphan_report(result, outputs["orphan_report"])

        outputs["orphan_summary"] = output_dir / "orphan_summary.csv"
        export_orphan_summary(result, outputs["orphan_summary"])

    # IEL reports
    if result.iel_candidates:
        outputs["iel_report"] = output_dir / "iel_candidates.csv"
        export_iel_report(result, outputs["iel_report"])

        outputs["iel_summary"] = output_dir / "iel_summary.csv"
        export_iel_summary(result, outputs["iel_summary"])

    # Harmonized summary (if harmonization was applied)
    if "cell_type_fine" in adata.obs.columns:
        outputs["harmonized_summary"] = output_dir / "harmonized_summary.csv"
        export_harmonized_summary(adata, outputs["harmonized_summary"])

    logger.info(f"All outputs saved to {output_dir}")

    return outputs
