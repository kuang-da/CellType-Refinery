"""Export functions for composition analysis results.

This module provides functions to export composition analysis results
to various formats (CSV, JSON).
"""

import json
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Optional

import pandas as pd

from .engine import CompositionResult, MultiColumnResult


def export_composition_by_sample(
    result: CompositionResult,
    output_dir: Path,
    prefix: str = "",
) -> Path:
    """Export per-sample composition.

    Parameters
    ----------
    result : CompositionResult
        Composition result
    output_dir : Path
        Output directory
    prefix : str
        File prefix

    Returns
    -------
    Path
        Output file path
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    filename = f"{prefix}composition_by_sample.csv" if prefix else "composition_by_sample.csv"
    path = output_dir / filename
    result.composition_by_sample.to_csv(path, index=False)
    return path


def export_composition_by_region(
    result: CompositionResult,
    output_dir: Path,
    prefix: str = "",
) -> Path:
    """Export regional composition aggregation.

    Parameters
    ----------
    result : CompositionResult
        Composition result
    output_dir : Path
        Output directory
    prefix : str
        File prefix

    Returns
    -------
    Path
        Output file path
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    filename = f"{prefix}composition_by_region.csv" if prefix else "composition_by_region.csv"
    path = output_dir / filename
    result.composition_by_region.to_csv(path, index=False)
    return path


def export_composition_by_donor(
    result: CompositionResult,
    output_dir: Path,
    prefix: str = "",
) -> Path:
    """Export per-donor composition aggregation.

    Parameters
    ----------
    result : CompositionResult
        Composition result
    output_dir : Path
        Output directory
    prefix : str
        File prefix

    Returns
    -------
    Path
        Output file path
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    filename = f"{prefix}composition_by_donor.csv" if prefix else "composition_by_donor.csv"
    path = output_dir / filename
    result.composition_by_donor.to_csv(path, index=False)
    return path


def export_composition_global(
    result: CompositionResult,
    output_dir: Path,
    prefix: str = "",
) -> Path:
    """Export global composition summary.

    Parameters
    ----------
    result : CompositionResult
        Composition result
    output_dir : Path
        Output directory
    prefix : str
        File prefix

    Returns
    -------
    Path
        Output file path
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    filename = f"{prefix}composition_global.csv" if prefix else "composition_global.csv"
    path = output_dir / filename
    result.composition_global.to_csv(path, index=False)
    return path


def export_composition_wide(
    result: CompositionResult,
    output_dir: Path,
    prefix: str = "",
) -> Path:
    """Export wide-format composition matrix.

    Parameters
    ----------
    result : CompositionResult
        Composition result
    output_dir : Path
        Output directory
    prefix : str
        File prefix

    Returns
    -------
    Path
        Output file path
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    filename = f"{prefix}composition_wide.csv" if prefix else "composition_wide.csv"
    path = output_dir / filename
    result.composition_wide.to_csv(path)
    return path


def export_diversity_by_sample(
    result: CompositionResult,
    output_dir: Path,
    prefix: str = "",
) -> Path:
    """Export per-sample diversity metrics.

    Parameters
    ----------
    result : CompositionResult
        Composition result
    output_dir : Path
        Output directory
    prefix : str
        File prefix

    Returns
    -------
    Path
        Output file path
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    filename = f"{prefix}diversity_by_sample.csv" if prefix else "diversity_by_sample.csv"
    path = output_dir / filename
    result.diversity_by_sample.to_csv(path, index=False)
    return path


def export_diversity_summary(
    result: CompositionResult,
    output_dir: Path,
    prefix: str = "",
) -> Path:
    """Export diversity summary.

    Parameters
    ----------
    result : CompositionResult
        Composition result
    output_dir : Path
        Output directory
    prefix : str
        File prefix

    Returns
    -------
    Path
        Output file path
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    filename = f"{prefix}diversity_summary.csv" if prefix else "diversity_summary.csv"
    path = output_dir / filename
    result.diversity_summary.to_csv(path, index=False)
    return path


def export_biology_by_sample(
    result: CompositionResult,
    output_dir: Path,
    prefix: str = "",
) -> Optional[Path]:
    """Export per-sample biology metrics.

    Parameters
    ----------
    result : CompositionResult
        Composition result
    output_dir : Path
        Output directory
    prefix : str
        File prefix

    Returns
    -------
    Path or None
        Output file path or None if no biology data
    """
    if result.biology_by_sample is None:
        return None

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    filename = f"{prefix}biology_by_sample.csv" if prefix else "biology_by_sample.csv"
    path = output_dir / filename
    result.biology_by_sample.to_csv(path, index=False)
    return path


def export_biology_by_region(
    result: CompositionResult,
    output_dir: Path,
    prefix: str = "",
) -> Optional[Path]:
    """Export per-region biology metrics.

    Parameters
    ----------
    result : CompositionResult
        Composition result
    output_dir : Path
        Output directory
    prefix : str
        File prefix

    Returns
    -------
    Path or None
        Output file path or None if no biology data
    """
    if result.biology_by_region is None:
        return None

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    filename = f"{prefix}biology_by_region.csv" if prefix else "biology_by_region.csv"
    path = output_dir / filename
    result.biology_by_region.to_csv(path, index=False)
    return path


def export_enrichment(
    result: CompositionResult,
    output_dir: Path,
    prefix: str = "",
) -> Optional[Path]:
    """Export regional enrichment results.

    Parameters
    ----------
    result : CompositionResult
        Composition result
    output_dir : Path
        Output directory
    prefix : str
        File prefix

    Returns
    -------
    Path or None
        Output file path or None if no enrichment data
    """
    if result.enrichment is None or len(result.enrichment) == 0:
        return None

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    filename = f"{prefix}regional_enrichment.csv" if prefix else "regional_enrichment.csv"
    path = output_dir / filename
    result.enrichment.to_csv(path, index=False)
    return path


def export_enrichment_summary_region(
    result: CompositionResult,
    output_dir: Path,
    prefix: str = "",
) -> Optional[Path]:
    """Export enrichment summary by region.

    Parameters
    ----------
    result : CompositionResult
        Composition result
    output_dir : Path
        Output directory
    prefix : str
        File prefix

    Returns
    -------
    Path or None
        Output file path or None if no enrichment data
    """
    if result.enrichment_summary_region is None:
        return None

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    filename = f"{prefix}enrichment_summary_region.csv" if prefix else "enrichment_summary_region.csv"
    path = output_dir / filename
    result.enrichment_summary_region.to_csv(path, index=False)
    return path


def export_enrichment_summary_cell_type(
    result: CompositionResult,
    output_dir: Path,
    prefix: str = "",
) -> Optional[Path]:
    """Export enrichment summary by cell type.

    Parameters
    ----------
    result : CompositionResult
        Composition result
    output_dir : Path
        Output directory
    prefix : str
        File prefix

    Returns
    -------
    Path or None
        Output file path or None if no enrichment data
    """
    if result.enrichment_summary_cell_type is None:
        return None

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    filename = f"{prefix}enrichment_summary_cell_type.csv" if prefix else "enrichment_summary_cell_type.csv"
    path = output_dir / filename
    result.enrichment_summary_cell_type.to_csv(path, index=False)
    return path


def export_provenance(
    result: CompositionResult,
    output_dir: Path,
    prefix: str = "",
) -> Path:
    """Export provenance information.

    Parameters
    ----------
    result : CompositionResult
        Composition result
    output_dir : Path
        Output directory
    prefix : str
        File prefix

    Returns
    -------
    Path
        Output file path
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    filename = f"{prefix}provenance.json" if prefix else "provenance.json"
    path = output_dir / filename

    with open(path, "w") as f:
        json.dump(result.provenance, f, indent=2, default=str)

    return path


def export_summary_json(
    result: CompositionResult,
    output_dir: Path,
    prefix: str = "",
) -> Path:
    """Export structured summary for downstream tools.

    Parameters
    ----------
    result : CompositionResult
        Composition result
    output_dir : Path
        Output directory
    prefix : str
        File prefix

    Returns
    -------
    Path
        Output file path
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    filename = f"{prefix}composition_summary.json" if prefix else "composition_summary.json"
    path = output_dir / filename

    # Build summary structure
    summary = {
        "metadata": {
            "timestamp": datetime.now().isoformat(),
            "cell_type_column": result.cell_type_column,
            "n_cells": result.provenance.get("n_cells", 0),
            "n_samples": result.provenance.get("n_samples", 0),
            "n_cell_types": result.provenance.get("n_cell_types", 0),
        },
        "cell_types": [],
        "diversity": {},
        "enrichment_summary": {},
    }

    # Add cell type summary
    if len(result.composition_global) > 0:
        for _, row in result.composition_global.iterrows():
            # Support both old and new column names for backward compatibility
            global_pct = row.get("global_pct", row.get("total_proportion", 0))
            mean_pct = row.get("mean_pct", row.get("mean_proportion", 0))
            summary["cell_types"].append({
                "name": row["cell_type"],
                "total_count": int(row["total_count"]),
                "global_pct": float(global_pct),
                "mean_pct": float(mean_pct),
                "cv": float(row["cv"]) if pd.notna(row["cv"]) else None,
                "n_samples_present": int(row["n_samples_present"]),
            })

    # Add diversity summary
    if len(result.diversity_summary) > 0:
        for _, row in result.diversity_summary.iterrows():
            metric = row["metric"]
            summary["diversity"][metric] = {
                "mean": float(row["mean"]),
                "std": float(row["std"]) if pd.notna(row["std"]) else None,
                "median": float(row["median"]),
            }

    # Add enrichment summary
    if result.enrichment_summary_region is not None and len(result.enrichment_summary_region) > 0:
        for _, row in result.enrichment_summary_region.iterrows():
            summary["enrichment_summary"][row["region"]] = {
                "n_enriched": int(row["n_enriched"]),
                "n_depleted": int(row["n_depleted"]),
            }

    with open(path, "w") as f:
        json.dump(summary, f, indent=2)

    return path


def export_all(
    result: CompositionResult,
    output_dir: Path,
    prefix: str = "",
) -> Dict[str, Path]:
    """Export all composition results.

    Parameters
    ----------
    result : CompositionResult
        Composition result
    output_dir : Path
        Output directory
    prefix : str
        File prefix

    Returns
    -------
    Dict[str, Path]
        Mapping of output name to file path
    """
    output_dir = Path(output_dir)
    paths = {}

    paths["composition_by_sample"] = export_composition_by_sample(result, output_dir, prefix)
    paths["composition_global"] = export_composition_global(result, output_dir, prefix)
    paths["composition_wide"] = export_composition_wide(result, output_dir, prefix)
    paths["diversity_by_sample"] = export_diversity_by_sample(result, output_dir, prefix)
    paths["diversity_summary"] = export_diversity_summary(result, output_dir, prefix)
    paths["provenance"] = export_provenance(result, output_dir, prefix)
    paths["summary_json"] = export_summary_json(result, output_dir, prefix)

    # Optional outputs
    if len(result.composition_by_region) > 0:
        paths["composition_by_region"] = export_composition_by_region(result, output_dir, prefix)

    if len(result.composition_by_donor) > 0:
        paths["composition_by_donor"] = export_composition_by_donor(result, output_dir, prefix)

    biology_sample = export_biology_by_sample(result, output_dir, prefix)
    if biology_sample:
        paths["biology_by_sample"] = biology_sample

    biology_region = export_biology_by_region(result, output_dir, prefix)
    if biology_region:
        paths["biology_by_region"] = biology_region

    enrichment = export_enrichment(result, output_dir, prefix)
    if enrichment:
        paths["enrichment"] = enrichment

    enrichment_region = export_enrichment_summary_region(result, output_dir, prefix)
    if enrichment_region:
        paths["enrichment_summary_region"] = enrichment_region

    enrichment_ct = export_enrichment_summary_cell_type(result, output_dir, prefix)
    if enrichment_ct:
        paths["enrichment_summary_cell_type"] = enrichment_ct

    return paths


def export_all_multi(
    result: MultiColumnResult,
    config: Any,  # Accept for API compatibility (may not use all fields)
    input_path: Path,
    output_dir: Path,
) -> Dict[str, Dict[str, Path]]:
    """Export results for multiple cell type columns.

    Creates subdirectory per column for organized output structure.

    Parameters
    ----------
    result : MultiColumnResult
        Multi-column result
    config : CompositionConfig
        Configuration used (for provenance)
    input_path : Path
        Path to input file (for provenance)
    output_dir : Path
        Output directory

    Returns
    -------
    Dict[str, Dict[str, Path]]
        Mapping of column name to output paths
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    all_paths = {}

    # Export per-column outputs into subdirectories (matching reference pattern)
    for col_name, col_result in result.results.items():
        col_dir = output_dir / col_name
        col_dir.mkdir(parents=True, exist_ok=True)
        all_paths[col_name] = export_all(col_result, col_dir, prefix="")

    # Build comparison metrics across columns (required for dashboard)
    comparison: Dict[str, Dict[str, Any]] = {
        "n_cell_types": {},
        "mean_shannon_entropy": {},
        "mean_simpson_index": {},
        "mean_evenness": {},
        "n_significant_enrichments": {},
    }

    for col, res in result.results.items():
        # n_cell_types from provenance or composition_global
        comparison["n_cell_types"][col] = res.provenance.get(
            "n_cell_types",
            len(res.composition_global) if res.composition_global is not None else 0
        )

        # Diversity metrics from diversity_summary DataFrame
        if res.diversity_summary is not None and len(res.diversity_summary) > 0:
            ds = res.diversity_summary
            shannon_row = ds[ds["metric"] == "shannon_entropy"]
            simpson_row = ds[ds["metric"] == "simpson_index"]
            evenness_row = ds[ds["metric"] == "evenness"]

            comparison["mean_shannon_entropy"][col] = round(
                float(shannon_row["mean"].iloc[0]) if len(shannon_row) > 0 else 0, 3
            )
            comparison["mean_simpson_index"][col] = round(
                float(simpson_row["mean"].iloc[0]) if len(simpson_row) > 0 else 0, 3
            )
            comparison["mean_evenness"][col] = round(
                float(evenness_row["mean"].iloc[0]) if len(evenness_row) > 0 else 0, 3
            )
        else:
            comparison["mean_shannon_entropy"][col] = 0
            comparison["mean_simpson_index"][col] = 0
            comparison["mean_evenness"][col] = 0

        # n_significant_enrichments from enrichment DataFrame
        if res.enrichment is not None and len(res.enrichment) > 0:
            comparison["n_significant_enrichments"][col] = int(
                res.enrichment["significant"].sum() if "significant" in res.enrichment.columns else 0
            )
        else:
            comparison["n_significant_enrichments"][col] = 0

    # Export multi-column summary at root
    summary = {
        "module": "celltype_refinery.core.composition",
        "version": "1.1.0",
        "timestamp": datetime.now().isoformat(),
        "mode": "multi_column",
        "columns_processed": result.columns_analyzed,
        "columns_analyzed": result.columns_analyzed,  # Backward compatibility
        "columns_skipped": getattr(result, "columns_skipped", []),
        "input_path": str(input_path),
        "comparison": comparison,
        "per_column_summary": {
            col: {
                "n_cells": res.provenance.get("n_cells", 0),
                "n_cell_types": res.provenance.get("n_cell_types", 0),
                "n_samples": res.provenance.get("n_samples", 0),
            }
            for col, res in result.results.items()
        },
        "total_execution_time_seconds": sum(
            res.provenance.get("duration_seconds", 0) for res in result.results.values()
        ),
    }
    summary_path = output_dir / "multi_column_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2, default=str)

    all_paths["_summary"] = {"multi_column_summary": summary_path}

    return all_paths
