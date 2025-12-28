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
            summary["cell_types"].append({
                "name": row["cell_type"],
                "total_count": int(row["total_count"]),
                "total_proportion": float(row["total_proportion"]),
                "mean_proportion": float(row["mean_proportion"]),
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
    output_dir: Path,
) -> Dict[str, Dict[str, Path]]:
    """Export results for multiple cell type columns.

    Parameters
    ----------
    result : MultiColumnResult
        Multi-column result
    output_dir : Path
        Output directory

    Returns
    -------
    Dict[str, Dict[str, Path]]
        Mapping of column name to output paths
    """
    output_dir = Path(output_dir)
    all_paths = {}

    for col_name, col_result in result.results.items():
        # Use column name as prefix
        prefix = f"{col_name}_" if len(result.columns_analyzed) > 1 else ""
        all_paths[col_name] = export_all(col_result, output_dir, prefix)

    # Export multi-column provenance
    prov_path = output_dir / "multi_column_provenance.json"
    with open(prov_path, "w") as f:
        json.dump(result.provenance, f, indent=2, default=str)

    return all_paths
