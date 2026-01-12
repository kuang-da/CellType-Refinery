"""Preprocessing module CLI runner.

Enables running preprocessing stages as:
    python -m celltype_refinery.core.preprocessing --stage A --metadata <path> --output <dir>
    python -m celltype_refinery.core.preprocessing --stage B --input <path> --output <dir>
    python -m celltype_refinery.core.preprocessing --stage C --input <path> --output <dir>
    python -m celltype_refinery.core.preprocessing --stage D --input <path> --output <dir>
    python -m celltype_refinery.core.preprocessing --stage E --input <path> --output <dir>

Usage Examples:
    # Run Stage A on fallopian tube data
    python -m celltype_refinery.core.preprocessing --stage A \
        --metadata data/fallopian_tube/metadata.csv \
        --output output/fallopian_tube/stage_a

    # Run Stage B (cell QC) on Stage A output
    python -m celltype_refinery.core.preprocessing --stage B \
        --input output/fallopian_tube/stage_a/metadata_verified.csv \
        --output output/fallopian_tube/stage_b

    # Run Stage C (normalization) on Stage B filtered matrices
    python -m celltype_refinery.core.preprocessing --stage C \
        --input output/fallopian_tube/stage_b/filtered/matrices \
        --output output/fallopian_tube/stage_c

    # Run Stage D (cross-sample alignment) on Stage C export
    python -m celltype_refinery.core.preprocessing --stage D \
        --input output/fallopian_tube/stage_c/export \
        --output output/fallopian_tube/stage_d \
        --variant clip-Q3__log1p

    # Run Stage E (batch correction) on Stage D aligned data
    python -m celltype_refinery.core.preprocessing --stage E \
        --input output/fallopian_tube/stage_d/aligned \
        --metadata output/fallopian_tube/stage_a/metadata_verified.csv \
        --channels-dir data/fallopian_tube/channels \
        --output output/fallopian_tube/stage_e

    # With verbose logging and log file
    python -m celltype_refinery.core.preprocessing --stage E \
        --input output/fallopian_tube/stage_d/aligned \
        --output output/fallopian_tube/stage_e \
        --log-dir logs/fallopian_tube \
        --verbose
"""

import argparse
import logging
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd
import scipy.sparse as sp

import re

from .config import LoaderConfig, QCConfig, NormalizationConfig, AlignmentConfig, BatchCorrectionConfig, PreprocessingConfig
from .loader import DataLoader, LoadResult
from .qc import CellQC, QCResult, REASON_COLUMNS
from .normalization import Normalizer, NormalizationResult, TRANSFORMS
from .alignment import CrossSampleAligner, AlignmentResult, AlignmentParams
from .batch import BatchCorrector, BatchCorrectionResult, BatchDiagnostics



# Canonical FT region order (anatomical: distal to proximal, from ovary toward uterus)
FT_REGION_ORDER: List[str] = ["fimbriae", "infundibulum", "ampulla", "isthmus"]


def sort_regions_canonical(regions: List[str], region_order: Optional[List[str]] = None) -> List[str]:
    """Sort regions by canonical anatomical order.

    Parameters
    ----------
    regions : List[str]
        Region names to sort
    region_order : Optional[List[str]]
        Custom region order (default: FT_REGION_ORDER)

    Returns
    -------
    List[str]
        Sorted region names
    """
    if region_order is None:
        region_order = FT_REGION_ORDER

    order_map = {r.lower(): i for i, r in enumerate(region_order)}

    def get_order(r: str) -> int:
        return order_map.get(r.lower(), len(region_order))

    return sorted(regions, key=get_order)


def setup_logging(
    verbose: bool = False,
    log_dir: Optional[Path] = None,
    log_filename: str = "stage_a.log",
) -> logging.Logger:
    """Configure logging for Stage A.

    Parameters
    ----------
    verbose : bool
        Enable debug-level logging
    log_dir : Optional[Path]
        Directory for log file (if None, console only)
    log_filename : str
        Log file name

    Returns
    -------
    logging.Logger
        Configured logger
    """
    level = logging.DEBUG if verbose else logging.INFO
    log_format = "%(asctime)s [%(levelname)s] %(message)s"
    date_format = "%Y-%m-%d %H:%M:%S"

    # Create logger
    logger = logging.getLogger("celltype_refinery.preprocessing")
    logger.setLevel(level)

    # Clear existing handlers
    logger.handlers.clear()

    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(level)
    console_handler.setFormatter(logging.Formatter(log_format, date_format))
    logger.addHandler(console_handler)

    # File handler (if log_dir provided)
    if log_dir is not None:
        log_dir.mkdir(parents=True, exist_ok=True)
        log_path = log_dir / log_filename
        file_handler = logging.FileHandler(log_path, mode="w")
        file_handler.setLevel(level)
        file_handler.setFormatter(logging.Formatter(log_format, date_format))
        logger.addHandler(file_handler)
        logger.info(f"Log file: {log_path}")

    return logger


def parse_issue_count(issues: List[str], prefix: str) -> int:
    """Extract count from issue string like 'prefix:count'.

    Parameters
    ----------
    issues : List[str]
        List of issue strings
    prefix : str
        Prefix to match (e.g., 'duplicate_ids_matrix')

    Returns
    -------
    int
        Count if found, else 0
    """
    for issue in issues:
        if issue.startswith(prefix + ":"):
            try:
                return int(issue.split(":")[1])
            except (IndexError, ValueError):
                return 0
    return 0


def load_result_to_row(
    result: LoadResult,
    registry_row: pd.Series,
) -> Dict[str, Any]:
    """Convert LoadResult to metadata_verified row format.

    Parameters
    ----------
    result : LoadResult
        Validation result from DataLoader
    registry_row : pd.Series
        Original registry row

    Returns
    -------
    Dict[str, Any]
        Row for metadata_verified.csv
    """
    # Start with original metadata columns
    row = registry_row.to_dict()

    # Add validation metrics
    row["matrix_rows"] = result.n_cells
    row["matrix_markers"] = len(result.markers)

    # Get metadata row count
    metadata_rows = 0
    if result.cell_metadata is not None:
        metadata_rows = len(result.cell_metadata)
    row["metadata_rows"] = metadata_rows

    # Parse issue counts
    row["duplicate_ids_matrix"] = parse_issue_count(result.issues, "duplicate_ids_matrix")
    row["duplicate_ids_metadata"] = parse_issue_count(result.issues, "duplicate_ids_metadata")
    row["missing_in_matrix"] = parse_issue_count(result.issues, "missing_in_matrix")
    row["missing_in_metadata"] = parse_issue_count(result.issues, "missing_in_metadata")

    # Status and issues
    row["status"] = result.status
    row["issues"] = ";".join(result.issues) if result.issues else ""

    return row


def compute_marker_coverage(
    cell_matrix: pd.DataFrame,
    markers: Set[str],
    cell_id_col: str = "cell_mask_id",
) -> Dict[str, float]:
    """Compute fraction of positive cells per marker.

    Parameters
    ----------
    cell_matrix : pd.DataFrame
        Cell-by-marker matrix
    markers : Set[str]
        Marker column names
    cell_id_col : str
        Cell ID column name to exclude

    Returns
    -------
    Dict[str, float]
        Marker -> fraction of positive cells
    """
    if cell_matrix is None or len(cell_matrix) == 0:
        return {}

    marker_cols = [c for c in cell_matrix.columns if c in markers]
    if not marker_cols:
        return {}

    n_cells = len(cell_matrix)
    coverage = (cell_matrix[marker_cols] > 0).sum() / n_cells
    return coverage.to_dict()


# =============================================================================
# Visualization Functions
# =============================================================================


def build_antibody_coverage_heatmap(
    coverage_records: List[Dict[str, Any]],
    output_dir: Path,
    dpi: int = 200,
) -> Optional[Path]:
    """Generate heatmap of antibody coverage across samples.

    Parameters
    ----------
    coverage_records : List[Dict]
        List of {sample_id, marker, coverage} records
    output_dir : Path
        Output directory for figures
    dpi : int
        Figure resolution

    Returns
    -------
    Optional[Path]
        Path to generated figure, or None if skipped
    """
    if not coverage_records:
        return None

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logging.warning("matplotlib not available, skipping heatmap")
        return None

    coverage_df = pd.DataFrame(coverage_records)
    if coverage_df.empty:
        return None

    # Pivot to sample × marker matrix
    pivot = coverage_df.pivot(
        index="sample_id", columns="marker", values="coverage"
    ).fillna(0)
    pivot = pivot.sort_index()

    # Create figure
    n_samples, n_markers = pivot.shape
    fig_width = max(12, n_markers * 0.3)
    fig_height = max(8, n_samples * 0.25)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    im = ax.imshow(pivot.values, aspect="auto", cmap="viridis", vmin=0, vmax=1)

    # Labels
    ax.set_xticks(range(n_markers))
    ax.set_xticklabels(pivot.columns, rotation=90, ha="center", fontsize=8)
    ax.set_yticks(range(n_samples))
    ax.set_yticklabels(pivot.index, fontsize=8)

    ax.set_xlabel("Marker")
    ax.set_ylabel("Sample")
    ax.set_title("Antibody Coverage (fraction of positive cells)")

    # Colorbar
    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("Fraction")

    fig.tight_layout()

    # Save
    output_path = output_dir / "stage_a_antibody_coverage.png"
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return output_path


def build_cell_count_region_plot(
    verified_df: pd.DataFrame,
    output_dir: Path,
    dpi: int = 200,
) -> Optional[Path]:
    """Generate boxplot of cell counts per region.

    Parameters
    ----------
    verified_df : pd.DataFrame
        Verified metadata with matrix_rows and region columns
    output_dir : Path
        Output directory for figures
    dpi : int
        Figure resolution

    Returns
    -------
    Optional[Path]
        Path to generated figure, or None if skipped
    """
    if verified_df.empty or "region" not in verified_df.columns:
        return None

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logging.warning("matplotlib not available, skipping boxplot")
        return None

    counts_df = verified_df.dropna(subset=["matrix_rows"]).copy()
    if counts_df.empty:
        return None

    counts_df["matrix_rows"] = pd.to_numeric(counts_df["matrix_rows"], errors="coerce")
    counts_df = counts_df.dropna(subset=["matrix_rows"])
    if counts_df.empty:
        return None

    available_regions = list(counts_df["region"].dropna().unique())
    if not available_regions:
        return None

    # Use canonical anatomical order (distal to proximal)
    regions = sort_regions_canonical(available_regions)

    data = [
        counts_df.loc[counts_df["region"] == region, "matrix_rows"].astype(float).values
        for region in regions
    ]

    # Create figure
    fig, ax = plt.subplots(figsize=(max(5, len(regions) * 1.2), 4))

    bp = ax.boxplot(data, tick_labels=regions, vert=True, patch_artist=True)

    # Color the boxes using region-specific colors
    region_colors = {
        "fimbriae": "#9b59b6",      # Purple
        "infundibulum": "#2ecc71",  # Green
        "ampulla": "#3498db",       # Blue
        "isthmus": "#e74c3c",       # Red
    }
    for patch, region in zip(bp["boxes"], regions):
        color = region_colors.get(region.lower(), "#95a5a6")
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    ax.set_title("Cell count per region")
    ax.set_xlabel("Region")
    ax.set_ylabel("Cells per sample")
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: format(int(x), ",")))

    fig.tight_layout()

    # Save
    output_path = output_dir / "stage_a_cell_count_region.png"
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return output_path


def build_donor_region_barplot(
    verified_df: pd.DataFrame,
    output_dir: Path,
    dpi: int = 200,
) -> Optional[Path]:
    """Generate grouped bar chart of cell counts by donor and region.

    Parameters
    ----------
    verified_df : pd.DataFrame
        Verified metadata with matrix_rows, donor, and region columns
    output_dir : Path
        Output directory for figures
    dpi : int
        Figure resolution

    Returns
    -------
    Optional[Path]
        Path to generated figure, or None if skipped
    """
    if verified_df.empty or not {"donor", "region"}.issubset(verified_df.columns):
        return None

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logging.warning("matplotlib not available, skipping bar plot")
        return None

    counts_df = verified_df.dropna(subset=["matrix_rows", "donor", "region"]).copy()
    if counts_df.empty:
        return None

    counts_df["matrix_rows"] = pd.to_numeric(counts_df["matrix_rows"], errors="coerce")
    counts_df = counts_df.dropna(subset=["matrix_rows"])
    if counts_df.empty:
        return None

    counts_df["donor"] = counts_df["donor"].astype(str)
    counts_df["region"] = counts_df["region"].astype(str)

    donors = sorted(counts_df["donor"].unique())
    available_regions = list(counts_df["region"].unique())
    if not donors or not available_regions:
        return None

    # Use canonical anatomical order (distal to proximal)
    regions = sort_regions_canonical(available_regions)

    # Pivot: donor × region
    pivot = (
        counts_df.groupby(["donor", "region"], dropna=False)["matrix_rows"]
        .sum()
        .unstack(fill_value=0)
        .reindex(index=donors, columns=regions, fill_value=0)
    )
    if pivot.empty:
        return None

    donors = list(pivot.index)
    regions = list(pivot.columns)

    # Region-specific colors (matching boxplot)
    region_colors = {
        "fimbriae": "#9b59b6",      # Purple
        "infundibulum": "#2ecc71",  # Green
        "ampulla": "#3498db",       # Blue
        "isthmus": "#e74c3c",       # Red
    }

    # Create figure
    fig_width = max(6, len(donors) * 1.4)
    fig, ax = plt.subplots(figsize=(fig_width, 5))

    x_positions = np.arange(len(donors))
    bar_width = 0.8 / max(1, len(regions))

    for idx, region in enumerate(regions):
        offset = (idx - (len(regions) - 1) / 2) * bar_width
        heights = pivot[region].to_numpy()
        color = region_colors.get(region.lower(), "#95a5a6")
        ax.bar(
            x_positions + offset,
            heights,
            bar_width,
            label=region,
            color=color,
        )

    ax.set_title("Cell counts per donor and region")
    ax.set_xlabel("Donor")
    ax.set_ylabel("Cells per sample")
    ax.set_xticks(x_positions)
    ax.set_xticklabels(donors, rotation=30, ha="right")
    ax.legend(title="Region", loc="upper right")
    ax.set_ylim(bottom=0)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: format(int(x), ",")))

    fig.tight_layout()

    # Save
    output_path = output_dir / "stage_a_donor_region_cell_counts.png"
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return output_path


def generate_figures(
    verified_df: pd.DataFrame,
    coverage_records: List[Dict[str, Any]],
    output_dir: Path,
    dpi: int = 200,
    logger: Optional[logging.Logger] = None,
) -> List[Path]:
    """Generate all Stage A visualization figures.

    Parameters
    ----------
    verified_df : pd.DataFrame
        Verified metadata
    coverage_records : List[Dict]
        Marker coverage records
    output_dir : Path
        Output directory for figures
    dpi : int
        Figure resolution
    logger : Optional[Logger]
        Logger instance

    Returns
    -------
    List[Path]
        Paths to generated figures
    """
    logger = logger or logging.getLogger(__name__)
    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)

    generated = []

    # 1. Antibody coverage heatmap
    path = build_antibody_coverage_heatmap(coverage_records, figure_dir, dpi)
    if path:
        logger.info(f"Generated: {path.name}")
        generated.append(path)

    # 2. Cell count by region boxplot
    path = build_cell_count_region_plot(verified_df, figure_dir, dpi)
    if path:
        logger.info(f"Generated: {path.name}")
        generated.append(path)

    # 3. Donor × region bar chart
    path = build_donor_region_barplot(verified_df, figure_dir, dpi)
    if path:
        logger.info(f"Generated: {path.name}")
        generated.append(path)

    return generated


# =============================================================================
# Stage B Visualization Functions
# =============================================================================


# Region-specific colors (consistent with Stage A)
FT_REGION_COLORS: Dict[str, str] = {
    "fimbriae": "#9b59b6",      # Purple
    "infundibulum": "#2ecc71",  # Green
    "ampulla": "#3498db",       # Blue
    "isthmus": "#e74c3c",       # Red
}


def build_area_distribution_summary(
    area_data: Dict[str, np.ndarray],
    output_dir: Path,
    title: str = "Cell Area Distribution",
    dpi: int = 200,
) -> Optional[Path]:
    """Generate overlaid KDE plots of area distributions.

    Parameters
    ----------
    area_data : Dict[str, np.ndarray]
        Sample ID -> area values
    output_dir : Path
        Output directory
    title : str
        Plot title
    dpi : int
        Figure resolution

    Returns
    -------
    Optional[Path]
        Path to generated figure
    """
    if not area_data:
        return None

    try:
        import matplotlib.pyplot as plt
        from scipy.stats import gaussian_kde
    except ImportError:
        return None

    fig, ax = plt.subplots(figsize=(10, 6))

    # Sample a few distributions for overlay
    samples_to_plot = list(area_data.keys())[:8]  # Limit overlay count

    for sample_id in samples_to_plot:
        values = area_data[sample_id]
        if len(values) < 10:
            continue
        # Filter positive values for KDE
        values = values[values > 0]
        if len(values) < 10:
            continue
        try:
            kde = gaussian_kde(values)
            x_range = np.linspace(values.min(), np.percentile(values, 99), 200)
            ax.plot(x_range, kde(x_range), alpha=0.6, label=sample_id)
        except Exception:
            continue

    ax.set_xlabel("Cell Area")
    ax.set_ylabel("Density")
    ax.set_title(title)
    if len(samples_to_plot) <= 10:
        ax.legend(fontsize=8)

    fig.tight_layout()
    output_path = output_dir / "stage_b_area_distribution_summary.png"
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return output_path


def build_area_by_region_plot(
    area_df: pd.DataFrame,
    output_dir: Path,
    suffix: str = "",
    dpi: int = 200,
) -> List[Path]:
    """Generate area distribution plots by region.

    Parameters
    ----------
    area_df : pd.DataFrame
        DataFrame with cell_area, region, donor columns
    output_dir : Path
        Output directory
    suffix : str
        Filename suffix (e.g., "_post_qc")
    dpi : int
        Figure resolution

    Returns
    -------
    List[Path]
        Paths to generated figures
    """
    if area_df.empty or "region" not in area_df.columns:
        return []

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return []

    generated = []
    available_regions = sort_regions_canonical(list(area_df["region"].dropna().unique()))

    for region in available_regions:
        region_data = area_df[area_df["region"] == region].copy()
        if region_data.empty:
            continue

        fig, ax = plt.subplots(figsize=(8, 5))

        # Box plot by donor if available
        if "donor" in region_data.columns:
            donors = sorted(region_data["donor"].dropna().unique())
            data = [
                region_data.loc[region_data["donor"] == d, "cell_area"].dropna().values
                for d in donors
            ]
            bp = ax.boxplot(data, tick_labels=donors, patch_artist=True)
            color = FT_REGION_COLORS.get(region.lower(), "#95a5a6")
            for patch in bp["boxes"]:
                patch.set_facecolor(color)
                patch.set_alpha(0.7)
            ax.set_xlabel("Donor")
        else:
            ax.hist(region_data["cell_area"].dropna(), bins=50, alpha=0.7,
                    color=FT_REGION_COLORS.get(region.lower(), "#95a5a6"))
            ax.set_xlabel("Cell Area")

        ax.set_ylabel("Count" if "donor" not in region_data.columns else "Cell Area")
        ax.set_title(f"Cell Area Distribution - {region.title()}{suffix}")

        fig.tight_layout()
        output_path = output_dir / f"stage_b_area_distribution_region_{region}{suffix}.png"
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        generated.append(output_path)

    return generated


def build_area_vs_intensity_scatter(
    scatter_df: pd.DataFrame,
    output_dir: Path,
    area_threshold: Optional[float] = None,
    intensity_threshold: Optional[float] = None,
    dpi: int = 200,
) -> Optional[Path]:
    """Generate scatter plot of cell area vs total intensity.

    Parameters
    ----------
    scatter_df : pd.DataFrame
        DataFrame with cell_area, total_intensity, removed columns
    output_dir : Path
        Output directory
    area_threshold : Optional[float]
        High area threshold line
    intensity_threshold : Optional[float]
        Low intensity threshold line
    dpi : int
        Figure resolution

    Returns
    -------
    Optional[Path]
        Path to generated figure
    """
    if scatter_df.empty:
        return None

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return None

    fig, ax = plt.subplots(figsize=(10, 8))

    # Separate kept and removed cells
    kept = scatter_df[~scatter_df["removed"]]
    removed = scatter_df[scatter_df["removed"]]

    ax.scatter(kept["cell_area"], kept["total_intensity"],
               alpha=0.3, s=1, c="#3498db", label="Kept")
    if not removed.empty:
        ax.scatter(removed["cell_area"], removed["total_intensity"],
                   alpha=0.5, s=2, c="#e74c3c", label="Removed")

    # Threshold lines
    if area_threshold is not None and not np.isnan(area_threshold):
        ax.axvline(area_threshold, color="#e67e22", linestyle="--",
                   label=f"High area threshold ({area_threshold:.0f})")

    if intensity_threshold is not None and not np.isnan(intensity_threshold):
        ax.axhline(intensity_threshold, color="#9b59b6", linestyle="--",
                   label=f"Low intensity threshold ({intensity_threshold:.0f})")

    ax.set_xlabel("Cell Area")
    ax.set_ylabel("Total Fluorescence Intensity")
    ax.set_title("Cell Area vs Total Intensity")
    ax.legend(markerscale=4)

    fig.tight_layout()
    output_path = output_dir / "stage_b_area_vs_intensity.png"
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return output_path


def build_qc_removal_summary(
    summary_df: pd.DataFrame,
    output_dir: Path,
    dpi: int = 200,
) -> Optional[Path]:
    """Generate bar chart of QC removal counts by reason.

    Parameters
    ----------
    summary_df : pd.DataFrame
        QC summary per sample
    output_dir : Path
        Output directory
    dpi : int
        Figure resolution

    Returns
    -------
    Optional[Path]
        Path to generated figure
    """
    if summary_df.empty:
        return None

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return None

    # Aggregate by reason
    reason_totals = {}
    for reason in REASON_COLUMNS:
        col = f"removed_{reason}"
        if col in summary_df.columns:
            reason_totals[reason] = int(summary_df[col].fillna(0).sum())

    if not reason_totals or sum(reason_totals.values()) == 0:
        return None

    labels = list(reason_totals.keys())
    values = [reason_totals[r] for r in labels]
    pretty_labels = [r.replace("_", " ").title() for r in labels]

    fig, ax = plt.subplots(figsize=(10, 5))
    bars = ax.bar(range(len(labels)), values, color="#3498db", alpha=0.8)

    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(pretty_labels, rotation=45, ha="right")
    ax.set_ylabel("Cells Removed")
    ax.set_title("QC Removal Summary by Reason")

    # Add value labels on bars
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                f"{val:,}", ha="center", va="bottom", fontsize=9)

    fig.tight_layout()
    output_path = output_dir / "stage_b_qc_removal_summary.png"
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return output_path


def build_removal_by_reason_plots(
    summary_df: pd.DataFrame,
    output_dir: Path,
    dpi: int = 200,
) -> List[Path]:
    """Generate per-reason bar charts showing removal counts per sample.

    Parameters
    ----------
    summary_df : pd.DataFrame
        QC summary per sample
    output_dir : Path
        Output directory
    dpi : int
        Figure resolution

    Returns
    -------
    List[Path]
        Paths to generated figures
    """
    if summary_df.empty:
        return []

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return []

    generated = []
    sample_ids = summary_df["sample_id"].astype(str).tolist()
    totals = summary_df["cells_total"].fillna(0).astype(int).tolist()

    for reason in REASON_COLUMNS:
        col = f"removed_{reason}"
        if col not in summary_df.columns:
            continue

        counts = summary_df[col].fillna(0).astype(int)
        if counts.sum() == 0:
            continue

        fig, ax = plt.subplots(figsize=(max(8, len(sample_ids) * 0.3), 5))

        # Calculate fractions
        fractions = [c / t * 100 if t > 0 else 0 for c, t in zip(counts, totals)]

        bars = ax.bar(range(len(sample_ids)), counts, color="#3498db", alpha=0.8)

        ax.set_xticks(range(len(sample_ids)))
        ax.set_xticklabels(sample_ids, rotation=90, fontsize=7)
        ax.set_ylabel("Cells Removed")
        pretty_reason = reason.replace("_", " ").title()
        ax.set_title(f"Cells Removed Per Sample - {pretty_reason}")

        fig.tight_layout()
        output_path = output_dir / f"stage_b_removed_{reason}.png"
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        generated.append(output_path)

    return generated


def generate_stage_b_figures(
    summary_df: pd.DataFrame,
    area_data: Dict[str, np.ndarray],
    scatter_df: pd.DataFrame,
    area_metadata_df: pd.DataFrame,
    filtered_area_metadata_df: pd.DataFrame,
    output_dir: Path,
    area_threshold: Optional[float] = None,
    intensity_threshold: Optional[float] = None,
    dpi: int = 200,
    logger: Optional[logging.Logger] = None,
) -> List[Path]:
    """Generate all Stage B visualization figures.

    Parameters
    ----------
    summary_df : pd.DataFrame
        QC summary per sample
    area_data : Dict[str, np.ndarray]
        Sample ID -> area values
    scatter_df : pd.DataFrame
        Scatter plot data
    area_metadata_df : pd.DataFrame
        Pre-QC area metadata
    filtered_area_metadata_df : pd.DataFrame
        Post-QC area metadata
    output_dir : Path
        Output directory
    area_threshold : Optional[float]
        High area threshold
    intensity_threshold : Optional[float]
        Low intensity threshold
    dpi : int
        Figure resolution
    logger : Optional[Logger]
        Logger instance

    Returns
    -------
    List[Path]
        Paths to generated figures
    """
    logger = logger or logging.getLogger(__name__)
    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)

    generated = []

    # 1. Area distribution summary
    path = build_area_distribution_summary(area_data, figure_dir, dpi=dpi)
    if path:
        logger.info(f"Generated: {path.name}")
        generated.append(path)

    # 2. Area distributions by region (pre-QC)
    paths = build_area_by_region_plot(area_metadata_df, figure_dir, suffix="", dpi=dpi)
    for p in paths:
        logger.info(f"Generated: {p.name}")
    generated.extend(paths)

    # 3. Area distributions by region (post-QC)
    paths = build_area_by_region_plot(filtered_area_metadata_df, figure_dir, suffix="_post_qc", dpi=dpi)
    for p in paths:
        logger.info(f"Generated: {p.name}")
    generated.extend(paths)

    # 4. Area vs intensity scatter
    path = build_area_vs_intensity_scatter(
        scatter_df, figure_dir,
        area_threshold=area_threshold,
        intensity_threshold=intensity_threshold,
        dpi=dpi,
    )
    if path:
        logger.info(f"Generated: {path.name}")
        generated.append(path)

    # 5. QC removal summary
    path = build_qc_removal_summary(summary_df, figure_dir, dpi=dpi)
    if path:
        logger.info(f"Generated: {path.name}")
        generated.append(path)

    # 6. Per-reason removal plots
    paths = build_removal_by_reason_plots(summary_df, figure_dir, dpi=dpi)
    for p in paths:
        logger.info(f"Generated: {p.name}")
    generated.extend(paths)

    return generated


# =============================================================================
# Stage B Runner
# =============================================================================


def run_stage_b(
    input_path: Path,
    output_dir: Path,
    config: Optional[QCConfig] = None,
    loader_config: Optional[LoaderConfig] = None,
    verbose: bool = False,
    skip_figures: bool = False,
    dpi: int = 200,
    log_dir: Optional[Path] = None,
) -> pd.DataFrame:
    """Run Stage B: Cell-level quality control.

    Parameters
    ----------
    input_path : Path
        Path to Stage A metadata_verified.csv
    output_dir : Path
        Output directory
    config : Optional[QCConfig]
        QC configuration
    loader_config : Optional[LoaderConfig]
        Loader configuration
    verbose : bool
        Enable verbose logging
    skip_figures : bool
        Skip figure generation
    dpi : int
        Figure resolution
    log_dir : Optional[Path]
        Directory for log file

    Returns
    -------
    pd.DataFrame
        QC summary per sample
    """
    logger = setup_logging(verbose, log_dir=log_dir, log_filename="stage_b.log")
    logger.info("Stage B: Cell-Level Quality Control")
    logger.info(f"Input: {input_path}")
    logger.info(f"Output: {output_dir}")

    # Setup directories
    output_dir.mkdir(parents=True, exist_ok=True)
    filtered_dir = output_dir / "filtered"
    filtered_dir.mkdir(exist_ok=True)
    filtered_matrix_dir = filtered_dir / "matrices"
    filtered_matrix_dir.mkdir(exist_ok=True)
    filtered_metadata_dir = filtered_dir / "metadata"
    filtered_metadata_dir.mkdir(exist_ok=True)

    # Load Stage A output
    metadata_df = pd.read_csv(input_path)
    if "sample_id" not in metadata_df.columns and "exp_id" in metadata_df.columns:
        metadata_df["sample_id"] = metadata_df["exp_id"]
    metadata_df["sample_id"] = metadata_df["sample_id"].astype(str)

    n_samples = len(metadata_df)
    logger.info(f"Found {n_samples} samples")

    # Initialize processors
    config = config or QCConfig()
    loader_config = loader_config or LoaderConfig()
    loader = DataLoader(loader_config)
    qc = CellQC(config)

    # Tracking variables
    summary_records: List[Dict[str, Any]] = []
    all_removals: List[Dict[str, Any]] = []
    area_data: Dict[str, np.ndarray] = {}
    area_metadata_frames: List[pd.DataFrame] = []
    filtered_area_metadata_frames: List[pd.DataFrame] = []
    scatter_records: List[Dict[str, Any]] = []
    area_thresholds: List[float] = []
    intensity_thresholds: List[float] = []

    rng = np.random.default_rng(42)

    # Process each sample
    for idx, row in metadata_df.iterrows():
        sample_id = str(row.get("sample_id", row.get("exp_id", f"sample_{idx}")))
        matrix_path = Path(row["cell_matrix_path"])
        metadata_path = Path(row["cell_metadata_path"])

        logger.info(f"[{idx+1}/{n_samples}] Processing {sample_id}")

        try:
            # Load data
            cell_matrix = loader.load_cell_matrix(matrix_path)
            cell_metadata = loader.load_cell_metadata(metadata_path)

            # Run QC
            result = qc.filter_sample(
                cell_matrix=cell_matrix,
                cell_metadata=cell_metadata,
                sample_id=sample_id,
                cell_id_col=loader_config.cell_id_col,
            )

            # Get area values for plots
            area_col = config.area_col
            if result.filtered_metadata is not None:
                if area_col in cell_metadata.columns:
                    area_values = pd.to_numeric(cell_metadata[area_col], errors="coerce").dropna().values
                    area_data[sample_id] = area_values

                    # Pre-QC area metadata
                    area_frame = pd.DataFrame({
                        "cell_area": area_values,
                        "sample_id": sample_id,
                        "donor": row.get("donor"),
                        "region": row.get("region"),
                    })
                    area_metadata_frames.append(area_frame)

                    # Post-QC area metadata
                    filtered_area = pd.to_numeric(
                        result.filtered_metadata[area_col], errors="coerce"
                    ).dropna().values
                    filtered_frame = pd.DataFrame({
                        "cell_area": filtered_area,
                        "sample_id": sample_id,
                        "donor": row.get("donor"),
                        "region": row.get("region"),
                    })
                    filtered_area_metadata_frames.append(filtered_frame)

            # Scatter plot data (sample for performance)
            if result.filtered_matrix is not None:
                # Compute total intensity from matrix
                marker_cols = loader.get_marker_columns(cell_matrix)
                total_intensity = cell_matrix[marker_cols].astype(float).sum(axis=1)
                area_series = pd.to_numeric(cell_metadata[area_col], errors="coerce")

                # Identify removed cells
                kept_ids = set(result.filtered_metadata[loader_config.cell_id_col])
                all_ids = set(cell_metadata[loader_config.cell_id_col])
                removed_ids = all_ids - kept_ids

                cell_metadata_indexed = cell_metadata.set_index(loader_config.cell_id_col)
                removed_mask = cell_metadata_indexed.index.isin(removed_ids)

                scatter_df = pd.DataFrame({
                    "cell_area": area_series.values,
                    "total_intensity": total_intensity.values,
                    "removed": removed_mask,
                }).dropna(subset=["cell_area", "total_intensity"])

                # Sample for scatter plot
                if len(scatter_df) > 5000:
                    scatter_df = scatter_df.sample(5000, random_state=rng.integers(0, 1_000_000))

                for _, scatter_row in scatter_df.iterrows():
                    scatter_records.append({
                        "cell_area": float(scatter_row["cell_area"]),
                        "total_intensity": float(scatter_row["total_intensity"]),
                        "removed": bool(scatter_row["removed"]),
                    })

            # Track thresholds
            area_p99 = np.percentile(area_data.get(sample_id, [0]), 99) if sample_id in area_data else np.nan
            if not np.isnan(area_p99):
                area_thresholds.append(area_p99)

            # Save filtered data
            if result.filtered_matrix is not None:
                result.filtered_matrix.to_csv(
                    filtered_matrix_dir / f"{sample_id}_filtered.csv", index=False
                )
            if result.filtered_metadata is not None:
                result.filtered_metadata.to_csv(
                    filtered_metadata_dir / f"{sample_id}_filtered_metadata.csv", index=False
                )

            # Record summary
            summary_row = result.to_dict()
            summary_records.append(summary_row)

            # Record removals
            for removal in result.removal_records:
                all_removals.append(removal)

            logger.info(
                f"  {result.cells_total:,} cells -> {result.cells_removed:,} removed "
                f"({result.removal_fraction*100:.1f}%)"
                + (" [CAPPED]" if result.capped_by_max else "")
            )

        except Exception as e:
            logger.error(f"  Error processing {sample_id}: {e}")
            summary_records.append({
                "sample_id": sample_id,
                "cells_total": 0,
                "cells_removed": 0,
                "removal_fraction": 0.0,
                "capped_by_max": False,
                "error": str(e),
            })

    # Create summary DataFrame
    summary_df = pd.DataFrame(summary_records)

    # Ensure column order
    column_order = [
        "sample_id", "cells_total", "cells_removed", "removal_fraction", "capped_by_max",
    ] + [f"removed_{r}" for r in REASON_COLUMNS]
    column_order = [c for c in column_order if c in summary_df.columns]
    summary_df = summary_df[[c for c in column_order if c in summary_df.columns]]

    # Write outputs
    summary_df.to_csv(output_dir / "qc_summary_per_sample.csv", index=False)
    logger.info(f"Wrote qc_summary_per_sample.csv")

    removals_df = pd.DataFrame(all_removals)
    removals_df.to_csv(output_dir / "qc_filtered_cells.csv", index=False)
    logger.info(f"Wrote qc_filtered_cells.csv ({len(removals_df)} cells)")

    # Aggregate statistics
    total_cells = int(summary_df["cells_total"].fillna(0).sum())
    total_removed = int(summary_df["cells_removed"].fillna(0).sum())
    removal_pct = (total_removed / total_cells * 100) if total_cells > 0 else 0.0
    logger.info(f"Total: {total_cells:,} cells, {total_removed:,} removed ({removal_pct:.2f}%)")

    # Generate figures
    if not skip_figures:
        logger.info("Generating figures...")

        area_metadata_df = pd.concat(area_metadata_frames, ignore_index=True) if area_metadata_frames else pd.DataFrame()
        filtered_area_df = pd.concat(filtered_area_metadata_frames, ignore_index=True) if filtered_area_metadata_frames else pd.DataFrame()
        scatter_df = pd.DataFrame(scatter_records) if scatter_records else pd.DataFrame()

        global_area_threshold = float(np.nanmedian(area_thresholds)) if area_thresholds else np.nan
        global_intensity_threshold = np.nan  # TODO: compute from data

        figures = generate_stage_b_figures(
            summary_df=summary_df,
            area_data=area_data,
            scatter_df=scatter_df,
            area_metadata_df=area_metadata_df,
            filtered_area_metadata_df=filtered_area_df,
            output_dir=output_dir,
            area_threshold=global_area_threshold,
            intensity_threshold=global_intensity_threshold,
            dpi=dpi,
            logger=logger,
        )
        logger.info(f"Generated {len(figures)} figures")

    # Summary
    n_with_removal = int((summary_df["cells_removed"] > 0).sum())
    logger.info(f"Stage B completed: {n_with_removal}/{len(summary_df)} samples had cells removed")

    return summary_df


# =============================================================================
# Stage C Visualization Functions
# =============================================================================


def build_cv_summary_plot(
    cv_df: pd.DataFrame,
    output_dir: Path,
    dpi: int = 200,
) -> Optional[Path]:
    """Generate coefficient of variation summary heatmap.

    Parameters
    ----------
    cv_df : pd.DataFrame
        CV data with sample_id, marker, cv columns
    output_dir : Path
        Output directory
    dpi : int
        Figure resolution

    Returns
    -------
    Optional[Path]
        Path to generated figure
    """
    if cv_df.empty:
        return None

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return None

    # Pivot to sample x marker matrix
    pivot = cv_df.pivot(index="sample_id", columns="marker", values="cv")

    fig, ax = plt.subplots(figsize=(14, max(6, len(pivot) * 0.3)))
    im = ax.imshow(pivot.values, aspect="auto", cmap="viridis")
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index, fontsize=7)
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=90, fontsize=7)
    ax.set_title("Coefficient of Variation by Sample and Marker")
    plt.colorbar(im, ax=ax, label="CV")

    fig.tight_layout()
    output_path = output_dir / "stage_c_cv_heatmap.png"
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return output_path


def build_transform_comparison_plot(
    stats_df: pd.DataFrame,
    output_dir: Path,
    sample_id: str,
    dpi: int = 200,
) -> Optional[Path]:
    """Generate transform comparison plot for a sample.

    Parameters
    ----------
    stats_df : pd.DataFrame
        Statistics data
    output_dir : Path
        Output directory
    sample_id : str
        Sample ID for title
    dpi : int
        Figure resolution

    Returns
    -------
    Optional[Path]
        Path to generated figure
    """
    if stats_df.empty:
        return None

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return None

    # Create comparison plot
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))

    # Zero fraction
    if "zero_frac" in stats_df.columns:
        ax = axes[0]
        values = stats_df.groupby("marker")["zero_frac"].mean()
        ax.bar(range(len(values)), values.values, alpha=0.7)
        ax.set_xticks(range(len(values)))
        ax.set_xticklabels(values.index, rotation=90, fontsize=6)
        ax.set_ylabel("Zero Fraction")
        ax.set_title("Zero Fraction by Marker")

    # Dynamic range
    if "dynamic_range" in stats_df.columns:
        ax = axes[1]
        values = stats_df.groupby("marker")["dynamic_range"].mean()
        ax.bar(range(len(values)), values.values, alpha=0.7, color="green")
        ax.set_xticks(range(len(values)))
        ax.set_xticklabels(values.index, rotation=90, fontsize=6)
        ax.set_ylabel("Dynamic Range")
        ax.set_title("Dynamic Range by Marker")

    # CV
    if "cv" in stats_df.columns:
        ax = axes[2]
        values = stats_df.groupby("marker")["cv"].mean()
        ax.bar(range(len(values)), values.values, alpha=0.7, color="orange")
        ax.set_xticks(range(len(values)))
        ax.set_xticklabels(values.index, rotation=90, fontsize=6)
        ax.set_ylabel("CV")
        ax.set_title("Coefficient of Variation by Marker")

    fig.suptitle(f"Normalization Statistics - {sample_id}")
    fig.tight_layout()
    output_path = output_dir / f"stage_c_stats_{sample_id}.png"
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return output_path


def build_normalization_summary_plot(
    summary_df: pd.DataFrame,
    output_dir: Path,
    dpi: int = 200,
) -> Optional[Path]:
    """Generate normalization summary bar chart.

    Parameters
    ----------
    summary_df : pd.DataFrame
        Summary data with sample_id, n_cells, n_markers
    output_dir : Path
        Output directory
    dpi : int
        Figure resolution

    Returns
    -------
    Optional[Path]
        Path to generated figure
    """
    if summary_df.empty:
        return None

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return None

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Cell counts
    ax = axes[0]
    ax.bar(range(len(summary_df)), summary_df["n_cells"], alpha=0.7)
    ax.set_xticks(range(len(summary_df)))
    ax.set_xticklabels(summary_df["sample_id"], rotation=90, fontsize=7)
    ax.set_ylabel("Cell Count")
    ax.set_title("Cells per Sample")

    # Mean zero fraction
    if "mean_zero_frac" in summary_df.columns:
        ax = axes[1]
        ax.bar(range(len(summary_df)), summary_df["mean_zero_frac"], alpha=0.7, color="orange")
        ax.set_xticks(range(len(summary_df)))
        ax.set_xticklabels(summary_df["sample_id"], rotation=90, fontsize=7)
        ax.set_ylabel("Mean Zero Fraction")
        ax.set_title("Mean Zero Fraction per Sample")

    fig.tight_layout()
    output_path = output_dir / "stage_c_normalization_summary.png"
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return output_path


def generate_stage_c_figures(
    summary_df: pd.DataFrame,
    cv_df: pd.DataFrame,
    output_dir: Path,
    dpi: int = 200,
    logger: Optional[logging.Logger] = None,
) -> List[Path]:
    """Generate all Stage C visualization figures.

    Parameters
    ----------
    summary_df : pd.DataFrame
        Normalization summary
    cv_df : pd.DataFrame
        CV data
    output_dir : Path
        Output directory
    dpi : int
        Figure resolution
    logger : Optional[Logger]
        Logger instance

    Returns
    -------
    List[Path]
        Paths to generated figures
    """
    logger = logger or logging.getLogger(__name__)
    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)

    generated = []

    # 1. CV heatmap
    if not cv_df.empty:
        path = build_cv_summary_plot(cv_df, figure_dir, dpi=dpi)
        if path:
            logger.info(f"Generated: {path.name}")
            generated.append(path)

    # 2. Normalization summary
    path = build_normalization_summary_plot(summary_df, figure_dir, dpi=dpi)
    if path:
        logger.info(f"Generated: {path.name}")
        generated.append(path)

    return generated


# =============================================================================
# Stage C Runner
# =============================================================================


def run_stage_c(
    input_dir: Path,
    output_dir: Path,
    config: Optional[NormalizationConfig] = None,
    loader_config: Optional[LoaderConfig] = None,
    variant: Optional[str] = None,
    verbose: bool = False,
    skip_figures: bool = False,
    dpi: int = 200,
    log_dir: Optional[Path] = None,
) -> pd.DataFrame:
    """Run Stage C: Within-sample normalization.

    Parameters
    ----------
    input_dir : Path
        Path to Stage B filtered matrices directory
    output_dir : Path
        Output directory
    config : Optional[NormalizationConfig]
        Normalization configuration
    loader_config : Optional[LoaderConfig]
        Loader configuration
    variant : Optional[str]
        Export variant (e.g., "clip-Q3__log1p")
    verbose : bool
        Enable verbose logging
    skip_figures : bool
        Skip figure generation
    dpi : int
        Figure resolution
    log_dir : Optional[Path]
        Directory for log file

    Returns
    -------
    pd.DataFrame
        Normalization summary per sample
    """
    logger = setup_logging(verbose, log_dir=log_dir, log_filename="stage_c.log")
    logger.info("Stage C: Within-Sample Normalization")
    logger.info(f"Input: {input_dir}")
    logger.info(f"Output: {output_dir}")

    # Setup directories
    output_dir.mkdir(parents=True, exist_ok=True)
    export_dir = output_dir / "export"
    export_dir.mkdir(exist_ok=True)
    params_dir = output_dir / "params"
    params_dir.mkdir(exist_ok=True)

    # Initialize config
    config = config or NormalizationConfig()
    loader_config = loader_config or LoaderConfig()

    # Parse variant
    if variant:
        parts = variant.split("__")
        if len(parts) == 2:
            bg_part, transform_part = parts
            # Parse bg_part: clip-Q3 -> mode=clip, p=3
            if "-" in bg_part:
                bg_mode, q_part = bg_part.split("-", 1)
                if q_part.startswith("Q"):
                    bg_percentile = int(q_part[1:])
                else:
                    bg_percentile = config.default_p
            else:
                bg_mode = config.default_bg_mode
                bg_percentile = config.default_p
            transform = transform_part
        else:
            bg_mode = config.default_bg_mode
            bg_percentile = config.default_p
            transform = config.default_transform
    else:
        bg_mode = config.default_bg_mode
        bg_percentile = config.default_p
        transform = config.default_transform

    variant_label = f"{bg_mode}-Q{bg_percentile}__{transform}"
    logger.info(f"Variant: {variant_label}")

    # Find input files
    matrix_files = sorted(input_dir.glob("*_filtered.csv"))
    if not matrix_files:
        # Try without _filtered suffix
        matrix_files = sorted(input_dir.glob("*.csv"))

    n_samples = len(matrix_files)
    logger.info(f"Found {n_samples} samples")

    if n_samples == 0:
        logger.error("No input matrices found")
        return pd.DataFrame()

    # Initialize normalizer
    normalizer = Normalizer(config)
    cell_id_col = loader_config.cell_id_col
    exclude_cols = set(loader_config.exclude_from_markers)

    # Tracking variables
    summary_records: List[Dict[str, Any]] = []
    cv_records: List[Dict[str, Any]] = []
    all_quantiles: List[pd.DataFrame] = []

    # Process each sample
    for idx, matrix_path in enumerate(matrix_files):
        # Extract sample ID
        sample_id = matrix_path.stem
        if sample_id.endswith("_filtered"):
            sample_id = sample_id[:-9]

        logger.info(f"[{idx+1}/{n_samples}] Processing {sample_id}")

        try:
            # Load matrix
            cell_matrix = pd.read_csv(matrix_path)
            if cell_id_col not in cell_matrix.columns:
                logger.error(f"  Missing cell_id column: {cell_id_col}")
                continue

            cell_matrix[cell_id_col] = cell_matrix[cell_id_col].astype(str)

            # Get marker columns
            marker_cols = [
                c for c in cell_matrix.columns
                if c != cell_id_col and c not in exclude_cols
            ]
            n_cells = len(cell_matrix)
            n_markers = len(marker_cols)

            # Normalize sample
            result = normalizer.normalize_sample(
                cell_matrix=cell_matrix,
                sample_id=sample_id,
                cell_id_col=cell_id_col,
                bg_mode=bg_mode,
                bg_percentile=bg_percentile,
                transform=transform,
            )

            # Save normalized matrix (only marker columns + cell_id)
            if result.normalized_matrix is not None:
                # Only include cell_id_col and marker columns
                output_cols = [cell_id_col] + marker_cols
                output_df = result.normalized_matrix[output_cols]
                output_path = export_dir / f"{sample_id}__{variant_label}.csv"
                output_df.to_csv(output_path, index=False)

            # Save quantiles
            if result.quantiles is not None:
                result.quantiles["sample_id"] = sample_id
                all_quantiles.append(result.quantiles.reset_index())

            # Compute statistics
            zero_fracs = []
            cvs = []
            for marker in marker_cols:
                if result.normalized_matrix is not None:
                    values = result.normalized_matrix[marker].values
                    zero_frac = normalizer.zero_fraction(values)
                    cv = normalizer.compute_cv(values)
                    zero_fracs.append(zero_frac)
                    cvs.append(cv)

                    cv_records.append({
                        "sample_id": sample_id,
                        "marker": marker,
                        "cv": cv,
                        "zero_frac": zero_frac,
                    })

            mean_zero_frac = float(np.nanmean(zero_fracs)) if zero_fracs else 0.0
            mean_cv = float(np.nanmean(cvs)) if cvs else 0.0

            summary_records.append({
                "sample_id": sample_id,
                "n_cells": n_cells,
                "n_markers": n_markers,
                "variant": variant_label,
                "mean_zero_frac": mean_zero_frac,
                "mean_cv": mean_cv,
            })

            logger.info(f"  {n_cells:,} cells, {n_markers} markers -> normalized")

        except Exception as e:
            logger.error(f"  Error processing {sample_id}: {e}")
            summary_records.append({
                "sample_id": sample_id,
                "n_cells": 0,
                "n_markers": 0,
                "variant": variant_label,
                "mean_zero_frac": float("nan"),
                "mean_cv": float("nan"),
                "error": str(e),
            })

    # Create summary DataFrame
    summary_df = pd.DataFrame(summary_records)
    summary_df.to_csv(output_dir / "normalization_summary.csv", index=False)
    logger.info(f"Wrote normalization_summary.csv")

    # Save CV data
    cv_df = pd.DataFrame(cv_records)
    cv_df.to_csv(output_dir / "cv_per_marker.csv", index=False)
    logger.info(f"Wrote cv_per_marker.csv")

    # Save combined quantiles
    if all_quantiles:
        quantiles_df = pd.concat(all_quantiles, ignore_index=True)
        quantiles_df.to_csv(params_dir / "quantiles.csv", index=False)
        logger.info(f"Wrote params/quantiles.csv")

    # Aggregate statistics
    total_cells = int(summary_df["n_cells"].fillna(0).sum())
    mean_zero = float(summary_df["mean_zero_frac"].mean()) if not summary_df.empty else 0.0
    logger.info(f"Total: {total_cells:,} cells normalized, mean zero fraction: {mean_zero:.3f}")

    # Generate figures
    if not skip_figures:
        logger.info("Generating figures...")
        figures = generate_stage_c_figures(
            summary_df=summary_df,
            cv_df=cv_df,
            output_dir=output_dir,
            dpi=dpi,
            logger=logger,
        )
        logger.info(f"Generated {len(figures)} figures")

    logger.info(f"Stage C completed: {len(summary_df)} samples normalized")

    return summary_df


# =============================================================================
# Stage D Visualization Functions
# =============================================================================


def build_alignment_histogram(
    before_values: np.ndarray,
    after_values: np.ndarray,
    output_dir: Path,
    intensity_label: str = "Transformed intensity",
    dpi: int = 200,
) -> Optional[Path]:
    """Generate histogram overlay of before/after alignment.

    Parameters
    ----------
    before_values : np.ndarray
        Intensity values before alignment
    after_values : np.ndarray
        Intensity values after alignment
    output_dir : Path
        Output directory
    intensity_label : str
        X-axis label
    dpi : int
        Figure resolution

    Returns
    -------
    Optional[Path]
        Path to generated figure
    """
    if before_values.size == 0 or after_values.size == 0:
        return None

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return None

    fig, ax = plt.subplots(figsize=(10, 6))

    ax.hist(before_values, bins=60, alpha=0.5, label="Before alignment", color="#3498db")
    ax.hist(after_values, bins=60, alpha=0.5, label="After alignment", color="#e74c3c")

    ax.set_xlabel(intensity_label)
    ax.set_ylabel("Frequency")
    ax.set_title("Distribution before vs after alignment")
    ax.legend()

    fig.tight_layout()
    output_path = output_dir / "stage_d_density_overlay.png"
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return output_path


def build_alignment_improvement_plot(
    metrics_df: pd.DataFrame,
    output_dir: Path,
    dpi: int = 200,
) -> Optional[Path]:
    """Generate bar chart of KS improvement per marker.

    Parameters
    ----------
    metrics_df : pd.DataFrame
        Alignment metrics
    output_dir : Path
        Output directory
    dpi : int
        Figure resolution

    Returns
    -------
    Optional[Path]
        Path to generated figure
    """
    if metrics_df.empty:
        return None

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return None

    metrics_df = metrics_df.copy()
    metrics_df["ks_improvement"] = metrics_df["ks_before"] - metrics_df["ks_after"]
    improvement = metrics_df.groupby("marker")["ks_improvement"].mean().sort_values(ascending=False)

    fig, ax = plt.subplots(figsize=(max(10, len(improvement) * 0.3), 5))

    colors = ["#2ecc71" if v > 0 else "#e74c3c" for v in improvement.values]
    ax.bar(range(len(improvement)), improvement.values, color=colors, alpha=0.8)

    ax.set_xticks(range(len(improvement)))
    ax.set_xticklabels(improvement.index, rotation=90, fontsize=7)
    ax.set_xlabel("Marker")
    ax.set_ylabel("ΔKS (improvement)")
    ax.set_title("KS improvement after alignment")
    ax.axhline(0, color="gray", linestyle="--", linewidth=0.5)

    fig.tight_layout()
    output_path = output_dir / "stage_d_alignment_improvement.png"
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return output_path


def build_slope_heatmap(
    coeff_df: pd.DataFrame,
    output_dir: Path,
    dpi: int = 200,
) -> Optional[Path]:
    """Generate heatmap of alignment slopes.

    Parameters
    ----------
    coeff_df : pd.DataFrame
        Alignment coefficients
    output_dir : Path
        Output directory
    dpi : int
        Figure resolution

    Returns
    -------
    Optional[Path]
        Path to generated figure
    """
    if coeff_df.empty:
        return None

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return None

    slope_pivot = coeff_df.pivot(index="sample_id", columns="marker", values="slope").fillna(1.0)

    n_samples, n_markers = slope_pivot.shape
    fig_width = max(12, n_markers * 0.3)
    fig_height = max(6, n_samples * 0.25)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    im = ax.imshow(slope_pivot.values, aspect="auto", cmap="coolwarm", vmin=0.5, vmax=1.5)

    ax.set_xticks(range(n_markers))
    ax.set_xticklabels(slope_pivot.columns, rotation=90, fontsize=7)
    ax.set_yticks(range(n_samples))
    ax.set_yticklabels(slope_pivot.index, fontsize=7)
    ax.set_title("Alignment Slope Heatmap")

    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("Slope")

    fig.tight_layout()
    output_path = output_dir / "stage_d_slope_heatmap.png"
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return output_path


def build_positive_fraction_plot(
    positive_df: pd.DataFrame,
    output_dir: Path,
    dpi: int = 200,
) -> Optional[Path]:
    """Generate bar chart of positive fraction consistency.

    Parameters
    ----------
    positive_df : pd.DataFrame
        Positive fraction data
    output_dir : Path
        Output directory
    dpi : int
        Figure resolution

    Returns
    -------
    Optional[Path]
        Path to generated figure
    """
    if positive_df.empty:
        return None

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return None

    consistency = positive_df.groupby("marker")["positive_fraction"].std().sort_values()

    fig, ax = plt.subplots(figsize=(max(10, len(consistency) * 0.3), 5))

    ax.bar(range(len(consistency)), consistency.values, color="#3498db", alpha=0.8)

    ax.set_xticks(range(len(consistency)))
    ax.set_xticklabels(consistency.index, rotation=90, fontsize=7)
    ax.set_xlabel("Marker")
    ax.set_ylabel("Std Dev")
    ax.set_title("Positive fraction consistency (stdev)")

    fig.tight_layout()
    output_path = output_dir / "stage_d_positive_fraction.png"
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return output_path


def generate_stage_d_figures(
    coeff_df: pd.DataFrame,
    metrics_df: pd.DataFrame,
    positive_df: pd.DataFrame,
    before_values: np.ndarray,
    after_values: np.ndarray,
    output_dir: Path,
    intensity_label: str = "Transformed intensity",
    dpi: int = 200,
    logger: Optional[logging.Logger] = None,
) -> List[Path]:
    """Generate all Stage D visualization figures.

    Parameters
    ----------
    coeff_df : pd.DataFrame
        Alignment coefficients
    metrics_df : pd.DataFrame
        Alignment metrics
    positive_df : pd.DataFrame
        Positive fraction data
    before_values : np.ndarray
        Values before alignment
    after_values : np.ndarray
        Values after alignment
    output_dir : Path
        Output directory
    intensity_label : str
        Intensity axis label
    dpi : int
        Figure resolution
    logger : Optional[Logger]
        Logger instance

    Returns
    -------
    List[Path]
        Paths to generated figures
    """
    logger = logger or logging.getLogger(__name__)
    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)

    generated = []

    # 1. Before/after histogram
    path = build_alignment_histogram(before_values, after_values, figure_dir, intensity_label, dpi)
    if path:
        logger.info(f"Generated: {path.name}")
        generated.append(path)

    # 2. KS improvement
    path = build_alignment_improvement_plot(metrics_df, figure_dir, dpi)
    if path:
        logger.info(f"Generated: {path.name}")
        generated.append(path)

    # 3. Slope heatmap
    path = build_slope_heatmap(coeff_df, figure_dir, dpi)
    if path:
        logger.info(f"Generated: {path.name}")
        generated.append(path)

    # 4. Positive fraction consistency
    path = build_positive_fraction_plot(positive_df, figure_dir, dpi)
    if path:
        logger.info(f"Generated: {path.name}")
        generated.append(path)

    return generated


# =============================================================================
# Stage D Runner
# =============================================================================


def run_stage_d(
    input_dir: Path,
    output_dir: Path,
    config: Optional[AlignmentConfig] = None,
    loader_config: Optional[LoaderConfig] = None,
    variant: Optional[str] = None,
    verbose: bool = False,
    skip_figures: bool = False,
    dpi: int = 200,
    log_dir: Optional[Path] = None,
) -> pd.DataFrame:
    """Run Stage D: Cross-sample alignment.

    Parameters
    ----------
    input_dir : Path
        Path to Stage C export directory
    output_dir : Path
        Output directory
    config : Optional[AlignmentConfig]
        Alignment configuration
    loader_config : Optional[LoaderConfig]
        Loader configuration
    variant : Optional[str]
        Normalization variant (e.g., "clip-Q3__log1p")
    verbose : bool
        Enable verbose logging
    skip_figures : bool
        Skip figure generation
    dpi : int
        Figure resolution
    log_dir : Optional[Path]
        Directory for log file

    Returns
    -------
    pd.DataFrame
        Alignment coefficients summary
    """
    logger = setup_logging(verbose, log_dir=log_dir, log_filename="stage_d.log")
    logger.info("Stage D: Cross-Sample Alignment")
    logger.info(f"Input: {input_dir}")
    logger.info(f"Output: {output_dir}")

    # Setup directories
    output_dir.mkdir(parents=True, exist_ok=True)
    aligned_dir = output_dir / "aligned"
    aligned_dir.mkdir(exist_ok=True)
    params_dir = output_dir / "params"
    params_dir.mkdir(exist_ok=True)
    metrics_dir = output_dir / "metrics"
    metrics_dir.mkdir(exist_ok=True)

    # Initialize config
    config = config or AlignmentConfig()
    loader_config = loader_config or LoaderConfig()
    cell_id_col = loader_config.cell_id_col
    exclude_cols = set(loader_config.exclude_from_markers)

    # Default variant
    variant = variant or "clip-Q3__log1p"
    logger.info(f"Variant: {variant}")

    # Find input files
    pattern = f"*__{variant}.csv"
    matrix_files = sorted(input_dir.glob(pattern))

    if not matrix_files:
        logger.error(f"No matrices found matching pattern: {pattern}")
        return pd.DataFrame()

    n_samples = len(matrix_files)
    logger.info(f"Found {n_samples} samples")

    # Initialize aligner
    aligner = CrossSampleAligner(config)
    rng = np.random.default_rng(config.random_seed)

    # Step 1: Load all samples and compute global targets
    logger.info("Loading samples and computing global targets...")
    sample_matrices: Dict[str, pd.DataFrame] = {}

    for matrix_path in matrix_files:
        # Extract sample ID
        sample_id = matrix_path.stem
        if "__" in sample_id:
            sample_id = sample_id.split("__", 1)[0]

        cell_matrix = pd.read_csv(matrix_path)
        if cell_id_col in cell_matrix.columns:
            cell_matrix[cell_id_col] = cell_matrix[cell_id_col].astype(str)
        sample_matrices[sample_id] = cell_matrix

    # Compute global targets
    global_targets = aligner.compute_global_targets(sample_matrices, cell_id_col)
    n_markers = len(global_targets)
    logger.info(f"Computed global targets for {n_markers} markers")

    # Step 2: Align each sample
    logger.info("Aligning samples...")
    coeff_records: List[Dict[str, Any]] = []
    metrics_records: List[Dict[str, Any]] = []
    positive_records: List[Dict[str, Any]] = []
    before_samples: List[np.ndarray] = []
    after_samples: List[np.ndarray] = []

    metric_sample_limit = config.metric_sample_limit
    fig_sample_limit = 4000

    for idx, (sample_id, cell_matrix) in enumerate(sample_matrices.items()):
        logger.info(f"[{idx+1}/{n_samples}] Aligning {sample_id}")

        # Get marker columns
        marker_cols = [
            c for c in cell_matrix.columns
            if c != cell_id_col and c not in exclude_cols
        ]

        # Align sample
        result = aligner.align_sample(
            cell_matrix=cell_matrix,
            sample_id=sample_id,
            global_targets=global_targets,
            cell_id_col=cell_id_col,
        )

        # Save aligned matrix (only marker columns + cell_id)
        if result.aligned_matrix is not None:
            output_cols = [cell_id_col] + marker_cols
            output_df = result.aligned_matrix[output_cols]
            output_df.to_csv(aligned_dir / f"{sample_id}_aligned.csv", index=False)

        # Record coefficients
        for param in result.params:
            coeff_records.append({
                "sample_id": sample_id,
                "marker": param.marker,
                "slope": param.slope,
                "intercept": param.intercept,
                "local_q5": param.local_lower,
                "local_q95": param.local_upper,
                "target_q5": param.target_lower,
                "target_q95": param.target_upper,
            })

        # Compute metrics for each marker
        intensity_df = cell_matrix[marker_cols].apply(pd.to_numeric, errors="coerce")
        aligned_df = result.aligned_matrix[marker_cols] if result.aligned_matrix is not None else intensity_df

        for marker in marker_cols:
            before_vals = intensity_df[marker].dropna().to_numpy(dtype=float)
            after_vals = aligned_df[marker].dropna().to_numpy(dtype=float)

            # Get global values for this marker (concatenate from all samples)
            global_vals = np.concatenate([
                sample_matrices[sid][marker].dropna().to_numpy(dtype=float)
                for sid in sample_matrices if marker in sample_matrices[sid].columns
            ])

            # Sample for metrics
            def sample_values(vals, limit):
                clean = vals[np.isfinite(vals)]
                if clean.size <= limit:
                    return clean
                indices = rng.choice(clean.size, size=limit, replace=False)
                return clean[indices]

            before_sampled = sample_values(before_vals, metric_sample_limit)
            after_sampled = sample_values(after_vals, metric_sample_limit)
            global_sampled = sample_values(global_vals, metric_sample_limit)

            # Compute distances
            ks_before = aligner.ks_distance(before_sampled, global_sampled)
            ks_after = aligner.ks_distance(after_sampled, global_sampled)
            emd_before = aligner.earth_movers_distance(before_sampled, global_sampled)
            emd_after = aligner.earth_movers_distance(after_sampled, global_sampled)

            metrics_records.append({
                "sample_id": sample_id,
                "marker": marker,
                "ks_before": ks_before,
                "ks_after": ks_after,
                "emd_before": emd_before,
                "emd_after": emd_after,
            })

            # Collect samples for histogram
            before_samples.append(sample_values(before_vals, fig_sample_limit))
            after_samples.append(sample_values(after_vals, fig_sample_limit))

            # Positive fraction
            positive_records.append({
                "sample_id": sample_id,
                "marker": marker,
                "positive_fraction": float(np.mean(after_vals > 0)),
            })

        n_cells = len(cell_matrix)
        logger.info(f"  {n_cells:,} cells aligned")

    # Create output DataFrames
    coeff_df = pd.DataFrame(coeff_records)
    metrics_df = pd.DataFrame(metrics_records)
    positive_df = pd.DataFrame(positive_records)

    # Save outputs
    coeff_df.to_csv(params_dir / "alignment_coefficients.csv", index=False)
    logger.info("Wrote params/alignment_coefficients.csv")

    metrics_df.to_csv(metrics_dir / "alignment_quality.csv", index=False)
    logger.info("Wrote metrics/alignment_quality.csv")

    # Aggregate statistics
    total_cells = sum(len(df) for df in sample_matrices.values())
    mean_ks_improvement = float((metrics_df["ks_before"] - metrics_df["ks_after"]).mean()) if not metrics_df.empty else 0.0
    logger.info(f"Total: {total_cells:,} cells aligned, mean KS improvement: {mean_ks_improvement:.4f}")

    # Generate figures
    if not skip_figures:
        logger.info("Generating figures...")

        # Concatenate before/after values
        before_concat = np.concatenate(before_samples) if before_samples else np.array([])
        after_concat = np.concatenate(after_samples) if after_samples else np.array([])

        # Derive intensity label
        if "__" in variant:
            transform = variant.split("__", 1)[1].strip()
            intensity_label = f"{transform.replace('_', ' ')} intensity"
        else:
            intensity_label = "Transformed intensity"

        figures = generate_stage_d_figures(
            coeff_df=coeff_df,
            metrics_df=metrics_df,
            positive_df=positive_df,
            before_values=before_concat,
            after_values=after_concat,
            output_dir=output_dir,
            intensity_label=intensity_label,
            dpi=dpi,
            logger=logger,
        )
        logger.info(f"Generated {len(figures)} figures")

    logger.info(f"Stage D completed: {n_samples} samples aligned")

    return coeff_df


# =============================================================================
# Stage E Helper Classes and Functions
# =============================================================================


class ChannelMetadataCache:
    """Cache for channel metadata lookup.

    Resolves imaging_color and imaging_cycle for each marker/donor.
    """

    def __init__(self, directory: Path) -> None:
        self.directory = directory
        self._cache: Dict[str, Dict[str, Tuple[Optional[str], Optional[str]]]] = {}
        self._missing: Set[Tuple[str, str]] = set()

    @staticmethod
    def _canonicalize_marker(name: str) -> str:
        """Normalize marker name for matching."""
        return re.sub(r"[^0-9A-Za-z]+", "", str(name).lower())

    def _resolve_path(self, donor: str) -> Path:
        """Find channel metadata file for donor."""
        pattern = f"{donor}_*.csv"
        matches = sorted(self.directory.glob(pattern))
        if not matches:
            raise FileNotFoundError(
                f"Channel metadata not found for donor {donor}: expected {pattern}"
            )
        return matches[0]

    def _load_map(self, donor: str) -> Dict[str, Tuple[Optional[str], Optional[str]]]:
        """Load channel map for a donor."""
        path = self._resolve_path(donor)
        df = pd.read_csv(path)
        if "antibody_name" not in df.columns:
            raise ValueError(f"Channel metadata for {donor} missing antibody_name column")

        df["__canon"] = df["antibody_name"].map(self._canonicalize_marker)
        mapping: Dict[str, Tuple[Optional[str], Optional[str]]] = {}

        for _, row in df.iterrows():
            canon = row["__canon"]
            color = row.get("Imaging Color") if "Imaging Color" in df.columns else None
            cycle = str(row.get("Imaging Cycle")) if "Imaging Cycle" in df.columns else None
            mapping[canon] = (color, cycle)

        return mapping

    def get(self, donor: str) -> Dict[str, Tuple[Optional[str], Optional[str]]]:
        """Get channel map for donor (cached)."""
        donor_key = str(donor)
        if donor_key not in self._cache:
            self._cache[donor_key] = self._load_map(donor_key)
        return self._cache[donor_key]

    def lookup(
        self,
        donor: str,
        marker: str,
        logger: Optional[logging.Logger] = None,
    ) -> Tuple[Optional[str], Optional[str]]:
        """Look up imaging_color and imaging_cycle for a marker."""
        mapping = self.get(donor)
        canon = self._canonicalize_marker(marker)
        if canon in mapping:
            return mapping[canon]
        key = (donor, canon)
        if key not in self._missing:
            if logger:
                logger.warning("Channel metadata missing for donor=%s marker=%s", donor, marker)
            self._missing.add(key)
        return None, None


# =============================================================================
# Stage E Visualization Functions
# =============================================================================


def build_pca_embedding(
    matrix: pd.DataFrame,
    random_seed: int = 2026,
) -> pd.DataFrame:
    """Compute PCA embedding from intensity matrix.

    Parameters
    ----------
    matrix : pd.DataFrame
        Intensity matrix (cells x markers)
    random_seed : int
        Random seed

    Returns
    -------
    pd.DataFrame
        PCA coordinates with PC1, PC2
    """
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    if matrix.empty:
        return pd.DataFrame()

    scaler = StandardScaler()
    scaled = scaler.fit_transform(matrix.to_numpy(dtype=float))
    n_components = min(2, matrix.shape[1])
    pca = PCA(n_components=n_components, random_state=random_seed)
    coords = pca.fit_transform(scaled)
    columns = ["PC1", "PC2"][:n_components]
    return pd.DataFrame(coords, index=matrix.index, columns=columns)


def build_pca_plot(
    coords: pd.DataFrame,
    metadata: pd.DataFrame,
    output_path: Path,
    label_columns: List[str],
    title_prefix: str,
    dpi: int = 200,
) -> Optional[Path]:
    """Generate PCA embedding plot colored by metadata columns.

    Parameters
    ----------
    coords : pd.DataFrame
        PCA coordinates
    metadata : pd.DataFrame
        Cell metadata
    output_path : Path
        Output path
    label_columns : List[str]
        Columns to use for coloring
    title_prefix : str
        Title prefix
    dpi : int
        Figure resolution

    Returns
    -------
    Optional[Path]
        Path to generated figure
    """
    if coords.empty:
        return None

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return None

    import math

    n_plots = len(label_columns)
    n_cols = min(3, n_plots)
    n_rows = math.ceil(n_plots / n_cols)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 4 * n_rows))
    axes = np.atleast_1d(axes).flatten()

    joined = coords.join(metadata, how="left")

    for idx, (axis, label) in enumerate(zip(axes, label_columns)):
        if label not in joined.columns:
            continue
        categories = joined[label].fillna("NA").astype(str)
        for category in sorted(categories.unique()):
            mask = categories == category
            axis.scatter(
                joined.loc[mask, coords.columns[0]],
                joined.loc[mask, coords.columns[1] if coords.shape[1] > 1 else coords.columns[0]],
                s=12,
                alpha=0.7,
                label=category,
            )
        axis.set_xlabel(coords.columns[0])
        axis.set_ylabel(coords.columns[1] if coords.shape[1] > 1 else coords.columns[0])
        axis.set_title(f"{title_prefix} by {label}")
        if idx > 0:  # Hide legend for first plot (too many samples)
            axis.legend(fontsize=7, markerscale=0.7)

    for axis in axes[n_plots:]:
        axis.axis("off")

    fig.tight_layout()
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return output_path


def build_drift_scatter_plot(
    stats_df: pd.DataFrame,
    output_path: Path,
    title: str = "Donor-associated intensity variation per marker",
    dpi: int = 200,
) -> Optional[Path]:
    """Generate drift scatter plot (eta2 vs -log10(p)).

    Parameters
    ----------
    stats_df : pd.DataFrame
        Marker statistics with eta2_batch, p_batch, stability
    output_path : Path
        Output path
    title : str
        Plot title
    dpi : int
        Figure resolution

    Returns
    -------
    Optional[Path]
        Path to generated figure
    """
    if stats_df.empty:
        return None

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return None

    df = stats_df.copy()
    df["neglog10_p"] = -np.log10(df["p_batch"].replace(0, 1e-300))

    fig, ax = plt.subplots(figsize=(8, 6))

    stability_colors = {
        "drifted": "#E24A33",
        "observe": "#FBC15E",
        "stable": "#348ABD",
    }

    for stability, subset in df.groupby("stability"):
        color = stability_colors.get(stability, "#888888")
        ax.scatter(
            subset["eta2_batch"],
            subset["neglog10_p"],
            label=stability,
            color=color,
            alpha=0.7,
        )

    ax.axvline(0.5, linestyle="--", linewidth=1, color="gray")
    ax.axhline(-np.log10(0.01), linestyle="--", linewidth=1, color="gray")
    ax.set_xlabel("η² (batch)")
    ax.set_ylabel("-log10(p)")
    ax.set_title(title)
    ax.legend()

    fig.tight_layout()
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return output_path


def build_cv_bar_plot(
    stats_df: pd.DataFrame,
    output_path: Path,
    title: str = "Coefficient of variation per marker",
    dpi: int = 200,
) -> Optional[Path]:
    """Generate CV bar plot.

    Parameters
    ----------
    stats_df : pd.DataFrame
        Marker statistics with marker and cv columns
    output_path : Path
        Output path
    title : str
        Plot title
    dpi : int
        Figure resolution

    Returns
    -------
    Optional[Path]
        Path to generated figure
    """
    if stats_df.empty:
        return None

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return None

    n_markers = len(stats_df)
    fig, ax = plt.subplots(figsize=(max(10, n_markers * 0.3), 5))

    ax.bar(range(n_markers), stats_df["cv"].fillna(0).values, color="#3498db", alpha=0.8)

    ax.set_xticks(range(n_markers))
    ax.set_xticklabels(stats_df["marker"].values, rotation=90, fontsize=7)
    ax.set_xlabel("Marker")
    ax.set_ylabel("CV")
    ax.set_title(title)

    fig.tight_layout()
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return output_path


def build_delta_cv_plot(
    baseline_stats: pd.DataFrame,
    corrected_stats: pd.DataFrame,
    output_path: Path,
    dpi: int = 200,
) -> Optional[Path]:
    """Generate delta CV bar plot (before - after correction).

    Parameters
    ----------
    baseline_stats : pd.DataFrame
        Baseline marker statistics
    corrected_stats : pd.DataFrame
        Corrected marker statistics
    output_path : Path
        Output path
    dpi : int
        Figure resolution

    Returns
    -------
    Optional[Path]
        Path to generated figure
    """
    if baseline_stats.empty or corrected_stats.empty:
        return None

    try:
        import matplotlib.pyplot as plt
        from matplotlib.patches import Patch
    except ImportError:
        return None

    # Merge
    base = baseline_stats[["marker", "cv", "stability"]].rename(
        columns={"cv": "cv_before", "stability": "stability_before"}
    )
    corr = corrected_stats[["marker", "cv"]].rename(columns={"cv": "cv_after"})
    merged = base.merge(corr, on="marker", how="inner")

    if merged.empty:
        return None

    merged["delta_cv"] = merged["cv_before"] - merged["cv_after"]
    merged = merged.sort_values("delta_cv", ascending=False)

    n_markers = len(merged)
    fig, ax = plt.subplots(figsize=(max(10, n_markers * 0.3), 5))

    stability_colors = {
        "drifted": "#E24A33",
        "observe": "#FBC15E",
        "stable": "#348ABD",
    }

    colors = [stability_colors.get(s, "#888888") for s in merged["stability_before"]]
    ax.bar(range(n_markers), merged["delta_cv"].values, color=colors, alpha=0.8)
    ax.axhline(0, color="black", linestyle="--", linewidth=1)

    ax.set_xticks(range(n_markers))
    ax.set_xticklabels(merged["marker"].values, rotation=90, fontsize=7)
    ax.set_xlabel("Marker")
    ax.set_ylabel("ΔCV = CV_before - CV_after")
    ax.set_title("Change in CV per marker after batch correction")

    legend_elements = [
        Patch(facecolor=stability_colors["drifted"], label="drifted"),
        Patch(facecolor=stability_colors["observe"], label="observe"),
        Patch(facecolor=stability_colors["stable"], label="stable"),
    ]
    ax.legend(handles=legend_elements, title="Baseline stability", loc="upper right")

    fig.tight_layout()
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return output_path


def generate_stage_e_figures(
    baseline_stats: pd.DataFrame,
    corrected_stats: Optional[pd.DataFrame],
    pca_coords: pd.DataFrame,
    pca_corrected_coords: Optional[pd.DataFrame],
    embedding_metadata: pd.DataFrame,
    output_dir: Path,
    dpi: int = 200,
    logger: Optional[logging.Logger] = None,
) -> List[Path]:
    """Generate all Stage E visualization figures.

    Parameters
    ----------
    baseline_stats : pd.DataFrame
        Baseline marker statistics
    corrected_stats : Optional[pd.DataFrame]
        Corrected marker statistics
    pca_coords : pd.DataFrame
        PCA coordinates before correction
    pca_corrected_coords : Optional[pd.DataFrame]
        PCA coordinates after correction
    embedding_metadata : pd.DataFrame
        Cell metadata for embedding
    output_dir : Path
        Output directory
    dpi : int
        Figure resolution
    logger : Optional[Logger]
        Logger instance

    Returns
    -------
    List[Path]
        Paths to generated figures
    """
    logger = logger or logging.getLogger(__name__)
    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)

    generated = []
    label_columns = ["sample_id", "donor", "region"]

    # 1. PCA before correction
    if not pca_coords.empty and not embedding_metadata.empty:
        path = build_pca_plot(
            pca_coords, embedding_metadata, figure_dir / "stage_e_pca_batch.png",
            label_columns, "Aligned PCA", dpi
        )
        if path:
            logger.info(f"Generated: {path.name}")
            generated.append(path)

    # 2. PCA after correction
    if pca_corrected_coords is not None and not pca_corrected_coords.empty:
        path = build_pca_plot(
            pca_corrected_coords, embedding_metadata, figure_dir / "stage_e_pca_corrected.png",
            label_columns, "Corrected PCA", dpi
        )
        if path:
            logger.info(f"Generated: {path.name}")
            generated.append(path)

    # 3. Drift scatter (baseline)
    path = build_drift_scatter_plot(baseline_stats, figure_dir / "stage_e_drift_scatter.png", dpi=dpi)
    if path:
        logger.info(f"Generated: {path.name}")
        generated.append(path)

    # 4. CV bar (baseline)
    path = build_cv_bar_plot(baseline_stats, figure_dir / "stage_e_cv_baseline.png", dpi=dpi)
    if path:
        logger.info(f"Generated: {path.name}")
        generated.append(path)

    # 5. Delta CV
    if corrected_stats is not None:
        path = build_delta_cv_plot(baseline_stats, corrected_stats, figure_dir / "stage_e_delta_cv.png", dpi)
        if path:
            logger.info(f"Generated: {path.name}")
            generated.append(path)

    return generated


# =============================================================================
# Stage E Runner
# =============================================================================


def run_stage_e(
    input_dir: Path,
    output_dir: Path,
    metadata_path: Path,
    channels_dir: Optional[Path] = None,
    config: Optional[BatchCorrectionConfig] = None,
    loader_config: Optional[LoaderConfig] = None,
    correction_mode: str = "residual",
    verbose: bool = False,
    skip_figures: bool = False,
    dpi: int = 200,
    log_dir: Optional[Path] = None,
) -> pd.DataFrame:
    """Run Stage E: Batch effect correction.

    Parameters
    ----------
    input_dir : Path
        Path to Stage D aligned directory
    output_dir : Path
        Output directory
    metadata_path : Path
        Path to Stage A metadata_verified.csv
    channels_dir : Optional[Path]
        Path to channel metadata directory
    config : Optional[BatchCorrectionConfig]
        Batch correction configuration
    loader_config : Optional[LoaderConfig]
        Loader configuration
    correction_mode : str
        Correction mode: "residual" or "none"
    verbose : bool
        Enable verbose logging
    skip_figures : bool
        Skip figure generation
    dpi : int
        Figure resolution
    log_dir : Optional[Path]
        Directory for log file

    Returns
    -------
    pd.DataFrame
        Marker statistics summary
    """
    logger = setup_logging(verbose, log_dir=log_dir, log_filename="stage_e.log")
    logger.info("Stage E: Batch Effect Correction")
    logger.info(f"Input: {input_dir}")
    logger.info(f"Output: {output_dir}")
    logger.info(f"Metadata: {metadata_path}")
    logger.info(f"Correction mode: {correction_mode}")

    # Setup directories
    output_dir.mkdir(parents=True, exist_ok=True)
    corrected_dir = output_dir / "corrected"
    corrected_dir.mkdir(exist_ok=True)

    # Initialize config
    config = config or BatchCorrectionConfig()
    config.correction_mode = correction_mode
    loader_config = loader_config or LoaderConfig()
    cell_id_col = loader_config.cell_id_col
    exclude_cols = set(loader_config.exclude_from_markers)

    # Load metadata
    metadata_df = pd.read_csv(metadata_path)
    if "sample_id" not in metadata_df.columns:
        raise ValueError("Metadata must have sample_id column")
    metadata_df["sample_id"] = metadata_df["sample_id"].astype(str)

    # Find aligned files
    aligned_files = sorted(input_dir.glob("*_aligned.csv"))
    if not aligned_files:
        logger.error(f"No aligned files found in {input_dir}")
        return pd.DataFrame()

    n_samples = len(aligned_files)
    logger.info(f"Found {n_samples} aligned samples")

    # Initialize channel cache if available
    channel_cache = None
    if channels_dir and channels_dir.exists():
        channel_cache = ChannelMetadataCache(channels_dir)
        logger.info(f"Using channel metadata from: {channels_dir}")

    # Initialize corrector
    if channels_dir:
        config.additional_batch_factors = ["imaging_color", "imaging_cycle"]
    corrector = BatchCorrector(config)
    rng = np.random.default_rng(config.random_seed)

    # Step 1: Load all samples and build summary
    logger.info("Loading samples and building summary...")
    sample_matrices: Dict[str, pd.DataFrame] = {}
    sample_metadata: Dict[str, Dict[str, Any]] = {}
    sample_counts: Dict[str, int] = {}
    summary_records: List[Dict[str, Any]] = []

    # For embeddings
    embedding_cells_per_sample = config.embedding_cells_per_sample
    embedding_rows: List[pd.DataFrame] = []
    embedding_meta_rows: List[pd.DataFrame] = []

    for aligned_path in aligned_files:
        sample_id = aligned_path.stem.replace("_aligned", "")

        # Load matrix
        cell_matrix = pd.read_csv(aligned_path)
        if cell_id_col in cell_matrix.columns:
            cell_matrix[cell_id_col] = cell_matrix[cell_id_col].astype(str)
        sample_matrices[sample_id] = cell_matrix
        sample_counts[sample_id] = len(cell_matrix)

        # Get metadata
        meta_row = metadata_df[metadata_df["sample_id"] == sample_id]
        if meta_row.empty:
            logger.warning(f"No metadata for sample {sample_id}")
            sample_metadata[sample_id] = {}
        else:
            sample_metadata[sample_id] = meta_row.iloc[0].to_dict()

        donor = sample_metadata[sample_id].get("donor")
        region = sample_metadata[sample_id].get("region")

        # Get marker columns
        marker_cols = [
            c for c in cell_matrix.columns
            if c != cell_id_col and c not in exclude_cols
        ]
        numeric = cell_matrix[marker_cols].apply(pd.to_numeric, errors="coerce")
        means = numeric.mean(axis=0)
        counts = numeric.notna().sum(axis=0)

        # Build summary records
        for marker in marker_cols:
            record = {
                "sample_id": sample_id,
                "marker": marker,
                "intensity": float(means.get(marker, np.nan)),
                "n_cells": int(counts.get(marker, 0)),
                config.batch_col: donor,
                config.biovar_col: region,
            }

            # Add channel metadata if available
            if channel_cache and donor:
                try:
                    color, cycle = channel_cache.lookup(donor, marker, logger)
                    record["imaging_color"] = color if pd.notna(color) else "missing_color"
                    record["imaging_cycle"] = str(cycle) if pd.notna(cycle) else "missing_cycle"
                except FileNotFoundError:
                    record["imaging_color"] = "missing_color"
                    record["imaging_cycle"] = "missing_cycle"

            summary_records.append(record)

        # Sample cells for embedding
        n_cells = len(cell_matrix)
        if embedding_cells_per_sample > 0 and n_cells > 0:
            take = min(embedding_cells_per_sample, n_cells)
            indices = rng.choice(n_cells, size=take, replace=False)
            marker_values = numeric.iloc[indices].reset_index(drop=True)
            cell_ids = cell_matrix[cell_id_col].iloc[indices].tolist()
            global_ids = [f"{sample_id}__{cid}" for cid in cell_ids]
            marker_values.index = global_ids

            meta_rows = pd.DataFrame({
                "cell_id": global_ids,
                "sample_id": sample_id,
                "donor": donor,
                "region": region,
            }).set_index("cell_id")

            embedding_rows.append(marker_values)
            embedding_meta_rows.append(meta_rows)

    summary_df = pd.DataFrame(summary_records)

    # Build embedding matrix
    embedding_matrix = pd.concat(embedding_rows, axis=0) if embedding_rows else pd.DataFrame()
    embedding_metadata = pd.concat(embedding_meta_rows, axis=0) if embedding_meta_rows else pd.DataFrame()

    # Step 2: Compute batch shifts using linear model
    logger.info("Computing batch correction shifts...")
    shift_lookup, diagnostics = corrector.compute_batch_shifts(summary_df)

    # Save model summary
    if diagnostics.model_summary:
        model_summary_path = output_dir / "linear_model_summary.txt"
        model_summary_path.write_text(diagnostics.model_summary, encoding="utf-8")
        logger.info(f"Wrote {model_summary_path}")

    # Step 3: Apply corrections
    corrected_summary_records: List[Dict[str, Any]] = []

    if correction_mode == "residual" and shift_lookup:
        logger.info("Applying batch corrections...")

        for sample_id, cell_matrix in sample_matrices.items():
            result = corrector.correct_sample(
                cell_matrix=cell_matrix,
                sample_id=sample_id,
                shift_lookup=shift_lookup,
                cell_id_col=cell_id_col,
            )

            # Save corrected matrix
            if result.corrected_matrix is not None:
                marker_cols = [c for c in result.corrected_matrix.columns if c != cell_id_col]
                output_df = result.corrected_matrix[[cell_id_col] + marker_cols]
                output_df.to_csv(corrected_dir / f"{sample_id}_corrected.csv", index=False)

                # Build corrected summary
                numeric = result.corrected_matrix[marker_cols].apply(pd.to_numeric, errors="coerce")
                means = numeric.mean(axis=0)
                counts = numeric.notna().sum(axis=0)

                for marker in marker_cols:
                    corrected_summary_records.append({
                        "sample_id": sample_id,
                        "marker": marker,
                        "intensity": float(means.get(marker, np.nan)),
                        "n_cells": int(counts.get(marker, 0)),
                        config.batch_col: sample_metadata[sample_id].get(config.batch_col),
                        config.biovar_col: sample_metadata[sample_id].get(config.biovar_col),
                    })

        logger.info(f"Wrote {len(sample_matrices)} corrected matrices")
    else:
        corrected_summary_records = [
            {
                "sample_id": r["sample_id"],
                "marker": r["marker"],
                "intensity": r["intensity"],
                "n_cells": r["n_cells"],
                config.batch_col: r.get(config.batch_col),
                config.biovar_col: r.get(config.biovar_col),
            }
            for r in summary_records
        ]

    corrected_summary_df = pd.DataFrame(corrected_summary_records)

    # Step 4: Compute marker statistics
    logger.info("Computing marker statistics...")
    baseline_stats = corrector.compute_marker_statistics(summary_df, sample_counts)
    baseline_stats.to_csv(output_dir / "batch_drift_summary.csv", index=False)
    logger.info("Wrote batch_drift_summary.csv")

    # Marker stability table
    stability_cols = ["marker", "eta2_batch", "p_batch", "stability"]
    if "eta2_biovar" in baseline_stats.columns:
        stability_cols.insert(3, "eta2_biovar")
        stability_cols.insert(4, "p_biovar")
    baseline_stats[stability_cols].to_csv(output_dir / "marker_stability_table.csv", index=False)
    logger.info("Wrote marker_stability_table.csv")

    corrected_stats = None
    if correction_mode == "residual":
        corrected_stats = corrector.compute_marker_statistics(corrected_summary_df, sample_counts)
        corrected_stats.to_csv(corrected_dir / "batch_drift_summary.csv", index=False)
        logger.info("Wrote corrected/batch_drift_summary.csv")

    # Step 5: Compute embeddings
    pca_coords = pd.DataFrame()
    pca_corrected_coords = None

    if not embedding_matrix.empty:
        logger.info("Computing PCA embeddings...")
        pca_coords = build_pca_embedding(embedding_matrix, config.random_seed)

        if correction_mode == "residual" and shift_lookup:
            # Apply shifts to embedding
            corrected_emb = embedding_matrix.copy()
            sample_ids = embedding_metadata["sample_id"]
            for marker in corrected_emb.columns:
                shifts = sample_ids.map(lambda sid: shift_lookup.get((sid, marker), 0.0)).to_numpy()
                corrected_emb[marker] = corrected_emb[marker].to_numpy(dtype=float) - shifts
            pca_corrected_coords = build_pca_embedding(corrected_emb, config.random_seed)

    # Summary
    total_cells = sum(sample_counts.values())
    n_drifted = (baseline_stats["stability"] == "drifted").sum() if "stability" in baseline_stats.columns else 0
    logger.info(f"Total: {total_cells:,} cells, {len(baseline_stats)} markers, {n_drifted} drifted")

    # Step 6: Generate figures
    if not skip_figures:
        logger.info("Generating figures...")
        figures = generate_stage_e_figures(
            baseline_stats=baseline_stats,
            corrected_stats=corrected_stats,
            pca_coords=pca_coords,
            pca_corrected_coords=pca_corrected_coords,
            embedding_metadata=embedding_metadata,
            output_dir=output_dir,
            dpi=dpi,
            logger=logger,
        )
        logger.info(f"Generated {len(figures)} figures")

    logger.info(f"Stage E completed: {n_samples} samples processed")

    return baseline_stats


# =============================================================================
# Stage F: Data Merge and Spatial Graph Construction
# =============================================================================
# TODO(FEAT-001): Add optional --marker-map validation after column sanitization
# Column names with "/" are sanitized to "_" (e.g., "Keratin 8/18" -> "Keratin 8_18")
# This can cause silent marker resolution failures in downstream annotation.
# See docs/FEATURE_REQUESTS.md for implementation details.


@dataclass
class MergeConfig:
    """Configuration for data merge stage.

    Attributes
    ----------
    cell_id_col : str
        Column name for cell IDs
    graph_mode : str
        Graph type: 'knn' or 'ball'
    k : int
        Number of neighbors for kNN graphs
    radius : float
        Radius for ball graphs (in spatial units)
    x_col : str
        Column name for x coordinate
    y_col : str
        Column name for y coordinate
    """

    cell_id_col: str = "cell_mask_id"
    graph_mode: str = "knn"
    k: int = 10
    radius: float = 20.0
    x_col: str = "x"
    y_col: str = "y"


def _sanitize_column_name(name: str) -> str:
    """Sanitize column name for AnnData compatibility."""
    sanitize_map = {"/": "_", "\\": "_"}
    for bad, replacement in sanitize_map.items():
        name = name.replace(bad, replacement)
    return name


def _sanitize_dataframe_columns(
    df: pd.DataFrame,
    exclude: Optional[set] = None,
    logger: Optional[logging.Logger] = None,
) -> pd.DataFrame:
    """Sanitize DataFrame column names for AnnData compatibility.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame
    exclude : set, optional
        Columns to exclude from sanitization
    logger : logging.Logger, optional
        Logger for reporting changes

    Returns
    -------
    pd.DataFrame
        DataFrame with sanitized column names
    """
    exclude = exclude or set()
    original_cols = list(df.columns)
    new_cols = []
    used = set()

    for col in original_cols:
        if col in exclude:
            new_col = col
        else:
            new_col = _sanitize_column_name(col)

        # Handle duplicates
        base = new_col
        count = 0
        while new_col in used:
            count += 1
            new_col = f"{base}__{count}"

        used.add(new_col)
        new_cols.append(new_col)

    df.columns = new_cols

    if logger:
        changes = {orig: new for orig, new in zip(original_cols, new_cols) if orig != new}
        if changes:
            logger.debug(f"Sanitized columns: {changes}")

    return df


def build_spatial_graph(
    coords: np.ndarray,
    mode: str = "knn",
    k: int = 10,
    radius: float = 20.0,
) -> sp.csr_matrix:
    """Build spatial neighborhood graph.

    Parameters
    ----------
    coords : np.ndarray
        Spatial coordinates (n_cells x 2)
    mode : str
        Graph type: 'knn' or 'ball'
    k : int
        Number of neighbors for kNN (including self)
    radius : float
        Radius for ball graphs

    Returns
    -------
    scipy.sparse.csr_matrix
        Sparse adjacency matrix
    """
    from sklearn.neighbors import NearestNeighbors, radius_neighbors_graph

    n_cells = coords.shape[0]
    if n_cells == 0:
        return sp.csr_matrix((0, 0))

    if mode == "knn":
        max_neighbors = n_cells - 1
        if max_neighbors < 1:
            return sp.csr_matrix((n_cells, n_cells))
        n_neighbors = min(max(k, 1), max_neighbors)
        nn = NearestNeighbors(n_neighbors=n_neighbors, metric="euclidean")
        nn.fit(coords)
        graph = nn.kneighbors_graph(mode="connectivity")
        # Remove self-loops
        graph = graph - sp.diags(graph.diagonal())
        return graph.tocsr()
    else:
        # Ball (radius) mode
        radius = max(radius, 1e-6)
        graph = radius_neighbors_graph(
            coords, radius=radius, mode="connectivity", include_self=False
        )
        return graph.tocsr()


def build_anndata(
    merged_tables: List[pd.DataFrame],
    raw_tables: List[Optional[pd.DataFrame]],
    corrected_tables: List[Optional[pd.DataFrame]],
    marker_cols: List[str],
    config: MergeConfig,
    logger: Optional[logging.Logger] = None,
) -> "anndata.AnnData":
    """Build AnnData object from merged sample data.

    Parameters
    ----------
    merged_tables : List[pd.DataFrame]
        Aligned data merged with cell metadata
    raw_tables : List[pd.DataFrame]
        Stage C normalized (not aligned) data
    corrected_tables : List[pd.DataFrame]
        Stage E batch-corrected data
    marker_cols : List[str]
        Marker column names
    config : MergeConfig
        Merge configuration
    logger : logging.Logger, optional
        Logger for progress messages

    Returns
    -------
    anndata.AnnData
        Merged AnnData object with layers
    """
    try:
        import anndata as ad
    except ImportError as exc:
        raise RuntimeError("anndata is required for Stage F. Install with: pip install anndata") from exc

    def _log(msg: str) -> None:
        if logger:
            logger.info(msg)

    # Concatenate all samples
    _log(f"  [1/4] Concatenating {len(merged_tables)} aligned tables...")
    merged_df = pd.concat(merged_tables, ignore_index=True)
    expression = merged_df[marker_cols].to_numpy(dtype=np.float32)
    _log(f"  [1/4] Done: {expression.shape[0]:,} cells x {expression.shape[1]} markers")

    # Build obs DataFrame
    obs_cols = [config.cell_id_col, "sample_id", "donor", "region", "exp_id"]
    base_set = set(marker_cols + obs_cols)
    morphology_cols = [c for c in merged_df.columns if c not in base_set]
    obs_df = merged_df[obs_cols + morphology_cols].copy()
    obs_df.index = pd.Index([f"cell_{i}" for i in range(len(obs_df))], name="cell_id")

    # Build var DataFrame
    var_df = pd.DataFrame(index=marker_cols)

    # Create AnnData
    adata = ad.AnnData(X=expression, obs=obs_df, var=var_df)
    adata.layers["aligned"] = expression.copy()

    # Add raw layer (Stage C normalized but not aligned)
    _log("  [2/4] Building 'raw' layer (concat + reindex)...")
    if any(tbl is not None for tbl in raw_tables):
        raw_dfs = []
        for tbl in raw_tables:
            if tbl is not None:
                tbl = tbl.copy()
                tbl["_idx"] = tbl[config.cell_id_col]
                raw_dfs.append(tbl)
        if raw_dfs:
            raw_concat = pd.concat(raw_dfs, ignore_index=True)
            raw_concat = raw_concat.set_index("_idx")
            # Reindex to match adata ordering
            raw_matrix = raw_concat.reindex(adata.obs[config.cell_id_col]).loc[:, marker_cols].to_numpy(dtype=np.float32)
            adata.layers["raw"] = raw_matrix
    _log("  [2/4] 'raw' layer added")

    # Add batch-corrected layer (Stage E)
    _log("  [3/4] Building 'batchcorr' layer (concat + reindex)...")
    if any(tbl is not None for tbl in corrected_tables):
        corr_dfs = []
        for tbl in corrected_tables:
            if tbl is not None:
                tbl = tbl.copy()
                tbl["_idx"] = tbl[config.cell_id_col]
                corr_dfs.append(tbl)
        if corr_dfs:
            corr_concat = pd.concat(corr_dfs, ignore_index=True)
            corr_concat = corr_concat.set_index("_idx")
            corr_matrix = corr_concat.reindex(adata.obs[config.cell_id_col]).loc[:, marker_cols].to_numpy(dtype=np.float32)
            adata.layers["batchcorr"] = corr_matrix
    _log("  [3/4] 'batchcorr' layer added")

    # Add spatial coordinates to obsm
    _log("  [4/4] Adding spatial coordinates...")
    if {config.x_col, config.y_col}.issubset(adata.obs.columns):
        adata.obsm["spatial"] = adata.obs[[config.x_col, config.y_col]].to_numpy(dtype=np.float32)

    return adata


def generate_stage_f_figures(
    degree_map: Dict[str, np.ndarray],
    all_coords: np.ndarray,
    output_dir: Path,
    dpi: int = 200,
    logger: Optional[logging.Logger] = None,
) -> List[Path]:
    """Generate Stage F visualization figures.

    Parameters
    ----------
    degree_map : Dict[str, np.ndarray]
        Sample ID -> degree array mapping
    all_coords : np.ndarray
        All spatial coordinates stacked
    output_dir : Path
        Output directory for figures
    dpi : int
        Figure resolution
    logger : logging.Logger, optional
        Logger

    Returns
    -------
    List[Path]
        List of generated figure paths
    """
    import matplotlib.pyplot as plt

    figures_dir = output_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    generated = []

    # 1. Neighbor degree distribution
    fig, ax = plt.subplots(figsize=(10, 6))

    # Combine all degrees for histogram
    all_degrees = np.concatenate([d for d in degree_map.values() if d.size > 0])
    if all_degrees.size > 0:
        ax.hist(all_degrees, bins=30, alpha=0.7, color="steelblue", edgecolor="white")
        ax.axvline(np.mean(all_degrees), color="red", linestyle="--", label=f"Mean: {np.mean(all_degrees):.1f}")
        ax.set_xlabel("Degree (number of neighbors)")
        ax.set_ylabel("Cell count")
        ax.set_title("Spatial Graph Degree Distribution")
        ax.legend()

    path = figures_dir / "stage_f_neighbor_degree.png"
    fig.tight_layout()
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    generated.append(path)
    if logger:
        logger.info(f"Generated: {path.name}")

    # 2. Spatial density map
    if all_coords.size > 0:
        fig, ax = plt.subplots(figsize=(8, 7))
        h = ax.hist2d(all_coords[:, 0], all_coords[:, 1], bins=50, cmap="viridis")
        fig.colorbar(h[3], ax=ax, label="Cell count")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_title("Spatial Cell Density Map")

        path = figures_dir / "stage_f_density_map.png"
        fig.tight_layout()
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        generated.append(path)
        if logger:
            logger.info(f"Generated: {path.name}")

    # 3. Per-sample cell counts
    sample_counts = {sid: len(deg) for sid, deg in degree_map.items()}
    if sample_counts:
        fig, ax = plt.subplots(figsize=(12, 6))
        samples = list(sample_counts.keys())
        counts = [sample_counts[s] for s in samples]

        bars = ax.bar(range(len(samples)), counts, color="steelblue", edgecolor="white")
        ax.set_xticks(range(len(samples)))
        ax.set_xticklabels(samples, rotation=45, ha="right", fontsize=8)
        ax.set_xlabel("Sample")
        ax.set_ylabel("Cell count")
        ax.set_title("Cell Counts per Sample")

        # Add value labels
        for bar, count in zip(bars, counts):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
                    f"{count:,}", ha="center", va="bottom", fontsize=7)

        path = figures_dir / "stage_f_cell_counts.png"
        fig.tight_layout()
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        generated.append(path)
        if logger:
            logger.info(f"Generated: {path.name}")

    return generated


def run_stage_f(
    aligned_dir: Path,
    raw_dir: Path,
    corrected_dir: Path,
    filtered_meta_dir: Path,
    metadata_path: Path,
    output_dir: Path,
    config: Optional[MergeConfig] = None,
    loader_config: Optional[LoaderConfig] = None,
    variant: Optional[str] = None,
    use_corrected: bool = True,
    verbose: bool = False,
    skip_figures: bool = False,
    dpi: int = 200,
    log_dir: Optional[Path] = None,
) -> "anndata.AnnData":
    """Run Stage F: Data merge and spatial graph construction.

    Merges all samples into a single AnnData object with multiple layers:
    - X / aligned: Stage D aligned data
    - raw: Stage C normalized (not aligned) data
    - batchcorr: Stage E batch-corrected data

    Also builds spatial neighborhood graphs per sample.

    Parameters
    ----------
    aligned_dir : Path
        Directory with Stage D aligned CSVs
    raw_dir : Path
        Directory with Stage C export CSVs
    corrected_dir : Path
        Directory with Stage E corrected CSVs
    filtered_meta_dir : Path
        Directory with Stage B filtered metadata CSVs
    metadata_path : Path
        Path to Stage A metadata_verified.csv
    output_dir : Path
        Output directory
    config : Optional[MergeConfig]
        Merge configuration
    loader_config : Optional[LoaderConfig]
        Loader configuration
    variant : Optional[str]
        Stage C export variant (e.g., clip-Q3__log1p)
    use_corrected : bool
        Whether to include batch-corrected layer
    verbose : bool
        Enable verbose logging
    skip_figures : bool
        Skip figure generation
    dpi : int
        Figure resolution
    log_dir : Optional[Path]
        Directory for log file

    Returns
    -------
    anndata.AnnData
        Merged AnnData object
    """
    # Setup logging
    logger = setup_logging(verbose, log_dir=log_dir, log_filename="stage_f.log")
    logger.info("Stage F: Data Merge and Spatial Graph Construction")
    logger.info(f"Aligned dir: {aligned_dir}")
    logger.info(f"Raw dir: {raw_dir}")
    logger.info(f"Corrected dir: {corrected_dir}")
    logger.info(f"Filtered meta dir: {filtered_meta_dir}")
    logger.info(f"Metadata: {metadata_path}")
    logger.info(f"Output: {output_dir}")
    logger.info(f"Use corrected: {use_corrected}")

    # Setup
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    graphs_dir = output_dir / "graphs"
    graphs_dir.mkdir(parents=True, exist_ok=True)

    config = config or MergeConfig()
    loader_config = loader_config or LoaderConfig()
    cell_id_col = config.cell_id_col

    # Parse variant for raw layer
    if variant is None:
        variant = "clip-Q3__log1p"
    logger.info(f"Raw layer variant: {variant}")

    # Load metadata
    meta_df = pd.read_csv(metadata_path)
    sample_ids = sorted(meta_df["sample_id"].astype(str).unique())
    n_samples = len(sample_ids)
    logger.info(f"Found {n_samples} samples")

    # Process samples
    merged_tables: List[pd.DataFrame] = []
    raw_tables: List[Optional[pd.DataFrame]] = []
    corrected_tables: List[Optional[pd.DataFrame]] = []
    marker_cols: Optional[List[str]] = None

    for i, sample_id in enumerate(sample_ids, 1):
        logger.info(f"[{i}/{n_samples}] Loading {sample_id}")

        # Get metadata row
        meta_row = meta_df[meta_df["sample_id"].astype(str) == sample_id].iloc[0]

        # Load aligned matrix (required)
        aligned_path = Path(aligned_dir) / f"{sample_id}_aligned.csv"
        if not aligned_path.exists():
            raise FileNotFoundError(f"Aligned matrix not found: {aligned_path}")

        aligned_df = pd.read_csv(aligned_path)
        aligned_df[cell_id_col] = aligned_df[cell_id_col].astype(str)
        aligned_df = _sanitize_dataframe_columns(aligned_df, exclude={cell_id_col}, logger=logger)
        markers = [c for c in aligned_df.columns if c != cell_id_col]

        if marker_cols is None:
            marker_cols = markers

        # Load raw matrix (Stage C export)
        raw_df = None
        raw_path = Path(raw_dir) / f"{sample_id}__{variant}.csv"
        if raw_path.exists():
            raw_df = pd.read_csv(raw_path)
            raw_df[cell_id_col] = raw_df[cell_id_col].astype(str)
            raw_df = _sanitize_dataframe_columns(raw_df, exclude={cell_id_col}, logger=logger)
        else:
            logger.warning(f"Raw matrix not found: {raw_path}")

        # Load corrected matrix (Stage E)
        corrected_df = None
        if use_corrected:
            corrected_path = Path(corrected_dir) / f"{sample_id}_corrected.csv"
            if corrected_path.exists():
                corrected_df = pd.read_csv(corrected_path)
                corrected_df[cell_id_col] = corrected_df[cell_id_col].astype(str)
                corrected_df = _sanitize_dataframe_columns(corrected_df, exclude={cell_id_col}, logger=logger)
            else:
                logger.warning(f"Corrected matrix not found: {corrected_path}")

        # Load cell metadata (from Stage B filtered)
        cell_meta_path = Path(filtered_meta_dir) / f"{sample_id}_filtered_metadata.csv"
        if not cell_meta_path.exists():
            # Fallback to original metadata
            orig_meta_path = meta_row.get("cell_metadata_path")
            if orig_meta_path and Path(orig_meta_path).exists():
                cell_meta_path = Path(orig_meta_path)
                logger.warning(f"Using original metadata: {cell_meta_path}")
            else:
                raise FileNotFoundError(f"Cell metadata not found for {sample_id}")

        cell_meta = pd.read_csv(cell_meta_path)
        cell_meta[cell_id_col] = cell_meta[cell_id_col].astype(str)
        cell_meta = _sanitize_dataframe_columns(cell_meta, exclude={cell_id_col}, logger=logger)

        # Merge aligned with cell metadata
        merged = aligned_df.merge(cell_meta, on=cell_id_col, how="left", suffixes=("", "_meta"))
        merged["sample_id"] = sample_id
        merged["donor"] = meta_row.get("donor")
        merged["region"] = meta_row.get("region")
        merged["exp_id"] = meta_row.get("exp_id")

        merged_tables.append(merged)
        raw_tables.append(raw_df)
        corrected_tables.append(corrected_df)

        logger.debug(f"  {sample_id}: {len(merged)} cells")

    if marker_cols is None:
        raise RuntimeError("No markers discovered from aligned data")

    # Build AnnData
    logger.info("Building AnnData object...")
    adata = build_anndata(
        merged_tables=merged_tables,
        raw_tables=raw_tables,
        corrected_tables=corrected_tables,
        marker_cols=marker_cols,
        config=config,
        logger=logger,
    )
    logger.info(f"AnnData: {adata.n_obs} cells, {adata.n_vars} markers")
    logger.info(f"Layers: {list(adata.layers.keys())}")

    # Save AnnData
    h5ad_path = output_dir / "merged_data.h5ad"
    adata.write_h5ad(h5ad_path)
    logger.info(f"Wrote {h5ad_path}")

    # Build spatial graphs per sample
    logger.info("Building spatial graphs...")
    degree_map: Dict[str, np.ndarray] = {}
    all_coords_list: List[np.ndarray] = []

    x_col = config.x_col
    y_col = config.y_col
    obs_df = adata.obs

    for sample_id in sample_ids:
        mask = obs_df["sample_id"] == sample_id

        if {x_col, y_col}.issubset(obs_df.columns):
            coords = obs_df.loc[mask, [x_col, y_col]].to_numpy(dtype=np.float32)
        else:
            coords = np.zeros((mask.sum(), 2), dtype=np.float32)

        all_coords_list.append(coords)

        # Build graph
        graph = build_spatial_graph(
            coords=coords,
            mode=config.graph_mode,
            k=config.k,
            radius=config.radius,
        )

        # Compute degrees
        degrees = np.asarray(graph.sum(axis=1)).ravel()
        degree_map[sample_id] = degrees

        # Save graph
        graph_path = graphs_dir / f"{sample_id}_neighbors.npz"
        sp.save_npz(graph_path, graph)

        logger.debug(f"  {sample_id}: {len(coords)} cells, mean degree {degrees.mean():.1f}")

    logger.info(f"Built {len(degree_map)} spatial graphs")

    # Generate figures
    if not skip_figures:
        logger.info("Generating figures...")
        all_coords = np.vstack(all_coords_list) if all_coords_list else np.array([])
        figures = generate_stage_f_figures(
            degree_map=degree_map,
            all_coords=all_coords,
            output_dir=output_dir,
            dpi=dpi,
            logger=logger,
        )
        logger.info(f"Generated {len(figures)} figures")

    total_cells = adata.n_obs
    logger.info(f"Stage F completed: {total_cells:,} cells merged, {len(marker_cols)} markers")

    return adata


# =============================================================================
# Main Stage A Runner
# =============================================================================


def run_stage_a(
    metadata_path: Path,
    output_dir: Path,
    config: Optional[LoaderConfig] = None,
    verbose: bool = False,
    skip_figures: bool = False,
    dpi: int = 200,
    log_dir: Optional[Path] = None,
) -> pd.DataFrame:
    """Run Stage A: Data loading and validation.

    Parameters
    ----------
    metadata_path : Path
        Path to metadata registry CSV
    output_dir : Path
        Output directory
    config : Optional[LoaderConfig]
        Loader configuration
    verbose : bool
        Enable verbose logging
    skip_figures : bool
        Skip figure generation
    dpi : int
        Figure resolution
    log_dir : Optional[Path]
        Directory for log file (if None, console only)

    Returns
    -------
    pd.DataFrame
        Verified metadata with validation results
    """
    logger = setup_logging(verbose, log_dir=log_dir)
    logger.info("Stage A: Data Loading and Validation")
    logger.info(f"Metadata: {metadata_path}")
    logger.info(f"Output: {output_dir}")

    # Setup
    output_dir.mkdir(parents=True, exist_ok=True)
    config = config or LoaderConfig()
    loader = DataLoader(config)

    # Load registry
    registry = loader.load_metadata_registry(metadata_path)
    n_samples = len(registry)
    logger.info(f"Found {n_samples} samples in registry")

    # Track reference markers for consistency check
    reference_markers = None
    results_rows = []
    coverage_records: List[Dict[str, Any]] = []

    # Process each sample
    for idx, row in registry.iterrows():
        sample_id = row.get("exp_id", row.get("sample_id", f"sample_{idx}"))
        matrix_path = Path(row["cell_matrix_path"])
        metadata_row_path = Path(row["cell_metadata_path"])

        logger.info(f"[{idx+1}/{n_samples}] Processing {sample_id}")

        # Load and validate
        result = loader.load_sample(
            matrix_path=matrix_path,
            metadata_path=metadata_row_path,
            sample_id=sample_id,
            reference_markers=reference_markers,
        )

        # Use first sample's markers as reference
        if reference_markers is None and result.markers:
            reference_markers = result.markers
            logger.debug(f"Reference markers set: {len(reference_markers)} markers")

        # Compute marker coverage for visualization
        if result.cell_matrix is not None and result.markers:
            coverage = compute_marker_coverage(
                result.cell_matrix,
                result.markers,
                config.cell_id_col,
            )
            for marker, fraction in coverage.items():
                coverage_records.append({
                    "sample_id": sample_id,
                    "marker": marker,
                    "coverage": float(fraction),
                })

        # Convert to output row
        out_row = load_result_to_row(result, row)
        results_rows.append(out_row)

        # Log status
        status_str = "OK" if result.status == "OK" else f"CHECK: {';'.join(result.issues)}"
        logger.info(f"  {result.n_cells:,} cells, {len(result.markers)} markers -> {status_str}")

    # Create output DataFrame
    output_df = pd.DataFrame(results_rows)

    # Ensure column order matches reference format
    column_order = [
        "exp_id", "sample_id", "donor", "region",
        "image_path", "cell_matrix_path", "cell_metadata_path", "Notes",
        "matrix_rows", "matrix_markers", "metadata_rows",
        "duplicate_ids_matrix", "duplicate_ids_metadata",
        "missing_in_matrix", "missing_in_metadata",
        "status", "issues",
    ]
    # Only include columns that exist
    column_order = [c for c in column_order if c in output_df.columns]
    output_df = output_df[column_order]

    # Write output
    output_path = output_dir / "metadata_verified.csv"
    output_df.to_csv(output_path, index=False)
    logger.info(f"Wrote {output_path}")

    # Generate figures
    if not skip_figures:
        logger.info("Generating figures...")
        figures = generate_figures(
            output_df,
            coverage_records,
            output_dir,
            dpi=dpi,
            logger=logger,
        )
        logger.info(f"Generated {len(figures)} figures")

    # Summary
    n_ok = (output_df["status"] == "OK").sum()
    n_check = (output_df["status"] == "CHECK").sum()
    logger.info(f"Summary: {n_ok} OK, {n_check} CHECK")

    return output_df


def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="CellType-Refinery Preprocessing Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Stage A: Data Loading and Validation
  python -m celltype_refinery.core.preprocessing --stage A \\
      --metadata data/fallopian_tube/metadata.csv \\
      --output output/fallopian_tube/stage_a

  # Stage B: Cell-Level Quality Control
  python -m celltype_refinery.core.preprocessing --stage B \\
      --input output/fallopian_tube/stage_a/metadata_verified.csv \\
      --output output/fallopian_tube/stage_b

  # Stage C: Within-Sample Normalization
  python -m celltype_refinery.core.preprocessing --stage C \\
      --input output/fallopian_tube/stage_b/filtered/matrices \\
      --output output/fallopian_tube/stage_c

  # Stage D: Cross-Sample Alignment
  python -m celltype_refinery.core.preprocessing --stage D \\
      --input output/fallopian_tube/stage_c/export \\
      --output output/fallopian_tube/stage_d \\
      --variant clip-Q3__log1p

  # Stage E: Batch Effect Removal
  python -m celltype_refinery.core.preprocessing --stage E \\
      --input output/fallopian_tube/stage_d/aligned \\
      --output output/fallopian_tube/stage_e \\
      --metadata output/fallopian_tube/stage_a/metadata_verified.csv \\
      --channels-dir data/fallopian_tube/channels \\
      --variant clip-Q3__log1p

  # Stage F: Data Merge and Spatial Graphs
  python -m celltype_refinery.core.preprocessing --stage F \\
      --aligned-dir output/fallopian_tube/stage_d/aligned \\
      --raw-dir output/fallopian_tube/stage_c/export \\
      --corrected-dir output/fallopian_tube/stage_e/corrected \\
      --filtered-meta-dir output/fallopian_tube/stage_b/filtered/metadata \\
      --metadata output/fallopian_tube/stage_a/metadata_verified.csv \\
      --output output/fallopian_tube/stage_f \\
      --variant clip-Q3__log1p

        """,
    )
    parser.add_argument(
        "--stage", "-s",
        type=str,
        choices=["A", "B", "C", "D", "E", "F"],
        default="A",
        help="Stage to run: A (data loading), B (cell QC), C (normalization), D (alignment), E (batch correction), F (merge)",
    )
    parser.add_argument(
        "--metadata", "-m",
        type=Path,
        default=None,
        help="Path to metadata registry CSV (Stage A input)",
    )
    parser.add_argument(
        "--input", "-i",
        type=Path,
        default=None,
        help="Path to input file/directory (B: metadata_verified.csv, C: filtered matrices, D: Stage C export, E: Stage D aligned)",
    )
    parser.add_argument(
        "--output", "-o",
        type=Path,
        required=True,
        help="Output directory",
    )
    parser.add_argument(
        "--config", "-c",
        type=Path,
        default=None,
        help="Path to preprocessing config YAML (optional)",
    )
    parser.add_argument(
        "--variant",
        type=str,
        default=None,
        help="Normalization variant for Stage C/D/E (e.g., clip-Q3__log1p)",
    )
    parser.add_argument(
        "--channels-dir",
        type=Path,
        default=None,
        help="Directory with channel metadata CSV files (Stage E, format: {donor}_{exp_id}.csv)",
    )
    parser.add_argument(
        "--correction-mode",
        type=str,
        choices=["residual", "none"],
        default="residual",
        help="Batch correction mode (Stage E): 'residual' applies correction, 'none' skips (default: residual)",
    )
    # Stage F arguments
    parser.add_argument(
        "--aligned-dir",
        type=Path,
        default=None,
        help="Directory with Stage D aligned CSVs (Stage F)",
    )
    parser.add_argument(
        "--raw-dir",
        type=Path,
        default=None,
        help="Directory with Stage C export CSVs (Stage F)",
    )
    parser.add_argument(
        "--corrected-dir",
        type=Path,
        default=None,
        help="Directory with Stage E corrected CSVs (Stage F)",
    )
    parser.add_argument(
        "--filtered-meta-dir",
        type=Path,
        default=None,
        help="Directory with Stage B filtered metadata CSVs (Stage F)",
    )
    parser.add_argument(
        "--graph-mode",
        type=str,
        choices=["knn", "ball"],
        default="knn",
        help="Spatial graph type: 'knn' (k-nearest neighbors) or 'ball' (radius-based) (Stage F)",
    )
    parser.add_argument(
        "--k",
        type=int,
        default=10,
        help="Number of neighbors for kNN graphs (Stage F, default: 10)",
    )
    parser.add_argument(
        "--radius",
        type=float,
        default=20.0,
        help="Radius for ball graphs in spatial units (Stage F, default: 20.0)",
    )
    parser.add_argument(
        "--use-corrected",
        action="store_true",
        default=True,
        help="Include batch-corrected layer from Stage E (Stage F, default: True)",
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable verbose logging",
    )
    parser.add_argument(
        "--skip-figures",
        action="store_true",
        help="Skip figure generation",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=200,
        help="Figure resolution (default: 200)",
    )
    parser.add_argument(
        "--log-dir", "-l",
        type=Path,
        default=None,
        help="Directory for log file (default: console only)",
    )

    args = parser.parse_args()

    # Load config if provided
    full_config = None
    if args.config:
        full_config = PreprocessingConfig.from_yaml(args.config)

    try:
        if args.stage == "A":
            # Stage A: Data Loading
            if args.metadata is None:
                parser.error("Stage A requires --metadata argument")

            loader_config = full_config.loader if full_config else None
            run_stage_a(
                metadata_path=args.metadata,
                output_dir=args.output,
                config=loader_config,
                verbose=args.verbose,
                skip_figures=args.skip_figures,
                dpi=args.dpi,
                log_dir=args.log_dir,
            )

        elif args.stage == "B":
            # Stage B: Cell QC
            if args.input is None:
                parser.error("Stage B requires --input argument (path to Stage A metadata_verified.csv)")

            qc_config = full_config.qc if full_config else None
            loader_config = full_config.loader if full_config else None
            run_stage_b(
                input_path=args.input,
                output_dir=args.output,
                config=qc_config,
                loader_config=loader_config,
                verbose=args.verbose,
                skip_figures=args.skip_figures,
                dpi=args.dpi,
                log_dir=args.log_dir,
            )

        elif args.stage == "C":
            # Stage C: Normalization
            if args.input is None:
                parser.error("Stage C requires --input argument (path to Stage B filtered matrices directory)")

            norm_config = full_config.normalization if full_config else None
            loader_config = full_config.loader if full_config else None
            run_stage_c(
                input_dir=args.input,
                output_dir=args.output,
                config=norm_config,
                loader_config=loader_config,
                variant=args.variant,
                verbose=args.verbose,
                skip_figures=args.skip_figures,
                dpi=args.dpi,
                log_dir=args.log_dir,
            )

        elif args.stage == "D":
            # Stage D: Cross-sample alignment
            if args.input is None:
                parser.error("Stage D requires --input argument (path to Stage C export directory)")

            align_config = full_config.alignment if full_config else None
            loader_config = full_config.loader if full_config else None
            run_stage_d(
                input_dir=args.input,
                output_dir=args.output,
                config=align_config,
                loader_config=loader_config,
                variant=args.variant,
                verbose=args.verbose,
                skip_figures=args.skip_figures,
                dpi=args.dpi,
                log_dir=args.log_dir,
            )

        elif args.stage == "E":
            # Stage E: Batch effect removal
            if args.input is None:
                parser.error("Stage E requires --input argument (path to Stage D aligned directory)")
            if args.metadata is None:
                parser.error("Stage E requires --metadata argument (path to Stage A metadata_verified.csv)")

            batch_config = full_config.batch_correction if full_config else None
            loader_config = full_config.loader if full_config else None
            run_stage_e(
                input_dir=args.input,
                output_dir=args.output,
                metadata_path=args.metadata,
                channels_dir=args.channels_dir,
                config=batch_config,
                loader_config=loader_config,
                correction_mode=args.correction_mode,
                verbose=args.verbose,
                skip_figures=args.skip_figures,
                dpi=args.dpi,
                log_dir=args.log_dir,
            )

        elif args.stage == "F":
            # Stage F: Data merge and spatial graphs
            if args.aligned_dir is None:
                parser.error("Stage F requires --aligned-dir argument (path to Stage D aligned directory)")
            if args.raw_dir is None:
                parser.error("Stage F requires --raw-dir argument (path to Stage C export directory)")
            if args.corrected_dir is None:
                parser.error("Stage F requires --corrected-dir argument (path to Stage E corrected directory)")
            if args.filtered_meta_dir is None:
                parser.error("Stage F requires --filtered-meta-dir argument (path to Stage B filtered metadata directory)")
            if args.metadata is None:
                parser.error("Stage F requires --metadata argument (path to Stage A metadata_verified.csv)")

            merge_config = MergeConfig(
                graph_mode=args.graph_mode,
                k=args.k,
                radius=args.radius,
            )
            loader_config = full_config.loader if full_config else None
            run_stage_f(
                aligned_dir=args.aligned_dir,
                raw_dir=args.raw_dir,
                corrected_dir=args.corrected_dir,
                filtered_meta_dir=args.filtered_meta_dir,
                metadata_path=args.metadata,
                output_dir=args.output,
                config=merge_config,
                loader_config=loader_config,
                variant=args.variant,
                use_corrected=args.use_corrected,
                verbose=args.verbose,
                skip_figures=args.skip_figures,
                dpi=args.dpi,
                log_dir=args.log_dir,
            )

    except Exception as e:
        logging.error(f"Stage {args.stage} failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
