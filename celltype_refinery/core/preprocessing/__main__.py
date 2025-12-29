"""Preprocessing module CLI runner.

Enables running preprocessing stages as:
    python -m celltype_refinery.core.preprocessing --stage A --metadata <path> --output <dir>
    python -m celltype_refinery.core.preprocessing --stage B --input <path> --output <dir>

Usage Examples:
    # Run Stage A on fallopian tube data
    python -m celltype_refinery.core.preprocessing --stage A \
        --metadata data/fallopian_tube/metadata.csv \
        --output output/fallopian_tube/stage_a

    # Run Stage B (cell QC) on Stage A output
    python -m celltype_refinery.core.preprocessing --stage B \
        --input output/fallopian_tube/stage_a/metadata_verified.csv \
        --output output/fallopian_tube/stage_b

    # With verbose logging and log file
    python -m celltype_refinery.core.preprocessing --stage B \
        --input output/fallopian_tube/stage_a/metadata_verified.csv \
        --output output/fallopian_tube/stage_b \
        --log-dir logs/fallopian_tube \
        --verbose
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd

from .config import LoaderConfig, QCConfig, PreprocessingConfig
from .loader import DataLoader, LoadResult
from .qc import CellQC, QCResult, REASON_COLUMNS


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

  # With config file
  python -m celltype_refinery.core.preprocessing --stage B \\
      --input output/fallopian_tube/stage_a/metadata_verified.csv \\
      --output output/fallopian_tube/stage_b \\
      --config configs/preprocessing.yaml
        """,
    )
    parser.add_argument(
        "--stage", "-s",
        type=str,
        choices=["A", "B"],
        default="A",
        help="Stage to run: A (data loading) or B (cell QC)",
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
        help="Path to input file (Stage B: metadata_verified.csv from Stage A)",
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

    except Exception as e:
        logging.error(f"Stage {args.stage} failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
