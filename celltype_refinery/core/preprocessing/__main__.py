"""Stage A data loading and validation runner.

Enables running the preprocessing module as:
    python -m celltype_refinery.core.preprocessing --metadata <path> --output <dir>

Usage Examples:
    # Run Stage A on fallopian tube data
    python -m celltype_refinery.core.preprocessing \
        --metadata data/fallopian_tube/metadata.csv \
        --output output/fallopian_tube/stage_a

    # With verbose logging
    python -m celltype_refinery.core.preprocessing \
        --metadata data/fallopian_tube/metadata.csv \
        --output output/fallopian_tube/stage_a \
        --verbose

    # Skip figure generation
    python -m celltype_refinery.core.preprocessing \
        --metadata data/fallopian_tube/metadata.csv \
        --output output/fallopian_tube/stage_a \
        --skip-figures
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd

from .config import LoaderConfig, PreprocessingConfig
from .loader import DataLoader, LoadResult


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
        description="CellType-Refinery Stage A: Data Loading and Validation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python -m celltype_refinery.core.preprocessing \\
      --metadata data/fallopian_tube/metadata.csv \\
      --output output/fallopian_tube/stage_a

  python -m celltype_refinery.core.preprocessing \\
      --metadata data/ovary/metadata.csv \\
      --output output/ovary/stage_a \\
      --config configs/preprocessing.yaml
        """,
    )
    parser.add_argument(
        "--metadata", "-m",
        type=Path,
        required=True,
        help="Path to metadata registry CSV",
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
    config = None
    if args.config:
        full_config = PreprocessingConfig.from_yaml(args.config)
        config = full_config.loader

    # Run Stage A
    try:
        run_stage_a(
            metadata_path=args.metadata,
            output_dir=args.output,
            config=config,
            verbose=args.verbose,
            skip_figures=args.skip_figures,
            dpi=args.dpi,
            log_dir=args.log_dir,
        )
    except Exception as e:
        logging.error(f"Stage A failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
