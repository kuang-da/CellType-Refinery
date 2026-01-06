#!/usr/bin/env python3
"""
CLI entry point for marker map visualization.

Usage:
    python -m celltype_refinery.viz.markers --marker-map data/markers.json --out output/viz
    python -m celltype_refinery.viz.markers --marker-map data/markers.json --out output/viz --formats tree sunburst
"""

import argparse
import json
import logging
import sys
from pathlib import Path

from . import generate_all_figures

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def main() -> int:
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Generate marker hierarchy visualizations (tree, sunburst, network, table, ascii)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate all visualizations
  python -m celltype_refinery.viz.markers \\
      --marker-map data/FT_cell_type_markers_v9.json \\
      --out output/marker_viz

  # Generate specific formats only
  python -m celltype_refinery.viz.markers \\
      --marker-map data/markers.json \\
      --out output/viz \\
      --formats tree sunburst

  # Custom prefix for output files
  python -m celltype_refinery.viz.markers \\
      --marker-map data/markers.json \\
      --out output/viz \\
      --prefix my_markers
        """,
    )
    parser.add_argument(
        "--marker-map", "-m",
        type=str,
        required=True,
        help="Path to marker map JSON file",
    )
    parser.add_argument(
        "--out", "-o",
        type=str,
        required=True,
        help="Output directory",
    )
    parser.add_argument(
        "--prefix", "-p",
        type=str,
        default=None,
        help="Filename prefix (default: derived from input filename)",
    )
    parser.add_argument(
        "--formats", "-f",
        nargs="+",
        choices=["tree", "sunburst", "network", "ascii", "table"],
        default=["tree", "sunburst", "network", "ascii", "table"],
        help="Visualization formats to generate (default: all)",
    )
    parser.add_argument(
        "--log-level",
        type=str,
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
        help="Logging level (default: INFO)",
    )

    args = parser.parse_args()

    # Set log level
    logging.getLogger().setLevel(getattr(logging, args.log_level))

    # Resolve paths
    marker_path = Path(args.marker_map)
    if not marker_path.exists():
        logger.error(f"Marker file not found: {marker_path}")
        return 1

    output_dir = Path(args.out)

    # Derive prefix from filename if not provided
    prefix = args.prefix
    if prefix is None:
        prefix = marker_path.stem  # e.g., "FT_cell_type_markers_v9"

    # Load marker map
    logger.info(f"Loading marker map from: {marker_path}")
    try:
        with open(marker_path, "r") as f:
            marker_map = json.load(f)
    except json.JSONDecodeError as e:
        logger.error(f"Invalid JSON in marker file: {e}")
        return 1

    logger.info(f"Generating visualizations: {args.formats}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Prefix: {prefix}")

    # Generate visualizations
    results = generate_all_figures(
        marker_map,
        output_dir,
        prefix=prefix,
        formats=tuple(args.formats),
    )

    # Print results
    print("\n" + "=" * 60)
    print("  Marker Visualization Complete")
    print("=" * 60)
    print(f"\nOutput directory: {output_dir}")
    print("\nGenerated files:")
    for fmt, path in results.items():
        if path is not None:
            print(f"  {fmt}: {path}")
        else:
            print(f"  {fmt}: FAILED")

    # Return non-zero if any failed
    failed = [fmt for fmt, path in results.items() if path is None]
    if failed:
        logger.warning(f"Some visualizations failed: {failed}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
