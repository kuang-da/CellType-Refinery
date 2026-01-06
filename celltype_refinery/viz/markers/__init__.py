"""Marker map visualization module.

Provides visualizations for hierarchical marker maps used in cell-type annotation:
- ASCII tree representation for console output
- Tree diagram (matplotlib PNG)
- Summary table (matplotlib PNG)
- Network graph (matplotlib PNG)
- Sunburst chart (interactive Plotly HTML)

Usage
-----
>>> from celltype_refinery.viz.markers import generate_all_figures
>>> import json
>>> with open("marker_map.json") as f:
...     marker_map = json.load(f)
>>> results = generate_all_figures(
...     marker_map,
...     output_dir="output/marker_viz",
...     prefix="FT_v9",
... )

CLI Usage
---------
python -m celltype_refinery.viz.markers --marker-map data/markers.json --out output/viz
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, Optional, Sequence

from .style import (
    # Constants
    CATEGORY_COLORS,
    MARKER_COLORS,
    # Data structures
    MarkerNode,
    # Utilities
    PathLike,
    format_marker_list,
    get_category_color,
    parse_marker_hierarchy,
)

logger = logging.getLogger(__name__)


def generate_all_figures(
    marker_map: Dict[str, Any],
    output_dir: PathLike,
    *,
    prefix: str = "marker_hierarchy",
    formats: Sequence[str] = ("tree", "sunburst", "network", "ascii", "table"),
) -> Dict[str, Optional[Path]]:
    """
    Generate all marker map visualizations at once.

    Parameters
    ----------
    marker_map : Dict[str, Any]
        Marker hierarchy loaded from JSON file.
    output_dir : PathLike
        Directory for output files.
    prefix : str
        Filename prefix.
    formats : Sequence[str]
        Which visualizations to generate: "tree", "sunburst", "network",
        "ascii", "table".

    Returns
    -------
    Dict[str, Optional[Path]]
        Mapping of format name to output path. None if generation failed.
    """
    from .hierarchy_ascii import print_marker_hierarchy_ascii
    from .hierarchy_network import plot_marker_network
    from .hierarchy_sunburst import plot_marker_hierarchy_sunburst
    from .hierarchy_tree import plot_marker_hierarchy_tree
    from .summary_table import plot_marker_summary_table

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    results: Dict[str, Optional[Path]] = {}

    if "ascii" in formats:
        try:
            ascii_output = print_marker_hierarchy_ascii(marker_map)
            ascii_path = output_dir / f"{prefix}_ascii.txt"
            ascii_path.write_text(ascii_output)
            results["ascii"] = ascii_path
            logger.info(f"Generated ASCII tree: {ascii_path}")
        except Exception as e:
            logger.warning(f"Failed to generate ASCII tree: {e}")
            results["ascii"] = None

    if "tree" in formats:
        try:
            tree_path = output_dir / f"{prefix}_tree.png"
            results["tree"] = plot_marker_hierarchy_tree(marker_map, tree_path)
            logger.info(f"Generated tree diagram: {tree_path}")
        except Exception as e:
            logger.warning(f"Failed to generate tree diagram: {e}")
            results["tree"] = None

    if "table" in formats:
        try:
            table_path = output_dir / f"{prefix}_table.png"
            results["table"] = plot_marker_summary_table(marker_map, table_path)
            logger.info(f"Generated summary table: {table_path}")
        except Exception as e:
            logger.warning(f"Failed to generate summary table: {e}")
            results["table"] = None

    if "network" in formats:
        try:
            network_path = output_dir / f"{prefix}_network.png"
            results["network"] = plot_marker_network(marker_map, network_path)
            logger.info(f"Generated network graph: {network_path}")
        except Exception as e:
            logger.warning(f"Failed to generate network graph: {e}")
            results["network"] = None

    if "sunburst" in formats:
        try:
            sunburst_path = output_dir / f"{prefix}_sunburst.html"
            results["sunburst"] = plot_marker_hierarchy_sunburst(marker_map, sunburst_path)
            logger.info(f"Generated sunburst chart: {sunburst_path}")
        except Exception as e:
            logger.warning(f"Failed to generate sunburst chart: {e}")
            results["sunburst"] = None

    return results


# Re-export individual plot functions for direct access
from .hierarchy_ascii import print_marker_hierarchy_ascii
from .hierarchy_network import plot_marker_network
from .hierarchy_sunburst import plot_marker_hierarchy_sunburst
from .hierarchy_tree import plot_marker_hierarchy_tree
from .summary_table import plot_marker_summary_table

__all__ = [
    # Main orchestrator
    "generate_all_figures",
    # Individual plot functions
    "print_marker_hierarchy_ascii",
    "plot_marker_hierarchy_tree",
    "plot_marker_summary_table",
    "plot_marker_network",
    "plot_marker_hierarchy_sunburst",
    # Constants
    "MARKER_COLORS",
    "CATEGORY_COLORS",
    # Data structures
    "MarkerNode",
    # Utilities
    "parse_marker_hierarchy",
    "format_marker_list",
    "get_category_color",
]
