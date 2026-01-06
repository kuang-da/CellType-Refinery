"""Visualization module for CellType-Refinery.

Provides plotting functions for composition, spatial, and marker visualizations.
Requires optional 'viz' dependencies: matplotlib, plotly.

Submodules
----------
markers
    Marker map hierarchy visualizations (tree, sunburst, network, etc.)
"""

from . import markers

__all__ = ["markers"]
