"""Command-line interface for CellType-Refinery.

Provides CLI commands for running pipeline stages.

Example Usage
-------------
    # From command line:
    celltype-refinery --help
    celltype-refinery cluster --input data.h5ad --out out/
    celltype-refinery annotate --input out/clustered.h5ad --marker-map markers.json
    celltype-refinery pipeline --config pipeline.yaml
"""

__version__ = "1.0.0"

from .main import cli, main

__all__ = [
    "__version__",
    "cli",
    "main",
]
