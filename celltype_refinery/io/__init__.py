"""I/O utilities for CellType-Refinery.

Provides logging, CSV I/O, and data loading utilities.
"""

from .logging import get_logger, get_timestamped_log_path, log_json, log_yaml
from .csv import (
    ensure_output_dir,
    load_cell_matrix,
    load_cell_metadata,
    load_metadata_registry,
    write_dataframe,
    unique_marker_columns,
    load_marker_map_from_source,
)

__all__ = [
    # Logging
    "get_logger",
    "get_timestamped_log_path",
    "log_json",
    "log_yaml",
    # CSV I/O
    "ensure_output_dir",
    "load_cell_matrix",
    "load_cell_metadata",
    "load_metadata_registry",
    "write_dataframe",
    "unique_marker_columns",
    "load_marker_map_from_source",
]
