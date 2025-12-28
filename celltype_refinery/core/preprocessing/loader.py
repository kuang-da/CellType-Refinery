"""Data loader for preprocessing pipeline (Stage A).

Handles loading cell matrices, cell metadata, and metadata registry
with validation and consistency checks.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd

from .config import LoaderConfig


@dataclass
class LoadResult:
    """Result from loading a single sample.

    Attributes
    ----------
    sample_id : str
        Sample identifier
    cell_matrix : pd.DataFrame
        Cell-by-marker matrix
    cell_metadata : pd.DataFrame
        Cell metadata with coordinates
    markers : Set[str]
        Marker names
    n_cells : int
        Number of cells
    issues : List[str]
        List of any issues found
    status : str
        'OK' or 'CHECK'
    """

    sample_id: str
    cell_matrix: Optional[pd.DataFrame] = None
    cell_metadata: Optional[pd.DataFrame] = None
    markers: Set[str] = field(default_factory=set)
    n_cells: int = 0
    issues: List[str] = field(default_factory=list)
    status: str = "OK"

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for reporting."""
        return {
            "sample_id": self.sample_id,
            "n_cells": self.n_cells,
            "n_markers": len(self.markers),
            "status": self.status,
            "issues": ";".join(self.issues) if self.issues else "",
        }


class DataLoader:
    """Data loader with validation.

    Parameters
    ----------
    config : LoaderConfig
        Loader configuration

    Example
    -------
    >>> from celltype_refinery.core.preprocessing import DataLoader, LoaderConfig
    >>> config = LoaderConfig(cell_id_col="cell_mask_id")
    >>> loader = DataLoader(config)
    >>> result = loader.load_sample(matrix_path, metadata_path, "sample_01")
    """

    def __init__(self, config: Optional[LoaderConfig] = None):
        self.config = config or LoaderConfig()

    def load_cell_matrix(self, path: Path) -> pd.DataFrame:
        """Load cell-by-marker intensity matrix.

        Parameters
        ----------
        path : Path
            Path to CSV file

        Returns
        -------
        pd.DataFrame
            Cell matrix with cell_id_col + marker columns
        """
        df = pd.read_csv(path)
        cell_id_col = self.config.cell_id_col

        if cell_id_col not in df.columns:
            raise ValueError(f"Cell ID column '{cell_id_col}' not found in {path}")

        df[cell_id_col] = df[cell_id_col].astype(str)
        return df

    def load_cell_metadata(self, path: Path) -> pd.DataFrame:
        """Load cell metadata with spatial coordinates.

        Parameters
        ----------
        path : Path
            Path to CSV file

        Returns
        -------
        pd.DataFrame
            Cell metadata
        """
        df = pd.read_csv(path)
        cell_id_col = self.config.cell_id_col

        if cell_id_col not in df.columns:
            raise ValueError(f"Cell ID column '{cell_id_col}' not found in {path}")

        df[cell_id_col] = df[cell_id_col].astype(str)
        return df

    def load_metadata_registry(self, path: Path) -> pd.DataFrame:
        """Load sample metadata registry.

        Parameters
        ----------
        path : Path
            Path to CSV file

        Returns
        -------
        pd.DataFrame
            Metadata registry with sample paths
        """
        df = pd.read_csv(path)
        sample_id_col = self.config.sample_id_col

        if sample_id_col not in df.columns:
            raise ValueError(f"Sample ID column '{sample_id_col}' not found in {path}")

        df[sample_id_col] = df[sample_id_col].astype(str)

        # Validate required columns
        missing = []
        for col in self.config.required_metadata_cols:
            if col not in df.columns:
                missing.append(col)

        if missing:
            raise ValueError(f"Required columns missing from metadata: {missing}")

        return df

    def _find_coordinate_cols(
        self, df: pd.DataFrame
    ) -> Optional[Tuple[str, str]]:
        """Find coordinate columns in dataframe.

        Returns
        -------
        Optional[Tuple[str, str]]
            (x_col, y_col) or None if not found
        """
        # Try primary coordinate columns
        x_col, y_col = self.config.coordinate_cols
        if x_col in df.columns and y_col in df.columns:
            return (x_col, y_col)

        # Try alternatives
        for alt_x, alt_y in self.config.coordinate_col_alternatives:
            if alt_x in df.columns and alt_y in df.columns:
                return (alt_x, alt_y)

        return None

    def load_sample(
        self,
        matrix_path: Path,
        metadata_path: Path,
        sample_id: str,
        reference_markers: Optional[Set[str]] = None,
    ) -> LoadResult:
        """Load and validate a single sample.

        Parameters
        ----------
        matrix_path : Path
            Path to cell matrix CSV
        metadata_path : Path
            Path to cell metadata CSV
        sample_id : str
            Sample identifier
        reference_markers : Set[str], optional
            Expected marker names for consistency check

        Returns
        -------
        LoadResult
            Loading result with data and validation status
        """
        result = LoadResult(sample_id=sample_id)
        cell_id_col = self.config.cell_id_col

        # Load matrix
        try:
            result.cell_matrix = self.load_cell_matrix(matrix_path)
            marker_cols = [
                c for c in result.cell_matrix.columns if c != cell_id_col
            ]
            result.markers = set(marker_cols)
            result.n_cells = len(result.cell_matrix)

            # Check for duplicates
            n_dups = result.cell_matrix[cell_id_col].duplicated().sum()
            if n_dups > 0:
                result.issues.append(f"duplicate_ids_matrix:{n_dups}")

        except FileNotFoundError:
            result.issues.append("cell_matrix_missing")
        except Exception as e:
            result.issues.append(f"cell_matrix_error:{e}")

        # Load metadata
        try:
            result.cell_metadata = self.load_cell_metadata(metadata_path)

            # Check for duplicates
            n_dups = result.cell_metadata[cell_id_col].duplicated().sum()
            if n_dups > 0:
                result.issues.append(f"duplicate_ids_metadata:{n_dups}")

            # Check coordinates
            coord_cols = self._find_coordinate_cols(result.cell_metadata)
            if coord_cols:
                x_col, y_col = coord_cols
                for axis, col in [("x", x_col), ("y", y_col)]:
                    min_val = pd.to_numeric(
                        result.cell_metadata[col], errors="coerce"
                    ).min()
                    if pd.notna(min_val) and float(min_val) < 0:
                        result.issues.append(f"{axis}_negative_min={min_val}")

        except FileNotFoundError:
            result.issues.append("cell_metadata_missing")
        except Exception as e:
            result.issues.append(f"cell_metadata_error:{e}")

        # Cross-validate
        if result.cell_matrix is not None and result.cell_metadata is not None:
            matrix_ids = set(result.cell_matrix[cell_id_col])
            metadata_ids = set(result.cell_metadata[cell_id_col])

            missing_in_matrix = len(metadata_ids - matrix_ids)
            missing_in_metadata = len(matrix_ids - metadata_ids)

            if missing_in_matrix > 0:
                result.issues.append(f"missing_in_matrix:{missing_in_matrix}")
            if missing_in_metadata > 0:
                result.issues.append(f"missing_in_metadata:{missing_in_metadata}")

            if len(result.cell_matrix) != len(result.cell_metadata):
                result.issues.append("row_mismatch")

        # Check marker consistency
        if reference_markers is not None and result.markers:
            missing = sorted(reference_markers - result.markers)
            extra = sorted(result.markers - reference_markers)
            if missing:
                result.issues.append("missing_markers:" + ",".join(missing))
            if extra:
                result.issues.append("extra_markers:" + ",".join(extra))

        result.status = "OK" if not result.issues else "CHECK"
        return result

    def align_to_common_ids(
        self,
        cell_matrix: pd.DataFrame,
        cell_metadata: pd.DataFrame,
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Align matrix and metadata to common cell IDs.

        Parameters
        ----------
        cell_matrix : pd.DataFrame
            Cell matrix
        cell_metadata : pd.DataFrame
            Cell metadata

        Returns
        -------
        Tuple[pd.DataFrame, pd.DataFrame]
            Aligned (matrix, metadata)
        """
        cell_id_col = self.config.cell_id_col

        common_ids = sorted(
            set(cell_matrix[cell_id_col]) & set(cell_metadata[cell_id_col])
        )

        matrix_aligned = cell_matrix[
            cell_matrix[cell_id_col].isin(common_ids)
        ].copy()
        metadata_aligned = cell_metadata[
            cell_metadata[cell_id_col].isin(common_ids)
        ].copy()

        return matrix_aligned, metadata_aligned

    def get_marker_columns(self, df: pd.DataFrame) -> List[str]:
        """Get marker column names from dataframe.

        Parameters
        ----------
        df : pd.DataFrame
            Cell matrix

        Returns
        -------
        List[str]
            Marker column names
        """
        cell_id_col = self.config.cell_id_col
        return [c for c in df.columns if c != cell_id_col]
