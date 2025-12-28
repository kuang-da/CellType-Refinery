"""CSV I/O utilities for CellType-Refinery.

Provides functions for loading cell matrices, metadata, and marker maps.
Tissue-specific column validation is configurable rather than hardcoded.
"""

from __future__ import annotations

import json
import logging
import re
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Union

import pandas as pd

try:
    import requests
    HAS_REQUESTS = True
except ImportError:
    HAS_REQUESTS = False

logger = logging.getLogger(__name__)

PathLike = Union[str, Path]

# Default required columns (can be overridden via config)
DEFAULT_REQUIRED_METADATA_COLUMNS = [
    "sample_id",
    "cell_matrix_path",
]


def ensure_output_dir(path: PathLike) -> Path:
    """Create the directory at path if it does not exist and return it.

    Parameters
    ----------
    path : PathLike
        Directory path to create.

    Returns
    -------
    Path
        The created/existing directory path.
    """
    directory = Path(path)
    directory.mkdir(parents=True, exist_ok=True)
    return directory


def _validate_metadata_columns(
    df: pd.DataFrame,
    required_columns: Optional[List[str]] = None,
) -> None:
    """Validate that required columns are present in metadata DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        Metadata DataFrame to validate.
    required_columns : List[str], optional
        List of required column names. If None, uses DEFAULT_REQUIRED_METADATA_COLUMNS.

    Raises
    ------
    ValueError
        If any required columns are missing.
    """
    if required_columns is None:
        required_columns = DEFAULT_REQUIRED_METADATA_COLUMNS
    missing = [col for col in required_columns if col not in df.columns]
    if missing:
        raise ValueError(f"Metadata table missing columns: {missing}")


def _resolve_path(value: str, base: Path) -> Path:
    """Resolve a path relative to a base directory."""
    candidate = Path(value)
    if not candidate.is_absolute():
        candidate = (base / candidate).resolve()
    return candidate


def load_metadata_registry(
    path: PathLike,
    required_columns: Optional[List[str]] = None,
    path_columns: Optional[List[str]] = None,
) -> pd.DataFrame:
    """Load metadata CSV and ensure required schema.

    Parameters
    ----------
    path : PathLike
        Path to metadata CSV file.
    required_columns : List[str], optional
        List of required column names for validation.
    path_columns : List[str], optional
        List of columns containing file paths to resolve relative to CSV location.
        Defaults to ["cell_matrix_path", "cell_metadata_path"] if present.

    Returns
    -------
    pd.DataFrame
        Loaded and validated metadata DataFrame.

    Raises
    ------
    FileNotFoundError
        If the metadata file does not exist.
    ValueError
        If required columns are missing.
    """
    csv_path = Path(path)
    if not csv_path.exists():
        raise FileNotFoundError(f"Metadata registry not found: {csv_path}")
    df = pd.read_csv(csv_path)
    _validate_metadata_columns(df, required_columns)

    # Resolve relative paths
    base_dir = csv_path.parent
    if path_columns is None:
        path_columns = ["cell_matrix_path", "cell_metadata_path", "image_path"]
    for col in path_columns:
        if col in df.columns:
            df[col] = df[col].apply(lambda p: str(_resolve_path(p, base_dir)))

    return df


def _validate_cell_table(df: pd.DataFrame, path: PathLike, id_column: str = "cell_mask_id") -> None:
    """Validate a cell table has the expected structure.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to validate.
    path : PathLike
        Path for error messages.
    id_column : str
        Expected name of the cell ID column.

    Raises
    ------
    ValueError
        If table is empty or missing expected ID column.
    """
    if df.empty:
        raise ValueError(f"Cell table {path} is empty")
    first_col = df.columns[0]
    if first_col != id_column:
        raise ValueError(
            f"Expected first column `{id_column}` in {path}, found `{first_col}`"
        )


def load_cell_matrix(
    path: PathLike,
    id_column: str = "cell_mask_id",
    drop_columns: Optional[List[str]] = None,
) -> pd.DataFrame:
    """Read a cell-by-marker matrix and apply light cleanup.

    Parameters
    ----------
    path : PathLike
        Path to cell matrix CSV file.
    id_column : str
        Name of the cell ID column (default: "cell_mask_id").
    drop_columns : List[str], optional
        Columns to drop from the matrix. Defaults to ["patch_id", "global_cell_id"].

    Returns
    -------
    pd.DataFrame
        Cell-by-marker matrix with ID column as string type.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    ValueError
        If the table is invalid.
    """
    csv_path = Path(path)
    if not csv_path.exists():
        raise FileNotFoundError(f"Cell matrix not found: {csv_path}")
    df = pd.read_csv(csv_path)
    _validate_cell_table(df, csv_path, id_column)

    if drop_columns is None:
        drop_columns = ["patch_id", "global_cell_id"]
    cols_to_drop = [c for c in drop_columns if c in df.columns]
    if cols_to_drop:
        df = df.drop(columns=cols_to_drop)

    df[id_column] = df[id_column].astype(str)
    return df


def load_cell_metadata(
    path: PathLike,
    id_column: str = "cell_mask_id",
) -> pd.DataFrame:
    """Read a cell-level metadata table.

    Parameters
    ----------
    path : PathLike
        Path to cell metadata CSV file.
    id_column : str
        Name of the cell ID column (default: "cell_mask_id").

    Returns
    -------
    pd.DataFrame
        Cell metadata DataFrame with ID column as string type.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    ValueError
        If the table is invalid.
    """
    csv_path = Path(path)
    if not csv_path.exists():
        raise FileNotFoundError(f"Cell metadata table not found: {csv_path}")
    df = pd.read_csv(csv_path)
    _validate_cell_table(df, csv_path, id_column)
    df[id_column] = df[id_column].astype(str)
    return df


def write_dataframe(df: pd.DataFrame, path: PathLike, *, index: bool = False) -> Path:
    """Write DataFrame to path ensuring the parent directory exists.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to write.
    path : PathLike
        Output path.
    index : bool
        Whether to write row index (default: False).

    Returns
    -------
    Path
        The output path.
    """
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=index)
    return output_path


def unique_marker_columns(df: pd.DataFrame, id_column: str = "cell_mask_id") -> Iterable[str]:
    """Return the marker/antibody columns from a cell-by-marker matrix.

    Parameters
    ----------
    df : pd.DataFrame
        Cell-by-marker matrix.
    id_column : str
        Name of the cell ID column to exclude.

    Returns
    -------
    Iterable[str]
        List of marker column names.
    """
    return [col for col in df.columns if col != id_column]


# ============================================================================
# Remote Data Loading (GitHub Gist)
# ============================================================================

GIST_URL_PATTERN = re.compile(r"^gist:([a-fA-F0-9]+)/(.+)$")
GITHUB_GIST_API = "https://api.github.com/gists/{gist_id}"
CACHE_METADATA_KEY = "_gist_cache_metadata"


def _fetch_gist_file(gist_id: str, filename: str) -> Dict[str, Any]:
    """Fetch a file from a GitHub gist.

    Parameters
    ----------
    gist_id : str
        The gist ID (e.g., "46154004f55e5743e43bda305c8f287e")
    filename : str
        The file name within the gist.

    Returns
    -------
    Dict[str, Any]
        Parsed JSON content as a dictionary.

    Raises
    ------
    ImportError
        If requests library is not installed.
    ValueError
        If the file is not found in the gist.
    """
    if not HAS_REQUESTS:
        raise ImportError(
            "The 'requests' library is required to fetch data from GitHub gists. "
            "Install it with: pip install requests"
        )

    api_url = GITHUB_GIST_API.format(gist_id=gist_id)
    logger.info("Fetching gist metadata from: %s", api_url)

    response = requests.get(api_url, timeout=30)
    response.raise_for_status()
    gist_data = response.json()

    files = gist_data.get("files", {})
    if filename not in files:
        available = list(files.keys())
        raise ValueError(
            f"File '{filename}' not found in gist {gist_id}. "
            f"Available files: {available}"
        )

    file_info = files[filename]
    raw_url = file_info.get("raw_url")
    if not raw_url:
        raise ValueError(f"No raw_url found for file '{filename}' in gist {gist_id}")

    logger.info("Fetching file content from: %s", raw_url)
    content_response = requests.get(raw_url, timeout=30)
    content_response.raise_for_status()

    return json.loads(content_response.text)


def _get_cache_path(gist_id: str, filename: str, cache_dir: Path) -> Path:
    """Generate a cache file path for a gist file."""
    safe_filename = re.sub(r"[^\w\-.]", "_", filename)
    cache_filename = f"gist_{gist_id[:8]}_{safe_filename}"
    return cache_dir / cache_filename


def _read_cached_marker_map(cache_path: Path) -> Optional[Dict[str, Any]]:
    """Read a cached marker map if it exists."""
    if not cache_path.exists():
        return None

    try:
        with open(cache_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        logger.info("Loaded cached marker map from: %s", cache_path)
        return data
    except (json.JSONDecodeError, IOError) as e:
        logger.warning("Failed to read cache file %s: %s", cache_path, e)
        return None


def _write_cached_marker_map(
    data: Dict[str, Any],
    cache_path: Path,
    gist_id: str,
    filename: str,
) -> None:
    """Write marker map to cache with metadata."""
    cached_data = data.copy()
    cached_data[CACHE_METADATA_KEY] = {
        "gist_id": gist_id,
        "filename": filename,
        "cached_at": datetime.utcnow().isoformat(),
        "source": f"gist:{gist_id}/{filename}",
    }

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with open(cache_path, "w", encoding="utf-8") as f:
        json.dump(cached_data, f, indent=2)
    logger.info("Cached marker map to: %s", cache_path)


def load_marker_map_from_source(
    source: str,
    cache_dir: Optional[PathLike] = None,
    refresh: bool = False,
) -> Dict[str, Any]:
    """Load a marker map from a local file or GitHub gist.

    Parameters
    ----------
    source : str
        Either a local file path OR a gist reference in the format
        "gist:<gist_id>/<filename>".
    cache_dir : PathLike, optional
        Directory to cache gist files. If None, uses current working directory.
    refresh : bool
        If True, re-download from gist even if a cached version exists.

    Returns
    -------
    Dict[str, Any]
        The marker map as a dictionary.

    Raises
    ------
    FileNotFoundError
        If the local file does not exist.
    ValueError
        If the gist reference is malformed or file not found.

    Examples
    --------
    >>> # Load from local file
    >>> marker_map = load_marker_map_from_source("markers.json")
    >>>
    >>> # Load from gist (caches locally)
    >>> marker_map = load_marker_map_from_source(
    ...     "gist:46154004f55e5743e43bda305c8f287e/markers.json"
    ... )
    """
    gist_match = GIST_URL_PATTERN.match(source)

    if gist_match:
        gist_id = gist_match.group(1)
        filename = gist_match.group(2)

        if cache_dir is None:
            cache_dir = Path(".")
        else:
            cache_dir = Path(cache_dir)

        cache_path = _get_cache_path(gist_id, filename, cache_dir)

        if not refresh:
            cached_data = _read_cached_marker_map(cache_path)
            if cached_data is not None:
                cached_data.pop(CACHE_METADATA_KEY, None)
                return cached_data

        logger.info("Fetching marker map from gist: %s/%s", gist_id, filename)
        data = _fetch_gist_file(gist_id, filename)
        _write_cached_marker_map(data, cache_path, gist_id, filename)

        return data

    else:
        local_path = Path(source)
        if not local_path.exists():
            raise FileNotFoundError(f"Marker map file not found: {local_path}")

        logger.info("Loading marker map from local file: %s", local_path)
        with open(local_path, "r", encoding="utf-8") as f:
            return json.load(f)
