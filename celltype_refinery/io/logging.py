"""Logging utilities for CellType-Refinery.

Provides timestamped file logging and structured log output (JSON, YAML).
"""

from __future__ import annotations

import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Any, Tuple, Union

import yaml

PathLike = Union[str, Path]


def get_timestamped_log_path(log_path: PathLike) -> Path:
    """Generate a timestamped log path from the base log path.

    Example: stage_h.log -> stage_h_20251209_080530.log

    Parameters
    ----------
    log_path : PathLike
        Base log file path.

    Returns
    -------
    Path
        Timestamped log path.
    """
    log_path = Path(log_path)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    stem = log_path.stem
    suffix = log_path.suffix or ".log"
    return log_path.parent / f"{stem}_{timestamp}{suffix}"


def get_logger(
    name: str,
    log_path: PathLike,
    level: int = logging.INFO,
    timestamped: bool = True,
) -> Tuple[logging.Logger, Path]:
    """Return a file logger configured for the given module.

    Parameters
    ----------
    name : str
        Logger name (typically module/stage name).
    log_path : PathLike
        Base path for log file.
    level : int
        Logging level (default: INFO).
    timestamped : bool
        If True, add timestamp to filename to preserve previous logs.
        If False, overwrite existing log file.

    Returns
    -------
    Tuple[logging.Logger, Path]
        Tuple of (logger, actual_log_path) where actual_log_path is the
        timestamped path if timestamped=True.
    """
    log_path = Path(log_path)

    if timestamped:
        actual_log_path = get_timestamped_log_path(log_path)
    else:
        actual_log_path = log_path
        actual_log_path.unlink(missing_ok=True)

    actual_log_path.parent.mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.propagate = False
    for handler in list(logger.handlers):
        logger.removeHandler(handler)
    handler = logging.FileHandler(actual_log_path, mode="a", encoding="utf-8")
    formatter = logging.Formatter(
        fmt="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return logger, actual_log_path


def _prepare_log_destination(log_path: PathLike) -> Path:
    """Ensure log destination directory exists."""
    path = Path(log_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def log_json(log_path: PathLike, record: dict[str, Any]) -> None:
    """Append a JSON line to log_path.

    Parameters
    ----------
    log_path : PathLike
        Path to log file.
    record : dict
        Dictionary to serialize as JSON.
    """
    path = _prepare_log_destination(log_path)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(record, default=str))
        handle.write("\n")


def log_yaml(
    log_path: PathLike,
    record: dict[str, Any],
    *,
    logger: logging.Logger | None = None,
) -> None:
    """Append a YAML document to log_path.

    Parameters
    ----------
    log_path : PathLike
        Path to log file.
    record : dict
        Dictionary to serialize as YAML.
    logger : logging.Logger, optional
        If provided, log to this logger instead of file.
    """
    yaml_text = yaml.safe_dump(record, sort_keys=False).rstrip("\n")
    message = f"{yaml_text}\n---"
    if logger is not None:
        logger.info("%s", message)
        return

    path = _prepare_log_destination(log_path)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(message)
        handle.write("\n")
