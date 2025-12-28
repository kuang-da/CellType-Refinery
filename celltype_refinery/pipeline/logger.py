"""Structured logging for pipeline execution."""

import logging
import sys
from datetime import datetime
from pathlib import Path
from typing import Optional


class ColoredFormatter(logging.Formatter):
    """Custom formatter with color support for console output."""

    def __init__(self, fmt: str, datefmt: str, colors: dict):
        super().__init__(fmt=fmt, datefmt=datefmt)
        self.colors = colors

    def format(self, record):
        levelname = record.levelname
        color = self.colors.get(levelname, self.colors["RESET"])
        reset = self.colors["RESET"]
        record.levelname = f"{color}{levelname}{reset}"
        return super().format(record)


class PipelineLogger:
    """Structured logging for pipeline execution.

    Provides file logging (detailed, persistent) and console logging
    (colored, user-friendly) with structured event logging for
    stage start/complete/error.

    Parameters
    ----------
    log_dir : str
        Directory for log files
    log_level : str
        Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
    log_name : str, optional
        Logger name. Default: "celltype_refinery"

    Attributes
    ----------
    log_dir : Path
        Directory for log files
    log_file : Path
        Path to the main pipeline log file
    logger : logging.Logger
        Python logger instance

    Example
    -------
    >>> logger = PipelineLogger("logs/", log_level="INFO")
    >>> logger.setup()
    >>> logger.log_stage_start("A", "Preprocessing")
    >>> logger.log_stage_complete("A", 45.2)
    """

    COLORS = {
        "DEBUG": "\033[0;36m",  # Cyan
        "INFO": "\033[0;34m",  # Blue
        "WARNING": "\033[1;33m",  # Yellow
        "ERROR": "\033[0;31m",  # Red
        "CRITICAL": "\033[1;31m",  # Bold Red
        "SUCCESS": "\033[0;32m",  # Green
        "RESET": "\033[0m",  # Reset
    }

    def __init__(
        self,
        log_dir: str,
        log_level: str = "INFO",
        log_name: str = "celltype_refinery",
    ):
        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.log_file = self.log_dir / f"pipeline_{timestamp}.log"

        self.log_level = getattr(logging, log_level.upper())
        self.logger = logging.getLogger(log_name)
        self.logger.setLevel(self.log_level)
        self.logger.handlers = []

    def setup(self) -> None:
        """Configure logging handlers.

        Sets up file handler (detailed logs) and console handler
        (colored output).
        """
        file_handler = logging.FileHandler(self.log_file, mode="w")
        file_handler.setLevel(self.log_level)
        file_handler.setFormatter(self._get_file_formatter())

        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(self.log_level)
        console_handler.setFormatter(self._get_console_formatter())

        self.logger.addHandler(file_handler)
        self.logger.addHandler(console_handler)

    def _get_file_formatter(self) -> logging.Formatter:
        """Get formatter for file logging (detailed, no colors)."""
        return logging.Formatter(
            fmt="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )

    def _get_console_formatter(self) -> logging.Formatter:
        """Get formatter for console logging (colored, concise)."""
        return ColoredFormatter(
            fmt="%(asctime)s - %(levelname)s - %(message)s",
            datefmt="%H:%M:%S",
            colors=self.COLORS,
        )

    def log_stage_start(self, stage_id: str, stage_name: str) -> None:
        """Log the start of a pipeline stage.

        Parameters
        ----------
        stage_id : str
            Stage identifier (e.g., "A", "B")
        stage_name : str
            Human-readable stage name
        """
        separator = "=" * 80
        self.logger.info(separator)
        self.logger.info(f"Starting Stage {stage_id}: {stage_name}")
        self.logger.info(separator)

    def log_stage_complete(self, stage_id: str, duration: float) -> None:
        """Log successful completion of a stage.

        Parameters
        ----------
        stage_id : str
            Stage identifier
        duration : float
            Execution time in seconds
        """
        duration_str = self.format_duration(duration)
        self.logger.info(f"Stage {stage_id} completed successfully in {duration_str}")

    def log_stage_error(self, stage_id: str, error: str) -> None:
        """Log a stage error.

        Parameters
        ----------
        stage_id : str
            Stage identifier
        error : str
            Error message or exception
        """
        self.logger.error(f"Stage {stage_id} failed: {error}")

    def log_info(self, message: str) -> None:
        """Log an info message."""
        self.logger.info(message)

    def log_warning(self, message: str) -> None:
        """Log a warning message."""
        self.logger.warning(message)

    def log_error(self, message: str) -> None:
        """Log an error message."""
        self.logger.error(message)

    def log_debug(self, message: str) -> None:
        """Log a debug message."""
        self.logger.debug(message)

    @staticmethod
    def format_duration(seconds: float) -> str:
        """Format duration in seconds to human-readable string.

        Parameters
        ----------
        seconds : float
            Duration in seconds

        Returns
        -------
        str
            Formatted string (e.g., "45.2s", "1m 23s", "2h 15m")
        """
        if seconds < 60:
            return f"{seconds:.1f}s"
        elif seconds < 3600:
            mins = int(seconds // 60)
            secs = int(seconds % 60)
            return f"{mins}m {secs}s"
        else:
            hours = int(seconds // 3600)
            mins = int((seconds % 3600) // 60)
            return f"{hours}h {mins}m"
