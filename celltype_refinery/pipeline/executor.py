"""Pipeline execution engine with checkpoint support."""

import json
import subprocess
import time
from datetime import datetime
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional

from .config import PipelineConfig
from .logger import PipelineLogger
from .stage import Stage


class PipelineExecutor:
    """Executes pipeline stages with validation, checkpointing, and error handling.

    Provides sequential stage execution with input/output validation,
    checkpoint/resume support, and detailed logging.

    Parameters
    ----------
    config : PipelineConfig
        Loaded PipelineConfig instance
    logger : PipelineLogger
        Initialized PipelineLogger instance
    state_file : str, optional
        Path to checkpoint state file

    Attributes
    ----------
    config : PipelineConfig
        Pipeline configuration
    logger : PipelineLogger
        Logger instance
    state_file : Path
        Path to checkpoint state file
    completed_stages : List[str]
        List of successfully completed stage IDs

    Example
    -------
    >>> config = PipelineConfig("pipeline.yaml")
    >>> config.load()
    >>> config.parse_stages()
    >>> logger = PipelineLogger("logs/")
    >>> logger.setup()
    >>> executor = PipelineExecutor(config, logger)
    >>> exit_code = executor.run()
    """

    def __init__(
        self,
        config: PipelineConfig,
        logger: PipelineLogger,
        state_file: Optional[str] = None,
    ):
        self.config = config
        self.logger = logger
        self.state_file = (
            Path(state_file) if state_file else Path(".pipeline_state.json")
        )
        self.completed_stages: List[str] = []

    def load_state(self) -> None:
        """Load checkpoint state from previous run.

        Reads the state file and loads the list of completed stages.
        If file doesn't exist, starts with empty state.
        """
        if not self.state_file.exists():
            self.logger.log_debug("No checkpoint file found, starting fresh")
            return

        try:
            with open(self.state_file, "r") as f:
                state = json.load(f)

            self.completed_stages = state.get("completed_stages", [])
            self.logger.log_info(
                f"Loaded checkpoint: {len(self.completed_stages)} stages completed"
            )

            if self.completed_stages:
                self.logger.log_info(f"Last completed: {self.completed_stages[-1]}")
        except Exception as e:
            self.logger.log_warning(f"Failed to load checkpoint: {e}")
            self.completed_stages = []

    def save_state(self) -> None:
        """Save current state to checkpoint file.

        Writes completed stages and metadata to JSON file.
        """
        state = {
            "completed_stages": self.completed_stages,
            "timestamp": datetime.now().isoformat(),
            "pipeline_version": self.config.raw_config.get("pipeline", {}).get(
                "version", "1.0"
            ),
        }

        self.state_file.parent.mkdir(parents=True, exist_ok=True)

        with open(self.state_file, "w") as f:
            json.dump(state, f, indent=2)

    def clear_state(self) -> None:
        """Clear checkpoint state (for fresh run)."""
        if self.state_file.exists():
            self.state_file.unlink()
            self.logger.log_info("Cleared checkpoint state")

        self.completed_stages = []

    def should_skip_stage(self, stage: Stage) -> bool:
        """Check if a stage should be skipped.

        For optional stages, checks if required conditions are met.

        Parameters
        ----------
        stage : Stage
            Stage object to check

        Returns
        -------
        bool
            True if stage should be skipped
        """
        if not stage.optional:
            return False

        # Check if config file exists for config-driven stages
        config_arg = stage.args.get("config")
        if config_arg:
            config_path = Path(config_arg)
            if not config_path.exists():
                self.logger.log_info(
                    f"Stage {stage.stage_id} config not found - skipping"
                )
                return True

        return False

    def execute_stage(self, stage: Stage, dry_run: bool = False) -> int:
        """Execute a single pipeline stage.

        Parameters
        ----------
        stage : Stage
            Stage object to execute
        dry_run : bool
            If True, only show what would be executed

        Returns
        -------
        int
            Exit code (0 = success, non-zero = failure)
        """
        if self.should_skip_stage(stage):
            self.logger.log_info(
                f"[SKIP] Stage {stage.stage_id} is optional and conditions not met"
            )
            self.completed_stages.append(stage.stage_id)
            self.save_state()
            return 0

        cmd = stage.get_command()

        if dry_run:
            self.logger.log_info(f"[DRY RUN] Would execute: {' '.join(cmd)}")
            return 0

        valid, errors = stage.validate_inputs()
        if not valid:
            self.logger.log_error(
                f"Input validation failed for stage {stage.stage_id}:"
            )
            for error in errors:
                self.logger.log_error(f"  - {error}")
            return 1

        self.logger.log_stage_start(stage.stage_id, stage.name)
        start_time = time.time()

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=False,
            )

            duration = time.time() - start_time

            if result.returncode != 0:
                self.logger.log_stage_error(
                    stage.stage_id, f"Exit code {result.returncode}"
                )
                self.logger.log_error(f"STDERR: {result.stderr[-1000:]}")
                return result.returncode

            valid, errors = stage.validate_outputs()
            if not valid:
                self.logger.log_error(
                    f"Output validation failed for stage {stage.stage_id}:"
                )
                for error in errors:
                    self.logger.log_error(f"  - {error}")
                return 1

            self.logger.log_stage_complete(stage.stage_id, duration)
            self.completed_stages.append(stage.stage_id)
            self.save_state()

            return 0

        except Exception as e:
            self.logger.log_stage_error(stage.stage_id, str(e))
            return 1

    def run(
        self,
        start_stage: Optional[str] = None,
        end_stage: Optional[str] = None,
        dry_run: bool = False,
        force: bool = False,
    ) -> int:
        """Execute pipeline from start_stage to end_stage.

        Parameters
        ----------
        start_stage : str, optional
            Stage ID to start from (default: first stage)
        end_stage : str, optional
            Stage ID to end at (default: last stage)
        dry_run : bool
            If True, show execution plan without running
        force : bool
            If True, ignore checkpoint and re-run all stages

        Returns
        -------
        int
            Exit code (0 = success, non-zero = failure)
        """
        order = self.config.get_execution_order()

        if start_stage:
            if start_stage not in order:
                self.logger.log_error(f"Start stage '{start_stage}' not found")
                return 1
            start_idx = order.index(start_stage)
            order = order[start_idx:]

        if end_stage:
            if end_stage not in order:
                self.logger.log_error(f"End stage '{end_stage}' not found")
                return 1
            end_idx = order.index(end_stage) + 1
            order = order[:end_idx]

        if force:
            self.clear_state()
        else:
            self.load_state()

        stage_names = " -> ".join(order)
        self.logger.log_info(f"Pipeline execution plan: {stage_names}")

        if dry_run:
            self.logger.log_info("DRY RUN MODE - No stages will be executed")

        for stage_id in order:
            if stage_id in self.completed_stages and not force:
                self.logger.log_info(f"[SKIP] Stage {stage_id} already completed")
                continue

            stage = self.config.stages[stage_id]
            exit_code = self.execute_stage(stage, dry_run)

            if exit_code != 0:
                self.logger.log_error(f"Pipeline failed at stage {stage_id}")
                return exit_code

        self.logger.log_info("Pipeline completed successfully")
        return 0

    def get_resume_stage(self) -> Optional[str]:
        """Get the next stage to resume from based on checkpoint.

        Returns
        -------
        Optional[str]
            Stage ID to resume from, or None if starting fresh or all complete
        """
        self.load_state()

        if not self.completed_stages:
            return None

        order = self.config.get_execution_order()
        last_completed = self.completed_stages[-1]

        if last_completed not in order:
            return None

        last_idx = order.index(last_completed)
        if last_idx + 1 < len(order):
            return order[last_idx + 1]

        return None  # All stages completed


class InMemoryExecutor:
    """Pipeline executor that runs stages in-memory without subprocess.

    Useful for testing and when all stages are Python functions.

    Parameters
    ----------
    logger : PipelineLogger, optional
        Logger instance

    Example
    -------
    >>> executor = InMemoryExecutor()
    >>> executor.register_stage("preprocess", preprocess_func)
    >>> executor.register_stage("cluster", cluster_func, depends_on=["preprocess"])
    >>> result = executor.run()
    """

    def __init__(self, logger: Optional[PipelineLogger] = None):
        self.logger = logger
        self.stages: Dict[str, Dict[str, Any]] = {}
        self.completed_stages: List[str] = []

    def register_stage(
        self,
        stage_id: str,
        func: Callable,
        depends_on: Optional[List[str]] = None,
        name: Optional[str] = None,
    ) -> None:
        """Register a stage function.

        Parameters
        ----------
        stage_id : str
            Stage identifier
        func : Callable
            Stage function to execute
        depends_on : List[str], optional
            List of stage IDs this stage depends on
        name : str, optional
            Human-readable stage name
        """
        self.stages[stage_id] = {
            "func": func,
            "depends_on": depends_on or [],
            "name": name or stage_id,
        }

    def _get_execution_order(self) -> List[str]:
        """Compute stage execution order via topological sort."""
        in_degree = {stage_id: 0 for stage_id in self.stages}

        for stage_id, stage in self.stages.items():
            for dep in stage["depends_on"]:
                in_degree[stage_id] += 1

        from collections import deque

        queue = deque([sid for sid, degree in in_degree.items() if degree == 0])
        order = []

        while queue:
            stage_id = queue.popleft()
            order.append(stage_id)

            for other_id, other_stage in self.stages.items():
                if stage_id in other_stage["depends_on"]:
                    in_degree[other_id] -= 1
                    if in_degree[other_id] == 0:
                        queue.append(other_id)

        if len(order) != len(self.stages):
            raise ValueError("Circular dependency detected")

        return order

    def run(self, **kwargs) -> Dict[str, Any]:
        """Execute all registered stages in order.

        Parameters
        ----------
        **kwargs
            Arguments passed to each stage function

        Returns
        -------
        Dict[str, Any]
            Map of stage_id to stage result
        """
        order = self._get_execution_order()
        results: Dict[str, Any] = {}

        for stage_id in order:
            stage = self.stages[stage_id]
            if self.logger:
                self.logger.log_stage_start(stage_id, stage["name"])

            start_time = time.time()
            try:
                result = stage["func"](**kwargs, stage_results=results)
                results[stage_id] = result
                self.completed_stages.append(stage_id)

                if self.logger:
                    duration = time.time() - start_time
                    self.logger.log_stage_complete(stage_id, duration)

            except Exception as e:
                if self.logger:
                    self.logger.log_stage_error(stage_id, str(e))
                raise

        return results
