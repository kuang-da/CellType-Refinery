"""Pipeline configuration loader and validator."""

import re
from collections import deque
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import yaml

from .stage import Stage


class PipelineConfig:
    """Loads and manages pipeline configuration from YAML files.

    Provides YAML-based configuration loading with template resolution,
    stage parsing, dependency validation, and execution order computation.

    Parameters
    ----------
    config_path : str
        Path to the YAML configuration file

    Attributes
    ----------
    config_path : Path
        Path to the configuration file
    raw_config : Dict[str, Any]
        Raw configuration dictionary loaded from YAML
    stages : Dict[str, Stage]
        Dictionary mapping stage_id to Stage objects
    global_settings : Dict[str, Any]
        Global configuration parameters

    Example
    -------
    >>> config = PipelineConfig("pipeline.yaml")
    >>> config.load()
    >>> config.parse_stages()
    >>> valid, errors = config.validate_dependencies()
    >>> order = config.get_execution_order()
    """

    def __init__(self, config_path: str):
        self.config_path = Path(config_path)
        self.raw_config: Dict[str, Any] = {}
        self.stages: Dict[str, Stage] = {}
        self.global_settings: Dict[str, Any] = {}

    def load(self) -> None:
        """Load YAML configuration file.

        Raises
        ------
        FileNotFoundError
            If config file doesn't exist
        yaml.YAMLError
            If YAML is malformed
        """
        if not self.config_path.exists():
            raise FileNotFoundError(f"Config file not found: {self.config_path}")

        with open(self.config_path, "r") as f:
            self.raw_config = yaml.safe_load(f)

        self.global_settings = {
            "pipeline": self.raw_config.get("pipeline", {}),
            "global": self.raw_config.get("global", {}),
        }

    def parse_stages(self) -> None:
        """Convert YAML stage definitions to Stage objects.

        Iterates through raw_config['stages'], resolves path templates,
        and creates Stage objects.

        Raises
        ------
        KeyError
            If required stage fields are missing
        """
        if "stages" not in self.raw_config:
            raise KeyError("No 'stages' section in configuration")

        for stage_id, stage_def in self.raw_config["stages"].items():
            if "script_module" not in stage_def:
                raise KeyError(
                    f"Stage '{stage_id}' missing required field 'script_module'"
                )

            resolved_inputs = {
                k: self.resolve_paths(v) for k, v in stage_def.get("inputs", {}).items()
            }
            resolved_outputs = {
                k: self.resolve_paths(v)
                for k, v in stage_def.get("outputs", {}).items()
            }
            resolved_args = {
                k: self.resolve_paths(v) if isinstance(v, str) else v
                for k, v in stage_def.get("args", {}).items()
            }
            resolved_required_files = [
                self.resolve_paths(path) for path in stage_def.get("required_files", [])
            ]

            stage = Stage(
                name=stage_def.get("name", stage_id),
                stage_id=stage_id,
                script_module=stage_def["script_module"],
                depends_on=stage_def.get("depends_on", []),
                inputs=resolved_inputs,
                outputs=resolved_outputs,
                args=resolved_args,
                required_files=resolved_required_files,
                optional=stage_def.get("optional", False),
            )

            self.stages[stage_id] = stage

    def resolve_paths(self, path_template: str) -> str:
        """Resolve path templates like {stages.A.outputs.metadata}.

        Templates can reference:
        - {global.param_name} - global parameters
        - {pipeline.param_name} - pipeline settings
        - {stages.STAGE_ID.outputs.output_name} - other stage outputs
        - {stages.STAGE_ID.inputs.input_name} - other stage inputs

        Parameters
        ----------
        path_template : str
            Path possibly containing {...} templates

        Returns
        -------
        str
            Resolved path with templates replaced by actual values
        """
        if "{" not in path_template:
            return path_template

        pattern = r"\{([^}]+)\}"

        def replace_template(match):
            ref = match.group(1)
            parts = ref.split(".")

            value = self.raw_config
            for part in parts:
                if isinstance(value, dict):
                    value = value.get(part, {})
                else:
                    return match.group(0)

                if isinstance(value, str):
                    break

            return str(value) if value and not isinstance(value, dict) else match.group(0)

        resolved = re.sub(pattern, replace_template, path_template)

        if resolved != path_template and "{" in resolved:
            return self.resolve_paths(resolved)

        return resolved

    def validate_dependencies(self) -> Tuple[bool, List[str]]:
        """Validate stage dependencies.

        Checks that all dependency stage_ids exist and there are no
        circular dependencies.

        Returns
        -------
        Tuple[bool, List[str]]
            (valid, errors) where valid is True if all dependencies are valid
        """
        errors = []

        for stage_id, stage in self.stages.items():
            for dep in stage.depends_on:
                if dep not in self.stages:
                    errors.append(
                        f"Stage '{stage_id}' depends on unknown stage '{dep}'"
                    )

        if len(errors) == 0:

            def has_cycle():
                visited = set()
                rec_stack = set()

                def visit(node):
                    if node not in self.stages:
                        return False

                    visited.add(node)
                    rec_stack.add(node)

                    for dep in self.stages[node].depends_on:
                        if dep not in self.stages:
                            continue

                        if dep not in visited:
                            if visit(dep):
                                return True
                        elif dep in rec_stack:
                            return True

                    rec_stack.remove(node)
                    return False

                for stage_id in self.stages:
                    if stage_id not in visited:
                        if visit(stage_id):
                            return True
                return False

            if has_cycle():
                errors.append("Circular dependency detected in stage dependencies")

        return (len(errors) == 0, errors)

    def get_execution_order(self) -> List[str]:
        """Compute stage execution order via topological sort.

        Uses Kahn's algorithm to compute a valid execution order
        respecting all dependencies.

        Returns
        -------
        List[str]
            List of stage_ids in execution order

        Raises
        ------
        ValueError
            If circular dependencies detected
        """
        in_degree = {stage_id: 0 for stage_id in self.stages}

        for stage_id, stage in self.stages.items():
            for dep in stage.depends_on:
                in_degree[stage_id] += 1

        queue = deque([sid for sid, degree in in_degree.items() if degree == 0])
        order = []

        while queue:
            stage_id = queue.popleft()
            order.append(stage_id)

            for other_id, other_stage in self.stages.items():
                if stage_id in other_stage.depends_on:
                    in_degree[other_id] -= 1
                    if in_degree[other_id] == 0:
                        queue.append(other_id)

        if len(order) != len(self.stages):
            raise ValueError(
                "Circular dependency detected - cannot compute execution order"
            )

        return order

    def get_stage(self, stage_id: str) -> Optional[Stage]:
        """Get a stage by its ID.

        Parameters
        ----------
        stage_id : str
            Stage identifier

        Returns
        -------
        Optional[Stage]
            Stage object, or None if not found
        """
        return self.stages.get(stage_id)

    def list_stages(self) -> List[str]:
        """List all stage IDs.

        Returns
        -------
        List[str]
            List of stage IDs
        """
        return list(self.stages.keys())

    def to_dict(self) -> Dict[str, Any]:
        """Convert configuration to dictionary.

        Returns
        -------
        Dict[str, Any]
            Configuration as dictionary
        """
        return {
            "pipeline": self.global_settings.get("pipeline", {}),
            "global": self.global_settings.get("global", {}),
            "stages": {
                stage_id: stage.to_dict()
                for stage_id, stage in self.stages.items()
            },
        }

    @classmethod
    def from_dict(cls, config_dict: Dict[str, Any]) -> "PipelineConfig":
        """Create PipelineConfig from dictionary.

        Parameters
        ----------
        config_dict : Dict[str, Any]
            Configuration dictionary

        Returns
        -------
        PipelineConfig
            Constructed config object
        """
        config = cls.__new__(cls)
        config.config_path = Path(".")
        config.raw_config = config_dict
        config.global_settings = {
            "pipeline": config_dict.get("pipeline", {}),
            "global": config_dict.get("global", {}),
        }
        config.stages = {}

        for stage_id, stage_def in config_dict.get("stages", {}).items():
            config.stages[stage_id] = Stage.from_dict(stage_def, stage_id)

        return config
