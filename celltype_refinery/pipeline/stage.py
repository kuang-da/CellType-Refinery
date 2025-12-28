"""Stage representation and validation for pipeline execution."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Tuple


@dataclass
class Stage:
    """Represents a single pipeline stage with inputs, outputs, and dependencies.

    Attributes
    ----------
    name : str
        Human-readable stage name (e.g., "Data Loading and Validation")
    stage_id : str
        Short identifier (e.g., "A", "B", "H")
    script_module : str
        Python module path (e.g., "celltype_refinery.stages.preprocessing")
    depends_on : List[str]
        List of stage IDs this stage depends on
    inputs : Dict[str, str]
        Dictionary mapping input names to file/directory paths
    outputs : Dict[str, str]
        Dictionary mapping output names to file/directory paths
    args : Dict[str, Any]
        Dictionary of command-line arguments for the stage
    required_files : List[str]
        List of files that must exist before stage can run
    optional : bool
        Whether this stage is optional and can be skipped

    Example
    -------
    >>> stage = Stage(
    ...     name="Preprocessing",
    ...     stage_id="A",
    ...     script_module="celltype_refinery.stages.preprocessing",
    ...     inputs={"raw_data": "/path/to/data.csv"},
    ...     outputs={"processed": "/path/to/output/"},
    ... )
    >>> valid, errors = stage.validate_inputs()
    """

    name: str
    stage_id: str
    script_module: str
    depends_on: List[str] = field(default_factory=list)
    inputs: Dict[str, str] = field(default_factory=dict)
    outputs: Dict[str, str] = field(default_factory=dict)
    args: Dict[str, Any] = field(default_factory=dict)
    required_files: List[str] = field(default_factory=list)
    optional: bool = False

    def validate_inputs(self) -> Tuple[bool, List[str]]:
        """Check if all input files/directories exist.

        Returns
        -------
        Tuple[bool, List[str]]
            (success, errors) where success is True if all inputs exist,
            errors contains list of missing file descriptions
        """
        errors = []

        for name, path in self.inputs.items():
            if not Path(path).exists():
                errors.append(f"Input '{name}' not found: {path}")

        for path in self.required_files:
            if not Path(path).exists():
                errors.append(f"Required file not found: {path}")

        return (len(errors) == 0, errors)

    def validate_outputs(self) -> Tuple[bool, List[str]]:
        """Check if all output files/directories exist after execution.

        Returns
        -------
        Tuple[bool, List[str]]
            (success, errors) where success is True if all outputs exist,
            errors contains list of missing file descriptions
        """
        errors = []

        for name, path in self.outputs.items():
            if not Path(path).exists():
                errors.append(f"Output '{name}' not found: {path}")

        return (len(errors) == 0, errors)

    def get_command(self) -> List[str]:
        """Build command-line arguments list for subprocess execution.

        Returns
        -------
        List[str]
            Command arguments suitable for subprocess.run()

        Example
        -------
        >>> stage.get_command()
        ['python', '-m', 'celltype_refinery.stages.preprocessing', '--input', '...']
        """
        cmd = ["python", "-m", self.script_module]

        for key, value in self.args.items():
            arg_name = key.replace("_", "-")

            if value is True:
                cmd.append(f"--{arg_name}")
            elif value is False or value is None:
                continue
            else:
                cmd.append(f"--{arg_name}")
                cmd.append(str(value))

        return cmd

    def to_dict(self) -> Dict[str, Any]:
        """Convert stage to dictionary for serialization.

        Returns
        -------
        Dict[str, Any]
            Stage configuration as dictionary
        """
        return {
            "name": self.name,
            "stage_id": self.stage_id,
            "script_module": self.script_module,
            "depends_on": self.depends_on,
            "inputs": self.inputs,
            "outputs": self.outputs,
            "args": self.args,
            "required_files": self.required_files,
            "optional": self.optional,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any], stage_id: str) -> "Stage":
        """Create Stage from dictionary.

        Parameters
        ----------
        data : Dict[str, Any]
            Stage configuration dictionary
        stage_id : str
            Stage identifier

        Returns
        -------
        Stage
            Constructed stage object
        """
        return cls(
            name=data.get("name", stage_id),
            stage_id=stage_id,
            script_module=data["script_module"],
            depends_on=data.get("depends_on", []),
            inputs=data.get("inputs", {}),
            outputs=data.get("outputs", {}),
            args=data.get("args", {}),
            required_files=data.get("required_files", []),
            optional=data.get("optional", False),
        )
