"""Unit tests for pipeline orchestration module."""

import pytest
import json
from pathlib import Path

from celltype_refinery.pipeline import (
    Stage,
    PipelineConfig,
    PipelineLogger,
    PipelineExecutor,
    InMemoryExecutor,
)


class TestStage:
    """Tests for Stage dataclass."""

    def test_create_stage(self):
        """Test creating a basic stage."""
        stage = Stage(
            name="Test Stage",
            stage_id="test",
            script_module="test.module",
        )
        assert stage.name == "Test Stage"
        assert stage.stage_id == "test"
        assert stage.script_module == "test.module"
        assert stage.depends_on == []
        assert stage.optional is False

    def test_stage_with_dependencies(self):
        """Test stage with dependencies."""
        stage = Stage(
            name="Dependent Stage",
            stage_id="B",
            script_module="test.module",
            depends_on=["A"],
        )
        assert "A" in stage.depends_on

    def test_validate_inputs_missing(self, tmp_path):
        """Test input validation with missing files."""
        stage = Stage(
            name="Test",
            stage_id="test",
            script_module="test.module",
            inputs={"data": str(tmp_path / "nonexistent.csv")},
        )
        valid, errors = stage.validate_inputs()
        assert not valid
        assert len(errors) == 1
        assert "not found" in errors[0]

    def test_validate_inputs_exists(self, tmp_path):
        """Test input validation with existing files."""
        data_file = tmp_path / "data.csv"
        data_file.write_text("col1,col2\n1,2")

        stage = Stage(
            name="Test",
            stage_id="test",
            script_module="test.module",
            inputs={"data": str(data_file)},
        )
        valid, errors = stage.validate_inputs()
        assert valid
        assert len(errors) == 0

    def test_get_command_basic(self):
        """Test command generation."""
        stage = Stage(
            name="Test",
            stage_id="test",
            script_module="test.module",
        )
        cmd = stage.get_command()
        assert cmd == ["python", "-m", "test.module"]

    def test_get_command_with_args(self):
        """Test command generation with arguments."""
        stage = Stage(
            name="Test",
            stage_id="test",
            script_module="test.module",
            args={"input_file": "/path/to/file", "verbose": True},
        )
        cmd = stage.get_command()
        assert "python" in cmd
        assert "-m" in cmd
        assert "test.module" in cmd
        assert "--input-file" in cmd
        assert "/path/to/file" in cmd
        assert "--verbose" in cmd

    def test_get_command_false_flag_skipped(self):
        """Test that False flags are skipped."""
        stage = Stage(
            name="Test",
            stage_id="test",
            script_module="test.module",
            args={"verbose": False, "debug": None},
        )
        cmd = stage.get_command()
        assert "--verbose" not in cmd
        assert "--debug" not in cmd

    def test_to_dict(self):
        """Test converting stage to dictionary."""
        stage = Stage(
            name="Test",
            stage_id="test",
            script_module="test.module",
            depends_on=["A"],
            args={"output": "/out"},
        )
        d = stage.to_dict()
        assert d["name"] == "Test"
        assert d["stage_id"] == "test"
        assert "A" in d["depends_on"]
        assert d["args"]["output"] == "/out"

    def test_from_dict(self):
        """Test creating stage from dictionary."""
        data = {
            "name": "From Dict",
            "script_module": "dict.module",
            "depends_on": ["X"],
            "args": {"key": "value"},
        }
        stage = Stage.from_dict(data, "from_dict")
        assert stage.stage_id == "from_dict"
        assert stage.name == "From Dict"
        assert stage.script_module == "dict.module"


class TestPipelineConfig:
    """Tests for PipelineConfig class."""

    @pytest.fixture
    def sample_config_file(self, tmp_path) -> Path:
        """Create sample pipeline config file."""
        config = {
            "pipeline": {"name": "Test Pipeline", "version": "1.0"},
            "global": {"output_dir": "/output"},
            "stages": {
                "A": {
                    "name": "Stage A",
                    "script_module": "test.stage_a",
                    "args": {"input": "{global.output_dir}/data.csv"},
                },
                "B": {
                    "name": "Stage B",
                    "script_module": "test.stage_b",
                    "depends_on": ["A"],
                },
            },
        }
        path = tmp_path / "pipeline.yaml"
        import yaml

        with open(path, "w") as f:
            yaml.dump(config, f)
        return path

    def test_init(self, tmp_path):
        """Test config initialization."""
        config = PipelineConfig(str(tmp_path / "config.yaml"))
        assert config.stages == {}

    def test_load(self, sample_config_file):
        """Test loading config from YAML."""
        config = PipelineConfig(str(sample_config_file))
        config.load()
        assert "pipeline" in config.global_settings
        assert config.global_settings["pipeline"]["name"] == "Test Pipeline"

    def test_parse_stages(self, sample_config_file):
        """Test parsing stages from config."""
        config = PipelineConfig(str(sample_config_file))
        config.load()
        config.parse_stages()

        assert "A" in config.stages
        assert "B" in config.stages
        assert config.stages["A"].name == "Stage A"
        assert "A" in config.stages["B"].depends_on

    def test_resolve_paths(self, sample_config_file):
        """Test path template resolution."""
        config = PipelineConfig(str(sample_config_file))
        config.load()
        config.parse_stages()

        # Check that {global.output_dir} was resolved
        assert config.stages["A"].args["input"] == "/output/data.csv"

    def test_validate_dependencies_valid(self, sample_config_file):
        """Test dependency validation with valid deps."""
        config = PipelineConfig(str(sample_config_file))
        config.load()
        config.parse_stages()

        valid, errors = config.validate_dependencies()
        assert valid
        assert len(errors) == 0

    def test_validate_dependencies_missing(self, tmp_path):
        """Test dependency validation with missing deps."""
        config_data = {
            "stages": {
                "A": {
                    "script_module": "test.a",
                    "depends_on": ["nonexistent"],
                },
            },
        }
        path = tmp_path / "config.yaml"
        import yaml

        with open(path, "w") as f:
            yaml.dump(config_data, f)

        config = PipelineConfig(str(path))
        config.load()
        config.parse_stages()

        valid, errors = config.validate_dependencies()
        assert not valid
        assert any("unknown stage" in e for e in errors)

    def test_get_execution_order(self, sample_config_file):
        """Test computing execution order."""
        config = PipelineConfig(str(sample_config_file))
        config.load()
        config.parse_stages()

        order = config.get_execution_order()
        assert order.index("A") < order.index("B")


class TestPipelineLogger:
    """Tests for PipelineLogger class."""

    def test_init(self, tmp_path):
        """Test logger initialization."""
        logger = PipelineLogger(str(tmp_path / "logs"))
        assert logger.log_dir.exists()

    def test_setup(self, tmp_path):
        """Test logger setup."""
        logger = PipelineLogger(str(tmp_path / "logs"))
        logger.setup()
        assert len(logger.logger.handlers) == 2

    def test_format_duration_seconds(self):
        """Test duration formatting for seconds."""
        assert PipelineLogger.format_duration(45.2) == "45.2s"

    def test_format_duration_minutes(self):
        """Test duration formatting for minutes."""
        result = PipelineLogger.format_duration(125)
        assert "m" in result
        assert "s" in result

    def test_format_duration_hours(self):
        """Test duration formatting for hours."""
        result = PipelineLogger.format_duration(7300)
        assert "h" in result
        assert "m" in result


class TestPipelineExecutor:
    """Tests for PipelineExecutor class."""

    @pytest.fixture
    def mock_config(self, tmp_path) -> PipelineConfig:
        """Create mock config for executor tests."""
        config_data = {
            "pipeline": {"version": "1.0"},
            "stages": {
                "A": {"name": "Test A", "script_module": "test.a"},
            },
        }
        path = tmp_path / "config.yaml"
        import yaml

        with open(path, "w") as f:
            yaml.dump(config_data, f)

        config = PipelineConfig(str(path))
        config.load()
        config.parse_stages()
        return config

    @pytest.fixture
    def mock_logger(self, tmp_path) -> PipelineLogger:
        """Create mock logger for tests."""
        logger = PipelineLogger(str(tmp_path / "logs"))
        logger.setup()
        return logger

    def test_init(self, mock_config, mock_logger, tmp_path):
        """Test executor initialization."""
        executor = PipelineExecutor(
            mock_config,
            mock_logger,
            state_file=str(tmp_path / "state.json"),
        )
        assert executor.completed_stages == []

    def test_save_and_load_state(self, mock_config, mock_logger, tmp_path):
        """Test state save and load."""
        state_file = tmp_path / "state.json"
        executor = PipelineExecutor(
            mock_config, mock_logger, state_file=str(state_file)
        )

        executor.completed_stages = ["A", "B"]
        executor.save_state()

        assert state_file.exists()

        # Create new executor and load state
        executor2 = PipelineExecutor(
            mock_config, mock_logger, state_file=str(state_file)
        )
        executor2.load_state()
        assert executor2.completed_stages == ["A", "B"]

    def test_clear_state(self, mock_config, mock_logger, tmp_path):
        """Test clearing state."""
        state_file = tmp_path / "state.json"
        executor = PipelineExecutor(
            mock_config, mock_logger, state_file=str(state_file)
        )

        executor.completed_stages = ["A"]
        executor.save_state()
        assert state_file.exists()

        executor.clear_state()
        assert not state_file.exists()
        assert executor.completed_stages == []


class TestInMemoryExecutor:
    """Tests for InMemoryExecutor class."""

    def test_register_stage(self):
        """Test registering stages."""
        executor = InMemoryExecutor()
        executor.register_stage("A", lambda **k: "result_a")
        assert "A" in executor.stages

    def test_register_with_deps(self):
        """Test registering stage with dependencies."""
        executor = InMemoryExecutor()
        executor.register_stage("A", lambda **k: None)
        executor.register_stage("B", lambda **k: None, depends_on=["A"])
        assert executor.stages["B"]["depends_on"] == ["A"]

    def test_run_simple(self):
        """Test running simple pipeline."""
        executor = InMemoryExecutor()
        executor.register_stage("A", lambda **k: "a_result")
        executor.register_stage(
            "B", lambda stage_results, **k: stage_results["A"] + "_b"
        )

        results = executor.run()
        assert results["A"] == "a_result"
        assert results["B"] == "a_result_b"

    def test_execution_order(self):
        """Test execution order with dependencies."""
        executor = InMemoryExecutor()
        order = []

        executor.register_stage("C", lambda **k: order.append("C"), depends_on=["B"])
        executor.register_stage("B", lambda **k: order.append("B"), depends_on=["A"])
        executor.register_stage("A", lambda **k: order.append("A"))

        executor.run()
        assert order == ["A", "B", "C"]
