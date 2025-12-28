"""Pipeline orchestration module.

Provides YAML-based configuration and stage execution with
checkpoint support and dependency resolution.

Example Usage
-------------
>>> from celltype_refinery.pipeline import (
...     PipelineConfig,
...     PipelineExecutor,
...     PipelineLogger,
... )
>>> # Load configuration
>>> config = PipelineConfig("pipeline.yaml")
>>> config.load()
>>> config.parse_stages()
>>> # Setup logging
>>> logger = PipelineLogger("logs/")
>>> logger.setup()
>>> # Execute pipeline
>>> executor = PipelineExecutor(config, logger)
>>> exit_code = executor.run()
"""

__version__ = "1.0.0"

# Stage representation
from .stage import Stage

# Configuration
from .config import PipelineConfig

# Logging
from .logger import (
    ColoredFormatter,
    PipelineLogger,
)

# Execution
from .executor import (
    InMemoryExecutor,
    PipelineExecutor,
)

__all__ = [
    # Version
    "__version__",
    # Stage
    "Stage",
    # Config
    "PipelineConfig",
    # Logging
    "ColoredFormatter",
    "PipelineLogger",
    # Execution
    "InMemoryExecutor",
    "PipelineExecutor",
]
