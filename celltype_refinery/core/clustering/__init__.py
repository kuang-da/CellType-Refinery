"""Clustering module for cell population identification (Stage H).

Provides Leiden clustering and targeted subclustering functionality
with optional GPU acceleration.

Pipeline Stages
---------------
- Stage H (Clustering): Leiden clustering, differential expression,
  marker scoring for cell-type annotation

Example Usage
-------------
>>> from celltype_refinery.core.clustering import (
...     ClusteringEngine, DERunner, StageHConfig,
... )
>>> # Configure and run clustering
>>> config = StageHConfig()
>>> engine = ClusteringEngine(config)
>>> result = engine.run_clustering(adata, cluster_key="leiden")
>>> # Run differential expression
>>> runner = DERunner(config)
>>> de_result = runner.run_de_tests(adata, cluster_key="leiden")
"""

__version__ = "1.0.0"

# Configuration classes
from .config import (
    ClusteringConfig,
    DEConfig,
    SubclusterConfig,
    ScoringConfig,
    StageHConfig,
)

# Clustering engine
from .engine import (
    ClusteringEngine,
    ClusteringResult,
    GPU_AVAILABLE,
)

# Differential expression
from .de import (
    DERunner,
    DEResult,
    default_de_runner,
)

__all__ = [
    # Version
    "__version__",
    # Config
    "ClusteringConfig",
    "DEConfig",
    "SubclusterConfig",
    "ScoringConfig",
    "StageHConfig",
    # Engine
    "ClusteringEngine",
    "ClusteringResult",
    "GPU_AVAILABLE",
    # DE
    "DERunner",
    "DEResult",
    "default_de_runner",
]
