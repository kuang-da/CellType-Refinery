"""Review module for cell-type annotation quality assessment.

This module provides rule-based flagging and validation of cell-type annotations
with support for tissue-specific expectations via YAML templates.

Example Usage
-------------
Basic review:

    >>> from celltype_refinery.core.review import ReviewEngine, TissueTemplate
    >>> template = TissueTemplate(Path("configs/tissues/fallopian_tube.yaml"))
    >>> engine = ReviewEngine(template=template)
    >>> result = engine.execute(
    ...     composition_dir=Path("out/composition/"),
    ...     spatial_dir=Path("out/spatial/"),
    ... )
    >>> print(f"Found {result.n_issues} issues ({result.n_critical} critical)")

With custom configuration:

    >>> from celltype_refinery.core.review import ReviewConfig
    >>> config = ReviewConfig.from_yaml("review_config.yaml")
    >>> engine = ReviewEngine(config=config)
    >>> result = engine.execute(composition_dir=Path("out/composition/"))

Export results:

    >>> from celltype_refinery.core.review import export_all
    >>> export_all(result, Path("out/review/"))
"""

__version__ = "1.0.0"

# Configuration
from .config import (
    ReviewConfig,
    RuleConfig,
    TissueTemplate,
    CellTypeExpectation,
    RegionalExpectation,
    BiologyRuleConfig,
    KnownArtifact,
)

# Aggregator
from .aggregator import MetricsAggregator

# Engine
from .engine import (
    ReviewEngine,
    ReviewResult,
    MultiColumnReviewResult,
)

# Rules
from .rules.base import BaseRule, FlaggedIssue
from .rules.registry import RuleRegistry

# Export
from .export import (
    export_json,
    export_csv,
    export_markdown,
    export_provenance,
    export_all,
)

__all__ = [
    # Version
    "__version__",
    # Config
    "ReviewConfig",
    "RuleConfig",
    "TissueTemplate",
    "CellTypeExpectation",
    "RegionalExpectation",
    "BiologyRuleConfig",
    "KnownArtifact",
    # Aggregator
    "MetricsAggregator",
    # Engine
    "ReviewEngine",
    "ReviewResult",
    "MultiColumnReviewResult",
    # Rules
    "BaseRule",
    "FlaggedIssue",
    "RuleRegistry",
    # Export
    "export_json",
    "export_csv",
    "export_markdown",
    "export_provenance",
    "export_all",
]
