"""Rule implementations for the Review module.

This package contains flagging rules organized by category:
- proportion: Cell type proportion anomalies
- spatial: Spatial distribution anomalies
- biology: Tissue-specific biology rules (config-driven)
- quality: Annotation quality metrics
"""

from .base import BaseRule, FlaggedIssue
from .registry import RuleRegistry

# Import rule modules to trigger registration
from . import proportion
from . import spatial
from . import biology
from . import quality

__all__ = [
    "BaseRule",
    "FlaggedIssue",
    "RuleRegistry",
]
