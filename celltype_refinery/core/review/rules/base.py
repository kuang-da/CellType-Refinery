"""Base classes for flagging rules.

Provides:
- FlaggedIssue: Dataclass for flagged issues
- BaseRule: Abstract base class for all rules
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any, Dict, List, Literal, Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from ..aggregator import MetricsAggregator
    from ..config import TissueTemplate, RuleConfig


@dataclass
class FlaggedIssue:
    """A flagged issue from a review rule.

    Attributes
    ----------
    cell_type : str
        Cell type involved
    rule_id : str
        Rule that triggered the flag
    severity : str
        Severity level (critical, warning, note)
    description : str
        Human-readable description
    evidence : str
        Supporting evidence/data
    recommendation : str
        Suggested action
    region : str
        Region involved (if applicable)
    value : float
        Observed value
    expected : str
        Expected value/range
    column : str
        Cell type column (for multi-column)
    metadata : Dict
        Additional metadata
    """

    cell_type: str
    rule_id: str
    severity: Literal["critical", "warning", "note"]
    description: str
    evidence: str
    recommendation: Optional[str] = None
    region: Optional[str] = None
    value: Optional[float] = None
    expected: Optional[str] = None
    column: Optional[str] = None
    metadata: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            "cell_type": self.cell_type,
            "rule_id": self.rule_id,
            "severity": self.severity,
            "description": self.description,
            "evidence": self.evidence,
            "recommendation": self.recommendation,
            "region": self.region,
            "value": self.value,
            "expected": self.expected,
            "column": self.column,
            **self.metadata,
        }


class BaseRule(ABC):
    """Abstract base class for flagging rules.

    All rules must implement:
    - rule_id: Unique identifier
    - category: Rule category (proportion, spatial, biology, quality)
    - check(): Main check method

    Attributes
    ----------
    rule_id : str
        Unique identifier for this rule
    category : str
        Rule category
    default_severity : str
        Default severity if not overridden
    description : str
        Human-readable description
    """

    rule_id: str = "BASE_RULE"
    category: str = "general"
    default_severity: str = "warning"
    description: str = "Base rule"

    def __init__(
        self,
        template: "TissueTemplate",
        config: Optional["RuleConfig"] = None,
    ):
        """Initialize rule with template and optional config.

        Parameters
        ----------
        template : TissueTemplate
            Tissue-specific expectations
        config : RuleConfig, optional
            Rule-specific configuration overrides
        """
        self.template = template
        self.config = config
        self._severity = self._resolve_severity()

    def _resolve_severity(self) -> str:
        """Resolve severity from config or default."""
        if self.config and self.config.severity:
            return self.config.severity
        return self.default_severity

    def get_threshold(self, name: str, default: Any = None) -> Any:
        """Get threshold value from config or default.

        Parameters
        ----------
        name : str
            Threshold name
        default : Any
            Default value if not in config

        Returns
        -------
        Any
            Threshold value
        """
        if self.config and name in self.config.thresholds:
            return self.config.thresholds[name]
        return default

    def create_issue(
        self,
        cell_type: str,
        description: str,
        evidence: str,
        recommendation: Optional[str] = None,
        severity: Optional[str] = None,
        region: Optional[str] = None,
        value: Optional[float] = None,
        expected: Optional[str] = None,
        **metadata,
    ) -> FlaggedIssue:
        """Create a FlaggedIssue with rule context.

        Parameters
        ----------
        cell_type : str
            Cell type involved
        description : str
            Issue description
        evidence : str
            Supporting evidence
        recommendation : str, optional
            Suggested action
        severity : str, optional
            Override default severity
        region : str, optional
            Region involved
        value : float, optional
            Observed value
        expected : str, optional
            Expected value/range
        **metadata
            Additional metadata

        Returns
        -------
        FlaggedIssue
            Created issue
        """
        return FlaggedIssue(
            cell_type=cell_type,
            rule_id=self.rule_id,
            severity=severity or self._severity,
            description=description,
            evidence=evidence,
            recommendation=recommendation,
            region=region,
            value=value,
            expected=expected,
            metadata=metadata,
        )

    @abstractmethod
    def check(self, aggregator: "MetricsAggregator") -> List[FlaggedIssue]:
        """Run the rule check.

        Parameters
        ----------
        aggregator : MetricsAggregator
            Data aggregator with loaded metrics

        Returns
        -------
        List[FlaggedIssue]
            List of flagged issues (empty if none)
        """
        pass

    def is_applicable(self, aggregator: "MetricsAggregator") -> bool:
        """Check if rule is applicable to current data.

        Override in subclasses for conditional rules.

        Parameters
        ----------
        aggregator : MetricsAggregator
            Data aggregator

        Returns
        -------
        bool
            True if rule should run
        """
        return True
