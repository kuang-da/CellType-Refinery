"""Configuration classes for the Review module.

Provides:
- RuleConfig: Per-rule configuration with thresholds
- TissueTemplate: Tissue-specific cell type expectations
- ReviewConfig: Overall review configuration
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import yaml


@dataclass
class CellTypeExpectation:
    """Expected characteristics for a cell type.

    Attributes
    ----------
    name : str
        Cell type name
    expected_pct : Tuple[float, float]
        Expected global percentage range (min, max)
    expected_regions : Dict[str, Tuple[float, float]]
        Expected percentage range per region
    expected_moran : Tuple[float, float]
        Expected Moran's I range (min, max)
    canonical_markers : List[str]
        Expected positive markers
    excluded_markers : List[str]
        Expected negative markers
    notes : str
        Additional notes
    """

    name: str
    expected_pct: Optional[Tuple[float, float]] = None
    expected_regions: Dict[str, Tuple[float, float]] = field(default_factory=dict)
    expected_moran: Optional[Tuple[float, float]] = None
    canonical_markers: List[str] = field(default_factory=list)
    excluded_markers: List[str] = field(default_factory=list)
    notes: str = ""


@dataclass
class RegionalExpectation:
    """Expected characteristics for a region.

    Attributes
    ----------
    name : str
        Region name
    dominant_types : List[str]
        Cell types expected to be dominant
    rare_types : List[str]
        Cell types expected to be rare/absent
    biology_metrics : Dict[str, Tuple[float, float]]
        Expected biology metric ranges
    """

    name: str
    dominant_types: List[str] = field(default_factory=list)
    rare_types: List[str] = field(default_factory=list)
    biology_metrics: Dict[str, Tuple[float, float]] = field(default_factory=dict)


@dataclass
class KnownArtifact:
    """Known problematic pattern to flag.

    Attributes
    ----------
    name : str
        Artifact name/identifier
    description : str
        Description of the artifact
    detection_method : str
        How to detect (sample_concentrated, extreme_clustering, etc.)
    threshold : float
        Detection threshold
    recommendation : str
        What to do about it
    """

    name: str
    description: str = ""
    detection_method: str = ""
    threshold: float = 0.0
    recommendation: str = ""


@dataclass
class BiologyRuleConfig:
    """Configuration for biology-based rules.

    These are tissue-specific rules that vary by tissue type.

    Attributes
    ----------
    enabled : bool
        Whether biology rules are enabled
    rules : List[Dict[str, Any]]
        List of biology rule definitions
    """

    enabled: bool = True
    rules: List[Dict[str, Any]] = field(default_factory=list)


class TissueTemplate:
    """Tissue-specific cell type expectations loaded from YAML template.

    Parameters
    ----------
    path : Path, optional
        Path to YAML template file

    Example template structure:
    ```yaml
    tissue: "template"
    version: "1.0"

    cell_types:
      Epithelium:
        expected_pct: [20, 60]
        expected_regions:
          fimbriae: [30, 70]
          isthmus: [5, 20]
        expected_moran: [0.2, 0.5]
        canonical_markers: ["Pan-CK", "E-cadherin"]

    regions:
      fimbriae:
        dominant_types: ["Epithelium", "Ciliated Epithelium"]
        rare_types: ["Smooth Muscle Cells"]

    biology_rules:
      enabled: true
      rules:
        - rule_id: "ES_RATIO"
          description: "Epithelial:Stromal ratio check"
          regional_thresholds:
            fimbriae: [0.5, 2.5]
            isthmus: [0.05, 0.3]

    quality_thresholds:
      max_unassigned_pct: 10.0
      max_hybrid_pct: 5.0
      max_low_confidence_pct: 20.0

    known_artifacts:
      - name: "sample-concentrated"
        description: "Cell type concentrated in single sample"
        detection_method: "sample_concentration"
        threshold: 0.5
    ```
    """

    def __init__(self, path: Optional[Path] = None):
        self.path = path
        self.tissue = "unknown"
        self.version = "1.0"
        self._cell_types: Dict[str, CellTypeExpectation] = {}
        self._regions: Dict[str, RegionalExpectation] = {}
        self._biology_rules = BiologyRuleConfig()
        self._quality_thresholds: Dict[str, float] = {}
        self._known_artifacts: List[KnownArtifact] = []

        if path is not None:
            self._load(path)

    def _load(self, path: Path) -> None:
        """Load template from YAML file."""
        with open(path) as f:
            data = yaml.safe_load(f)

        self.tissue = data.get("tissue", "unknown")
        self.version = data.get("version", "1.0")

        # Load cell type expectations
        for name, ct_data in data.get("cell_types", {}).items():
            self._cell_types[name] = CellTypeExpectation(
                name=name,
                expected_pct=tuple(ct_data["expected_pct"])
                if "expected_pct" in ct_data
                else None,
                expected_regions={
                    r: tuple(v) for r, v in ct_data.get("expected_regions", {}).items()
                },
                expected_moran=tuple(ct_data["expected_moran"])
                if "expected_moran" in ct_data
                else None,
                canonical_markers=ct_data.get("canonical_markers", []),
                excluded_markers=ct_data.get("excluded_markers", []),
                notes=ct_data.get("notes", ""),
            )

        # Load regional expectations
        for name, reg_data in data.get("regions", {}).items():
            self._regions[name] = RegionalExpectation(
                name=name,
                dominant_types=reg_data.get("dominant_types", []),
                rare_types=reg_data.get("rare_types", []),
                biology_metrics={
                    k: tuple(v)
                    for k, v in reg_data.get("biology_metrics", {}).items()
                },
            )

        # Load biology rules
        bio_data = data.get("biology_rules", {})
        self._biology_rules = BiologyRuleConfig(
            enabled=bio_data.get("enabled", True),
            rules=bio_data.get("rules", []),
        )

        # Load quality thresholds
        self._quality_thresholds = data.get("quality_thresholds", {})

        # Load known artifacts
        for art_data in data.get("known_artifacts", []):
            self._known_artifacts.append(
                KnownArtifact(
                    name=art_data.get("name", ""),
                    description=art_data.get("description", ""),
                    detection_method=art_data.get("detection_method", ""),
                    threshold=art_data.get("threshold", 0.0),
                    recommendation=art_data.get("recommendation", ""),
                )
            )

    def get_cell_type(self, name: str) -> Optional[CellTypeExpectation]:
        """Get cell type expectation by name."""
        return self._cell_types.get(name)

    def get_expected_pct_range(
        self, cell_type: str
    ) -> Optional[Tuple[float, float]]:
        """Get expected global percentage range for cell type."""
        ct = self._cell_types.get(cell_type)
        return ct.expected_pct if ct else None

    def get_expected_moran_range(
        self, cell_type: str
    ) -> Optional[Tuple[float, float]]:
        """Get expected Moran's I range for cell type."""
        ct = self._cell_types.get(cell_type)
        return ct.expected_moran if ct else None

    def get_region(self, name: str) -> Optional[RegionalExpectation]:
        """Get regional expectation by name."""
        return self._regions.get(name)

    def get_biology_rules(self) -> BiologyRuleConfig:
        """Get biology rule configuration."""
        return self._biology_rules

    def get_quality_threshold(self, name: str, default: float = 0.0) -> float:
        """Get quality threshold by name."""
        return self._quality_thresholds.get(name, default)

    def get_known_artifacts(self) -> List[Dict[str, Any]]:
        """Get known artifacts as list of dicts."""
        return [
            {
                "name": a.name,
                "description": a.description,
                "detection_method": a.detection_method,
                "threshold": a.threshold,
                "recommendation": a.recommendation,
            }
            for a in self._known_artifacts
        ]

    def get_all_cell_types(self) -> List[str]:
        """Get all defined cell type names."""
        return list(self._cell_types.keys())

    def get_all_regions(self) -> List[str]:
        """Get all defined region names."""
        return list(self._regions.keys())

    @classmethod
    def empty(cls) -> "TissueTemplate":
        """Create empty template with no expectations."""
        return cls(path=None)


@dataclass
class RuleConfig:
    """Configuration for a single flagging rule.

    Attributes
    ----------
    rule_id : str
        Rule identifier
    enabled : bool
        Whether rule is active
    severity : str
        Default severity (critical, warning, note)
    thresholds : Dict[str, Any]
        Rule-specific thresholds
    """

    rule_id: str
    enabled: bool = True
    severity: str = "warning"
    thresholds: Dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "RuleConfig":
        """Create from dictionary."""
        return cls(
            rule_id=data.get("rule_id", ""),
            enabled=data.get("enabled", True),
            severity=data.get("severity", "warning"),
            thresholds=data.get("thresholds", {}),
        )


@dataclass
class ReviewConfig:
    """Overall configuration for the Review module.

    Attributes
    ----------
    rules : Dict[str, RuleConfig]
        Per-rule configurations
    template_path : Path
        Path to tissue template
    cell_type_cols : List[str]
        Cell type columns to review
    skip_rules : List[str]
        Rules to skip
    """

    rules: Dict[str, RuleConfig] = field(default_factory=dict)
    template_path: Optional[Path] = None
    cell_type_cols: List[str] = field(default_factory=list)
    skip_rules: List[str] = field(default_factory=list)

    @classmethod
    def from_yaml(cls, path: Path) -> "ReviewConfig":
        """Load configuration from YAML file."""
        with open(path) as f:
            data = yaml.safe_load(f)

        rules = {}
        for rule_data in data.get("rules", []):
            rule = RuleConfig.from_dict(rule_data)
            rules[rule.rule_id] = rule

        template_path = None
        if "template_path" in data:
            template_path = Path(data["template_path"])

        return cls(
            rules=rules,
            template_path=template_path,
            cell_type_cols=data.get("cell_type_cols", []),
            skip_rules=data.get("skip_rules", []),
        )

    @classmethod
    def default(cls) -> "ReviewConfig":
        """Create default configuration with all rules enabled."""
        return cls()

    def get_rule_config(self, rule_id: str) -> Optional[RuleConfig]:
        """Get configuration for specific rule."""
        return self.rules.get(rule_id)

    def is_rule_enabled(self, rule_id: str) -> bool:
        """Check if rule is enabled."""
        if rule_id in self.skip_rules:
            return False
        config = self.rules.get(rule_id)
        return config.enabled if config else True

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            "rules": {
                rid: {
                    "rule_id": r.rule_id,
                    "enabled": r.enabled,
                    "severity": r.severity,
                    "thresholds": r.thresholds,
                }
                for rid, r in self.rules.items()
            },
            "template_path": str(self.template_path) if self.template_path else None,
            "cell_type_cols": self.cell_type_cols,
            "skip_rules": self.skip_rules,
        }
