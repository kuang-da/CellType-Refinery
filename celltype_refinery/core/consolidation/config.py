"""Configuration dataclasses for consolidation module.

This module defines configuration structures for:
- Manual overrides (cluster_id -> final_label)
- Relabel rules (from_label -> to_label)
- Orphan rescue settings
- IEL (intraepithelial immune) rescue settings
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml


@dataclass
class OverrideRule:
    """Manual override for a specific cluster.

    Attributes
    ----------
    cluster_id : str
        The cluster ID to override
    final_label : str
        The final label to assign
    reason : str
        Reason for override (for audit trail)
    """

    cluster_id: str
    final_label: str
    reason: str = ""


@dataclass
class RelabelRule:
    """Global relabel rule (applies to all matching labels).

    Attributes
    ----------
    from_label : str
        Original label to match
    to_label : str
        New label to assign
    reason : str
        Reason for relabeling (for audit trail)
    """

    from_label: str
    to_label: str
    reason: str = ""


@dataclass
class OrphanRule:
    """Custom rule for a specific orphan subtype.

    Attributes
    ----------
    subtype : str
        The subtype name (e.g., "Lymphatic Endothelium")
    action : str
        Action to take: "rescue", "keep_unassigned", "expert_review"
    min_score : float, optional
        Minimum subtype score to consider (overrides global threshold)
    suffix : str, optional
        Custom suffix for this subtype (overrides global suffix)
    flag : str
        Optional flag for tracking
    """

    subtype: str
    action: str = "rescue"
    min_score: Optional[float] = None
    suffix: Optional[str] = None
    flag: str = ""


@dataclass
class OrphanRescueConfig:
    """Configuration for orphan rescue.

    Attributes
    ----------
    enabled : bool
        Whether orphan rescue is enabled
    subtype_threshold : float
        Minimum subtype score to consider as orphan candidate
    root_fail_threshold : float
        Maximum root score to consider as "failed" root gate
    suffix : str
        Suffix to add to rescued labels (default: "(orphan)")
    rules : List[OrphanRule]
        Custom rules for specific subtypes
    """

    enabled: bool = True
    subtype_threshold: float = 0.8
    root_fail_threshold: float = 0.5
    suffix: str = "(orphan)"
    rules: List[OrphanRule] = field(default_factory=list)


@dataclass
class IELRescueConfig:
    """Configuration for IEL (intraepithelial immune) rescue.

    Attributes
    ----------
    enabled : bool
        Whether IEL rescue is enabled
    cd45_min_pos_frac : float
        Minimum CD45 positive fraction to consider for rescue
    lymphoid_score_threshold : float
        Minimum lymphoid score to classify as IEL
    myeloid_score_threshold : float
        Minimum myeloid score to classify as tissue-resident macrophage
    suffix : str
        Suffix to add to rescued labels (empty by default)
    """

    enabled: bool = False
    cd45_min_pos_frac: float = 0.30
    lymphoid_score_threshold: float = 1.0
    myeloid_score_threshold: float = 1.0
    suffix: str = ""


@dataclass
class ConsolidationConfig:
    """Main configuration for consolidation engine.

    Attributes
    ----------
    version : str
        Config version for tracking
    description : str
        Optional description of this config
    overrides : List[OverrideRule]
        Manual cluster overrides
    relabel_rules : List[RelabelRule]
        Global relabel rules
    orphan_rescue : OrphanRescueConfig
        Orphan rescue configuration
    iel_rescue : IELRescueConfig
        IEL rescue configuration
    """

    version: str = "1.0"
    description: str = ""
    overrides: List[OverrideRule] = field(default_factory=list)
    relabel_rules: List[RelabelRule] = field(default_factory=list)
    orphan_rescue: OrphanRescueConfig = field(default_factory=OrphanRescueConfig)
    iel_rescue: IELRescueConfig = field(default_factory=IELRescueConfig)

    def get_override_map(self) -> Dict[str, OverrideRule]:
        """Get mapping from cluster_id to override rule."""
        return {o.cluster_id: o for o in self.overrides}

    def get_relabel_map(self) -> Dict[str, RelabelRule]:
        """Get mapping from from_label to relabel rule."""
        return {r.from_label: r for r in self.relabel_rules}

    @classmethod
    def from_yaml(cls, path: Path) -> "ConsolidationConfig":
        """Load configuration from YAML file.

        Parameters
        ----------
        path : Path
            Path to YAML config file

        Returns
        -------
        ConsolidationConfig
            Loaded configuration
        """
        path = Path(path)
        with open(path) as f:
            data = yaml.safe_load(f)

        if data is None:
            return cls()

        # Parse overrides
        overrides = []
        for o in data.get("overrides", []):
            overrides.append(
                OverrideRule(
                    cluster_id=str(o["cluster_id"]),
                    final_label=o["final_label"],
                    reason=o.get("reason", ""),
                )
            )

        # Parse relabel rules
        relabel_rules = []
        for r in data.get("relabel_rules", []):
            relabel_rules.append(
                RelabelRule(
                    from_label=r["from_label"],
                    to_label=r["to_label"],
                    reason=r.get("reason", ""),
                )
            )

        # Parse orphan rescue config
        orphan_data = data.get("orphan_rescue", {})
        orphan_rules = []
        for r in orphan_data.get("rules", []):
            orphan_rules.append(
                OrphanRule(
                    subtype=r["subtype"],
                    action=r.get("action", "rescue"),
                    min_score=r.get("min_score"),
                    suffix=r.get("suffix"),
                    flag=r.get("flag", ""),
                )
            )
        orphan_rescue = OrphanRescueConfig(
            enabled=orphan_data.get("enabled", True),
            subtype_threshold=orphan_data.get("subtype_threshold", 0.8),
            root_fail_threshold=orphan_data.get("root_fail_threshold", 0.5),
            suffix=orphan_data.get("suffix", "(orphan)"),
            rules=orphan_rules,
        )

        # Parse IEL rescue config
        iel_data = data.get("iel_rescue", {})
        iel_rescue = IELRescueConfig(
            enabled=iel_data.get("enabled", False),
            cd45_min_pos_frac=iel_data.get("cd45_min_pos_frac", 0.30),
            lymphoid_score_threshold=iel_data.get("lymphoid_score_threshold", 1.0),
            myeloid_score_threshold=iel_data.get("myeloid_score_threshold", 1.0),
            suffix=iel_data.get("suffix", ""),
        )

        return cls(
            version=data.get("version", "1.0"),
            description=data.get("description", ""),
            overrides=overrides,
            relabel_rules=relabel_rules,
            orphan_rescue=orphan_rescue,
            iel_rescue=iel_rescue,
        )

    def to_yaml(self, path: Path) -> None:
        """Save configuration to YAML file.

        Parameters
        ----------
        path : Path
            Output path for YAML file
        """
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)

        data: Dict[str, Any] = {
            "version": self.version,
            "description": self.description,
        }

        # Overrides
        if self.overrides:
            data["overrides"] = [
                {
                    "cluster_id": o.cluster_id,
                    "final_label": o.final_label,
                    "reason": o.reason,
                }
                for o in self.overrides
            ]

        # Relabel rules
        if self.relabel_rules:
            data["relabel_rules"] = [
                {
                    "from_label": r.from_label,
                    "to_label": r.to_label,
                    "reason": r.reason,
                }
                for r in self.relabel_rules
            ]

        # Orphan rescue
        orphan_data: Dict[str, Any] = {
            "enabled": self.orphan_rescue.enabled,
            "subtype_threshold": self.orphan_rescue.subtype_threshold,
            "root_fail_threshold": self.orphan_rescue.root_fail_threshold,
            "suffix": self.orphan_rescue.suffix,
        }
        if self.orphan_rescue.rules:
            orphan_data["rules"] = [
                {
                    "subtype": r.subtype,
                    "action": r.action,
                    "min_score": r.min_score,
                    "suffix": r.suffix,
                    "flag": r.flag,
                }
                for r in self.orphan_rescue.rules
            ]
        data["orphan_rescue"] = orphan_data

        # IEL rescue
        data["iel_rescue"] = {
            "enabled": self.iel_rescue.enabled,
            "cd45_min_pos_frac": self.iel_rescue.cd45_min_pos_frac,
            "lymphoid_score_threshold": self.iel_rescue.lymphoid_score_threshold,
            "myeloid_score_threshold": self.iel_rescue.myeloid_score_threshold,
            "suffix": self.iel_rescue.suffix,
        }

        with open(path, "w") as f:
            yaml.dump(data, f, default_flow_style=False, sort_keys=False)
