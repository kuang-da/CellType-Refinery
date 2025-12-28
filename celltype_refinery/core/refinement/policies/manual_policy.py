"""
ManualPolicy: YAML configuration-based refinement for cell-type annotations.

This policy parses YAML configuration files and generates RefinePlans
for manual expert curation. Supports:
- overrides: Direct cell-type label assignments
- subcluster: Targeted subclustering with custom parameters
- merge: Combining multiple clusters into one label
- rescoring: Optional rescoring after modifications
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Set

import yaml

try:
    import scanpy as sc
except ImportError:
    sc = None

from ..plan import (
    RefinePlan,
    SubclusterOp,
    MergeOp,
    OverrideOp,
    RescoreOp,
)


@dataclass
class ManualPolicyConfig:
    """Configuration parsed from YAML file."""
    version: str = "1.0"
    iteration: int = 1
    overrides: List[Dict[str, Any]] = field(default_factory=list)
    subcluster: List[Dict[str, Any]] = field(default_factory=list)
    merge: List[Dict[str, Any]] = field(default_factory=list)
    rescoring: Dict[str, Any] = field(default_factory=dict)


class ManualPolicy:
    """Generates RefinePlan from YAML configuration.

    This policy parses a YAML configuration file and converts it to
    a RefinePlan with operations for manual expert curation.

    YAML Format:
    ```yaml
    version: 1.0
    iteration: 1

    overrides:
      - cluster_id: "3"
        cell_type: "Ciliated_Epithelium"
        reason: "High FOXJ1 expression"

    subcluster:
      - cluster_id: "5"
        resolution: 0.3
        focus_markers: [CD68, CD163, CD206]
        reason: "M1/M2 macrophage heterogeneity"

    merge:
      - source_clusters: ["17", "18"]
        target_label: "Rare_Immune"
        reason: "Low cell count"

    rescoring:
      mode: "smart"
      recompute_de: false
    ```

    Example:
        >>> policy = ManualPolicy.from_yaml(Path("curation.yaml"))
        >>> plan = policy.generate_plan(adata)
        >>> print(plan.summary())
        {'override': 2, 'subcluster': 1, 'merge': 1, 'rescore': 1}
    """

    def __init__(
        self,
        config: ManualPolicyConfig,
        config_path: Optional[Path] = None,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize ManualPolicy with parsed configuration.

        Parameters
        ----------
        config : ManualPolicyConfig
            Parsed configuration
        config_path : Path, optional
            Path to source config file (for provenance)
        logger : logging.Logger, optional
            Logger for tracking operations
        """
        self.config = config
        self.config_path = config_path
        self.logger = logger or logging.getLogger(__name__)

    @classmethod
    def from_yaml(cls, config_path: Path, logger: Optional[logging.Logger] = None) -> "ManualPolicy":
        """Create ManualPolicy from YAML configuration file.

        Parameters
        ----------
        config_path : Path
            Path to YAML config file
        logger : logging.Logger, optional
            Logger for tracking operations

        Returns
        -------
        ManualPolicy
            Policy instance with parsed configuration

        Raises
        ------
        FileNotFoundError
            If config file doesn't exist
        yaml.YAMLError
            If config is invalid YAML
        """
        _logger = logger or logging.getLogger(__name__)
        config_path = Path(config_path)

        if not config_path.exists():
            raise FileNotFoundError(
                f"Config file not found: {config_path}\n"
                "Use ManualPolicy.generate_template() to create a template."
            )

        _logger.info("Loading curation config: %s", config_path)

        with open(config_path, "r") as f:
            raw_config = yaml.safe_load(f) or {}

        config = ManualPolicyConfig(
            version=raw_config.get("version", "1.0"),
            iteration=raw_config.get("iteration", 1),
            overrides=raw_config.get("overrides", []),
            subcluster=raw_config.get("subcluster", []),
            merge=raw_config.get("merge", []),
            rescoring=raw_config.get("rescoring", {}),
        )

        _logger.info("Config loaded successfully")
        _logger.info("  Version: %s", config.version)
        _logger.info("  Iteration: %d", config.iteration)
        _logger.info("  Overrides: %d", len(config.overrides))
        _logger.info("  Subcluster rules: %d", len(config.subcluster))
        _logger.info("  Merge rules: %d", len(config.merge))

        return cls(config=config, config_path=config_path, logger=_logger)

    @classmethod
    def from_dict(cls, config_dict: Dict[str, Any], logger: Optional[logging.Logger] = None) -> "ManualPolicy":
        """Create ManualPolicy from dictionary.

        Parameters
        ----------
        config_dict : Dict[str, Any]
            Configuration dictionary
        logger : logging.Logger, optional
            Logger for tracking operations

        Returns
        -------
        ManualPolicy
            Policy instance
        """
        config = ManualPolicyConfig(
            version=config_dict.get("version", "1.0"),
            iteration=config_dict.get("iteration", 1),
            overrides=config_dict.get("overrides", []),
            subcluster=config_dict.get("subcluster", []),
            merge=config_dict.get("merge", []),
            rescoring=config_dict.get("rescoring", {}),
        )
        return cls(config=config, logger=logger)

    def generate_plan(self, adata: Optional["sc.AnnData"] = None) -> RefinePlan:
        """Generate RefinePlan from YAML configuration.

        Parameters
        ----------
        adata : sc.AnnData, optional
            AnnData object for validation (optional)

        Returns
        -------
        RefinePlan
            Plan with operations from YAML config
        """
        plan = RefinePlan(
            metadata={
                "policy": "ManualPolicy",
                "created_at": datetime.now().isoformat(),
                "config_version": self.config.version,
                "iteration": self.config.iteration,
                "config_path": str(self.config_path) if self.config_path else None,
            }
        )

        # Validate if adata provided
        if adata is not None:
            self._validate_config(adata)

        # Add override operations
        for rule in self.config.overrides:
            plan.add_override(
                cluster_id=str(rule["cluster_id"]),
                new_label=rule["cell_type"],
                reason=rule.get("reason", "Manual override"),
                source="manual",
            )

        # Add merge operations
        for rule in self.config.merge:
            plan.add_merge(
                source_clusters=[str(c) for c in rule["source_clusters"]],
                target_label=rule["target_label"],
                reason=rule.get("reason", "Manual merge"),
                source="manual",
            )

        # Add subcluster operations
        for rule in self.config.subcluster:
            plan.add_subcluster(
                cluster_id=str(rule["cluster_id"]),
                resolution=rule.get("resolution", 0.4),
                n_pcs=rule.get("n_pcs", 30),
                neighbors_k=rule.get("neighbors_k", 15),
                focus_markers=rule.get("focus_markers"),
                min_cells=rule.get("min_cells", 100),
                reason=rule.get("reason", "Manual subcluster"),
                source="manual",
            )

        # Add rescore operation if configured
        rescore_config = self.config.rescoring
        if rescore_config.get("mode", "none") != "none":
            plan.add_rescore(
                mode=rescore_config.get("mode", "smart"),
                recompute_de=rescore_config.get("recompute_de", False),
                reason="Manual rescoring after modifications",
            )

        # Sort operations in deterministic order
        plan.sort_operations()

        # Update metadata
        plan.metadata["n_overrides"] = len(self.config.overrides)
        plan.metadata["n_subclusters"] = len(self.config.subcluster)
        plan.metadata["n_merges"] = len(self.config.merge)

        return plan

    def _validate_config(self, adata: "sc.AnnData") -> None:
        """Validate configuration against AnnData object.

        Raises ValueError if config contains invalid cluster IDs.
        """
        self.logger.info("Validating config against input data...")

        # Get available cluster IDs (check both lvl0 and lvl1)
        available_lvl0 = set(adata.obs["cluster_lvl0"].astype(str).unique())
        available_lvl1 = set()
        if "cluster_lvl1" in adata.obs:
            available_lvl1 = set(adata.obs["cluster_lvl1"].dropna().astype(str).unique())
        all_available = available_lvl0 | available_lvl1

        # Validate override cluster IDs
        for i, rule in enumerate(self.config.overrides):
            if "cluster_id" not in rule:
                raise ValueError(f"Override rule {i} missing 'cluster_id' field")
            if "cell_type" not in rule:
                raise ValueError(f"Override rule {i} missing 'cell_type' field")

            cluster_id = str(rule["cluster_id"])
            if cluster_id not in all_available:
                raise ValueError(
                    f"Override rule {i}: cluster_id '{cluster_id}' not found in data. "
                    f"Available clusters: {sorted(all_available)}"
                )

        # Validate subcluster cluster IDs
        for i, rule in enumerate(self.config.subcluster):
            if "cluster_id" not in rule:
                raise ValueError(f"Subcluster rule {i} missing 'cluster_id' field")

            cluster_id = str(rule["cluster_id"])
            if cluster_id not in all_available:
                raise ValueError(
                    f"Subcluster rule {i}: cluster_id '{cluster_id}' not found in data. "
                    f"Available clusters: {sorted(all_available)}"
                )

            # Check min_cells if specified
            if "min_cells" in rule:
                n_cells = (adata.obs["cluster_lvl0"].astype(str) == cluster_id).sum()
                if "cluster_lvl1" in adata.obs:
                    n_cells = max(
                        n_cells,
                        (adata.obs["cluster_lvl1"].astype(str) == cluster_id).sum()
                    )
                if n_cells < rule["min_cells"]:
                    self.logger.warning(
                        "Subcluster rule %d: cluster %s has only %d cells (min_cells=%d)",
                        i, cluster_id, n_cells, rule["min_cells"],
                    )

        # Validate merge cluster IDs
        for i, rule in enumerate(self.config.merge):
            if "source_clusters" not in rule:
                raise ValueError(f"Merge rule {i} missing 'source_clusters' field")
            if "target_label" not in rule:
                raise ValueError(f"Merge rule {i} missing 'target_label' field")

            for cluster_id in rule["source_clusters"]:
                cluster_id = str(cluster_id)
                if cluster_id not in all_available:
                    raise ValueError(
                        f"Merge rule {i}: cluster_id '{cluster_id}' not found in data. "
                        f"Available clusters: {sorted(all_available)}"
                    )

        self.logger.info("Config validation passed")

    @staticmethod
    def generate_template(
        adata: Optional["sc.AnnData"] = None,
        output_path: Optional[Path] = None,
    ) -> str:
        """Generate a template YAML configuration file.

        Parameters
        ----------
        adata : sc.AnnData, optional
            If provided, template includes actual cluster IDs
        output_path : Path, optional
            If provided, write template to this file

        Returns
        -------
        str
            YAML template content
        """
        template_lines = [
            "# Manual curation configuration for cell-type refinement",
            "# Generated by ManualPolicy.generate_template()",
            "",
            "version: 1.0",
            "iteration: 1",
            "",
            "# Override rules: Directly assign cell types to clusters",
            "overrides:",
        ]

        if adata is not None and "cluster_lvl0" in adata.obs:
            clusters = sorted(adata.obs["cluster_lvl0"].astype(str).unique())[:3]
            for cid in clusters:
                template_lines.extend([
                    f"  - cluster_id: \"{cid}\"",
                    "    cell_type: \"Cell_Type_Name\"  # TODO: Set correct cell type",
                    "    reason: \"TODO: Add reasoning\"",
                    "",
                ])
        else:
            template_lines.extend([
                "  - cluster_id: \"0\"",
                "    cell_type: \"Cell_Type_Name\"",
                "    reason: \"Expert reasoning for this assignment\"",
                "",
            ])

        template_lines.extend([
            "# Subcluster rules: Re-cluster specific clusters with custom parameters",
            "subcluster:",
            "  - cluster_id: \"5\"",
            "    resolution: 0.3  # Lower = fewer subclusters",
            "    min_cells: 200",
            "    # focus_markers: [\"CD68\", \"CD163\"]  # Optional: restrict to specific markers",
            "    reason: \"Heterogeneous population needs finer resolution\"",
            "",
            "# Merge rules: Combine multiple clusters into one",
            "merge:",
            "  - source_clusters: [\"10\", \"11\"]",
            "    target_label: \"Merged_Population\"",
            "    reason: \"Similar populations with low cell counts\"",
            "",
            "# Rescoring options (optional)",
            "rescoring:",
            "  mode: \"none\"  # Options: none, smart, full",
            "  recompute_de: false",
        ])

        template = "\n".join(template_lines)

        if output_path is not None:
            output_path = Path(output_path)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            with open(output_path, "w") as f:
                f.write(template)

        return template

    def get_affected_clusters(self) -> Set[str]:
        """Get all cluster IDs affected by this policy's rules."""
        clusters = set()

        for rule in self.config.overrides:
            clusters.add(str(rule["cluster_id"]))

        for rule in self.config.subcluster:
            clusters.add(str(rule["cluster_id"]))

        for rule in self.config.merge:
            clusters.update(str(c) for c in rule["source_clusters"])

        return clusters

    def __repr__(self) -> str:
        return (f"ManualPolicy(overrides={len(self.config.overrides)}, "
                f"subclusters={len(self.config.subcluster)}, "
                f"merges={len(self.config.merge)})")
