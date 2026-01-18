"""Centralized organ-specific configuration.

This module provides a centralized configuration system for organ-specific
settings like region ordering and colors. All modules (composition, annotation,
preprocessing, spatial) should use this single source of truth.

Example
-------
>>> from celltype_refinery.config import get_organ_config
>>> config = get_organ_config("fallopian_tube")
>>> sorted_regions = config.sort_regions(["isthmus", "fimbriae"])
>>> print(sorted_regions)
['fimbriae', 'isthmus']
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import yaml


@dataclass
class OrganConfig:
    """Configuration for an organ type.

    Attributes
    ----------
    organ_name : str
        Canonical organ name (lowercase, underscores)
    region_order : List[str]
        Anatomical region ordering for visualizations
    region_colors : Dict[str, str]
        Hex colors for each region
    aliases : List[str]
        Alternative names for this organ
    patterns : Dict[str, List[str]]
        Cell type classification patterns (optional)
    """

    organ_name: str
    region_order: List[str] = field(default_factory=list)
    region_colors: Dict[str, str] = field(default_factory=dict)
    aliases: List[str] = field(default_factory=list)
    patterns: Dict[str, List[str]] = field(default_factory=dict)

    def sort_regions(self, regions: Iterable[str]) -> List[str]:
        """Sort regions by canonical anatomical order.

        Parameters
        ----------
        regions : Iterable[str]
            Regions to sort

        Returns
        -------
        List[str]
            Sorted regions. Unknown regions are placed at the end alphabetically.
        """
        if not self.region_order:
            # No canonical order defined, return alphabetical
            return sorted(regions, key=lambda r: str(r).lower())

        # Build order map (case-insensitive)
        order_map = {r.lower(): i for i, r in enumerate(self.region_order)}

        def get_order(region: str) -> tuple:
            r_lower = str(region).lower()
            if r_lower in order_map:
                return (0, order_map[r_lower], r_lower)
            # Unknown regions: sort alphabetically after known regions
            return (1, 0, r_lower)

        return sorted(regions, key=get_order)

    def get_region_color(self, region: str, default: str = "#808080") -> str:
        """Get color for a region (case-insensitive).

        Parameters
        ----------
        region : str
            Region name
        default : str
            Default color if region not found

        Returns
        -------
        str
            Hex color code
        """
        if not self.region_colors:
            return default

        # Try exact match first
        if region in self.region_colors:
            return self.region_colors[region]

        # Try case-insensitive match
        region_lower = region.lower()
        for key, color in self.region_colors.items():
            if key.lower() == region_lower:
                return color

        return default

    def get_region_colors_ordered(self) -> Dict[str, str]:
        """Get region colors in canonical order.

        Returns
        -------
        Dict[str, str]
            Ordered dict of region -> color
        """
        if not self.region_order:
            return dict(self.region_colors)

        ordered = {}
        for region in self.region_order:
            color = self.get_region_color(region)
            if color != "#808080":  # Only include if has defined color
                ordered[region] = color
        return ordered

    @classmethod
    def from_yaml(cls, path: Path) -> "OrganConfig":
        """Load organ config from YAML file.

        Parameters
        ----------
        path : Path
            Path to YAML file

        Returns
        -------
        OrganConfig
            Loaded configuration
        """
        with open(path) as f:
            data = yaml.safe_load(f)

        return cls(
            organ_name=data.get("organ", ""),
            region_order=data.get("regions", {}).get("order", []),
            region_colors=data.get("regions", {}).get("colors", {}),
            aliases=data.get("aliases", []),
            patterns=data.get("patterns", {}),
        )


# =============================================================================
# Registry
# =============================================================================

# Global registry of organ configurations
ORGAN_CONFIG_REGISTRY: Dict[str, OrganConfig] = {}

# Aliases map organ aliases to canonical names
ORGAN_ALIASES: Dict[str, str] = {}

# Flag to track if builtin configs have been loaded
_BUILTINS_LOADED = False


def register_organ_config(config: OrganConfig) -> None:
    """Register an organ configuration.

    Parameters
    ----------
    config : OrganConfig
        Configuration to register
    """
    organ_name = config.organ_name.lower().replace(" ", "_").replace("-", "_")
    ORGAN_CONFIG_REGISTRY[organ_name] = config

    # Register aliases
    for alias in config.aliases:
        ORGAN_ALIASES[alias.lower()] = organ_name


def get_organ_config(organ: Optional[str] = None) -> OrganConfig:
    """Get organ configuration by name or alias.

    Parameters
    ----------
    organ : str, optional
        Organ name or alias. If None or "generic", returns generic config.

    Returns
    -------
    OrganConfig
        Organ configuration

    Raises
    ------
    ValueError
        If organ is not registered
    """
    _ensure_builtins_loaded()

    if organ is None:
        return _get_generic_config()

    # Normalize name
    organ_lower = organ.lower().replace(" ", "_").replace("-", "_")

    # Check for "generic"
    if organ_lower == "generic":
        return _get_generic_config()

    # Check aliases
    if organ_lower in ORGAN_ALIASES:
        organ_lower = ORGAN_ALIASES[organ_lower]

    # Look up in registry
    if organ_lower not in ORGAN_CONFIG_REGISTRY:
        available = list(ORGAN_CONFIG_REGISTRY.keys())
        alias_info = [f"{k} -> {v}" for k, v in ORGAN_ALIASES.items()]
        raise ValueError(
            f"Unknown organ: '{organ}'. "
            f"Available: {available}. "
            f"Aliases: {alias_info}"
        )

    return ORGAN_CONFIG_REGISTRY[organ_lower]


def list_available_organs() -> List[str]:
    """List all available organ names.

    Returns
    -------
    List[str]
        Registered organ names (canonical)
    """
    _ensure_builtins_loaded()
    return sorted(ORGAN_CONFIG_REGISTRY.keys())


def list_organ_aliases() -> Dict[str, str]:
    """List all organ aliases.

    Returns
    -------
    Dict[str, str]
        Map of alias -> canonical name
    """
    _ensure_builtins_loaded()
    return dict(ORGAN_ALIASES)


def _get_generic_config() -> OrganConfig:
    """Get generic config (no region ordering, alphabetical fallback)."""
    if "generic" not in ORGAN_CONFIG_REGISTRY:
        ORGAN_CONFIG_REGISTRY["generic"] = OrganConfig(
            organ_name="generic",
            region_order=[],
            region_colors={},
            aliases=[],
            patterns={},
        )
    return ORGAN_CONFIG_REGISTRY["generic"]


def _ensure_builtins_loaded() -> None:
    """Ensure builtin configs are loaded."""
    global _BUILTINS_LOADED
    if _BUILTINS_LOADED:
        return

    _load_builtin_configs()
    _BUILTINS_LOADED = True


def _load_builtin_configs() -> None:
    """Load builtin organ configs from configs/tissues/*.yaml."""
    # Find configs directory relative to this module
    # Expected location: CellType-Refinery/configs/tissues/
    module_dir = Path(__file__).parent
    package_root = module_dir.parent.parent  # celltype_refinery -> CellType-Refinery
    configs_dir = package_root / "configs" / "tissues"

    if not configs_dir.exists():
        return

    for yaml_path in configs_dir.glob("*.yaml"):
        # Skip template files
        if "template" in yaml_path.name.lower():
            continue

        try:
            config = OrganConfig.from_yaml(yaml_path)
            if config.organ_name:
                register_organ_config(config)
        except Exception as e:
            # Log warning but don't fail
            import warnings
            warnings.warn(f"Failed to load organ config from {yaml_path}: {e}")


