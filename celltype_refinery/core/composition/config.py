"""Configuration for composition analysis.

This module provides tissue-agnostic configuration for cell-type composition
analysis. Patterns and region ordering are loaded from YAML configuration
rather than being hardcoded.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List

import yaml


@dataclass
class PatternConfig:
    """Pattern-based cell type classification.

    Patterns are used to classify cell types into biological categories
    (e.g., epithelial, stromal, immune). Each pattern is a list of
    substrings that match cell type labels.

    Attributes
    ----------
    epithelial : List[str]
        Patterns matching epithelial cell types
    stromal : List[str]
        Patterns matching stromal cell types
    smooth_muscle : List[str]
        Patterns matching smooth muscle cell types
    immune : List[str]
        Patterns matching immune cell types
    endothelial : List[str]
        Patterns matching endothelial cell types
    ciliated : List[str]
        Patterns matching ciliated epithelial cell types (FT-specific)
    secretory : List[str]
        Patterns matching secretory epithelial cell types (FT-specific)
    custom : Dict[str, List[str]]
        Additional custom pattern categories
    """

    epithelial: List[str] = field(default_factory=lambda: [
        "Epithelium", "Epithelial", "Ciliated", "Secretory"
    ])
    stromal: List[str] = field(default_factory=lambda: [
        "Stromal", "Stroma", "Fibroblast"
    ])
    smooth_muscle: List[str] = field(default_factory=lambda: [
        "Smooth", "SMA", "muscle"
    ])
    immune: List[str] = field(default_factory=lambda: [
        "CD45", "Immune", "Myeloid", "Lymphocyte", "Macrophage", "T_cell", "B_cell"
    ])
    endothelial: List[str] = field(default_factory=lambda: [
        "Endothel", "Vascular", "Lymphatic"
    ])
    ciliated: List[str] = field(default_factory=lambda: [
        "Ciliated", "FOXJ1"
    ])
    secretory: List[str] = field(default_factory=lambda: [
        "Secretory", "PAX8", "Glandular", "Peg"
    ])
    custom: Dict[str, List[str]] = field(default_factory=dict)

    def get_all_categories(self) -> Dict[str, List[str]]:
        """Get all pattern categories including custom ones.

        Returns
        -------
        Dict[str, List[str]]
            Mapping of category name to pattern list
        """
        categories = {
            "epithelial": self.epithelial,
            "stromal": self.stromal,
            "smooth_muscle": self.smooth_muscle,
            "immune": self.immune,
            "endothelial": self.endothelial,
            "ciliated": self.ciliated,
            "secretory": self.secretory,
        }
        categories.update(self.custom)
        return categories

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "PatternConfig":
        """Create from dictionary.

        Parameters
        ----------
        data : Dict[str, Any]
            Dictionary with pattern configuration

        Returns
        -------
        PatternConfig
            Configuration instance
        """
        standard_keys = {"epithelial", "stromal", "smooth_muscle", "immune", "endothelial", "ciliated", "secretory"}
        standard = {k: v for k, v in data.items() if k in standard_keys}
        custom = {k: v for k, v in data.items() if k not in standard_keys and k != "custom"}

        if "custom" in data:
            custom.update(data["custom"])

        return cls(**standard, custom=custom)


@dataclass
class EnrichmentConfig:
    """Configuration for regional enrichment tests.

    Attributes
    ----------
    min_samples_per_region : int
        Minimum samples per region for valid test
    correction_method : str
        Multiple testing correction method (bonferroni, fdr_bh, holm)
    alpha : float
        Significance threshold
    min_cells_per_type : int
        Minimum cells per type to include in test
    """

    min_samples_per_region: int = 3
    correction_method: str = "fdr_bh"
    alpha: float = 0.05
    min_cells_per_type: int = 10

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "EnrichmentConfig":
        """Create from dictionary."""
        return cls(**{k: v for k, v in data.items() if k in cls.__dataclass_fields__})


@dataclass
class DiversityConfig:
    """Configuration for diversity metrics.

    Attributes
    ----------
    min_cells_per_sample : int
        Minimum cells per sample for valid diversity calculation
    rare_threshold : float
        Proportion threshold for rare cell type classification
    """

    min_cells_per_sample: int = 100
    rare_threshold: float = 0.01

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "DiversityConfig":
        """Create from dictionary."""
        return cls(**{k: v for k, v in data.items() if k in cls.__dataclass_fields__})


@dataclass
class VisualizationConfig:
    """Configuration for composition visualizations.

    Attributes
    ----------
    dpi : int
        Figure resolution
    top_n : int
        Number of top cell types to show in plots
    figsize_composition : tuple
        Figure size for composition plots (width, height)
    figsize_heatmap : tuple
        Figure size for heatmaps (width, height)
    cmap : str
        Colormap for heatmaps
    """

    dpi: int = 200
    top_n: int = 15
    figsize_composition: tuple = (12, 8)
    figsize_heatmap: tuple = (14, 10)
    cmap: str = "RdBu_r"

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "VisualizationConfig":
        """Create from dictionary."""
        config_data = {}
        for k, v in data.items():
            if k in cls.__dataclass_fields__:
                if k.startswith("figsize") and isinstance(v, list):
                    config_data[k] = tuple(v)
                else:
                    config_data[k] = v
        return cls(**config_data)


@dataclass
class CompositionConfig:
    """Main configuration for composition analysis.

    Attributes
    ----------
    patterns : PatternConfig
        Cell type pattern classification
    enrichment : EnrichmentConfig
        Regional enrichment test settings
    diversity : DiversityConfig
        Diversity metric settings
    visualization : VisualizationConfig
        Visualization settings
    region_column : str
        Column name for region in adata.obs
    sample_column : str
        Column name for sample in adata.obs
    donor_column : str
        Column name for donor in adata.obs
    region_order : List[str]
        Ordered list of regions for consistent output
    region_colors : Dict[str, str]
        Hex colors for regions (used in visualizations)
    cell_type_columns : List[str]
        Cell type columns to analyze
    skip_enrichment : bool
        Skip enrichment tests
    skip_biology : bool
        Skip biology metrics
    skip_viz : bool
        Skip visualizations
    """

    patterns: PatternConfig = field(default_factory=PatternConfig)
    enrichment: EnrichmentConfig = field(default_factory=EnrichmentConfig)
    diversity: DiversityConfig = field(default_factory=DiversityConfig)
    visualization: VisualizationConfig = field(default_factory=VisualizationConfig)
    region_column: str = "region"
    sample_column: str = "sample_id"
    donor_column: str = "donor"
    region_order: List[str] = field(default_factory=list)
    region_colors: Dict[str, str] = field(default_factory=dict)
    cell_type_columns: List[str] = field(default_factory=lambda: ["cell_type"])
    skip_enrichment: bool = False
    skip_biology: bool = False
    skip_viz: bool = False

    @classmethod
    def from_yaml(cls, path: Path) -> "CompositionConfig":
        """Load configuration from YAML file.

        Parameters
        ----------
        path : Path
            Path to YAML configuration file

        Returns
        -------
        CompositionConfig
            Configuration instance
        """
        with open(path) as f:
            data = yaml.safe_load(f)
        return cls.from_dict(data)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "CompositionConfig":
        """Create from dictionary.

        Parameters
        ----------
        data : Dict[str, Any]
            Configuration dictionary

        Returns
        -------
        CompositionConfig
            Configuration instance
        """
        patterns = PatternConfig.from_dict(data.get("patterns", {}))
        enrichment = EnrichmentConfig.from_dict(data.get("enrichment", {}))
        diversity = DiversityConfig.from_dict(data.get("diversity", {}))
        visualization = VisualizationConfig.from_dict(data.get("visualization", {}))

        return cls(
            patterns=patterns,
            enrichment=enrichment,
            diversity=diversity,
            visualization=visualization,
            region_column=data.get("region_column", "region"),
            sample_column=data.get("sample_column", "sample_id"),
            donor_column=data.get("donor_column", "donor"),
            region_order=data.get("region_order", []),
            region_colors=data.get("region_colors", {}),
            cell_type_columns=data.get("cell_type_columns", ["cell_type"]),
            skip_enrichment=data.get("skip_enrichment", False),
            skip_biology=data.get("skip_biology", False),
            skip_viz=data.get("skip_viz", False),
        )

    @classmethod
    def default(cls) -> "CompositionConfig":
        """Create default configuration.

        Returns
        -------
        CompositionConfig
            Default configuration instance
        """
        return cls()

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary.

        Returns
        -------
        Dict[str, Any]
            Configuration as dictionary
        """
        return {
            "patterns": {
                "epithelial": self.patterns.epithelial,
                "stromal": self.patterns.stromal,
                "smooth_muscle": self.patterns.smooth_muscle,
                "immune": self.patterns.immune,
                "endothelial": self.patterns.endothelial,
                "custom": self.patterns.custom,
            },
            "enrichment": {
                "min_samples_per_region": self.enrichment.min_samples_per_region,
                "correction_method": self.enrichment.correction_method,
                "alpha": self.enrichment.alpha,
                "min_cells_per_type": self.enrichment.min_cells_per_type,
            },
            "diversity": {
                "min_cells_per_sample": self.diversity.min_cells_per_sample,
                "rare_threshold": self.diversity.rare_threshold,
            },
            "visualization": {
                "dpi": self.visualization.dpi,
                "top_n": self.visualization.top_n,
                "figsize_composition": list(self.visualization.figsize_composition),
                "figsize_heatmap": list(self.visualization.figsize_heatmap),
                "cmap": self.visualization.cmap,
            },
            "region_column": self.region_column,
            "sample_column": self.sample_column,
            "donor_column": self.donor_column,
            "region_order": self.region_order,
            "region_colors": self.region_colors,
            "cell_type_columns": self.cell_type_columns,
            "skip_enrichment": self.skip_enrichment,
            "skip_biology": self.skip_biology,
            "skip_viz": self.skip_viz,
        }

    def to_yaml(self, path: Path) -> None:
        """Save configuration to YAML file.

        Parameters
        ----------
        path : Path
            Output path for YAML file
        """
        with open(path, "w") as f:
            yaml.dump(self.to_dict(), f, default_flow_style=False, sort_keys=False)

    def apply_organ_config(self, organ_config: "OrganConfig") -> None:
        """Apply organ-specific configuration.

        Updates region_order and region_colors from the organ configuration.
        Only updates fields that are empty (doesn't override explicit settings).

        Parameters
        ----------
        organ_config : OrganConfig
            Organ configuration (from celltype_refinery.config)
        """
        # Only apply if current values are empty (don't override explicit config)
        if not self.region_order and organ_config.region_order:
            self.region_order = organ_config.region_order.copy()
        if not self.region_colors and organ_config.region_colors:
            self.region_colors = organ_config.region_colors.copy()


# Type hint for OrganConfig (avoid circular import)
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from celltype_refinery.config import OrganConfig
