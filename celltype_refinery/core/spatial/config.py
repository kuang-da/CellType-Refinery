"""Configuration for spatial analysis.

Provides dataclasses for configuring:
- Permutation test settings (parallelization, number of permutations)
- Multiple testing correction (FDR, Bonferroni, Holm)
- Visualization settings
- Column mappings
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
import logging

import yaml

logger = logging.getLogger(__name__)


@dataclass
class PermutationConfig:
    """Configuration for permutation tests.

    Attributes
    ----------
    n_permutations : int
        Number of permutations for null model
    seed : int
        Random seed for reproducibility
    n_threads : int
        Number of threads for parallel permutations
    batch_size : int
        Batch size for joblib parallel processing
    large_sample_threshold : int
        Samples with more cells than this use parallel permutations
    """

    n_permutations: int = 100
    seed: int = 42
    n_threads: int = 8
    batch_size: int = 16
    large_sample_threshold: int = 50000

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "PermutationConfig":
        """Create PermutationConfig from dictionary."""
        return cls(
            n_permutations=data.get("n_permutations", 100),
            seed=data.get("seed", 42),
            n_threads=data.get("n_threads", 8),
            batch_size=data.get("batch_size", 16),
            large_sample_threshold=data.get("large_sample_threshold", 50000),
        )


@dataclass
class CorrectionConfig:
    """Configuration for multiple testing correction.

    Attributes
    ----------
    method : str
        Correction method: "fdr_bh", "bonferroni", "holm", "none"
    alpha : float
        Significance threshold
    """

    method: str = "fdr_bh"
    alpha: float = 0.05

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "CorrectionConfig":
        """Create CorrectionConfig from dictionary."""
        return cls(
            method=data.get("method", "fdr_bh"),
            alpha=data.get("alpha", 0.05),
        )


@dataclass
class VisualizationConfig:
    """Configuration for visualization settings."""

    dpi: int = 200
    figsize_heatmap: Tuple[int, int] = (12, 10)
    figsize_network: Tuple[int, int] = (12, 12)
    figsize_bar: Tuple[int, int] = (12, 6)
    cmap_enrichment: str = "RdBu_r"
    z_score_vmin: float = -5.0
    z_score_vmax: float = 5.0
    interaction_threshold: float = 0.5

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "VisualizationConfig":
        """Create VisualizationConfig from dictionary."""
        figsize_heatmap = data.get("figsize_heatmap", [12, 10])
        figsize_network = data.get("figsize_network", [12, 12])
        figsize_bar = data.get("figsize_bar", [12, 6])
        if isinstance(figsize_heatmap, list):
            figsize_heatmap = tuple(figsize_heatmap)
        if isinstance(figsize_network, list):
            figsize_network = tuple(figsize_network)
        if isinstance(figsize_bar, list):
            figsize_bar = tuple(figsize_bar)
        return cls(
            dpi=data.get("dpi", 200),
            figsize_heatmap=figsize_heatmap,
            figsize_network=figsize_network,
            figsize_bar=figsize_bar,
            cmap_enrichment=data.get("cmap_enrichment", "RdBu_r"),
            z_score_vmin=data.get("z_score_vmin", -5.0),
            z_score_vmax=data.get("z_score_vmax", 5.0),
            interaction_threshold=data.get("interaction_threshold", 0.5),
        )


@dataclass
class SpatialConfig:
    """Main configuration for spatial analysis.

    Attributes
    ----------
    version : str
        Configuration version
    description : str
        Optional description
    cell_type_col : str
        Default column for cell type labels
    cell_type_cols : List[str]
        Columns for multi-column analysis
    sample_col : str
        Column for sample IDs
    region_col : str
        Column for anatomical regions
    graphs_dir : Path, optional
        Directory containing spatial graphs
    min_cells : int
        Minimum cells per sample to include
    min_cells_per_type : int
        Minimum cells per cell type to include
    permutation : PermutationConfig
        Permutation test settings
    correction : CorrectionConfig
        Multiple testing correction settings
    visualization : VisualizationConfig
        Visualization settings
    """

    version: str = "1.0"
    description: str = ""

    # Column defaults
    cell_type_col: str = "cell_type"
    cell_type_cols: List[str] = field(default_factory=lambda: ["cell_type"])
    sample_col: str = "sample_id"
    region_col: str = "region"

    # Graph settings
    graphs_dir: Optional[Path] = None

    # Thresholds
    min_cells: int = 100
    min_cells_per_type: int = 50

    # Sub-configs
    permutation: PermutationConfig = field(default_factory=PermutationConfig)
    correction: CorrectionConfig = field(default_factory=CorrectionConfig)
    visualization: VisualizationConfig = field(default_factory=VisualizationConfig)

    @classmethod
    def from_yaml(cls, path: Path) -> "SpatialConfig":
        """Load configuration from YAML file."""
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"Config file not found: {path}")

        with open(path) as f:
            data = yaml.safe_load(f)

        logger.info(f"Loaded config from {path}")

        # Parse nested configs
        permutation = PermutationConfig()
        if "permutation" in data:
            permutation = PermutationConfig.from_dict(data["permutation"])

        correction = CorrectionConfig()
        if "correction" in data or "enrichment" in data:
            corr_data = data.get("correction", data.get("enrichment", {}))
            correction = CorrectionConfig.from_dict(corr_data)

        visualization = VisualizationConfig()
        if "visualization" in data:
            visualization = VisualizationConfig.from_dict(data["visualization"])

        # Column mappings
        columns = data.get("columns", {})
        cell_type_cols = columns.get("cell_type_cols", ["cell_type"])

        # Thresholds
        thresholds = data.get("thresholds", {})

        # Graphs dir
        graphs_dir = data.get("graphs", {}).get("dir", None)
        if graphs_dir:
            graphs_dir = Path(graphs_dir)

        return cls(
            version=data.get("version", "1.0"),
            description=data.get("description", ""),
            permutation=permutation,
            correction=correction,
            visualization=visualization,
            cell_type_col=columns.get("cell_type", "cell_type"),
            cell_type_cols=cell_type_cols,
            sample_col=columns.get("sample", "sample_id"),
            region_col=columns.get("region", "region"),
            graphs_dir=graphs_dir,
            min_cells=thresholds.get("min_cells", 100),
            min_cells_per_type=thresholds.get("min_cells_per_type", 50),
        )

    @classmethod
    def default(cls) -> "SpatialConfig":
        """Return default configuration."""
        return cls()

    def to_dict(self) -> Dict[str, Any]:
        """Convert configuration to dictionary."""
        return {
            "version": self.version,
            "description": self.description,
            "columns": {
                "cell_type": self.cell_type_col,
                "cell_type_cols": self.cell_type_cols,
                "sample": self.sample_col,
                "region": self.region_col,
            },
            "graphs": {
                "dir": str(self.graphs_dir) if self.graphs_dir else None,
            },
            "thresholds": {
                "min_cells": self.min_cells,
                "min_cells_per_type": self.min_cells_per_type,
            },
            "permutation": {
                "n_permutations": self.permutation.n_permutations,
                "seed": self.permutation.seed,
                "n_threads": self.permutation.n_threads,
                "batch_size": self.permutation.batch_size,
                "large_sample_threshold": self.permutation.large_sample_threshold,
            },
            "correction": {
                "method": self.correction.method,
                "alpha": self.correction.alpha,
            },
            "visualization": {
                "dpi": self.visualization.dpi,
                "figsize_heatmap": list(self.visualization.figsize_heatmap),
                "figsize_network": list(self.visualization.figsize_network),
                "figsize_bar": list(self.visualization.figsize_bar),
                "cmap_enrichment": self.visualization.cmap_enrichment,
                "z_score_vmin": self.visualization.z_score_vmin,
                "z_score_vmax": self.visualization.z_score_vmax,
                "interaction_threshold": self.visualization.interaction_threshold,
            },
        }
