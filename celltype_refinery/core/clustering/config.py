"""Configuration classes for clustering module.

All clustering parameters are configurable for tissue-agnostic operation.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml


@dataclass
class ClusteringConfig:
    """Configuration for Leiden clustering.

    Attributes
    ----------
    n_pcs : int
        Number of principal components for neighbors
    neighbors_k : int
        k for the neighborhood graph
    resolution : float
        Leiden resolution for clustering
    scale_clip : float
        Value clipping during scaling
    random_seed : int
        Random seed for reproducibility
    use_gpu : bool
        Use GPU acceleration if available
    compute_umap : bool
        Compute UMAP embeddings for visualization
    min_marker_std : float
        Drop markers with std below this value
    """

    n_pcs: int = 30
    neighbors_k: int = 15
    resolution: float = 0.6
    scale_clip: float = 10.0
    random_seed: int = 1337
    use_gpu: bool = True
    compute_umap: bool = False
    min_marker_std: float = 1e-3


@dataclass
class DEConfig:
    """Configuration for differential expression analysis.

    Attributes
    ----------
    method : str
        DE method (wilcoxon, t-test, etc.)
    n_genes : int
        Number of top genes to return per cluster
    layer : str
        Layer to use for DE analysis
    tie_correct : bool
        Apply tie correction for Wilcoxon test
    de_bonus : float
        Score bonus for curated markers in DE top list
    top_frac : float
        Fraction of panel for top-K calculation
    min_k : int
        Minimum K for DE top-K
    commonness_alpha : float
        Exponent for marker commonness penalty
    """

    method: str = "wilcoxon"
    n_genes: int = 12
    layer: str = "batchcorr"
    tie_correct: bool = True
    de_bonus: float = 0.5
    top_frac: float = 0.2
    min_k: int = 3
    commonness_alpha: float = 0.5


@dataclass
class SubclusterConfig:
    """Configuration for subclustering.

    Attributes
    ----------
    min_cells : int
        Minimum cluster size for subclustering
    score_threshold : float
        Clusters below this score are considered for subclustering
    resolution : float
        Leiden resolution for subclustering
    max_iterations : int
        Maximum subclustering iterations
    """

    min_cells: int = 500
    score_threshold: float = 1.0
    resolution: float = 0.4
    max_iterations: int = 3


@dataclass
class ScoringConfig:
    """Configuration for marker scoring.

    Attributes
    ----------
    positive_quantile : float
        Quantile used to define marker positivity
    anti_marker_weight : float
        Weight for anti-marker penalty
    use_idf_weights : bool
        Use IDF-based marker weighting
    """

    positive_quantile: float = 0.75
    anti_marker_weight: float = 0.5
    use_idf_weights: bool = False


@dataclass
class StageHConfig:
    """Master configuration for Stage H clustering and annotation.

    Attributes
    ----------
    clustering : ClusteringConfig
        Clustering configuration
    de : DEConfig
        Differential expression configuration
    subcluster : SubclusterConfig
        Subclustering configuration
    scoring : ScoringConfig
        Marker scoring configuration
    technical_markers : List[str]
        Markers to exclude from analysis
    """

    clustering: ClusteringConfig = field(default_factory=ClusteringConfig)
    de: DEConfig = field(default_factory=DEConfig)
    subcluster: SubclusterConfig = field(default_factory=SubclusterConfig)
    scoring: ScoringConfig = field(default_factory=ScoringConfig)
    technical_markers: List[str] = field(
        default_factory=lambda: ["DAPI", "Collagen IV", "Beta-actin"]
    )

    @classmethod
    def from_yaml(cls, path: Path) -> "StageHConfig":
        """Load configuration from YAML file."""
        with open(path) as f:
            data = yaml.safe_load(f) or {}

        # Handle nested stage_h section
        if "stage_h" in data:
            data = data["stage_h"]

        return cls(
            clustering=ClusteringConfig(**data.get("clustering", {})),
            de=DEConfig(**data.get("de", {})),
            subcluster=SubclusterConfig(**data.get("subcluster", {})),
            scoring=ScoringConfig(**data.get("scoring", {})),
            technical_markers=data.get(
                "technical_markers", ["DAPI", "Collagen IV", "Beta-actin"]
            ),
        )

    @classmethod
    def default(cls) -> "StageHConfig":
        """Create default configuration."""
        return cls()

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            "clustering": {
                "n_pcs": self.clustering.n_pcs,
                "neighbors_k": self.clustering.neighbors_k,
                "resolution": self.clustering.resolution,
                "scale_clip": self.clustering.scale_clip,
                "random_seed": self.clustering.random_seed,
                "use_gpu": self.clustering.use_gpu,
                "compute_umap": self.clustering.compute_umap,
                "min_marker_std": self.clustering.min_marker_std,
            },
            "de": {
                "method": self.de.method,
                "n_genes": self.de.n_genes,
                "layer": self.de.layer,
                "tie_correct": self.de.tie_correct,
                "de_bonus": self.de.de_bonus,
            },
            "subcluster": {
                "min_cells": self.subcluster.min_cells,
                "score_threshold": self.subcluster.score_threshold,
                "resolution": self.subcluster.resolution,
            },
            "scoring": {
                "positive_quantile": self.scoring.positive_quantile,
                "anti_marker_weight": self.scoring.anti_marker_weight,
                "use_idf_weights": self.scoring.use_idf_weights,
            },
            "technical_markers": self.technical_markers,
        }
