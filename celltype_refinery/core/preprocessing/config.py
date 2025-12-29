"""Configuration classes for preprocessing stages.

All preprocessing parameters are configurable via YAML for
tissue-agnostic operation.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import yaml


@dataclass
class LoaderConfig:
    """Configuration for data loading (Stage A).

    Attributes
    ----------
    cell_id_col : str
        Column name for cell identifiers
    sample_id_col : str
        Column name for sample identifiers
    required_metadata_cols : List[str]
        Required metadata columns
    coordinate_cols : Tuple[str, str]
        Column names for x, y coordinates
    exclude_from_markers : List[str]
        Columns to exclude from marker list (in addition to cell_id_col)
    """

    cell_id_col: str = "cell_mask_id"
    sample_id_col: str = "sample_id"
    required_metadata_cols: List[str] = field(
        default_factory=lambda: ["sample_id", "cell_matrix_path", "cell_metadata_path"]
    )
    coordinate_cols: Tuple[str, str] = ("x", "y")
    coordinate_col_alternatives: List[Tuple[str, str]] = field(
        default_factory=lambda: [
            ("cell_centroid_x", "cell_centroid_y"),
            ("patch_x", "patch_y"),
        ]
    )
    exclude_from_markers: List[str] = field(
        default_factory=lambda: ["patch_id", "global_cell_id"]
    )


@dataclass
class QCConfig:
    """Configuration for cell QC (Stage B).

    Attributes
    ----------
    area_percentile_low : float
        Lower percentile for area filtering
    area_percentile_high : float
        Upper percentile for area filtering
    nucleus_ratio_min : float
        Minimum nucleus-to-cell area ratio
    nucleus_ratio_max : float
        Maximum nucleus-to-cell area ratio
    intensity_percentile_low : float
        Lower percentile for intensity filtering
    autofluorescence_percentile_high : float
        Upper percentile for autofluorescence filtering
    max_removal_fraction : float
        Maximum fraction of cells to remove per sample
    area_col : str
        Column name for cell area
    nucleus_area_col : str
        Column name for nucleus area
    exclude_from_markers : List[str]
        Columns to exclude from marker list (in addition to cell_id_col)
    """

    area_percentile_low: float = 1.0
    area_percentile_high: float = 99.0
    nucleus_ratio_min: float = 0.1
    nucleus_ratio_max: float = 0.9
    intensity_percentile_low: float = 1.0
    autofluorescence_percentile_high: float = 99.0
    max_removal_fraction: float = 0.15
    area_col: str = "cell_area"
    nucleus_area_col: str = "nucleus_area"
    exclude_from_markers: List[str] = field(
        default_factory=lambda: ["patch_id", "global_cell_id"]
    )


@dataclass
class NormalizationConfig:
    """Configuration for within-sample normalization (Stage C).

    Attributes
    ----------
    mode : str
        Normalization mode (explore or export)
    p_candidates : List[int]
        Percentile candidates for background estimation
    bg_modes : List[str]
        Background correction modes (clip, center)
    transforms : List[str]
        Variance stabilization transforms
    default_bg_mode : str
        Default background mode for export
    default_p : int
        Default percentile for export
    default_transform : str
        Default transform for export
    random_seed : int
        Random seed for reproducibility
    """

    mode: str = "export"
    p_candidates: List[int] = field(default_factory=lambda: [1, 3, 5, 7, 10])
    bg_modes: List[str] = field(default_factory=lambda: ["clip", "center"])
    transforms: List[str] = field(
        default_factory=lambda: ["raw", "log1p", "asinh_c5", "asinh_c10"]
    )
    default_bg_mode: str = "clip"
    default_p: int = 3
    default_transform: str = "log1p"
    random_seed: int = 42


@dataclass
class AlignmentConfig:
    """Configuration for cross-sample alignment (Stage D).

    Attributes
    ----------
    lower_percentile : float
        Lower percentile for alignment
    upper_percentile : float
        Upper percentile for alignment
    clip_to_targets : bool
        Whether to clip aligned values to target range
    metric_sample_limit : int
        Sample limit for distance metrics
    random_seed : int
        Random seed for reproducibility
    """

    lower_percentile: float = 5.0
    upper_percentile: float = 95.0
    clip_to_targets: bool = True
    metric_sample_limit: int = 5000
    random_seed: int = 2024


@dataclass
class BatchCorrectionConfig:
    """Configuration for batch effect correction (Stage E).

    Attributes
    ----------
    enabled : bool
        Whether to apply batch correction
    batch_col : str
        Column for batch/technical factor
    biovar_col : str
        Column for biological factor to preserve
    correction_mode : str
        Correction mode (residual or none)
    additional_batch_factors : List[str]
        Additional batch factors (e.g., imaging_color, imaging_cycle)
    embedding_cells_per_sample : int
        Cells per sample for embedding
    random_seed : int
        Random seed for reproducibility
    """

    enabled: bool = True
    batch_col: str = "donor"
    biovar_col: str = "region"
    correction_mode: str = "residual"
    additional_batch_factors: List[str] = field(default_factory=list)
    embedding_cells_per_sample: int = 1500
    random_seed: int = 2026


@dataclass
class MergeConfig:
    """Configuration for data merging (Stage F).

    Attributes
    ----------
    graph_mode : str
        Neighborhood graph mode (knn or ball)
    k : int
        Number of neighbors for kNN graph
    radius : float
        Radius for ball graph
    include_raw_layer : bool
        Whether to include raw layer
    include_aligned_layer : bool
        Whether to include aligned layer
    include_batchcorr_layer : bool
        Whether to include batch-corrected layer
    random_seed : int
        Random seed for reproducibility
    """

    graph_mode: str = "knn"
    k: int = 10
    radius: float = 20.0
    include_raw_layer: bool = True
    include_aligned_layer: bool = True
    include_batchcorr_layer: bool = True
    random_seed: int = 303


@dataclass
class PreprocessingConfig:
    """Master configuration for preprocessing pipeline.

    Attributes
    ----------
    loader : LoaderConfig
        Stage A configuration
    qc : QCConfig
        Stage B configuration
    normalization : NormalizationConfig
        Stage C configuration
    alignment : AlignmentConfig
        Stage D configuration
    batch_correction : BatchCorrectionConfig
        Stage E configuration
    merge : MergeConfig
        Stage F configuration
    skip_stages : List[str]
        Stages to skip
    """

    loader: LoaderConfig = field(default_factory=LoaderConfig)
    qc: QCConfig = field(default_factory=QCConfig)
    normalization: NormalizationConfig = field(default_factory=NormalizationConfig)
    alignment: AlignmentConfig = field(default_factory=AlignmentConfig)
    batch_correction: BatchCorrectionConfig = field(default_factory=BatchCorrectionConfig)
    merge: MergeConfig = field(default_factory=MergeConfig)
    skip_stages: List[str] = field(default_factory=list)

    @classmethod
    def from_yaml(cls, path: Path) -> "PreprocessingConfig":
        """Load configuration from YAML file."""
        with open(path) as f:
            data = yaml.safe_load(f) or {}

        # Handle nested preprocessing section
        if "preprocessing" in data:
            data = data["preprocessing"]

        return cls(
            loader=LoaderConfig(**data.get("loader", {})),
            qc=QCConfig(**data.get("qc", {})),
            normalization=NormalizationConfig(**data.get("normalization", {})),
            alignment=AlignmentConfig(**data.get("alignment", {})),
            batch_correction=BatchCorrectionConfig(**data.get("batch_correction", {})),
            merge=MergeConfig(**data.get("merge", {})),
            skip_stages=data.get("skip_stages", []),
        )

    @classmethod
    def default(cls) -> "PreprocessingConfig":
        """Create default configuration."""
        return cls()

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            "loader": {
                "cell_id_col": self.loader.cell_id_col,
                "sample_id_col": self.loader.sample_id_col,
                "required_metadata_cols": self.loader.required_metadata_cols,
                "coordinate_cols": self.loader.coordinate_cols,
            },
            "qc": {
                "area_percentile_low": self.qc.area_percentile_low,
                "area_percentile_high": self.qc.area_percentile_high,
                "nucleus_ratio_min": self.qc.nucleus_ratio_min,
                "nucleus_ratio_max": self.qc.nucleus_ratio_max,
                "intensity_percentile_low": self.qc.intensity_percentile_low,
                "autofluorescence_percentile_high": self.qc.autofluorescence_percentile_high,
                "max_removal_fraction": self.qc.max_removal_fraction,
            },
            "normalization": {
                "mode": self.normalization.mode,
                "default_bg_mode": self.normalization.default_bg_mode,
                "default_p": self.normalization.default_p,
                "default_transform": self.normalization.default_transform,
            },
            "alignment": {
                "lower_percentile": self.alignment.lower_percentile,
                "upper_percentile": self.alignment.upper_percentile,
                "clip_to_targets": self.alignment.clip_to_targets,
            },
            "batch_correction": {
                "enabled": self.batch_correction.enabled,
                "batch_col": self.batch_correction.batch_col,
                "biovar_col": self.batch_correction.biovar_col,
                "correction_mode": self.batch_correction.correction_mode,
            },
            "merge": {
                "graph_mode": self.merge.graph_mode,
                "k": self.merge.k,
                "radius": self.merge.radius,
            },
            "skip_stages": self.skip_stages,
        }
