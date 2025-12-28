"""Preprocessing module for data loading and quality control.

Provides functions for loading cell matrices, QC filtering,
normalization, alignment, batch correction, and data merging.

Pipeline Stages
---------------
- Stage A (Loader): Data loading and validation
- Stage B (QC): Cell-level quality control
- Stage C (Normalization): Background correction and variance stabilization
- Stage D (Alignment): Cross-sample percentile alignment
- Stage E (Batch): Batch effect correction
- Stage F (Merge): Data merging and spatial graph construction

Example Usage
-------------
>>> from celltype_refinery.core.preprocessing import (
...     DataLoader, LoaderConfig,
...     CellQC, QCConfig,
...     Normalizer, NormalizationConfig,
...     CrossSampleAligner, AlignmentConfig,
...     BatchCorrector, BatchCorrectionConfig,
...     DataMerger, MergeConfig,
... )
>>> # Load and validate data
>>> loader = DataLoader(LoaderConfig())
>>> result = loader.load_sample(matrix_path, metadata_path, "sample_01")
>>> # Quality control
>>> qc = CellQC(QCConfig())
>>> qc_result = qc.filter_sample(result.cell_matrix, result.cell_metadata, "sample_01")
"""

__version__ = "1.0.0"

# Configuration classes
from .config import (
    LoaderConfig,
    QCConfig,
    NormalizationConfig,
    AlignmentConfig,
    BatchCorrectionConfig,
    MergeConfig,
    PreprocessingConfig,
)

# Stage A: Data loading
from .loader import (
    DataLoader,
    LoadResult,
)

# Stage B: Cell QC
from .qc import (
    CellQC,
    QCResult,
    REASON_COLUMNS,
)

# Stage C: Normalization
from .normalization import (
    Normalizer,
    NormalizationResult,
    TransformSpec,
    TRANSFORMS,
)

# Stage D: Cross-sample alignment
from .alignment import (
    CrossSampleAligner,
    AlignmentResult,
    AlignmentParams,
)

# Stage E: Batch correction
from .batch import (
    BatchCorrector,
    BatchCorrectionResult,
    BatchDiagnostics,
)

# Stage F: Data merging
from .merge import (
    DataMerger,
    MergeResult,
)

__all__ = [
    # Version
    "__version__",
    # Config
    "LoaderConfig",
    "QCConfig",
    "NormalizationConfig",
    "AlignmentConfig",
    "BatchCorrectionConfig",
    "MergeConfig",
    "PreprocessingConfig",
    # Stage A: Loader
    "DataLoader",
    "LoadResult",
    # Stage B: QC
    "CellQC",
    "QCResult",
    "REASON_COLUMNS",
    # Stage C: Normalization
    "Normalizer",
    "NormalizationResult",
    "TransformSpec",
    "TRANSFORMS",
    # Stage D: Alignment
    "CrossSampleAligner",
    "AlignmentResult",
    "AlignmentParams",
    # Stage E: Batch correction
    "BatchCorrector",
    "BatchCorrectionResult",
    "BatchDiagnostics",
    # Stage F: Merge
    "DataMerger",
    "MergeResult",
]
