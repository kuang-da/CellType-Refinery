"""Unit tests for preprocessing module."""

import pytest
import numpy as np
import pandas as pd
from pathlib import Path

from celltype_refinery.core.preprocessing import (
    LoaderConfig,
    QCConfig,
    NormalizationConfig,
    AlignmentConfig,
    BatchCorrectionConfig,
    MergeConfig,
    PreprocessingConfig,
    DataLoader,
    CellQC,
    Normalizer,
    CrossSampleAligner,
)


class TestLoaderConfig:
    """Tests for LoaderConfig dataclass."""

    def test_default_values(self):
        """Test default configuration values."""
        config = LoaderConfig()
        assert config.cell_id_col == "cell_mask_id"
        assert config.sample_id_col == "sample_id"
        assert len(config.required_metadata_cols) > 0

    def test_custom_cell_id(self):
        """Test custom cell ID column."""
        config = LoaderConfig(cell_id_col="custom_id")
        assert config.cell_id_col == "custom_id"


class TestQCConfig:
    """Tests for QCConfig dataclass."""

    def test_default_percentiles(self):
        """Test default QC percentile values."""
        config = QCConfig()
        assert config.area_percentile_low == 1.0
        assert config.area_percentile_high == 99.0
        assert config.nucleus_ratio_min == 0.1
        assert config.nucleus_ratio_max == 0.9

    def test_max_removal_fraction(self):
        """Test max removal fraction default."""
        config = QCConfig()
        assert config.max_removal_fraction == 0.15


class TestNormalizationConfig:
    """Tests for NormalizationConfig dataclass."""

    def test_default_values(self):
        """Test default normalization config."""
        config = NormalizationConfig()
        assert config.default_bg_mode == "clip"
        assert config.default_transform == "log1p"
        assert "log1p" in config.transforms


class TestAlignmentConfig:
    """Tests for AlignmentConfig dataclass."""

    def test_default_percentiles(self):
        """Test default alignment percentiles."""
        config = AlignmentConfig()
        assert config.lower_percentile == 5.0
        assert config.upper_percentile == 95.0
        assert config.clip_to_targets is True


class TestDataLoader:
    """Tests for DataLoader class."""

    @pytest.fixture
    def sample_matrix(self, tmp_path) -> Path:
        """Create sample cell matrix file."""
        df = pd.DataFrame({
            "cell_mask_id": [1, 2, 3, 4, 5],
            "Marker1": [1.0, 2.0, 3.0, 4.0, 5.0],
            "Marker2": [5.0, 4.0, 3.0, 2.0, 1.0],
        })
        path = tmp_path / "matrix.csv"
        df.to_csv(path, index=False)
        return path

    @pytest.fixture
    def sample_metadata(self, tmp_path) -> Path:
        """Create sample cell metadata file."""
        df = pd.DataFrame({
            "cell_mask_id": [1, 2, 3, 4, 5],
            "x": [100, 200, 300, 400, 500],
            "y": [50, 60, 70, 80, 90],
            "cell_area": [100, 120, 110, 130, 105],
        })
        path = tmp_path / "metadata.csv"
        df.to_csv(path, index=False)
        return path

    def test_init_default_config(self):
        """Test loader initialization with default config."""
        loader = DataLoader()
        assert loader.config.cell_id_col == "cell_mask_id"

    def test_load_cell_matrix(self, sample_matrix):
        """Test loading cell matrix from file."""
        loader = DataLoader()
        df = loader.load_cell_matrix(sample_matrix)
        assert len(df) == 5
        assert "Marker1" in df.columns
        assert "cell_mask_id" in df.columns

    def test_load_cell_metadata(self, sample_metadata):
        """Test loading cell metadata from file."""
        loader = DataLoader()
        df = loader.load_cell_metadata(sample_metadata)
        assert len(df) == 5
        assert "x" in df.columns
        assert "y" in df.columns

    def test_load_sample(self, sample_matrix, sample_metadata):
        """Test loading complete sample."""
        loader = DataLoader()
        result = loader.load_sample(
            sample_matrix, sample_metadata, "test_sample"
        )
        assert result.sample_id == "test_sample"
        assert result.cell_matrix is not None
        assert result.cell_metadata is not None
        assert len(result.markers) == 2


class TestCellQC:
    """Tests for CellQC class."""

    @pytest.fixture
    def sample_data(self):
        """Create sample data for QC testing."""
        np.random.seed(42)
        n_cells = 100

        matrix = pd.DataFrame({
            "cell_mask_id": list(range(n_cells)),
            "Marker1": np.random.lognormal(0, 1, n_cells),
            "Marker2": np.random.lognormal(0, 1, n_cells),
        })

        metadata = pd.DataFrame({
            "cell_mask_id": list(range(n_cells)),
            "x": np.random.uniform(0, 1000, n_cells),
            "y": np.random.uniform(0, 1000, n_cells),
            "cell_area": np.random.uniform(50, 500, n_cells),
            "nucleus_area": np.random.uniform(10, 100, n_cells),
        })

        return matrix, metadata

    def test_init_default_config(self):
        """Test QC initialization with default config."""
        qc = CellQC()
        assert qc.config.max_removal_fraction == 0.15

    def test_filter_sample(self, sample_data):
        """Test QC filtering on sample."""
        matrix, metadata = sample_data
        qc = CellQC()
        result = qc.filter_sample(matrix, metadata, "test")

        assert result.sample_id == "test"
        assert result.filtered_matrix is not None
        assert result.cells_total > 0
        assert result.cells_total >= result.cells_removed


class TestNormalizer:
    """Tests for Normalizer class."""

    @pytest.fixture
    def sample_matrix(self):
        """Create sample matrix for normalization."""
        np.random.seed(42)
        return pd.DataFrame({
            "cell_mask_id": list(range(50)),
            "Marker1": np.random.lognormal(0, 2, 50),
            "Marker2": np.random.lognormal(0, 2, 50),
        })

    def test_init_default_config(self):
        """Test normalizer initialization."""
        norm = Normalizer()
        assert norm.config.default_transform == "log1p"

    def test_apply_log1p_transform(self, sample_matrix):
        """Test log1p transformation."""
        norm = Normalizer()
        values = sample_matrix["Marker1"].values
        result = norm.apply_transform(values, "log1p")
        assert np.all(np.isfinite(result))
        assert np.all(result >= 0)

    def test_apply_asinh_transform(self, sample_matrix):
        """Test asinh transformation."""
        norm = Normalizer()
        values = sample_matrix["Marker1"].values
        result = norm.apply_transform(values, "asinh_c5")
        assert np.all(np.isfinite(result))


class TestCrossSampleAligner:
    """Tests for CrossSampleAligner class."""

    @pytest.fixture
    def sample_matrices(self):
        """Create sample matrices for alignment."""
        np.random.seed(42)
        matrices = {}
        for i in range(3):
            matrices[f"sample_{i}"] = pd.DataFrame({
                "cell_mask_id": list(range(50)),
                "Marker1": np.random.lognormal(i * 0.5, 1, 50),
                "Marker2": np.random.lognormal(i * 0.3, 1, 50),
            })
        return matrices

    def test_init_default_config(self):
        """Test aligner initialization."""
        aligner = CrossSampleAligner()
        assert aligner.config.lower_percentile == 5.0
        assert aligner.config.upper_percentile == 95.0

    def test_compute_global_targets(self, sample_matrices):
        """Test computing global target percentiles."""
        aligner = CrossSampleAligner()
        targets = aligner.compute_global_targets(sample_matrices)

        assert "Marker1" in targets
        assert "Marker2" in targets
        assert "target_lower" in targets["Marker1"]
        assert "target_upper" in targets["Marker1"]

    def test_align_sample(self, sample_matrices):
        """Test aligning a single sample."""
        aligner = CrossSampleAligner()
        targets = aligner.compute_global_targets(sample_matrices)

        result = aligner.align_sample(
            sample_matrices["sample_0"],
            "sample_0",
            targets,
        )

        assert result.sample_id == "sample_0"
        assert result.aligned_matrix is not None
        assert len(result.params) > 0
