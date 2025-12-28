"""Unit tests for clustering module."""

import pytest
import numpy as np
import pandas as pd

from celltype_refinery.core.clustering import (
    ClusteringConfig,
    DEConfig,
    StageHConfig,
    ClusteringEngine,
    ClusteringResult,
)


class TestClusteringConfig:
    """Tests for ClusteringConfig dataclass."""

    def test_default_values(self):
        """Test default configuration values."""
        config = ClusteringConfig()
        assert config.n_pcs == 30
        assert config.neighbors_k == 15
        assert config.resolution == 0.6
        assert config.scale_clip == 10.0
        assert config.random_seed == 1337
        assert config.use_gpu is True
        assert config.compute_umap is False

    def test_custom_values(self):
        """Test custom configuration values."""
        config = ClusteringConfig(
            n_pcs=50,
            resolution=0.8,
            use_gpu=False,
        )
        assert config.n_pcs == 50
        assert config.resolution == 0.8
        assert config.use_gpu is False


class TestDEConfig:
    """Tests for DEConfig dataclass."""

    def test_default_values(self):
        """Test default DE configuration values."""
        config = DEConfig()
        assert config.method == "wilcoxon"
        assert config.n_genes == 12
        assert config.layer == "batchcorr"
        assert config.tie_correct is True

    def test_custom_method(self):
        """Test custom DE method."""
        config = DEConfig(method="t-test")
        assert config.method == "t-test"


class TestStageHConfig:
    """Tests for StageHConfig dataclass."""

    def test_default_values(self):
        """Test master config with default values."""
        config = StageHConfig()
        assert config.clustering.n_pcs == 30
        assert config.de.method == "wilcoxon"
        assert config.subcluster.min_cells == 500
        assert config.scoring.positive_quantile == 0.75
        assert "DAPI" in config.technical_markers

    def test_from_yaml(self, tmp_path):
        """Test loading config from YAML."""
        yaml_content = """
stage_h:
  clustering:
    n_pcs: 25
    resolution: 0.5
  de:
    method: t-test
    n_genes: 15
  technical_markers:
    - DAPI
    - Custom_Marker
"""
        yaml_file = tmp_path / "config.yaml"
        yaml_file.write_text(yaml_content)

        config = StageHConfig.from_yaml(yaml_file)
        assert config.clustering.n_pcs == 25
        assert config.clustering.resolution == 0.5
        assert config.de.method == "t-test"
        assert config.de.n_genes == 15
        assert "Custom_Marker" in config.technical_markers

    def test_to_dict(self):
        """Test converting config to dictionary."""
        config = StageHConfig()
        d = config.to_dict()
        assert "clustering" in d
        assert "de" in d
        assert "subcluster" in d
        assert d["clustering"]["n_pcs"] == 30


class TestClusteringEngine:
    """Tests for ClusteringEngine class."""

    @pytest.fixture
    def mock_adata(self):
        """Create mock AnnData for testing."""
        import anndata as ad

        np.random.seed(42)
        n_cells = 200
        n_markers = 15

        X = np.random.lognormal(0, 1, (n_cells, n_markers)).astype(np.float32)
        obs = pd.DataFrame({
            "sample_id": pd.Categorical(["S1"] * 100 + ["S2"] * 100),
            "donor": pd.Categorical(["D1"] * 100 + ["D2"] * 100),
        })
        obs.index = pd.Index([f"cell_{i}" for i in range(n_cells)])
        var = pd.DataFrame(index=[f"Marker_{i}" for i in range(n_markers)])

        return ad.AnnData(X=X, obs=obs, var=var)

    def test_init_default_config(self):
        """Test engine initialization with default config."""
        engine = ClusteringEngine()
        assert engine.config is not None
        assert engine.config.clustering.n_pcs == 30

    def test_init_custom_config(self):
        """Test engine initialization with custom config."""
        config = StageHConfig(clustering=ClusteringConfig(n_pcs=20))
        engine = ClusteringEngine(config)
        assert engine.config.clustering.n_pcs == 20

    def test_gpu_available_property(self):
        """Test GPU availability check."""
        engine = ClusteringEngine()
        # Should not raise, returns bool
        assert isinstance(engine.gpu_available, bool)

    def test_subset_adata_by_sample(self, mock_adata):
        """Test AnnData subsetting by sample ID."""
        engine = ClusteringEngine()
        subset = engine.subset_adata(mock_adata, sample_ids=["S1"])
        assert subset.n_obs == 100
        assert all(subset.obs["sample_id"] == "S1")

    def test_subset_adata_empty_raises(self, mock_adata):
        """Test that subsetting to zero cells raises error."""
        engine = ClusteringEngine()
        with pytest.raises(ValueError, match="removed all cells"):
            engine.subset_adata(mock_adata, sample_ids=["nonexistent"])

    def test_filter_low_variance_markers(self, mock_adata):
        """Test low variance marker filtering."""
        engine = ClusteringEngine()
        # Add constant marker
        mock_adata.X[:, 0] = 1.0

        filtered, dropped = engine.filter_low_variance_markers(mock_adata, min_std=0.1)
        assert "Marker_0" in dropped
        assert filtered.n_vars < mock_adata.n_vars

    def test_exclude_technical_markers(self, mock_adata):
        """Test technical marker exclusion."""
        engine = ClusteringEngine()
        # Rename a marker to DAPI
        mock_adata.var.index = ["DAPI"] + [f"Marker_{i}" for i in range(1, mock_adata.n_vars)]

        filtered, excluded = engine.exclude_technical_markers(mock_adata, ["DAPI"])
        assert "DAPI" in excluded
        assert "DAPI_intensity" in filtered.obs.columns
        assert "DAPI" not in filtered.var_names

    def test_select_layer(self, mock_adata):
        """Test layer selection."""
        engine = ClusteringEngine()
        mock_adata.layers["batchcorr"] = mock_adata.X.copy()

        layer = engine.select_layer(mock_adata, "batchcorr")
        assert layer == "batchcorr"
        assert "stage_h_input" in mock_adata.layers

    def test_select_layer_fallback(self, mock_adata):
        """Test layer fallback when requested layer not found."""
        engine = ClusteringEngine()
        layer = engine.select_layer(mock_adata, "nonexistent")
        assert layer == "X"

    def test_run_clustering_cpu(self, mock_adata):
        """Test clustering with CPU (no GPU)."""
        config = StageHConfig(clustering=ClusteringConfig(use_gpu=False))
        engine = ClusteringEngine(config)

        result = engine.run_clustering(mock_adata, use_gpu=False)

        assert isinstance(result, ClusteringResult)
        assert result.n_clusters > 0
        assert result.cluster_key == "cluster_lvl0"
        assert "cluster_lvl0" in mock_adata.obs.columns
        assert len(result.cluster_sizes) == result.n_clusters


class TestClusteringResult:
    """Tests for ClusteringResult dataclass."""

    def test_default_values(self):
        """Test default result values."""
        result = ClusteringResult()
        assert result.n_clusters == 0
        assert result.cluster_key == "cluster_lvl0"
        assert result.cluster_sizes == {}
        assert result.dropped_markers == []
        assert result.excluded_markers == []
