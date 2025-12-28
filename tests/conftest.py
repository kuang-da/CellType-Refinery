"""Pytest configuration and shared fixtures for CellType-Refinery tests."""

import sys
from pathlib import Path

import pytest
import numpy as np
import pandas as pd

# Add package to path
PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

# Import mock data generators
from tests.fixtures import (
    create_mock_adata,
    create_clustered_adata,
    create_annotated_adata,
    create_minimal_expression_matrix,
    create_minimal_cell_metadata,
)


# ============================================================================
# Mock Data Fixtures
# ============================================================================


@pytest.fixture
def minimal_obs() -> pd.DataFrame:
    """Create minimal cell observation DataFrame."""
    np.random.seed(42)
    n_cells = 100
    n_clusters = 5

    cluster_ids = np.repeat([str(i) for i in range(n_clusters)], n_cells // n_clusters)

    cell_type_map = {
        "0": "Epithelium",
        "1": "Immune Cells",
        "2": "Mesenchymal Cells",
        "3": "Endothelium",
        "4": "Unassigned",
    }

    return pd.DataFrame({
        "cluster_lvl0": pd.Categorical(cluster_ids),
        "cell_type_lvl0_auto": [cell_type_map.get(c, "Unassigned") for c in cluster_ids],
        "sample_id": np.repeat(["S1", "S2"], n_cells // 2),
        "donor": np.repeat(["D1", "D2"], n_cells // 2),
        "region": np.repeat(["region_A", "region_B"], n_cells // 2),
    })


@pytest.fixture
def minimal_marker_map() -> dict:
    """Create minimal marker map for testing."""
    return {
        "Epithelium": {
            "markers": ["Marker0", "Marker1", "Marker2"],
            "anti_markers": ["Marker10"],
            "children": {
                "Ciliated Epithelium": {
                    "markers": ["Marker0", "Marker3"],
                    "anti_markers": [],
                },
                "Secretory Epithelium": {
                    "markers": ["Marker1", "Marker4"],
                    "anti_markers": [],
                },
            },
        },
        "Immune Cells": {
            "markers": ["Marker5", "Marker6"],
            "anti_markers": ["Marker0"],
            "children": {
                "T Cells": {"markers": ["Marker5"], "anti_markers": []},
                "Macrophages": {"markers": ["Marker6"], "anti_markers": []},
            },
        },
        "Mesenchymal Cells": {
            "markers": ["Marker7", "Marker8"],
            "anti_markers": [],
            "children": {},
        },
        "Endothelium": {
            "markers": ["Marker9", "Marker10"],
            "anti_markers": [],
            "children": {},
        },
    }


@pytest.fixture
def minimal_cluster_annotations() -> pd.DataFrame:
    """Create minimal cluster annotations DataFrame."""
    return pd.DataFrame({
        "cluster_id": ["0", "1", "2", "3", "4"],
        "assigned_label": ["Epithelium", "Immune Cells", "Mesenchymal Cells",
                          "Endothelium", "Unassigned"],
        "assigned_score": [2.5, 1.8, 2.0, 1.5, 0.5],
        "n_cells": [100, 80, 90, 70, 60],
        "top2_gap": [0.8, 0.6, 0.7, 0.5, 0.2],
        "is_ambiguous_root": [False, False, False, False, False],
    })


@pytest.fixture
def tmp_output_dir(tmp_path: Path) -> Path:
    """Create a temporary output directory for tests."""
    output_dir = tmp_path / "output"
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


# ============================================================================
# Path Fixtures
# ============================================================================


@pytest.fixture
def ft_stage_h_output() -> Path:
    """Path to frozen Stage H output for equivalence tests.

    Returns None if the file doesn't exist (for CI environments).
    """
    path = Path("/workspaces/1-spatial_frs_analysis/ft/out/stage_h/coarse_clusters.h5ad")
    if path.exists():
        return path
    return None


@pytest.fixture
def ft_marker_map() -> Path:
    """Path to FT marker map for equivalence tests."""
    path = Path("/workspaces/1-spatial_frs_analysis/ft/data/FT_cell_type_markers_v9.json")
    if path.exists():
        return path
    return None


# ============================================================================
# AnnData Fixtures
# ============================================================================


@pytest.fixture
def mock_adata():
    """Create mock AnnData with 500 cells and 20 markers."""
    return create_mock_adata(n_cells=500, n_markers=20, n_clusters=5)


@pytest.fixture
def small_adata():
    """Create small AnnData for quick tests."""
    return create_mock_adata(n_cells=100, n_markers=10, n_clusters=3)


@pytest.fixture
def clustered_adata():
    """Create AnnData that looks like Stage H output."""
    return create_clustered_adata()


@pytest.fixture
def annotated_adata():
    """Create AnnData with cell type annotations."""
    return create_annotated_adata()


# ============================================================================
# Configuration Fixtures
# ============================================================================


@pytest.fixture
def sample_pipeline_config(tmp_path) -> Path:
    """Create sample pipeline configuration file."""
    import yaml

    config = {
        "pipeline": {
            "name": "Test Pipeline",
            "version": "1.0",
        },
        "global": {
            "output_dir": str(tmp_path / "output"),
        },
        "stages": {
            "preprocess": {
                "name": "Preprocessing",
                "script_module": "celltype_refinery.stages.preprocess",
                "inputs": {"raw": str(tmp_path / "raw")},
                "outputs": {"processed": str(tmp_path / "processed")},
            },
            "cluster": {
                "name": "Clustering",
                "script_module": "celltype_refinery.stages.cluster",
                "depends_on": ["preprocess"],
            },
        },
    }

    path = tmp_path / "pipeline.yaml"
    with open(path, "w") as f:
        yaml.dump(config, f)

    return path


@pytest.fixture
def sample_tissue_config(tmp_path) -> Path:
    """Create sample tissue configuration file."""
    import yaml

    config = {
        "version": "1.0",
        "tissue": "test_tissue",
        "gating": {
            "min_coverage": {0: 0.5, 1: 0.4, 2: 0.3},
            "min_pos_frac": {0: 0.3, 1: 0.2, 2: 0.15},
        },
        "patterns": {
            "epithelial": ["Epithelium", "Epithelial"],
            "stromal": ["Stromal", "Stroma"],
            "immune": ["Immune", "CD45"],
        },
    }

    path = tmp_path / "tissue.yaml"
    with open(path, "w") as f:
        yaml.dump(config, f)

    return path
