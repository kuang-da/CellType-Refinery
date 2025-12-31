---
sidebar_position: 1
---

# Installation

CellType-Refinery can be installed from source using pip.

## Requirements

- Python >= 3.10
- pip

## From Source

Clone the repository and install in editable mode:

```bash
git clone https://github.com/kimpenn/CellType-Refinery.git
cd CellType-Refinery
pip install -e .
```

## Optional Dependencies

### Visualization Support

For generating plots and figures:

```bash
pip install -e ".[viz]"
```

This includes `matplotlib` and `seaborn`.

### GPU Acceleration

For RAPIDS/cuGraph-accelerated Leiden clustering:

```bash
pip install -e ".[gpu]"
```

This includes `rapids-singlecell` and `cupy-cuda12x`.

### Development Tools

For testing and code quality:

```bash
pip install -e ".[dev]"
```

This includes `pytest`, `pytest-cov`, `ruff`, and `mypy`.

### All Dependencies

Install everything:

```bash
pip install -e ".[all]"
```

## Verify Installation

After installation, verify CellType-Refinery is working:

```bash
# Check CLI is available
celltype-refinery --help

# Check Python import
python -c "from celltype_refinery.core.annotation import AnnotationEngine; print('OK')"
```

## Dependencies

### Core Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| numpy | \>= 1.24 | Numerical computing |
| pandas | \>= 2.0 | Data manipulation |
| scipy | \>= 1.10 | Statistical functions |
| scikit-learn | \>= 1.2 | Machine learning utilities |
| anndata | \>= 0.10 | AnnData objects |
| scanpy | \>= 1.9 | Single-cell analysis |
| leidenalg | \>= 0.10 | Leiden clustering |
| pyyaml | \>= 6.0 | Configuration files |
| click | \>= 8.0 | CLI framework |
| tqdm | \>= 4.65 | Progress bars |
| statsmodels | \>= 0.14 | Statistical models |

### Optional Dependencies

| Extra | Packages | Purpose |
|-------|----------|---------|
| `viz` | matplotlib, seaborn | Visualization |
| `gpu` | rapids-singlecell, cupy | GPU acceleration |
| `dev` | pytest, ruff, mypy | Development tools |
