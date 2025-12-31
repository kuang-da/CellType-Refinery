---
sidebar_position: 3
---

# Configuration

CellType-Refinery uses YAML configuration files for tissue-specific parameters and pipeline settings.

## Configuration Files

### Tissue Templates

Tissue-specific settings are defined in YAML templates:

```yaml
# configs/tissues/my_tissue.yaml
version: "1.0"
tissue: "my_tissue"

# Hierarchical gating parameters
gating:
  # Root marker requirements
  root_hard_requirements:
    "Immune Cells":
      marker: "CD45"
      min_pos_frac: 0.30
    "Epithelium":
      marker: "EpCAM"
      min_pos_frac: 0.25

  # Default gating thresholds
  defaults:
    min_coverage: 0.3
    min_pos_frac: 0.2

# Pattern matching for cell type grouping
patterns:
  epithelial: ["Epithelium", "Epithelial", "Keratin"]
  immune: ["CD45", "Immune", "Myeloid", "Lymphoid"]
  stromal: ["Fibroblast", "Mesenchymal", "Stromal"]
  endothelial: ["Endothelium", "CD31", "VWF"]

# Regional configuration (optional)
regions:
  expected:
    - name: "proximal"
      cell_types: ["Epithelium", "Stromal"]
    - name: "distal"
      cell_types: ["Epithelium", "Immune Cells"]
```

### Pipeline Configuration

Full pipeline settings for automated runs:

```yaml
# configs/pipelines/full.yaml
version: "1.0"
name: "full_pipeline"

# Input/output paths
paths:
  input: "data/merged.h5ad"
  output: "output"
  marker_map: "configs/marker_maps/markers.json"
  tissue_config: "configs/tissues/my_tissue.yaml"

# Stage parameters
stages:
  clustering:
    enabled: true
    resolution: 0.6
    n_pcs: 30
    neighbors_k: 15
    use_gpu: true

  annotation:
    enabled: true
    layer: "batchcorr"

  refinement:
    enabled: true
    mode: "auto"  # auto, manual, or hybrid
    execute: true
    score_threshold: 1.0
    min_cells: 500

  analysis:
    composition: true
    spatial: true
    review: true
```

## Environment Variables

Override configuration via environment variables:

| Variable | Description | Default |
|----------|-------------|---------|
| `CTR_CONFIG_DIR` | Configuration directory | `./configs` |
| `CTR_OUTPUT_DIR` | Default output directory | `./output` |
| `CTR_LOG_LEVEL` | Logging verbosity | `INFO` |
| `CTR_USE_GPU` | Enable GPU acceleration | `false` |

## CLI Configuration

Override config file settings via CLI flags:

```bash
celltype-refinery annotate \
  --input data.h5ad \
  --marker-map markers.json \
  --layer batchcorr \
  --min-coverage 0.3 \
  --min-pos-frac 0.2 \
  --out output/
```

## Configuration Precedence

Settings are applied in order (later overrides earlier):

1. Package defaults
2. Tissue template (`configs/tissues/*.yaml`)
3. Pipeline config (`configs/pipelines/*.yaml`)
4. Environment variables
5. CLI arguments

## Validation

Validate configuration before running:

```bash
celltype-refinery validate-config --config pipeline.yaml
```

Or in Python:

```python
from celltype_refinery.pipeline import PipelineConfig

config = PipelineConfig.load("pipeline.yaml")
errors = config.validate()
if errors:
    for error in errors:
        print(f"Error: {error}")
```

## Next Steps

- [First Annotation](first-annotation) - Complete walkthrough
- [Marker Maps](../configuration/marker-maps) - Marker map reference
- [Tissue Templates](../configuration/tissue-templates) - Template reference
