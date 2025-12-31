---
sidebar_position: 2
---

# Quickstart

Get your first cell-type annotation running in 5 minutes.

## Prerequisites

- CellType-Refinery installed ([Installation guide](installation))
- A clustered AnnData file (`.h5ad`) with Leiden clusters
- A marker map JSON file

## Python API

### Basic Annotation

```python
from celltype_refinery.core.annotation import AnnotationEngine
import scanpy as sc

# Load your clustered data
adata = sc.read_h5ad("clustered_data.h5ad")

# Create annotation engine with marker map
engine = AnnotationEngine(marker_map_path="markers.json")

# Run annotation
result = engine.run(adata, output_dir="output/annotation")

# Check results
print(f"Annotated {result.n_cells} cells")
print(f"Found {result.n_types} cell types")
```

### With Refinement

```python
from celltype_refinery.core.annotation import AnnotationEngine
from celltype_refinery.core.refinement import RefinementEngine

# Annotate
engine = AnnotationEngine(marker_map_path="markers.json")
result = engine.run(adata, output_dir="output/annotation")

# Refine with automatic policy
refiner = RefinementEngine()
refiner.load(adata)
plan = refiner.create_auto_plan(score_threshold=1.0)
refiner.execute(plan)

# Save refined results
adata.write_h5ad("output/refined.h5ad")
```

## CLI

### Full Pipeline

Run the complete annotation pipeline from config:

```bash
celltype-refinery pipeline --config configs/pipelines/full.yaml
```

### Step-by-Step

Run individual stages for more control:

```bash
# 1. Cluster the merged data
celltype-refinery cluster \
  --input merged.h5ad \
  --out output/clustered

# 2. Annotate with marker map
celltype-refinery annotate \
  --input output/clustered/clustered.h5ad \
  --marker-map markers.json \
  --out output/annotated

# 3. Refine annotations (auto mode)
celltype-refinery refine \
  --input output/annotated/annotated.h5ad \
  --auto \
  --execute \
  --out output/refined

# 4. Run composition and spatial analysis
celltype-refinery analyze \
  --input output/refined/refined.h5ad \
  --out output/analysis
```

## Marker Map Format

Marker maps are JSON files defining the cell-type hierarchy:

```json
{
  "_marker_map_metadata": {
    "version": "1.0",
    "tissue": "my_tissue"
  },
  "Epithelium": {
    "markers": ["EpCAM", "Cytokeratin"],
    "anti_markers": ["CD45"],
    "subtypes": {
      "Ciliated": {
        "markers": ["FOXJ1", "Acetylated Tubulin"]
      },
      "Secretory": {
        "markers": ["MUC5B", "PAX8"]
      }
    }
  },
  "Immune Cells": {
    "markers": ["CD45"],
    "anti_markers": ["EpCAM"],
    "subtypes": {
      "Myeloids": {
        "markers": ["CD68", "CD14"]
      },
      "Lymphoids": {
        "markers": ["CD3", "CD19", "CD56"]
      }
    }
  }
}
```

## Output Files

After annotation, you'll find:

| File | Description |
|------|-------------|
| `annotated.h5ad` | AnnData with cell-type labels in `.obs` |
| `cluster_annotations.csv` | Per-cluster annotation summary |
| `marker_scores.csv` | Scoring matrix for all clusters |
| `mapping_table.csv` | Cluster ID to cell-type mapping |

## Next Steps

- [Configuration Guide](configuration) - Customize parameters
- [First Annotation Tutorial](first-annotation) - Detailed walkthrough
- [Core Workflows](../core-workflows/workflow-overview) - Learn usage patterns
