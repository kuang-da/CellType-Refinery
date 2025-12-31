---
sidebar_position: 4
---

# First Annotation Tutorial

A complete walkthrough of annotating your first dataset with CellType-Refinery.

## Overview

In this tutorial, you'll:

1. Prepare your input data
2. Create a marker map
3. Run clustering
4. Annotate cell types
5. Review and refine results

## Prerequisites

- CellType-Refinery installed
- A merged AnnData file with normalized expression data
- Knowledge of expected cell types in your tissue

## Step 1: Prepare Input Data

Your input should be an AnnData (`.h5ad`) file with:

- **X matrix**: Normalized expression values
- **obs**: Cell metadata (sample_id, region, etc.)
- **var**: Marker/antibody names

```python
import scanpy as sc

# Load your data
adata = sc.read_h5ad("merged_data.h5ad")

# Check structure
print(f"Cells: {adata.n_obs}")
print(f"Markers: {adata.n_vars}")
print(f"Layers: {list(adata.layers.keys())}")
```

## Step 2: Create Marker Map

Define your cell-type hierarchy in JSON:

```json
{
  "_marker_map_metadata": {
    "version": "1.0",
    "tissue": "example_tissue",
    "description": "Example marker map"
  },
  "Epithelium": {
    "markers": ["EpCAM", "Pan-Cytokeratin"],
    "anti_markers": ["CD45", "Vimentin"],
    "subtypes": {
      "Type A": {
        "markers": ["Marker1", "Marker2"]
      },
      "Type B": {
        "markers": ["Marker3", "Marker4"]
      }
    }
  },
  "Immune Cells": {
    "markers": ["CD45"],
    "anti_markers": ["EpCAM"],
    "subtypes": {
      "T Cells": {
        "markers": ["CD3"]
      },
      "Macrophages": {
        "markers": ["CD68", "CD163"]
      }
    }
  },
  "Stromal": {
    "markers": ["Vimentin"],
    "anti_markers": ["CD45", "EpCAM"]
  }
}
```

Save as `marker_map.json`.

## Step 3: Run Clustering

Cluster cells using Leiden algorithm:

```bash
celltype-refinery cluster \
  --input merged_data.h5ad \
  --resolution 0.6 \
  --n-pcs 30 \
  --out output/clustering
```

Or in Python:

```python
from celltype_refinery.core.clustering import ClusteringEngine

engine = ClusteringEngine(resolution=0.6, n_pcs=30)
result = engine.run(adata, output_dir="output/clustering")

print(f"Found {result.n_clusters} clusters")
```

## Step 4: Annotate Cell Types

Run annotation with your marker map:

```bash
celltype-refinery annotate \
  --input output/clustering/clustered.h5ad \
  --marker-map marker_map.json \
  --out output/annotation
```

Or in Python:

```python
from celltype_refinery.core.annotation import AnnotationEngine

engine = AnnotationEngine(marker_map_path="marker_map.json")
result = engine.run(adata, output_dir="output/annotation")

# Review results
print(result.summary())
```

## Step 5: Review Results

Check the annotation outputs:

```python
import pandas as pd

# Load cluster annotations
annotations = pd.read_csv("output/annotation/cluster_annotations.csv")
print(annotations[["cluster_id", "assigned_label", "score", "confidence"]])

# Check cell-type distribution
print(adata.obs["cell_type_curated"].value_counts())
```

### Understanding Scores

| Score Range | Confidence | Interpretation |
|-------------|------------|----------------|
| \> 2.0 | HIGH | Strong marker expression |
| 1.0 - 2.0 | MEDIUM | Moderate expression |
| 0.5 - 1.0 | LOW | Weak expression |
| \< 0.5 | VERY_LOW | Consider refinement |

## Step 6: Refine Annotations

If some clusters have low confidence, refine them:

```bash
# Diagnostic mode first (no changes)
celltype-refinery refine \
  --input output/annotation/annotated.h5ad \
  --auto \
  --out output/refinement

# Review diagnostic report
cat output/refinement/diagnostic_report.csv

# Execute refinement
celltype-refinery refine \
  --input output/annotation/annotated.h5ad \
  --auto \
  --execute \
  --out output/refinement
```

## Step 7: Run Analysis

Generate composition and spatial analysis:

```bash
celltype-refinery analyze \
  --input output/refinement/refined.h5ad \
  --out output/analysis
```

## Output Structure

After completing all steps:

```
output/
├── clustering/
│   ├── clustered.h5ad
│   └── cluster_stats.csv
├── annotation/
│   ├── annotated.h5ad
│   ├── cluster_annotations.csv
│   ├── marker_scores.csv
│   └── mapping_table.csv
├── refinement/
│   ├── refined.h5ad
│   ├── diagnostic_report.csv
│   └── curation_log.json
└── analysis/
    ├── composition/
    ├── spatial/
    └── review/
```

## Troubleshooting

### High Unassigned Rate

If many cells are unassigned:

1. Check marker names match between data and marker map
2. Lower gating thresholds (`--min-coverage 0.2`)
3. Review marker expression with `sc.pl.dotplot()`

### Incorrect Annotations

If cell types seem wrong:

1. Check marker specificity in your data
2. Add anti-markers to distinguish similar types
3. Use manual overrides for known misassignments

### Low Confidence Scores

If scores are consistently low:

1. Verify data normalization
2. Check for batch effects
3. Consider fewer, more specific markers

## Next Steps

- [Core Workflows](../core-workflows/workflow-overview) - Learn advanced patterns
- [Refinement Guide](../modules/refinement/overview) - Deep dive into refinement
- [Configuration Reference](../configuration/marker-maps) - Full parameter docs
