---
sidebar_position: 1
---

# Clustering Overview

The clustering module performs Leiden clustering and differential expression analysis to identify distinct cell populations in your data.

:::tip Pipeline Context
Clustering is part of Stage H. For the complete pipeline, see [Annotation Pipeline](/docs/methodology/annotation-pipeline.md).
:::

## Pipeline Architecture

```text
┌─────────────────────────────────────────────────────────────────────────┐
│                         CLUSTERING PIPELINE                             │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                         │
│  Expression Matrix (n_cells × n_markers)                                │
│         │                                                               │
│         ▼                                                               │
│  ┌─────────────────────────────────────────────────────────────────┐   │
│  │ 1. SCALE                                                        │   │
│  │    • Zero-center each marker                                    │   │
│  │    • Clip values at ±10 standard deviations                     │   │
│  └─────────────────────────────────────────────────────────────────┘   │
│         │                                                               │
│         ▼                                                               │
│  ┌─────────────────────────────────────────────────────────────────┐   │
│  │ 2. PCA                                                          │   │
│  │    • n_components = 30 (default)                                │   │
│  │    • Stores in adata.obsm["X_pca"]                              │   │
│  └─────────────────────────────────────────────────────────────────┘   │
│         │                                                               │
│         ▼                                                               │
│  ┌─────────────────────────────────────────────────────────────────┐   │
│  │ 3. NEIGHBORS                                                    │   │
│  │    • k = 15 (default)                                           │   │
│  │    • Uses PCA coordinates                                       │   │
│  └─────────────────────────────────────────────────────────────────┘   │
│         │                                                               │
│         ▼                                                               │
│  ┌─────────────────────────────────────────────────────────────────┐   │
│  │ 4. LEIDEN                                                       │   │
│  │    • resolution = 0.6 (default)                                 │   │
│  │    • GPU: rapids_singlecell / CPU: igraph                       │   │
│  │    • Stores in adata.obs["cluster_lvl0"]                        │   │
│  └─────────────────────────────────────────────────────────────────┘   │
│         │                                                               │
│         ▼                                                               │
│  ┌─────────────────────────────────────────────────────────────────┐   │
│  │ 5. UMAP (optional)                                              │   │
│  │    • 2D embedding for visualization                             │   │
│  └─────────────────────────────────────────────────────────────────┘   │
│                                                                         │
└─────────────────────────────────────────────────────────────────────────┘
```

## Preprocessing

Before clustering begins, several preprocessing steps prepare the data:

### Layer Selection

The pipeline automatically selects the best available data layer in this priority order:

1. **batchcorr** - Batch-corrected expression values (preferred)
2. **aligned** - Aligned expression values
3. **X** - Raw expression matrix (fallback)

### Technical Marker Exclusion

Technical markers that should not influence clustering are automatically moved to `.obs`:

- **DAPI** - Nuclear stain
- **Collagen IV** - Structural marker
- **Beta-actin** - Housekeeping gene

These markers remain available for visualization but do not affect cluster assignments.

### Low-Variance Filtering

Markers with very low variance (std < 1e-3) are excluded from clustering to prevent numerical instability and improve cluster quality.

:::warning GPU Non-Determinism
GPU-accelerated Leiden produces different cluster counts across runs (~5 clusters).
For reproducibility, use `--no-gpu` or run clustering once and use `--annotation-only` for iterations.
:::

## Key Parameters

| Parameter | Default | Description | When to Adjust |
|-----------|---------|-------------|----------------|
| `resolution` | 0.6 | Leiden clustering resolution | Increase for more clusters, decrease for fewer |
| `n_pcs` | 30 | Number of principal components | Lower if markers < 30; increase for complex datasets |
| `neighbors_k` | 15 | k for k-NN graph construction | Increase for smoother clusters; decrease to preserve rare populations |
| `use_gpu` | false | Enable GPU acceleration | Set true for large datasets (>100k cells) |
| `min_cells` | 10 | Minimum cells per cluster | Increase to filter noise clusters |
| `random_state` | 42 | Random seed for reproducibility | Change to test stability |

## GPU Acceleration

The clustering module supports GPU acceleration via RAPIDS and cuGraph for significantly faster processing on large datasets.

### When GPU is Available

- **rapids_singlecell** handles PCA and neighbor computation
- **cuGraph** performs GPU-accelerated Leiden clustering
- Typical speedup: 10-50x for datasets >100k cells

### Automatic Fallback

If GPU is unavailable or `--no-gpu` is specified:

- **scanpy** handles PCA and neighbor computation
- **igraph** performs CPU-based Leiden clustering

### Hardware Requirements

- NVIDIA GPU with CUDA support
- RAPIDS libraries installed (`rapids-singlecell`, `cugraph`)

## CLI Usage

### Basic Usage

```bash
celltype-refinery cluster \
  --input merged.h5ad \
  --resolution 0.6 \
  --n-pcs 30 \
  --out output/clustered
```

### High-Resolution Clustering

```bash
celltype-refinery cluster \
  --input merged.h5ad \
  --resolution 1.2 \
  --neighbors-k 10 \
  --out output/high_res
```

### GPU-Accelerated Processing

```bash
celltype-refinery cluster \
  --input large_dataset.h5ad \
  --use-gpu \
  --resolution 0.6 \
  --out output/gpu_clustered
```

### Reproducible CPU-Only Run

```bash
celltype-refinery cluster \
  --input merged.h5ad \
  --no-gpu \
  --random-state 42 \
  --out output/reproducible
```

### Annotation-Only Mode (Skip Reclustering)

```bash
celltype-refinery cluster \
  --input already_clustered.h5ad \
  --annotation-only \
  --out output/reannotated
```

## See Also

- [Parameters Reference](./parameters.md) - Complete parameter documentation
- [Tuning Guide](./tuning-guide.md) - Best practices for parameter tuning
- [Annotation Pipeline](/docs/methodology/annotation-pipeline.md) - Full pipeline overview
- [Marker Interpretation](/docs/reference/marker-interpretation.md) - Understanding marker expression
