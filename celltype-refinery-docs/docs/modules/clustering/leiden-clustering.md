---
sidebar_position: 2
---

# Leiden Clustering

Community detection using the Leiden algorithm.

```mermaid
flowchart LR
    subgraph Input
        A[merged.h5ad<br/>batchcorr layer]
    end

    subgraph Pipeline["Clustering Pipeline"]
        B[Scale Data] --> C[PCA<br/>n_pcs=30]
        C --> D[k-NN Graph<br/>k=15]
        D --> E[Leiden<br/>resolution=0.6]
    end

    subgraph Output
        F[cluster_lvl0<br/>assignments]
    end

    A --> B
    E --> F

    style E fill:#FFD700
```

## Algorithm

1. Compute PCA embeddings
2. Build k-NN graph
3. Run Leiden clustering
4. Assign cluster labels

## GPU Acceleration

GPU-accelerated Leiden via RAPIDS is non-deterministic. Cluster counts may vary between runs.

## CLI

```bash
celltype-refinery cluster leiden \
  --input merged.h5ad \
  --resolution 0.6 \
  --use-gpu \
  --out output/
```
