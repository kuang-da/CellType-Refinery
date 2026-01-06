---
sidebar_position: 7
---

# Data Merging

Merge all samples into a unified AnnData object.

```mermaid
flowchart TD
    subgraph Sources["Source Data"]
        S1[Normalized<br/>Stage C]
        S2[Aligned<br/>Stage D]
        S3[Corrected<br/>Stage E]
        S4[Metadata<br/>Stage A]
    end

    subgraph AnnData["merged_data.h5ad"]
        direction TB
        subgraph Layers
            L1[raw layer]
            L2[aligned layer]
            L3[batchcorr layer]
        end
        subgraph Meta["Cell Metadata"]
            M1[obs: cell info]
            M2[var: markers]
        end
        subgraph Graphs["Spatial Graphs"]
            G1[k-NN Graph]
            G2[Radius Graph]
        end
    end

    S1 --> L1
    S2 --> L2
    S3 --> L3
    S4 --> M1
    S4 --> M2
    M1 --> G1
    M1 --> G2

    style AnnData fill:#E6E6FA
```

## Layers Created

| Layer | Source | Description |
|-------|--------|-------------|
| `raw` | Normalization | Normalized, not aligned |
| `aligned` | Alignment | Cross-sample aligned |
| `batchcorr` | Batch Correction | Batch-corrected |

## Spatial Graphs

Builds k-NN or radius-based spatial graphs for each sample.

## CLI

```bash
celltype-refinery preprocess merge \
  --input corrected/ \
  --graph-mode knn \
  --out output/merged.h5ad
```
