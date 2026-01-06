---
sidebar_position: 4
---

# Within-Sample Normalization

Normalize expression values within each sample.

```mermaid
flowchart LR
    subgraph Input
        A[Raw Intensities]
    end

    subgraph Step1["Step 1: Background Correction"]
        B[Quantile Clipping]
        B1[Q3 Clip]
        B2[Q5 Clip]
    end

    subgraph Step2["Step 2: Variance Stabilization"]
        C1[log1p]
        C2[asinh c=5]
        C3[asinh c=10]
    end

    subgraph Output
        D[Normalized Matrix]
    end

    A --> B
    B --> B1
    B --> B2
    B1 --> C1
    B1 --> C2
    B2 --> C2
    B2 --> C3
    C1 --> D
    C2 --> D
    C3 --> D
```

**Default variant**: `clip-Q3__log1p` (Q3 clipping + log1p transform)

## Methods

### Background Correction

- Quantile-based clipping (Q3, Q5)
- Centering

### Variance Stabilization

- `asinh` transform (cofactor 5 or 10)
- `log1p` transform

## Default Variant

`clip-Q3__log1p`: Q3 quantile clipping + log1p transform

## CLI

```bash
celltype-refinery preprocess normalize \
  --input filtered/ \
  --method clip-Q3__log1p \
  --out output/normalized
```
