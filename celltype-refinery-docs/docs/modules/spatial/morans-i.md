---
sidebar_position: 3
---

# Moran's I

Moran's I measures spatial autocorrelation—whether cells of the same type cluster together or are dispersed in space.

## Computation Flow

```mermaid
flowchart LR
    subgraph Input
        A[Cell positions] --> C
        B[Cell type labels] --> C
    end

    subgraph Computation
        C[Build spatial graph] --> D[Binary cell-type indicator]
        D --> E[Compute local Moran's I]
        E --> F[Aggregate to global I]
    end

    subgraph Output
        F --> G[Moran's I per cell type]
        F --> H[Z-score / p-value]
    end
```

## Interpretation Decision Tree

```mermaid
flowchart TB
    A[Moran's I Value] --> B{I > 0.3?}
    B -->|Yes| C[Strong Clustering]
    B -->|No| D{I > 0?}
    D -->|Yes| E[Weak Clustering]
    D -->|No| F{I ≈ 0?}
    F -->|Yes| G[Random Distribution]
    F -->|No| H[Dispersion]

    C --> I[Expected for: Epithelium, Tumor nests]
    E --> J[Expected for: Fibroblasts, Endothelium]
    G --> K[Unexpected: May indicate annotation issues]
    H --> L[Expected for: Rare immune cells, Mobile populations]
```

## Interpretation

| Value | Meaning | Biological Example |
|-------|---------|-------------------|
| \> 0.5 | Strong clustering | Epithelial sheets, tumor nests |
| 0.1–0.5 | Moderate clustering | Fibroblast networks |
| ≈ 0 | Random | Well-mixed immune infiltrate |
| \< 0 | Dispersion | Patrolling immune cells |

## Statistical Significance

The z-score indicates whether the observed Moran's I is significantly different from random:

- **|z| > 2.58**: Significant at p < 0.01
- **|z| > 1.96**: Significant at p < 0.05
- **|z| < 1.96**: Not significant
