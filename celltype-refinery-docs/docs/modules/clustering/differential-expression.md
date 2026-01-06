---
sidebar_position: 3
---

# Differential Expression

Identify cluster-specific markers via DE analysis.

```mermaid
flowchart LR
    subgraph Input
        A[Clustered<br/>AnnData]
    end

    subgraph Analysis["DE Analysis"]
        B[For each cluster]
        C[Cluster vs Rest<br/>Wilcoxon test]
        D[Rank by<br/>log fold change]
        E[Select top N<br/>genes]
    end

    subgraph Output
        F[de_genes.csv<br/>Top markers]
        G[de_stats.csv<br/>Full statistics]
    end

    A --> B
    B --> C
    C --> D
    D --> E
    E --> F
    E --> G
```

## Method

Wilcoxon rank-sum test per cluster vs. rest.

## Output

- `de_genes.csv`: Top DE genes per cluster
- `de_stats.csv`: Full statistics

## CLI

```bash
celltype-refinery cluster de \
  --input clustered.h5ad \
  --n-genes 50 \
  --out output/
```
