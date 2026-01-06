---
sidebar_position: 1
---

# Composition Overview

The composition module computes cell-type statistics, diversity metrics, and tissue-specific biology metrics from annotated single-cell data.

## Module Architecture

```mermaid
flowchart LR
    subgraph Input
        A[consolidated.h5ad]
        B[config.yaml]
    end

    subgraph CompositionEngine
        A --> C[Load AnnData]
        B --> D[Load Config]
        C --> E[Aggregation]
        D --> E
        E --> F[Diversity]
        E --> G[Biology]
        E --> H[Enrichment]
    end

    subgraph Outputs
        F --> I[diversity_by_sample.csv]
        G --> J[ft_biology_by_region.csv]
        H --> K[regional_enrichment.csv]
        E --> L[composition_by_*.csv]
    end
```

## Analysis Pipeline

```mermaid
flowchart TB
    subgraph Aggregation
        A1[Count cells per type] --> A2[Group by sample]
        A2 --> A3[Group by region]
        A3 --> A4[Group by donor]
        A4 --> A5[Compute proportions]
    end

    subgraph Diversity
        A5 --> D1[Shannon Entropy]
        A5 --> D2[Simpson Index]
        A5 --> D3[Pielou Evenness]
    end

    subgraph Biology
        A5 --> B1[E:S Ratio]
        A5 --> B2[Ciliated:Secretory]
        A5 --> B3[Immune Infiltration]
    end

    subgraph Enrichment
        A3 --> E1[Mann-Whitney U]
        E1 --> E2[Multiple Testing Correction]
        E2 --> E3[Fold Change + p-value]
    end
```

## Features

- **Composition by sample/region/donor**: Hierarchical aggregation of cell-type counts and proportions
- **Diversity metrics**: Shannon entropy, Simpson index, Pielou's evenness per sample
- **FT biology metrics**: Epithelial:stromal ratio, ciliated:secretory ratio, immune infiltration
- **Regional enrichment tests**: Mann-Whitney U tests with Bonferroni/FDR correction
