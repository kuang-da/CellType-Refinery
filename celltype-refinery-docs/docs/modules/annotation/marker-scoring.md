---
sidebar_position: 4
---

# Marker Scoring

Calculate marker expression scores per cluster.

```mermaid
flowchart LR
    subgraph Input
        A[Cluster<br/>Expression]
        B[Marker<br/>Definition]
    end

    subgraph Components["Score Components"]
        C[Mean Enrichment<br/>avg expression]
        D[Positive Fraction<br/>% expressing]
        E[DE Bonus<br/>+0.5 if top DE]
        F[Anti-Penalty<br/>-anti expression]
    end

    subgraph Formula["Final Score"]
        G["score = enrichment<br/>+ positive<br/>+ de_bonus<br/>- anti_penalty"]
    end

    A --> C
    A --> D
    A --> E
    B --> F
    C --> G
    D --> G
    E --> G
    F --> G

    style G fill:#FFD700
```

## Components

- **Mean Enrichment**: Average marker expression
- **Positive Fraction**: Cells expressing marker
- **DE Bonus**: Markers in top DE genes
- **Anti-Penalty**: Anti-marker expression

## CLI

```bash
celltype-refinery annotate score \
  --input clustered.h5ad \
  --marker-map markers.json \
  --out output/scores.csv
```
