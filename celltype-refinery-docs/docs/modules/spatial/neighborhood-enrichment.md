---
sidebar_position: 2
---

# Neighborhood Enrichment

Which cell types co-locate more than expected?

```mermaid
flowchart LR
    subgraph Input
        A[Spatial<br/>Graph]
        B[Cell Type<br/>Labels]
    end

    subgraph Analysis
        C[Count neighbors<br/>per type pair]
        D[Permute labels<br/>N times]
        E[Compare observed<br/>vs expected]
    end

    subgraph Output
        F["Z-score Matrix<br/>Type A × Type B"]
    end

    A --> C
    B --> C
    C --> D
    D --> E
    E --> F

    subgraph Interpretation
        G["+Z = Co-localize"]
        H["-Z = Avoid"]
    end

    F --> G
    F --> H

    style G fill:#90EE90
    style H fill:#FFB6C1
```

## Method

Permutation-based Z-score for cell-type pairs.

## Output

Z-score matrix: positive = co-localization, negative = avoidance.
