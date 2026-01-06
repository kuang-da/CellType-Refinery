---
sidebar_position: 2
---

# Diversity Metrics

Quantify cell-type diversity.

```mermaid
flowchart LR
    subgraph Input
        A[Cell Type<br/>Proportions]
    end

    subgraph Metrics["Diversity Metrics"]
        B["Shannon Entropy<br/>-Σ pᵢ log pᵢ"]
        C["Simpson Index<br/>1 - Σ pᵢ²"]
        D["Evenness<br/>H / H_max"]
    end

    subgraph Interpretation
        E[Higher = More<br/>Diverse]
        F[Lower = More<br/>Dominant]
        G[1.0 = Perfectly<br/>Even]
    end

    A --> B
    A --> C
    A --> D
    B --> E
    C --> F
    D --> G
```

## Metrics

| Metric | Formula | Interpretation |
|--------|---------|----------------|
| Shannon | -Σ pᵢ log pᵢ | Overall diversity |
| Simpson | 1 - Σ pᵢ² | Dominance inverse |
| Evenness | H / H_max | Distribution equality |
