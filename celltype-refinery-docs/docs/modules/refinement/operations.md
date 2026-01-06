---
sidebar_position: 3
---

# Refinement Operations

Available operations in refinement plans.

```mermaid
flowchart LR
    subgraph Override
        O1[Cluster] --> O2[Direct<br/>Label]
    end

    subgraph Merge
        M1[Cluster A]
        M2[Cluster B]
        M1 --> M3[Combined<br/>Label]
        M2 --> M3
    end

    subgraph Subcluster
        S1[Cluster] --> S2[Re-cluster<br/>at 0.3]
        S2 --> S3[Sub-labels<br/>0:1, 0:2]
    end

    subgraph Relabel
        R1[Old Label] --> R2[New Label]
    end

    subgraph Rescore
        RS1[Cluster] --> RS2[Recompute<br/>Scores]
    end
```

## Operation Types

| Type | Description |
|------|-------------|
| Override | Direct label assignment |
| Merge | Combine clusters |
| Subcluster | Re-cluster at finer resolution |
| Relabel | Rename without re-clustering |
| Rescore | Recompute marker scores |
