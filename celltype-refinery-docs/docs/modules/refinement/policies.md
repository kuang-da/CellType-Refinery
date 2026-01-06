---
sidebar_position: 2
---

# Refinement Policies

Policies generate refinement plans.

```mermaid
flowchart TD
    subgraph Policies["Policy Selection"]
        direction TB
        AUTO["AutoPolicy<br/>Automatic"]
        MANUAL["ManualPolicy<br/>YAML Config"]
    end

    subgraph AutoCriteria["Auto Criteria"]
        A1[Score < threshold]
        A2[Cluster size > min]
        A3[High DE heterogeneity]
    end

    subgraph ManualOps["Manual Operations"]
        M1[Overrides]
        M2[Merges]
        M3[Subclusters]
        M4[Relabels]
    end

    subgraph Output
        PLAN[Refinement<br/>Plan]
    end

    AUTO --> A1
    AUTO --> A2
    AUTO --> A3
    MANUAL --> M1
    MANUAL --> M2
    MANUAL --> M3
    MANUAL --> M4

    A1 --> PLAN
    A2 --> PLAN
    A3 --> PLAN
    M1 --> PLAN
    M2 --> PLAN
    M3 --> PLAN
    M4 --> PLAN

    style AUTO fill:#ADD8E6
    style MANUAL fill:#FFD700
```

## AutoPolicy

Automatic candidate selection based on:
- Score thresholds
- Cluster size
- DE heterogeneity

## ManualPolicy

Parse YAML configuration files for:
- Overrides
- Merges
- Subclusters
- Relabels
