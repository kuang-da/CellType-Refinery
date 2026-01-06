---
sidebar_position: 4
---

# Diagnostic Mode

Preview refinement recommendations without executing.

```mermaid
flowchart TD
    A[Input:<br/>annotated.h5ad] --> B[Analyze<br/>Clusters]

    B --> C{Score<br/>High?}
    C -->|Yes| SKIP[SKIP<br/>Good annotation]
    C -->|No| D{Large<br/>cluster?}

    D -->|Yes| SUB[SUBCLUSTER<br/>Split heterogeneous]
    D -->|No| REL[RELABEL<br/>Update label]

    subgraph Output["diagnostic_report.csv"]
        SKIP
        SUB
        REL
    end

    style SKIP fill:#90EE90
    style SUB fill:#FFD700
    style REL fill:#ADD8E6
```

## Output

- `diagnostic_report.csv`: Per-cluster recommendations
- Recommendations: SUBCLUSTER, RELABEL, SKIP

## CLI

```bash
celltype-refinery refine \
  --input annotated.h5ad \
  --auto \
  --out output/diagnostic
```
