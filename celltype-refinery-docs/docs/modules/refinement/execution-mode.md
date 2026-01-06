---
sidebar_position: 5
---

# Execution Mode

Apply refinement plans to the data.

```mermaid
flowchart LR
    subgraph Input
        A[annotated.h5ad]
        B[Refinement Plan]
    end

    subgraph Engine["Refinement Engine"]
        C[Apply Operations<br/>in order]
    end

    subgraph Output
        D[refined.h5ad<br/>Updated labels]
        E[curation_log.json<br/>Operation history]
    end

    A --> C
    B --> C
    C --> D
    C --> E

    style Engine fill:#E6E6FA
```

## CLI

```bash
celltype-refinery refine \
  --input annotated.h5ad \
  --auto \
  --execute \
  --out output/refined
```

## Output

- `refined.h5ad`: AnnData with refined labels
- `curation_log.json`: Operation history
