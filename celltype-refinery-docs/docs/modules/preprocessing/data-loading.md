---
sidebar_position: 2
---

# Data Loading

Load and validate input data files.

```mermaid
flowchart LR
    subgraph Input["Input Files"]
        A[Metadata CSV]
        B[Cell Matrices]
        C[Spatial Coords]
    end

    subgraph Validation["Validation Checks"]
        D{Columns Match?}
        E{No Duplicates?}
        F{Coords Valid?}
    end

    subgraph Output["Output"]
        G[Verified Metadata]
        H[Validated Matrices]
    end

    A --> D
    B --> D
    D -->|Yes| E
    D -->|No| X[Error: Column Mismatch]
    E -->|Yes| F
    E -->|No| Y[Error: Duplicate IDs]
    C --> F
    F -->|Yes| G
    F -->|Yes| H
    F -->|No| Z[Error: Invalid Coords]
```

## Supported Formats

- CSV cell matrices
- TSV metadata files
- H5AD AnnData objects

## Validation Checks

- Column consistency across samples
- Duplicate cell IDs
- Missing metadata fields
- Coordinate anomalies

## CLI

```bash
celltype-refinery preprocess load \
  --metadata metadata.csv \
  --out output/loaded
```
