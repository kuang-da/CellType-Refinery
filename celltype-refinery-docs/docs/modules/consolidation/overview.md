---
sidebar_position: 1
---

# Consolidation Overview

The consolidation module finalizes cell-type annotations after iterative refinement, producing a single `cell_type_final` column with production-ready labels.

## Module Architecture

```mermaid
flowchart LR
    subgraph Inputs
        A[refined.h5ad] --> E
        B[diagnostic_report.csv] --> E
        C[marker_scores.csv] --> E
        D[config.yaml] --> E
    end

    subgraph Engine
        E[ConsolidationEngine] --> F{Orphan Rescue?}
        F -->|Yes| G[Detect Orphan Subtypes]
        F -->|No| H[Skip Rescue]
        G --> I{IEL Rescue?}
        H --> I
        I -->|Yes| J[Detect IEL Candidates]
        I -->|No| K[Skip IEL]
        J --> L[Harmonize Labels]
        K --> L
    end

    subgraph Outputs
        L --> M[consolidated.h5ad]
        L --> N[orphan_candidates.csv]
        L --> O[provenance.json]
    end
```

## Decision Flow

```mermaid
flowchart TB
    A[Cell from Refinement] --> B{Has curated label?}
    B -->|Yes| C[Keep curated label]
    B -->|No| D{Strong subtype markers?}
    D -->|Yes| E{Weak root markers?}
    D -->|No| F[Assign Unassigned]
    E -->|Yes| G[Orphan Candidate]
    E -->|No| C
    G --> H{Plausible lineage?}
    H -->|Yes| I[Rescue as orphan subtype]
    H -->|No| F
    C --> J[cell_type_final]
    I --> J
    F --> J
```

## Features

- **Single output column**: `cell_type_final` with clean, production-ready labels
- **Orphan rescue**: Detects cells with strong subtype markers but weak root markers
- **IEL rescue**: Recovers intraepithelial lymphocytes from epithelial clusters
- **Label harmonization**: Standardizes naming conventions across annotations
