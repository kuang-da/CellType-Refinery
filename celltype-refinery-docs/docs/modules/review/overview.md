---
sidebar_position: 1
---

# Review Overview

The review module performs quality control and validation of cell-type annotations using configurable flagging rules.

## Module Architecture

```mermaid
flowchart LR
    subgraph Input
        A[composition/] --> E
        B[spatial/] --> E
        C[stage_n/] --> E
        D[tissue_template.yaml] --> E
    end

    subgraph ReviewEngine
        E[MetricsAggregator] --> F[RuleRegistry]
        F --> G[Proportion Rules]
        F --> H[Spatial Rules]
        F --> I[Biology Rules]
        F --> J[Quality Rules]
    end

    subgraph Output
        G --> K[flagged_issues.csv]
        H --> K
        I --> K
        J --> K
        K --> L[review_report.md]
        K --> M[review_summary.json]
    end
```

## Flagging Decision Flow

```mermaid
flowchart TB
    A[Load Metrics] --> B[Execute Rule Registry]
    B --> C{For each rule}
    C --> D{Check condition}
    D -->|Violated| E[Create FlaggedIssue]
    D -->|OK| F[Skip]
    E --> G{Severity?}
    G -->|Critical| H[CRITICAL: Action Required]
    G -->|Warning| I[WARNING: Review Recommended]
    G -->|Note| J[NOTE: For Information]
    H --> K[Aggregate Issues]
    I --> K
    J --> K
    F --> C
    K --> L[Generate Report]
```

## Features

- **Configurable flagging rules**: 17 rules across proportion, spatial, biology, and quality categories
- **Multi-column comparison**: Compare annotations across `cell_type_phenocycler`, `cell_type_multiomics`, etc.
- **Issue tracking**: Structured output with severity levels (critical, warning, note)
- **Tissue templates**: Expected patterns defined in YAML for tissue-specific validation
