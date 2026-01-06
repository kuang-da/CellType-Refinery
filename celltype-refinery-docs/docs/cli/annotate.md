---
sidebar_position: 4
---

# annotate

Run marker-based annotation.

```mermaid
flowchart LR
    subgraph Input
        A[clustered.h5ad]
        B[markers.json]
    end

    subgraph Process
        C[Load Markers] --> D[Score Clusters]
        D --> E[Hierarchical<br/>Gating]
    end

    subgraph Output
        F[annotated.h5ad]
        G[marker_scores.csv]
    end

    A --> C
    B --> C
    E --> F
    E --> G
```

## Usage

```bash
celltype-refinery annotate [OPTIONS]
```

## Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--input` | PATH | - | Input H5AD file |
| `--marker-map` | PATH | - | Marker map JSON |
| `--layer` | TEXT | batchcorr | Expression layer |
| `--min-coverage` | FLOAT | 0.3 | Minimum coverage |
| `--min-pos-frac` | FLOAT | 0.2 | Minimum positive fraction |
| `--out` | PATH | - | Output directory |
