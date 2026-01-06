---
sidebar_position: 2
---

# preprocess

Run preprocessing stages.

```mermaid
flowchart LR
    subgraph Stages["--stage options"]
        A[load] --> B[qc]
        B --> C[normalize]
        C --> D[align]
        D --> E[batch]
        E --> F[merge]
    end

    subgraph Output
        O[Preprocessed<br/>Data]
    end

    F --> O
```

## Usage

```bash
celltype-refinery preprocess [OPTIONS]
```

## Options

| Option | Type | Description |
|--------|------|-------------|
| `--input` | PATH | Input data directory |
| `--config` | PATH | Configuration file |
| `--out` | PATH | Output directory |
| `--stage` | TEXT | Specific stage to run |
