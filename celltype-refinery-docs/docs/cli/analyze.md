---
sidebar_position: 7
---

# analyze

Run composition and spatial analysis.

```mermaid
flowchart LR
    subgraph Input
        A[consolidated.h5ad]
    end

    subgraph Flags["Analysis Flags"]
        F1["--composition<br/>Cell type stats"]
        F2["--spatial<br/>Neighborhood analysis"]
        F3["--review<br/>QC flagging"]
    end

    subgraph Output
        O1[composition/]
        O2[spatial/]
        O3[review/]
    end

    A --> F1
    A --> F2
    A --> F3
    F1 --> O1
    F2 --> O2
    F3 --> O3
```

## Usage

```bash
celltype-refinery analyze [OPTIONS]
```

## Options

| Option | Type | Description |
|--------|------|-------------|
| `--input` | PATH | Input H5AD file |
| `--composition` | FLAG | Run composition |
| `--spatial` | FLAG | Run spatial |
| `--review` | FLAG | Run review |
| `--out` | PATH | Output directory |
