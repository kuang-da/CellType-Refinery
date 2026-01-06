---
sidebar_position: 6
---

# consolidate

Run final consolidation.

```mermaid
flowchart LR
    subgraph Input
        A[refined.h5ad]
        B[--config yaml]
    end

    subgraph Features["Optional Features"]
        F1[--enable-orphan-rescue]
        F2[Manual overrides]
        F3[Label harmonization]
    end

    subgraph Output
        O[consolidated.h5ad<br/>cell_type_final]
    end

    A --> Features
    B --> Features
    Features --> O
```

## Usage

```bash
celltype-refinery consolidate [OPTIONS]
```

## Options

| Option | Type | Description |
|--------|------|-------------|
| `--input` | PATH | Input H5AD file |
| `--config` | PATH | Consolidation config |
| `--enable-orphan-rescue` | FLAG | Enable orphan rescue |
| `--out` | PATH | Output directory |
