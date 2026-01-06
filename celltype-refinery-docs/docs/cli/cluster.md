---
sidebar_position: 3
---

# cluster

Run Leiden clustering.

```mermaid
flowchart LR
    subgraph Input
        A[merged.h5ad]
    end

    subgraph Params["Key Parameters"]
        P1[--resolution 0.6]
        P2[--n-pcs 30]
        P3[--neighbors-k 15]
        P4[--use-gpu]
    end

    subgraph Output
        O[clustered.h5ad<br/>+ cluster_lvl0]
    end

    A --> Params
    Params --> O
```

## Usage

```bash
celltype-refinery cluster [OPTIONS]
```

## Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--input` | PATH | - | Input H5AD file |
| `--resolution` | FLOAT | 0.6 | Leiden resolution |
| `--n-pcs` | INT | 30 | Principal components |
| `--neighbors-k` | INT | 15 | k for k-NN |
| `--use-gpu` | FLAG | False | Use GPU |
| `--out` | PATH | - | Output directory |
