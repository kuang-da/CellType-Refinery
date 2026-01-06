---
sidebar_position: 1
---

# CLI Overview

Command-line interface for CellType-Refinery.

```mermaid
flowchart LR
    subgraph CLI["celltype-refinery"]
        A[preprocess] --> B[cluster]
        B --> C[annotate]
        C --> D[refine]
        D --> E[consolidate]
        E --> F[analyze]
    end

    subgraph Alt["Alternative"]
        P[pipeline<br/>--config yaml]
    end

    P -.-> A
    P -.-> F
```

## Installation

The CLI is installed with the package:

```bash
pip install -e .
celltype-refinery --help
```

## Commands

| Command | Description |
|---------|-------------|
| [preprocess](preprocess) | Run preprocessing stages |
| [cluster](cluster) | Leiden clustering |
| [annotate](annotate) | Marker-based annotation |
| [refine](refine) | Iterative refinement |
| [consolidate](consolidate) | Final consolidation |
| [analyze](analyze) | Composition and spatial |
| [pipeline](pipeline) | Full pipeline from config |

## Global Options

```bash
celltype-refinery --help
celltype-refinery --version
celltype-refinery --verbose COMMAND
celltype-refinery --debug COMMAND
```
