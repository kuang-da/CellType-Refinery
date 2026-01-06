---
sidebar_position: 5
---

# refine

Run iterative refinement.

```mermaid
flowchart TD
    subgraph Modes["Execution Modes"]
        M1["--auto<br/>Automatic policy"]
        M2["--config yaml<br/>Manual policy"]
        M3["--auto --config<br/>Hybrid"]
    end

    subgraph Execute
        D{--execute?}
        D -->|No| DIAG[Diagnostic<br/>Report Only]
        D -->|Yes| RUN[Apply<br/>Refinements]
    end

    M1 --> D
    M2 --> D
    M3 --> D

    style DIAG fill:#FFFACD
    style RUN fill:#90EE90
```

## Usage

```bash
celltype-refinery refine [OPTIONS]
```

## Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--input` | PATH | - | Input H5AD file |
| `--auto` | FLAG | False | Enable auto policy |
| `--config` | PATH | - | Manual config YAML |
| `--execute` | FLAG | False | Execute refinement |
| `--score-threshold` | FLOAT | 1.0 | Score threshold |
| `--min-cells` | INT | 500 | Minimum cells |
| `--out` | PATH | - | Output directory |

## Examples

```bash
# Diagnostic mode
celltype-refinery refine --input annotated.h5ad --auto --out output/

# Execute mode
celltype-refinery refine --input annotated.h5ad --auto --execute --out output/

# Hybrid mode
celltype-refinery refine --input annotated.h5ad --auto --config curation.yaml --execute --out output/
```
