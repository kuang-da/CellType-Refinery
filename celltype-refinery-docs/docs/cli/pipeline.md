---
sidebar_position: 8
---

# pipeline

Run full pipeline from configuration.

```mermaid
flowchart LR
    subgraph Input
        C[pipeline.yaml]
        M[--mode auto/manual/hybrid]
    end

    subgraph Pipeline["Pipeline Runner"]
        P[Parse Config] --> V[Validate]
        V --> E[Execute Stages]
    end

    subgraph Output
        O[Full Pipeline<br/>Output]
    end

    C --> P
    M --> P
    E --> O
```

## Usage

```bash
celltype-refinery pipeline [OPTIONS]
```

## Options

| Option | Type | Description |
|--------|------|-------------|
| `--config` | PATH | Pipeline configuration YAML |
| `--mode` | TEXT | Execution mode (auto, manual, hybrid) |
| `--dry-run` | FLAG | Show plan without executing |

## Example

```bash
celltype-refinery pipeline --config pipeline.yaml --mode auto
```
