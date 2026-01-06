---
sidebar_position: 3
---

# Pipeline YAML

Configuration for full pipeline runs.

```mermaid
flowchart LR
    subgraph Config["pipeline.yaml"]
        V[version]
        N[name]
        P[paths]
        S[stages]
    end

    subgraph Paths
        I[input]
        O[output]
        M[marker_map]
    end

    subgraph Stages["Stage Toggles"]
        C[clustering<br/>enabled: true]
        A[annotation<br/>enabled: true]
        R[refinement<br/>enabled: true]
    end

    P --> I
    P --> O
    P --> M
    S --> C
    S --> A
    S --> R

    style Config fill:#E6E6FA
```

## Example

```yaml
version: "1.0"
name: "my_pipeline"

paths:
  input: "data/merged.h5ad"
  output: "output"
  marker_map: "markers.json"

stages:
  clustering:
    enabled: true
    resolution: 0.6
  annotation:
    enabled: true
  refinement:
    enabled: true
    mode: "auto"
```
