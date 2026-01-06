---
sidebar_position: 3
---

# Cell Quality Control

Remove low-quality cells based on morphology and intensity metrics.

```mermaid
flowchart TD
    A[Input Cell] --> B{Area OK?<br/>50-5000 px}
    B -->|No| FAIL[Remove Cell]
    B -->|Yes| C{Nucleus Ratio?<br/>0.1-0.9}
    C -->|No| FAIL
    C -->|Yes| D{Total Intensity?<br/>>100}
    D -->|No| FAIL
    D -->|Yes| E{Autofluorescence?<br/><500}
    E -->|No| FAIL
    E -->|Yes| PASS[Keep Cell]

    style PASS fill:#90EE90
    style FAIL fill:#FFB6C1
```

## QC Criteria

| Metric | Default Threshold | Description |
|--------|-------------------|-------------|
| Cell area | 50-5000 px | Remove debris and doublets |
| Nucleus ratio | 0.1-0.9 | Nuclear/cell area ratio |
| Total intensity | \>100 | Minimum signal |
| Autofluorescence | \<500 | Maximum background |

## CLI

```bash
celltype-refinery preprocess qc \
  --input loaded/ \
  --max-removal 0.15 \
  --out output/filtered
```
