---
sidebar_position: 5
---

# Cross-Sample Alignment

Align marker intensities across samples.

```mermaid
flowchart LR
    subgraph Samples["Per-Sample Data"]
        S1[Sample 1]
        S2[Sample 2]
        S3[Sample N]
    end

    subgraph Reference["Global Reference"]
        R1[Compute Q5<br/>per marker]
        R2[Compute Q95<br/>per marker]
    end

    subgraph Mapping["Percentile Mapping"]
        M[Linear Scale<br/>Q5 → 0<br/>Q95 → 1]
    end

    subgraph Output["Aligned Output"]
        O[Unified<br/>Intensity Scale]
    end

    S1 --> R1
    S2 --> R1
    S3 --> R1
    S1 --> R2
    S2 --> R2
    S3 --> R2
    R1 --> M
    R2 --> M
    M --> O
```

## Method

Global percentile mapping using Q5 and Q95 reference values.

## CLI

```bash
celltype-refinery preprocess align \
  --input normalized/ \
  --out output/aligned
```
