---
sidebar_position: 6
---

# Batch Effect Correction

Remove technical variation while preserving biological effects.

## Technical Effects Modeled

- Donor effects
- Imaging color
- Imaging cycle

## Biological Effects Preserved

- Marker effects
- Region effects
- Marker-region interactions

## CLI

```bash
celltype-refinery preprocess batch \
  --input aligned/ \
  --mode residual \
  --out output/corrected
```
