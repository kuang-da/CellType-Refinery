---
sidebar_position: 6
---

# Focus Controls

Focus controls restrict refinement to a subset of clusters, enabling iterative targeted refinement.

## Overview

Two focus mechanisms are available:
- `--focus-labels`: Match clusters by their assigned label
- `--focus-clusters`: Match clusters by explicit ID

## --focus-labels

Match clusters by assigned label (exact match or path prefix):

```bash
celltype-refinery refine \
  --input stage_h/coarse_clusters.h5ad \
  --auto --execute \
  --focus-labels "Immune Cells,Epithelium"
```

### Matching Behavior

```text
--focus-labels "Immune Cells" matches:
  ✓ "Immune Cells"
  ✓ "Immune Cells / Myeloids"
  ✓ "Immune Cells / Myeloids / Macrophage"
  ✗ "Epithelium"
  ✗ "Epithelium / Ciliated"
```

## --focus-clusters

Match clusters by explicit ID:

```bash
celltype-refinery refine \
  --input stage_h/coarse_clusters.h5ad \
  --auto --execute \
  --focus-clusters "3,5,7,12"
```

Only clusters 3, 5, 7, and 12 are eligible for refinement.

## Combined Usage (Intersection)

When both are specified, the **intersection** is used:

```bash
--focus-labels "Immune Cells" --focus-clusters "3,5,7,10"
```

Only clusters that are BOTH:
1. Assigned to "Immune Cells" (or descendant), AND
2. In the explicit cluster set {3, 5, 7, 10}

## Iterative Workflow Example

```text
Iteration 1: Refine Immune cells only
  --focus-labels "Immune Cells"

Iteration 2: Refine Epithelium only
  --focus-labels "Epithelium"

Iteration 3: Refine specific problematic clusters
  --focus-clusters "15,23,31"
```

## Validation

- Operations targeting out-of-scope clusters raise `ValueError`
- Invalid cluster IDs in `--focus-clusters` are warned (not errors)
- ManualPolicy operations are validated against focus controls

## See Also

- [Iterative Refinement](/docs/modules/refinement/iterative-refinement.md)
- [Refinement Decision Logic](/docs/methodology/refinement-decision-logic.md)
