---
sidebar_position: 7
---

# Iterative Refinement

Stage I supports iterative refinement chains until convergence.

## Refinement Chain

```text
H → I (iteration 1) → I (iteration 2) → I (iteration 3) → ...
```

Each iteration:
1. Uses previous output's `cluster_lvl1` as new `cluster_lvl0`
2. Selects new candidates based on updated scores
3. Further refines clusters that still need improvement
4. Tracks full lineage in provenance

## Convergence

Typically converges in 2-3 iterations when no more candidates meet criteria:
- All clusters above score threshold
- No heterogeneous parent signals remain
- All clusters below min_cells threshold

## Workflow Example

### Iteration 1: Initial Refinement

```bash
celltype-refinery refine \
  --input stage_h/coarse_clusters.h5ad \
  --auto --execute \
  --out stage_i_v1
```

### Iteration 2: Further Refinement

```bash
celltype-refinery refine \
  --input stage_i_v1/refined.h5ad \
  --auto --execute \
  --out stage_i_v2
```

### Iteration 3: Targeted Manual Fixes

```bash
celltype-refinery refine \
  --input stage_i_v2/refined.h5ad \
  --config manual_fixes.yaml --execute \
  --out stage_i_final
```

## Cluster ID Evolution

Subclustering creates hierarchical IDs:

```text
Original cluster "3" subclustered:
  → "3:0", "3:1", "3:2"

Further subclustering "3:1":
  → "3:1:0", "3:1:1"
```

## Provenance Tracking

Each iteration stores provenance in `adata.uns["stage_refine_iteration_N"]`:

```python
{
  "parent_stage": "I",
  "parent_file": "stage_i_v1/refined.h5ad",
  "iteration": 2,
  "plan_summary": {"subcluster": 5, "relabel": 3},
  "clusters_modified": ["3", "7", "12"],
  "subclusters_created": ["3:0", "3:1", "7:0", "7:1"],
  "n_cells_modified": 15420,
  "timestamp": "2025-01-12T14:30:00"
}
```

## Best Practices

1. **Start with diagnostic mode** to review recommendations
2. **Use focus controls** for targeted refinement per cell type category
3. **Save each iteration** for reproducibility
4. **Monitor cluster counts** - should stabilize by iteration 3
5. **Check provenance** to understand refinement history

## See Also

- [Focus Controls](/docs/modules/refinement/focus-controls.md)
- [Refinement Decision Logic](/docs/methodology/refinement-decision-logic.md)
- [Diagnostic Mode](/docs/modules/refinement/diagnostic-mode.md)
