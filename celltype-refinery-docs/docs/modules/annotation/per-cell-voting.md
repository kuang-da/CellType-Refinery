---
sidebar_position: 6
---

# Per-Cell Voting

Per-cell voting provides evidence for downstream refinement when cluster-level annotation is ambiguous.

:::warning Critical Concept
Per-cell voting is **EVIDENCE ONLY**. It does NOT override cluster-level labels.
The assigned_label always reflects the cluster-level decision.
:::

## When Voting Occurs

Voting is triggered when the gap between top candidates is below the `min_gap` threshold:
- Level 0: gap < 0.50
- Level 1: gap < 0.30
- Level 2+: gap < 0.20

## Algorithm

```text
┌─────────────────────────────────────────────────────────────────────────┐
│  PER-CELL VOTING                                                        │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                         │
│  For each cell in cluster:                                              │
│                                                                         │
│    For each candidate:                                                  │
│      cell_score[candidate] = mean(Z-scores for candidate's markers)     │
│                                                                         │
│    Sort candidates by cell_score                                        │
│    cell_gap = best_score - runner_up_score                              │
│                                                                         │
│    IF cell_gap < per_cell_gap_threshold (0.1):                          │
│      vote = "Uncertain"                                                 │
│    ELSE:                                                                │
│      vote = best_candidate                                              │
│                                                                         │
│  Aggregate votes into composition dictionary                            │
│                                                                         │
└─────────────────────────────────────────────────────────────────────────┘
```

## Output Format

The composition is stored as a JSON object:

```json
{
  "Ciliated Epithelium": 0.55,
  "Glandular Epithelium": 0.40,
  "Uncertain": 0.05,
  "_evidence_only": true,
  "_top_by_score": "Ciliated Epithelium",
  "_runner_up_by_score": "Glandular Epithelium"
}
```

## Important: Evidence vs Override

| Aspect | Cluster-Level | Per-Cell Voting |
|--------|---------------|-----------------|
| **Purpose** | Assign labels | Collect evidence |
| **Authority** | Authoritative | Advisory |
| **Modifies labels** | Yes | NO |
| **Used by** | All downstream | Stage I refinement |

## Usage in Stage I

Per-cell voting composition helps Stage I decide:
- **Homogeneous (e.g., 95% one type)**: Cluster is well-assigned
- **Heterogeneous (e.g., 55%/40% split)**: Cluster may need subclustering
- **High Uncertain (e.g., 30%)**: Markers may need improvement

## See Also

- [Hierarchical Gating Algorithm](/docs/methodology/hierarchical-gating-algorithm.md) - When voting is triggered
- [Refinement Decision Logic](/docs/methodology/refinement-decision-logic.md) - How voting informs refinement
