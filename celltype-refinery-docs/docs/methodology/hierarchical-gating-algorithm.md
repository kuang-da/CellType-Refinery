---
sidebar_position: 4
---

# Hierarchical Gating Algorithm

:::info Prerequisites
This page explains the assignment algorithm. For scoring details, see [Marker Scoring Algorithm](./marker-scoring-algorithm.md) first.
:::

## Overview: Top-Down Traversal

The hierarchical gating algorithm is the core decision engine that assigns cell type labels to clusters. It traverses a marker hierarchy from roots to leaves, making competitive decisions at each level based on marker expression evidence.

**Key principles:**

- **Start at root categories** (Epithelium, Immune, Mesenchymal, etc.)
- **Descend through the hierarchy** following the strongest signal
- **Only siblings compete at each level** - unrelated branches never compete directly
- **Stop when uncertain** - ambiguous cases preserve the most specific confident label

The algorithm prioritizes accuracy over specificity. It is better to assign a cluster as "Immune" with high confidence than to force a guess between "T Cell" and "B Cell" when the evidence is insufficient.

### Marker Hierarchy Structure

The hierarchy defines parent-child relationships between cell type categories. Each node has associated positive markers (and optionally anti-markers) that define its expression signature.

```text
                    ┌─────────┐
         ┌──────────┤  ROOT   ├──────────┐
         ▼          └─────────┘          ▼
    ┌──────────┐                   ┌───────────┐
    │Epithelium│                   │ Immune    │  ... (other roots)
    └────┬─────┘                   │  Cells    │
         │                         └─────┬─────┘
    ┌────┴────┐                    ┌─────┴─────┐
    ▼         ▼                    ▼           ▼
 Ciliated  Glandular           Myeloids    Lymphoids
                                   │
                              ┌────┴────┐
                              ▼         ▼
                         Granulocytes  Monocytes
```

In this example:
- **Epithelium** and **Immune Cells** are root nodes - they compete at level 0
- **Ciliated** and **Glandular** are siblings under Epithelium - they compete at level 1
- **Myeloids** and **Lymphoids** are siblings under Immune Cells - they compete at level 1
- **Granulocytes** and **Monocytes** are siblings under Myeloids - they compete at level 2

A cluster assigned to "Myeloids" will never be compared against "Ciliated" because they exist in different branches of the hierarchy.

---

## Gate Check Sequence

Before a node can be considered as a candidate for assignment, it must pass a series of gate checks. These gates filter out nodes that lack sufficient evidence, ensuring only well-supported candidates compete for selection.

```text
┌─────────────────────────────────────────────────────────────────────────┐
│  GATE CHECKS (in order)                                                 │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                         │
│  Level 0 only:                                                          │
│  ╔═══════════════════════════════════════════════════════════════════╗  │
│  ║ CHECK 1: root_hard_requirements                                   ║  │
│  ║   IF marker pos_frac < required → REJECT                          ║  │
│  ║   IF marker enrichment < required → REJECT                        ║  │
│  ╠═══════════════════════════════════════════════════════════════════╣  │
│  ║ CHECK 2: root_veto_markers                                        ║  │
│  ║   IF veto marker pos_frac > max_allowed → REJECT                  ║  │
│  ╚═══════════════════════════════════════════════════════════════════╝  │
│                                                                         │
│  All levels:                                                            │
│  ╔═══════════════════════════════════════════════════════════════════╗  │
│  ║ CHECK 3: coverage                                                 ║  │
│  ║   IF coverage < min_coverage[level] → REJECT                      ║  │
│  ╠═══════════════════════════════════════════════════════════════════╣  │
│  ║ CHECK 4: positive fraction (OR gate)                              ║  │
│  ║   pos_ok = mean_positive_fraction >= min_pos_frac[level]          ║  │
│  ║   frac_ok = frac_markers_on >= min_frac_markers_on[level]         ║  │
│  ║   IF NOT pos_ok AND NOT frac_ok → REJECT                          ║  │
│  ╠═══════════════════════════════════════════════════════════════════╣  │
│  ║ CHECK 5: enrichment                                               ║  │
│  ║   IF mean_enrichment < min_enrichment[level] → REJECT             ║  │
│  ╠═══════════════════════════════════════════════════════════════════╣  │
│  ║ CHECK 6: anti-marker conflict                                     ║  │
│  ║   IF anti_penalty > anti_penalty_hard_gate (1.0) → REJECT         ║  │
│  ╚═══════════════════════════════════════════════════════════════════╝  │
│                                                                         │
│  All checks passed → ACCEPT                                             │
│                                                                         │
└─────────────────────────────────────────────────────────────────────────┘
```

### Gate Check Details

**CHECK 1: Root Hard Requirements** (Level 0 only)

Some root categories have mandatory markers that must be expressed. For example, Immune cells typically require CD45 expression. If a cluster lacks CD45, it cannot be assigned to Immune regardless of other marker evidence.

- **pos_frac threshold**: Minimum fraction of cells expressing the required marker
- **enrichment threshold**: Minimum log-fold enrichment vs. global background

**CHECK 2: Root Veto Markers** (Level 0 only)

Certain markers actively exclude assignment to a root. For example, high Pan-Cytokeratin expression (an epithelial marker) should veto Immune assignment. This prevents lineage misclassification even when some markers are shared.

**CHECK 3: Coverage**

Coverage measures what fraction of the node's expected markers are actually present in the panel. Low coverage means insufficient evidence to make a decision. The threshold relaxes at deeper levels because specialized subtypes often have fewer defining markers.

**CHECK 4: Positive Fraction (OR Gate)**

This is a dual-condition gate - either condition passing is sufficient:
- **pos_ok**: The average positive fraction across markers exceeds the threshold
- **frac_ok**: A sufficient fraction of individual markers are "on" (above their own thresholds)

The OR logic allows assignment when either the overall signal is strong OR when multiple markers show consistent (if weaker) positivity.

**CHECK 5: Enrichment**

Mean enrichment must exceed the level-specific threshold. Negative thresholds at deeper levels allow assignment even when expression is slightly below global average, which can occur for rare cell types.

**CHECK 6: Anti-Marker Conflict**

Anti-markers are markers that should be OFF for a given cell type. If anti-markers show strong expression (anti_penalty > 1.0), the node is rejected. This prevents assigning a cluster to a cell type when contradictory evidence is present.

---

## Base Gating Parameters

The gating thresholds vary by hierarchy level. Deeper levels have relaxed thresholds because:
- Specialized subtypes often have fewer and weaker distinguishing markers
- Once the lineage is established, less evidence is needed for subtype decisions
- The hierarchy structure itself provides contextual constraint

| Level | min_coverage | min_pos_frac | min_enrichment | min_gap | min_frac_markers_on |
|-------|--------------|--------------|----------------|---------|---------------------|
| 0 (Root) | 0.50 | 0.30 | 0.00 | 0.50 | 0.40 |
| 1 | 0.40 | 0.20 | -0.50 | 0.30 | 0.30 |
| 2 | 0.30 | 0.15 | -1.00 | 0.20 | 0.25 |
| 3+ | 0.30 | 0.15 | -1.00 | 0.20 | 0.25 |

### Parameter Interpretation

- **min_coverage**: Minimum fraction of expected markers present in panel
- **min_pos_frac**: Minimum mean positive fraction across markers
- **min_enrichment**: Minimum mean log-enrichment (negative values allow below-average expression)
- **min_gap**: Minimum score difference between winner and runner-up (for confident selection)
- **min_frac_markers_on**: Minimum fraction of markers that are individually "on"

These parameters can be adjusted via the configuration system. See the [Tuning Guide](./tuning-guide.md) for optimization strategies.

---

## Assignment Algorithm

The algorithm proceeds in a top-down fashion, first selecting the best root category, then descending through the hierarchy until reaching a leaf node or encountering uncertainty.

```text
┌─────────────────────────────────────────────────────────────────────────┐
│  HIERARCHICAL ASSIGNMENT ALGORITHM                                      │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                         │
│  For each cluster:                                                      │
│                                                                         │
│  1. IDENTIFY ROOTS                                                      │
│     • Extract labels that are never children (Epithelium, Immune, ...)  │
│                                                                         │
│  2. EVALUATE ROOTS                                                      │
│     • Check root_hard_requirements (e.g., CD45 for Immune)              │
│     • Check root_veto_markers (e.g., Pan-CK blocks Immune)              │
│     • Check standard gates (coverage, pos_frac, enrichment)             │
│     • Collect passing roots                                             │
│                                                                         │
│  3. SELECT BEST ROOT                                                    │
│     • Sort passing roots by score                                       │
│     • If margin < root_gap_threshold (0.25): ambiguous_root             │
│     • Otherwise: proceed with best root                                 │
│                                                                         │
│  4. DESCEND HIERARCHY                                                   │
│     • Get children of current node                                      │
│     • Check gates for each child (level-appropriate thresholds)         │
│     • Sort passing children by score                                    │
│     • If no children pass: stop_reason = "no_child_passed"              │
│     • If gap < min_gap[level]: stop_reason = "ambiguous_siblings"       │
│     • Otherwise: select best child, continue descent                    │
│                                                                         │
│  5. RECORD RESULT                                                       │
│     • assigned_label = final node reached                               │
│     • confidence = min_margin_along_path                                │
│     • stop_reason = reason for stopping                                 │
│                                                                         │
└─────────────────────────────────────────────────────────────────────────┘
```

### Algorithm Pseudocode

The following pseudocode illustrates the core assignment logic:

```python
def assign_label(cluster, hierarchy):
    # Step 1: Find best root
    passing_roots = [r for r in roots if passes_gate(cluster, r)]
    if not passing_roots:
        return "Unassigned", "no_root_passed"

    best_root = max(passing_roots, key=score)
    if gap_to_runner_up < root_gap_threshold:
        return best_root, "ambiguous_root"

    # Step 2: Descend through hierarchy
    current = best_root
    while has_children(current):
        children = get_children(current)
        passing = [c for c in children if passes_gate(cluster, c)]

        if not passing:
            return current, "no_child_passed"

        best_child = max(passing, key=score)
        gap = best_child.score - runner_up.score

        if gap < min_gap[level]:
            return current, "ambiguous_siblings"

        current = best_child

    return current, "leaf_reached"
```

### Key Decision Points

1. **Root Selection**: The algorithm first identifies which root category (lineage) the cluster belongs to. This is the most critical decision because it constrains all downstream choices.

2. **Gate Filtering**: At each level, only nodes that pass all gate checks are considered as candidates. This prevents weak evidence from influencing decisions.

3. **Score-Based Ranking**: Among passing candidates, nodes are ranked by their composite score (see [Marker Scoring Algorithm](./marker-scoring-algorithm.md)).

4. **Gap-Based Confidence**: The algorithm checks if the margin between the best and second-best candidate exceeds the level-specific threshold. If not, it stops and reports ambiguity.

5. **Descent Termination**: The algorithm stops when it reaches a leaf node, encounters ambiguity, or when no children pass gates.

---

## Stop Reasons Explained

The algorithm records why it stopped descending, providing diagnostic information for result interpretation and downstream refinement.

| Stop Reason | Meaning | Assigned Label | What to Do |
|-------------|---------|----------------|------------|
| `leaf_reached` | Best case - reached deepest level | Deepest node | None needed |
| `ambiguous_siblings` | Top 2 children too close | Parent node | Consider refinement in Stage I |
| `no_child_passed` | All children failed gates | Current node | Check marker coverage |
| `no_root_passed` | No root category fit | "Unassigned" | Check panel, lower thresholds |
| `ambiguous_root` | Top 2 roots too close | "Root1~Root2" | Check root markers |

### Detailed Stop Reason Analysis

**leaf_reached**

This is the ideal outcome. The algorithm descended all the way to a leaf node with confident decisions at every level. The assigned label is as specific as the hierarchy allows.

**ambiguous_siblings**

The two top-scoring children at some level had scores too close together (gap < min_gap). Rather than making an uncertain call, the algorithm assigns the parent label. This preserves accuracy at the cost of specificity.

For example, if "T Cell" and "B Cell" both score highly under "Lymphoid" with insufficient gap, the cluster is assigned "Lymphoid" with stop_reason="ambiguous_siblings".

**no_child_passed**

All children of the current node failed gate checks. This often indicates:
- Missing markers in the panel for the expected subtypes
- A true intermediate or transitional cell state
- Thresholds that are too stringent for the data

The cluster is assigned to the current (parent) node.

**no_root_passed**

No root category passed all gate checks. This results in an "Unassigned" label and typically indicates:
- Clusters of low-quality cells or doublets
- Cell types not represented in the marker hierarchy
- Panel design issues (missing key lineage markers)

**ambiguous_root**

Two or more root categories scored too closely. The label is formatted as "Root1~Root2" to indicate the ambiguity. This can happen with:
- Doublets containing cells from multiple lineages
- Rare hybrid cell types
- Panel limitations preventing clear lineage discrimination

---

## Confidence Calculation

Each assignment includes a confidence score that reflects the minimum margin encountered along the assignment path.

```text
confidence = min_margin_along_path

Where margin at each level = winner_score - runner_up_score

Special cases:
  • No competition (single candidate): margin = infinity → sentinel 1e6
  • no_root_passed: confidence = 0.0
```

### Confidence Interpretation

- **High confidence (> 1.0)**: Clear winner at every level; robust assignment
- **Moderate confidence (0.3 - 1.0)**: Generally reliable but some levels had closer competition
- **Low confidence (< 0.3)**: Assignment is tentative; consider for manual review

The confidence score is the **minimum** margin along the path, not the average. This means a single close call anywhere in the hierarchy will lower the overall confidence, which is appropriate since uncertainty at any level propagates to the final assignment.

### Margin Calculation Example

Consider a path: Root (Immune) -> Lymphoid -> T Cell

- Root level: Immune scores 2.5, Epithelium scores 1.0 -> margin = 1.5
- Level 1: Lymphoid scores 1.8, Myeloid scores 0.9 -> margin = 0.9
- Level 2: T Cell scores 1.2, B Cell scores 0.8 -> margin = 0.4

**confidence = min(1.5, 0.9, 0.4) = 0.4**

The close competition between T Cell and B Cell dominates the confidence, even though root selection was highly confident.

---

## Per-Cell Voting (Evidence Collection)

When siblings are too close (ambiguous_siblings), the algorithm performs per-cell voting to gather additional evidence. This information is stored for potential use in downstream refinement stages.

```text
╔═══════════════════════════════════════════════════════════════════════╗
║  CRITICAL: Voting is EVIDENCE ONLY                                    ║
║                                                                       ║
║  • assigned_label = PARENT (not any child)                            ║
║  • stop_reason = "ambiguous_siblings" (not "voted")                   ║
║  • Composition stored for audit and Stage I refinement hints          ║
╚═══════════════════════════════════════════════════════════════════════╝
```

### How Per-Cell Voting Works

When the cluster-level decision is ambiguous:

1. **Individual cell scoring**: Each cell in the cluster is scored against the competing siblings
2. **Vote assignment**: Each cell "votes" for its highest-scoring candidate
3. **Composition recording**: The vote distribution is stored (e.g., 60% T Cell, 40% B Cell)

### Important Constraints

Per-cell voting does **not** change the assignment. The label remains the parent node, and the stop_reason remains "ambiguous_siblings". The voting information is purely diagnostic.

This design prevents over-commitment to uncertain calls. If the cluster truly contains a mixture of cell types, this will be addressed in refinement stages where clusters can be split.

### Voting Output Format

The voting results are stored in the cluster metadata:

```python
{
    "cluster_id": 42,
    "assigned_label": "Lymphoid",
    "stop_reason": "ambiguous_siblings",
    "ambiguous_candidates": ["T Cell", "B Cell"],
    "cell_votes": {
        "T Cell": 0.62,
        "B Cell": 0.38
    },
    "recommendation": "consider_split"  # or "likely_homogeneous"
}
```

---

## Handling Edge Cases

### Missing Markers

When key markers are missing from the panel, the algorithm adapts:

1. **Coverage check fails**: If too few markers are present, the node is rejected
2. **Reduced discrimination**: Available markers may not distinguish between siblings
3. **Conservative assignment**: The algorithm stops at a higher level rather than guessing

### Single-Child Nodes

Some hierarchy nodes may have only one child. In these cases:
- If the child passes gates, it is automatically selected (no competition)
- The margin is set to the sentinel value (1e6) indicating "no competition"
- This does not reduce confidence since there was no ambiguity

### Negative Markers Only

Some cell types are defined primarily by what they do NOT express (e.g., CD4-/CD8- T cells). The algorithm handles this through:
- Anti-marker scoring contributing to the overall score
- Gate checks on positive fraction can be bypassed if sufficient anti-marker evidence exists
- Coverage calculations include negative markers in the denominator

---

## Performance Considerations

The hierarchical structure provides significant computational efficiency:

1. **Pruning**: Once a root is selected, all other branches are ignored
2. **Early stopping**: Ambiguous decisions prevent unnecessary computation on deeper levels
3. **Parallelization**: Different clusters can be processed independently

For large datasets, the algorithm scales linearly with the number of clusters, not with the size of the hierarchy.

---

## See Also

- [Marker Scoring Algorithm](./marker-scoring-algorithm.md) - How individual node scores are calculated
- [Refinement Decision Logic](./refinement-decision-logic.md) - What happens after initial assignment
- [Tuning Guide](./tuning-guide.md) - How to adjust parameters for your data
