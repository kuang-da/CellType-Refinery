---
sidebar_position: 5
---

# Refinement Decision Logic

:::info Prerequisites
This page explains Stage I refinement criteria. For initial annotation, see [Hierarchical Gating](./hierarchical-gating-algorithm.md) first.
:::

## Overview: The Refinement Question

Stage H produces initial annotations. Some clusters are confident, others need refinement.

The key question: **SUBCLUSTER (re-cluster) or RELABEL (instant reassignment)?**

This decision determines whether a cluster undergoes computational re-clustering to discover
hidden subpopulations, or receives a simple label update based on existing evidence. Making
the wrong choice wastes computational resources or misses biological signal.

```text
┌─────────────────────────────────────────────────────────────────┐
│                        POLICY-BASED ARCHITECTURE                │
│                                                                 │
│    ┌─────────────────┐         ┌─────────────────┐              │
│    │   AutoPolicy    │         │  ManualPolicy   │              │
│    │  (--auto flag)  │         │ (--config YAML) │              │
│    └────────┬────────┘         └────────┬────────┘              │
│             │                           │                       │
│             ▼                           ▼                       │
│        base_plan                   overlay_plan                 │
│             │                           │                       │
│             └───────────┬───────────────┘                       │
│                         ▼                                       │
│              ┌─────────────────────┐                            │
│              │    merge_plans()    │                            │
│              │  (overlay wins on   │                            │
│              │   conflicts)        │                            │
│              └──────────┬──────────┘                            │
│                         ▼                                       │
│                    merged_plan                                  │
│                         │                                       │
│                         ▼                                       │
│              ┌─────────────────────┐                            │
│              │  RefinementEngine   │                            │
│              │  (execute plan on   │                            │
│              │   AnnData object)   │                            │
│              └──────────┬──────────┘                            │
│                         ▼                                       │
│                      Output                                     │
└─────────────────────────────────────────────────────────────────┘
```

The architecture separates policy definition from execution. AutoPolicy provides algorithmic
recommendations while ManualPolicy allows expert overrides. Both produce plans that merge
before execution.

## AutoPolicy Selection Criteria

AutoPolicy evaluates each cluster against a series of criteria to determine the appropriate
refinement action. These criteria are evaluated in order, with the first matching criterion
determining the action.

### Criterion 1: Low Confidence

```text
IF score < score_threshold (default: 1.0)
AND n_cells >= min_cells (default: 500)
→ SUBCLUSTER
```

**Rationale**: Large, poorly-scoring clusters likely contain mixed populations. Subclustering
may reveal hidden subpopulations.

**Implementation Details**:
- The score threshold is configurable via `--score-threshold`
- Minimum cell count prevents fragmenting small clusters into noise
- Score is computed from the hierarchical gating algorithm output
- Clusters below this threshold show weak marker expression patterns

**Example Scenario**:
A cluster of 2,000 cells labeled "Immune" has a score of 0.6. This suggests the marker
signature is ambiguous. Subclustering may reveal distinct T-cell and B-cell populations
that were merged due to insufficient resolution.

### Criterion 2a: Homogeneous Parent

```text
IF assigned at parent level (e.g., "Epithelium" not "Ciliated Epithelium")
AND best_child_score > subtype_signal_threshold (default: 0.5)
AND (only_one_passing_child OR gap_to_runner_up >= heterogeneity_gap)
→ RELABEL to best child
```

**Rationale**: Clear subtype signal with no competition. No need to re-cluster - instant
relabel is sufficient.

**Key Parameters**:
- `subtype_signal_threshold`: Minimum score for a child to be considered viable (default: 0.5)
- `heterogeneity_gap`: Required score difference between top candidates (default: 0.3)

**Example Scenario**:
A cluster labeled "Epithelium" shows a child score of 0.85 for "Ciliated Epithelium" and
0.3 for "Secretory Epithelium". The gap of 0.55 exceeds the threshold of 0.3, so the
cluster is relabeled to "Ciliated Epithelium" without re-clustering.

### Criterion 2b: Heterogeneous Parent

```text
IF assigned at parent level
AND best_child_score > subtype_signal_threshold
AND multiple_passing_children
AND gap < heterogeneity_gap (default: 0.3)
→ SUBCLUSTER
```

**Rationale**: Competing child signals suggest a mixed population. Subclustering should
separate the different cell types.

**Detection Logic**:
1. Count children with scores above `subtype_signal_threshold`
2. If more than one child passes AND
3. The difference between top two is less than `heterogeneity_gap`
4. The cluster contains mixed populations requiring separation

**Example Scenario**:
A cluster labeled "T Cell" has child scores of 0.7 for "CD4+ T Cell" and 0.65 for
"CD8+ T Cell". Both exceed the threshold, and the gap of 0.05 is below 0.3. This
indicates a mix of CD4+ and CD8+ T cells that subclustering can separate.

### Criterion 3: Mixed Population (Extension)

```text
IF parent has high score
AND all children below subtype_signal_threshold
→ SUBCLUSTER
```

**Rationale**: Strong parent signal but weak children = novel subtype or noise that needs
investigation.

**Biological Interpretation**:
- The cells clearly belong to the parent category (strong marker expression)
- No known child subtype matches the expression pattern
- Possibilities include:
  - Novel subtype not in the marker map
  - Transitional state between subtypes
  - Technical noise requiring further investigation

**Example Scenario**:
A cluster labeled "Macrophage" has a parent score of 0.9 but all child scores (M1, M2,
tissue-resident) are below 0.4. This may represent a novel macrophage subtype or
activation state not captured in the current marker hierarchy.

### Criterion 4: Weak Leaf (Extension)

```text
IF at leaf level
AND score low
AND marker heterogeneity high
→ SUBCLUSTER
```

**Rationale**: Leaf with high variance may contain distinct subpopulations.

**Heterogeneity Metrics**:
- Coefficient of variation across marker genes
- Bimodal distribution detection in key markers
- Entropy of marker expression within cluster

**Example Scenario**:
A cluster labeled "Regulatory T Cell" (a leaf node) has a score of 0.5 and shows
bimodal expression of FOXP3. Subclustering may reveal true Tregs versus activated
conventional T cells.

## Decision Flowchart

The following flowchart summarizes the complete decision logic. Each cluster traverses
this tree to determine its refinement action.

```text
┌─────────────────────────────────────────────────────────────────┐
│                     REFINEMENT DECISION TREE                    │
└─────────────────────────────────────────────────────────────────┘

                         Cluster Assessment
                                 │
                 ┌───────────────┴───────────────┐
                 ▼                               ▼
         score < threshold?              At parent level?
                 │                               │
           ┌─────┴─────┐                  ┌──────┴──────┐
           YES         NO                 YES           NO
           │           │                   │            │
           ▼           ▼                   ▼            ▼
      n_cells ≥ 500?  SKIP         Has child signal?   SKIP
           │                               │
     ┌─────┴─────┐                  ┌──────┴──────┐
     YES         NO                 YES           NO
     │           │                   │            │
     ▼           ▼                   ▼            ▼
  SUBCLUSTER   SKIP           Gap ≥ threshold?   SKIP
                                     │
                              ┌──────┴──────┐
                              YES           NO
                              │             │
                              ▼             ▼
                           RELABEL     SUBCLUSTER
```

**Flowchart Interpretation**:
- Left branch: Low-confidence clusters evaluated for subclustering based on size
- Right branch: Parent-level clusters evaluated for relabeling vs subclustering
- SKIP: Cluster is already optimally annotated, no action needed

## Operation Types

Operations execute in deterministic order: `override → merge → subcluster → relabel → rescore`

This ordering ensures that:
1. Manual overrides take precedence over all automated decisions
2. Cluster merges happen before subclustering to avoid redundant computation
3. Subclustering precedes relabeling to generate new clusters first
4. Rescoring happens last to reflect all structural changes

| Operation | Description | Source |
|-----------|-------------|--------|
| **OVERRIDE** | Direct label assignment | ManualPolicy (YAML) |
| **MERGE** | Combine multiple clusters | ManualPolicy |
| **SUBCLUSTER** | Re-cluster with finer resolution | AutoPolicy or ManualPolicy |
| **RELABEL** | Instant label change (no re-clustering) | AutoPolicy |
| **RESCORE** | Recompute scores with updated markers | AutoPolicy |

### Operation Details

**OVERRIDE**: Assigns a specific label regardless of algorithmic recommendations. Useful
when domain expertise contradicts automated classification.

**MERGE**: Combines two or more clusters into a single cluster. Typically used when
over-clustering has split a homogeneous population.

**SUBCLUSTER**: Applies Leiden clustering at higher resolution to subdivide a cluster.
Parameters include resolution factor and minimum cluster size.

**RELABEL**: Changes the cluster label without re-clustering. Fast operation that
updates metadata only.

**RESCORE**: Recomputes hierarchical gating scores using potentially updated marker
definitions. Triggered when marker maps change between iterations.

## Diagnostic vs Execution Mode

Stage I operates in two distinct modes to support both exploration and production workflows.

```text
DIAGNOSTIC MODE (default, no --execute)
├─ Generates diagnostic_report.csv
├─ Shows recommendations (SUBCLUSTER, RELABEL, or SKIP)
└─ Does NOT modify data

EXECUTION MODE (with --execute)
├─ Builds execution plan from recommendations
├─ Executes all operations
├─ Exports refined.h5ad
└─ Stores provenance in AnnData
```

### Diagnostic Mode

The default mode generates a comprehensive report without modifying data. This allows
users to review recommendations before committing changes.

**Output Files**:
- `diagnostic_report.csv`: Per-cluster recommendations with scores
- `score_distribution.png`: Visualization of score distributions
- `decision_summary.json`: Machine-readable decision rationale

**Report Columns**:
- `cluster_id`: Original cluster identifier
- `current_label`: Label from Stage H
- `recommendation`: SUBCLUSTER, RELABEL, or SKIP
- `target_label`: Proposed new label (for RELABEL)
- `confidence`: Algorithmic confidence in recommendation
- `rationale`: Human-readable explanation

### Execution Mode

Execution mode applies the refinement plan to produce updated annotations.

**Provenance Tracking**:
All operations are logged in `adata.uns['refinement_provenance']`:
- Operation type and parameters
- Source cluster(s) and target label(s)
- Timestamp and iteration number
- Score changes before/after

## Iterative Refinement

Stage I supports iterative chains until convergence:

```text
H → I (iteration 1) → I (iteration 2) → I (iteration 3) → ...
```

Each iteration:
1. Uses previous output's cluster_lvl1 as new cluster_lvl0
2. Selects new candidates based on updated scores
3. Further refines clusters that still need improvement
4. Tracks full lineage in provenance

Typically converges in 2-3 iterations.

### Convergence Criteria

Iteration stops when any of the following conditions are met:
- No clusters meet refinement criteria (all SKIP)
- Maximum iteration count reached (default: 5)
- Score improvement falls below threshold (default: 0.05)
- Manual termination via configuration

### Lineage Tracking

Each cell maintains a complete lineage record:
```text
Cell X: Immune → T Cell → CD8+ T Cell → Effector CD8+ T Cell
        (H)       (I.1)     (I.2)         (I.3)
```

This enables:
- Tracing annotation provenance
- Comparing intermediate states
- Debugging unexpected classifications

## Configuration Reference

### AutoPolicy Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `score_threshold` | 1.0 | Minimum score to skip refinement |
| `min_cells` | 500 | Minimum cluster size for subclustering |
| `subtype_signal_threshold` | 0.5 | Minimum child score for relabeling |
| `heterogeneity_gap` | 0.3 | Required gap between top candidates |
| `max_iterations` | 5 | Maximum refinement iterations |

### ManualPolicy YAML Format

```yaml
refinements:
  - cluster: "Cluster_15"
    action: override
    label: "Plasma Cell"

  - cluster: "Cluster_8"
    action: subcluster
    resolution: 1.2

  - clusters: ["Cluster_3", "Cluster_7"]
    action: merge
    label: "Fibroblast"
```

## See Also

- [Hierarchical Gating Algorithm](./hierarchical-gating-algorithm.md) - Initial annotation
- [Tuning Guide](./tuning-guide.md) - Refinement parameters
