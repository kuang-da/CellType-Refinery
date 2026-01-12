---
sidebar_position: 1
---

# Cell Type Annotation Methodology

This document describes the core methodological contribution of CellType-Refinery: a hierarchical marker-based annotation system that combines the scalability of automated clustering with the biological rigor of expert-defined marker hierarchies. Unlike traditional approaches that rely on post-hoc manual annotation of unsupervised clusters, our method embeds domain knowledge directly into the assignment algorithm, producing interpretable and reproducible cell type labels at scale.

## The Problem

### Traditional Clustering Falls Short at Scale

The standard single-cell annotation workflow follows a familiar pattern: cluster cells using an unsupervised algorithm (typically Leiden or Louvain), then manually inspect each cluster to assign biological labels. This approach has served the field well for small datasets, but breaks down as datasets grow:

- **Manual bottleneck**: A dataset with 50-100 clusters requires significant expert time to annotate, and this work must be repeated whenever clustering parameters change.
- **Inconsistent labels**: Different annotators may assign different labels to similar clusters, and the same annotator may be inconsistent across sessions.
- **No audit trail**: When a cluster is labeled "CD8+ T cells," there is often no record of which markers were considered, which were present, or how close the call was.
- **Binary decisions**: Traditional annotation forces a single label per cluster, even when the evidence is genuinely ambiguous.

### Marker Ambiguity and Overlapping Cell Types

Biological complexity compounds these challenges:

- **Shared markers**: Many markers are expressed by multiple cell types. CD45 is expressed by all immune cells; EpCAM is expressed by multiple epithelial subtypes. A cluster expressing CD45 alone cannot be assigned to a specific immune lineage.
- **Continuous phenotypes**: Cell states exist on a continuum. A "transitional" cell may express markers from two canonical types, making assignment genuinely difficult.
- **Panel limitations**: Multiplexed imaging panels (20-60 markers) cannot capture the full transcriptome. A panel may lack the discriminating markers for certain cell types.
- **Technical noise**: Spillover, autofluorescence, and batch effects can obscure true biological signal.

These realities mean that any automated annotation system must handle ambiguity gracefully, rather than pretending it does not exist.

## Key Insight: Hierarchical Marker-Based Gating

The central insight of CellType-Refinery is that cell type annotation should mirror how experts actually think about cell identity: as a hierarchical decision tree where each branch point is defined by specific marker requirements.

Consider how a pathologist would identify a neutrophil:

1. First, confirm the cell is immune (CD45+, not epithelial)
2. Then, confirm it is myeloid (not lymphoid)
3. Then, confirm it is granulocytic (not monocytic)
4. Finally, distinguish neutrophil from eosinophil based on specific markers

This hierarchical reasoning has critical implications for algorithm design:

- **Order matters**: You cannot meaningfully compare a neutrophil to an epithelial cell. They exist at different positions in the hierarchy.
- **Parents gate children**: A cluster cannot be considered for "Neutrophil" unless it first passes the gates for "Immune," "Myeloid," and "Granulocyte."
- **Siblings compete**: At each level, only nodes sharing a parent are compared. Neutrophils compete with eosinophils, not with T cells.

This structure is not just computationally convenient; it reflects biological reality. Cell types are organized hierarchically, and annotation should respect that organization.

### Beyond Simple Clustering

Traditional approaches treat clustering and annotation as separate steps:

```
Cells → Cluster → Manual Review → Labels
```

CellType-Refinery integrates domain knowledge into the assignment process:

```
Cells → Cluster → Hierarchical Gating → Traceable Labels + Confidence Scores
```

The key differences:

| Traditional Approach | CellType-Refinery |
|---------------------|-------------------|
| Labels assigned post-hoc by expert | Labels assigned algorithmically using expert-defined rules |
| No confidence scores | Every assignment has a confidence score |
| No audit trail | Every decision step is logged |
| Binary labels only | Ambiguous cases flagged for review |
| Restart from scratch if markers change | Re-run annotation without re-clustering |

## Design Principles

The annotation algorithm is built on five core principles, each addressing a specific failure mode of traditional approaches:

### 1. Parents as Gates

A cluster must pass the parent's gate before comparing with siblings. This prevents biologically meaningless comparisons and ensures that each level of the hierarchy adds specificity.

```
Wrong: Compare Epithelium vs Immune vs Neutrophils directly
Right: Must pass Immune gate before comparing Neutrophils vs Eosinophils
```

For example, a cluster cannot be considered for "Ciliated Epithelium" unless it first demonstrates strong epithelial markers (EpCAM, E-cadherin). This prevents a cluster of immune cells with weak off-target EpCAM signal from being erroneously assigned an epithelial label.

### 2. Siblings-Only Competition

At each level, only compare nodes that share a parent. This maintains biological coherence and prevents the algorithm from making comparisons that an expert would never make.

```
Level 0: Epithelium vs Immune vs Mesenchymal vs Endothelium
Level 1 (under Immune): Lymphoids vs Myeloids
Level 2 (under Myeloids): Granulocytes vs Monocytes

Never compare: Neutrophils vs Ciliated Epithelium (different branches)
```

This principle ensures that the discriminating markers are always appropriate for the comparison being made.

### 3. Traceable Decisions

Every assignment step is logged with full context. For any cluster, you can trace exactly which candidates were considered, which passed or failed gates (and why), and how the final selection was made.

This traceability serves multiple purposes:
- **Debugging**: When an assignment seems wrong, you can inspect the decision trace to understand why
- **Quality control**: Reviewers can audit the algorithm's reasoning
- **Reproducibility**: The same inputs always produce the same outputs, with full documentation

### 4. Ambiguity Detection

Close calls are flagged for review rather than forced into a single label. The algorithm distinguishes between:

- **Confident assignments**: Clear winner with large margin over alternatives
- **Ambiguous assignments**: Two or more candidates within the minimum gap threshold

When ambiguity is detected, the algorithm:
- Assigns the parent label (not a child)
- Records the stop reason as "ambiguous_siblings"
- Collects per-cell voting evidence for downstream refinement

This approach acknowledges biological reality: some clusters genuinely contain mixed populations or transitional states.

### 5. Evidence vs Override

Per-cell voting provides evidence but never overrides cluster-level assignments. When siblings are too close to call at the cluster level, the algorithm computes per-cell scores for each candidate and records the vote composition.

This evidence is invaluable for refinement decisions:
- If 80% of cells vote for one type and 20% for another, subclustering may reveal two distinct populations
- If votes are evenly split, the cluster may represent a genuine transitional state

Critically, this voting does not override the cluster-level assignment. The cluster is labeled with the parent type, and the voting composition is stored as evidence for expert review or Stage I refinement.

## Pipeline Overview

The complete annotation pipeline consists of two major stages that work together iteratively:

```text
+-------------------------------------------------------------------+
|                    CELL TYPE ANNOTATION PIPELINE                  |
+-------------------------------------------------------------------+
|                                                                   |
|  STAGE H: Initial Clustering & Annotation                         |
|  ----------------------------------------                         |
|  1. Preprocessing (layer selection, marker resolution)            |
|  2. Clustering (PCA -> k-NN -> Leiden)                            |
|  3. Differential Expression (Wilcoxon rank-sum)                   |
|  4. Marker Scoring (enrichment + positive + DE - anti)            |
|  5. Hierarchical Assignment (top-down with gating)                |
|  6. Per-Cell Voting (evidence collection)                         |
|                                                                   |
|                           |                                       |
|                           v                                       |
|                                                                   |
|  STAGE I: Refinement                                              |
|  -------------------                                              |
|  - Diagnostic mode: identify clusters needing attention           |
|  - Subclustering: split heterogeneous clusters                    |
|  - Relabeling: instant reassignment when clear signal exists      |
|  - Iterative: H -> I -> I -> I until convergence                  |
|                                                                   |
+-------------------------------------------------------------------+
```

### Stage H: Initial Clustering and Annotation

Stage H performs the heavy lifting of clustering and initial annotation:

1. **Preprocessing**: Select the appropriate expression layer, resolve marker names against the panel, exclude technical markers (DAPI, housekeeping genes)

2. **Clustering**: Standard single-cell workflow (scale, PCA, k-NN graph, Leiden community detection) with optional GPU acceleration

3. **Differential Expression**: Wilcoxon rank-sum test identifies markers enriched in each cluster vs. all others, providing orthogonal evidence for annotation

4. **Marker Scoring**: For each (cluster, cell type) pair, compute a composite score incorporating expression enrichment, positive fraction, DE bonus, and anti-marker penalty

5. **Hierarchical Assignment**: Top-down traversal through the marker hierarchy, selecting the best-scoring child at each level until a leaf is reached or ambiguity is detected

6. **Per-Cell Voting**: When siblings are too close, collect per-cell scores to record vote composition as evidence

### Stage I: Refinement

Stage I addresses the clusters that Stage H could not confidently annotate:

- **Diagnostic Mode**: Identify clusters with low confidence, ambiguous assignments, or unusual marker patterns
- **Subclustering**: Re-cluster ambiguous clusters at higher resolution to separate mixed populations
- **Relabeling**: When new evidence (from subclustering or marker map updates) provides clear signal, instantly reassign labels
- **Iteration**: The refinement loop can be repeated until annotation quality converges

The iterative nature of this workflow acknowledges that perfect annotation rarely happens in one pass. Stage H provides a strong starting point; Stage I refines the difficult cases.

## Algorithm Summary

The hierarchical assignment algorithm can be summarized in pseudocode:

```
For each cluster:
    1. Evaluate all root cell types (Epithelium, Immune, etc.)
    2. Select the best-scoring root that passes the root gate
    3. While current node has children:
        a. Evaluate all children of current node
        b. Filter to those passing the child gate
        c. If no children pass: stop, assign current node
        d. If gap between top two children < threshold: stop, flag ambiguous
        e. Otherwise: select best child, continue descent
    4. Record final assignment with confidence and stop reason
```

The confidence score is the minimum margin along the entire path from root to final assignment. This captures the "weakest link" in the decision chain.

## What Makes This Different

CellType-Refinery is not the first tool to attempt automated cell type annotation. What distinguishes this approach:

1. **Hierarchy-aware**: Most tools treat cell types as a flat list. CellType-Refinery respects the biological hierarchy.

2. **Interpretable**: Every decision is logged and can be inspected. No black-box classifiers.

3. **Uncertainty-aware**: Ambiguity is detected and flagged, not hidden.

4. **Iterative refinement**: The two-stage workflow allows progressive improvement.

5. **Configurable**: Gating thresholds, marker definitions, and anti-marker rules can all be customized per tissue and per cell type.

## See Also

- [Annotation Pipeline](./annotation-pipeline.md) - Detailed pipeline flow with implementation notes
- [Marker Scoring Algorithm](./marker-scoring-algorithm.md) - The scoring formula and its components
- [Hierarchical Gating](./hierarchical-gating-algorithm.md) - Assignment logic and gate checking
- [Refinement Decision Logic](./refinement-decision-logic.md) - When to subcluster vs relabel
- [Tuning Guide](./tuning-guide.md) - Parameter adjustment for different tissues and panels
