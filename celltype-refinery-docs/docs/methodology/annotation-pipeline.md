---
sidebar_position: 2
---

# Annotation Pipeline

This page provides a complete overview of the cell type annotation pipeline, showing how data flows through Stage H (initial clustering and annotation) and Stage I (refinement).

## Complete Pipeline Flow

The annotation pipeline transforms raw merged expression data into refined cell type assignments through a multi-stage process. The following diagram shows the complete data flow from inputs through all processing stages to final outputs.

```text
                               INPUTS
                                 │
         ┌───────────────────────┼───────────────────────┐
         ▼                       ▼                       ▼
 ┌───────────────┐    ┌─────────────────────┐    ┌─────────────────┐
 │ merged_data   │    │ cell_type_markers   │    │   CLI Args      │
 │   .h5ad       │    │   .json             │    │ (resolution,    │
 │ (Stage F)     │    │ (marker hierarchy)  │    │  n_pcs, etc.)   │
 └───────┬───────┘    └──────────┬──────────┘    └────────┬────────┘
         │                       │                        │
         └───────────────────────┼────────────────────────┘
                                 │
                                 ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  1. PREPROCESSING                                                           │
│  ┌────────────────────────────────────────────────────────────────────────┐ │
│  │ - Select expression layer (batchcorr / aligned / raw)                  │ │
│  │ - Filter low-variance markers (std < 1e-3)                             │ │
│  │ - Exclude technical markers (DAPI, Collagen IV, Beta-actin) -> .obs    │ │
│  │ - Load and resolve marker sets from JSON hierarchy                     │ │
│  └────────────────────────────────────────────────────────────────────────┘ │
└─────────────────────────────────────────────────────────────────────────────┘
                                 │
                                 ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  2. CLUSTERING (GPU-accelerated if available)                               │
│  ┌────────────────────────────────────────────────────────────────────────┐ │
│  │ - Scale features (zero-center, clip at +/-10 std)                      │ │
│  │ - PCA dimensionality reduction (default: 30 components)                │ │
│  │ - k-NN graph construction (default: k=15)                              │ │
│  │ - Leiden community detection (default: resolution=0.6)                 │ │
│  │ - Optional: UMAP embedding for visualization                           │ │
│  └────────────────────────────────────────────────────────────────────────┘ │
└─────────────────────────────────────────────────────────────────────────────┘
                                 │
                                 ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  3. DIFFERENTIAL EXPRESSION                                                 │
│  ┌────────────────────────────────────────────────────────────────────────┐ │
│  │ - Wilcoxon rank-sum test with tie correction (cluster vs rest)         │ │
│  │ - Dynamic K: ceil(panel_size x de_top_frac) clamped to [de_min_k, 12]  │ │
│  │ - Rank-based lookup for continuous DE bonus scoring                    │ │
│  │ - Commonness penalty for markers shared across many cell types         │ │
│  └────────────────────────────────────────────────────────────────────────┘ │
└─────────────────────────────────────────────────────────────────────────────┘
                                 │
                                 ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  4. MARKER SCORING (per cluster x cell type)                                │
│  ┌────────────────────────────────────────────────────────────────────────┐ │
│  │                                                                        │ │
│  │  score = mean_enrichment + mean_positive + de_component - anti_penalty │ │
│  │          -----------------------------------------------------         │ │
│  │               |              |           |            |                │ │
│  │               v              v           v            v                │ │
│  │         z-score vs     % cells      bonus if    penalty for           │ │
│  │         global median  above Q75    marker is   anti-markers           │ │
│  │                                     DE hit      being expressed        │ │
│  │                                                                        │ │
│  │  Additional metrics: coverage, frac_markers_on, n_markers_on           │ │
│  └────────────────────────────────────────────────────────────────────────┘ │
└─────────────────────────────────────────────────────────────────────────────┘
                                 │
                                 ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  5. HIERARCHICAL ASSIGNMENT (top-down traversal)                            │
│  ┌────────────────────────────────────────────────────────────────────────┐ │
│  │                                                                        │ │
│  │              Marker Hierarchy Example:                                 │ │
│  │                                                                        │ │
│  │                    ┌─────────┐                                         │ │
│  │         ┌──────────┤  ROOT   ├──────────┐                              │ │
│  │         ▼          └─────────┘          ▼                              │ │
│  │    ┌──────────┐                   ┌───────────┐                        │ │
│  │    │Epithelium│                   │ Immune    │  ... (other roots)     │ │
│  │    └────┬─────┘                   │  Cells    │                        │ │
│  │         │                         └─────┬─────┘                        │ │
│  │    ┌────┴────┐                    ┌─────┴─────┐                        │ │
│  │    ▼         ▼                    ▼           ▼                        │ │
│  │ Ciliated  Glandular           Myeloids    Lymphoids                    │ │
│  │                                   │                                    │ │
│  │                              ┌────┴────┐                               │ │
│  │                              ▼         ▼                               │ │
│  │                         Granulocytes  Monocytes                        │ │
│  │                                                                        │ │
│  │  Algorithm:                                                            │ │
│  │  1. Find best passing ROOT (score > threshold, coverage > min)         │ │
│  │  2. Descend: compare SIBLINGS only at each level                       │ │
│  │  3. Stop when: leaf reached, no child passes gate, or ambiguous        │ │
│  │                                                                        │ │
│  └────────────────────────────────────────────────────────────────────────┘ │
└─────────────────────────────────────────────────────────────────────────────┘
                                 │
                                 ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  6. PER-CELL VOTING (Evidence Collection)                                   │
│  ┌────────────────────────────────────────────────────────────────────────┐ │
│  │ - Triggered when sibling gap < min_gap threshold                       │ │
│  │ - Computes per-cell Z-score for each candidate's markers               │ │
│  │ - Records vote composition as EVIDENCE ONLY (does NOT override labels) │ │
│  │ - Stored for downstream audit and Stage I refinement                   │ │
│  └────────────────────────────────────────────────────────────────────────┘ │
└─────────────────────────────────────────────────────────────────────────────┘
                                 │
                                 ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  7. QC & DEBUG OUTPUTS                                                      │
│  ┌────────────────────────────────────────────────────────────────────────┐ │
│  │ - Cluster size distribution plot                                       │ │
│  │ - Confidence distribution plot                                         │ │
│  │ - Depth distribution plot                                              │ │
│  │ - Stop reason summary plot                                             │ │
│  │ - Near-miss report (close-call clusters)                               │ │
│  │ - Marker evidence table (--expand-markers)                             │ │
│  │ - Audit cards HTML (--generate-audit-cards)                            │ │
│  └────────────────────────────────────────────────────────────────────────┘ │
└─────────────────────────────────────────────────────────────────────────────┘
                                 │
                                 ▼
                              OUTPUTS
                                 │
      ┌──────────────┬───────────┼───────────┬──────────────┐
      ▼              ▼           ▼           ▼              ▼
 ┌─────────┐  ┌───────────┐ ┌─────────┐ ┌─────────┐  ┌───────────┐
 │coarse_  │  │cluster_   │ │marker_  │ │decision_│  │figures/   │
 │clusters │  │annotations│ │scores   │ │steps    │  │ *.png     │
 │.h5ad    │  │.csv       │ │.csv     │ │.csv     │  │           │
 └─────────┘  └───────────┘ └─────────┘ └─────────┘  └───────────┘
      │              │           │           │              │
      ▼              ▼           ▼           ▼              ▼
  AnnData with    Summary     Full       Hierarchical   QC plots
  cluster_lvl0    per         scoring    decision       (confidence,
  assignments     cluster     matrix     trace          depth, etc.)
```

---

## Stage Details

### 1. Preprocessing

The preprocessing stage prepares expression data for clustering and annotation.

#### Layer Selection Priority

Stage H operates on one expression layer from the merged AnnData, selected in order of preference:

```text
Priority: batchcorr > aligned > X (raw)
```

The selected layer is copied to `adata.layers["stage_h_input"]` and pinned to `adata.X` for downstream processing.

#### Technical Marker Exclusion

Technical markers are excluded from clustering features but preserved for QC:

```python
TECHNICAL_MARKERS = ["DAPI", "Collagen IV", "Beta-actin"]
```

These are moved to `.obs` columns named `{marker}_intensity` and are not used in clustering or scoring calculations.

#### Low-Variance Filtering

Markers with standard deviation < 1e-3 are dropped to avoid numerical issues in scaling and PCA. This prevents constant or near-constant markers from causing division-by-zero errors or dominating principal components.

#### Marker Resolution from JSON

Marker names from the JSON hierarchy are resolved against the AnnData panel:

| Marker Map Entry | AnnData var_names | Result |
|------------------|-------------------|--------|
| "CD45" | "CD45" | Resolved |
| "Pan-Cytokeratin" | "Pan-CK" | Alias matched |
| "CD3epsilon" | (not found) | Missing - excluded |

Missing markers are logged and excluded from scoring for that cell type.

---

### 2. Clustering

The clustering stage groups cells based on expression similarity using the following pipeline.

#### Clustering Pipeline

```text
Expression Matrix (n_cells x n_markers)
         │
         ▼
┌─────────────────────────────────────────────────────────────────┐
│ 1. SCALE                                                        │
│    - Zero-center each marker                                    │
│    - Clip values at +/-10 standard deviations                   │
│    - Prevents outlier domination in PCA                         │
└─────────────────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────────────────┐
│ 2. PCA                                                          │
│    - n_components = 30 (default)                                │
│    - SVD solver: arpack                                         │
│    - Stores in adata.obsm["X_pca"]                              │
└─────────────────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────────────────┐
│ 3. NEIGHBORS                                                    │
│    - k = 15 (default)                                           │
│    - Uses PCA coordinates                                       │
│    - Stores connectivity in adata.obsp                          │
└─────────────────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────────────────┐
│ 4. LEIDEN                                                       │
│    - resolution = 0.6 (default)                                 │
│    - GPU: rapids_singlecell (RAPIDS/cuGraph)                    │
│    - CPU: igraph backend, n_iterations=2                        │
│    - random_state = 1337                                        │
│    - Stores in adata.obs["cluster_lvl0"]                        │
└─────────────────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────────────────┐
│ 5. UMAP (optional, --compute-umap)                              │
│    - 2D embedding for visualization                             │
│    - Stores in adata.obsm["X_umap"]                             │
└─────────────────────────────────────────────────────────────────┘
```

#### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--n-pcs` | 30 | Number of PCA components |
| `--neighbors-k` | 15 | k-NN neighbors for graph construction |
| `--resolution` | 0.6 | Leiden resolution (higher = more clusters) |

#### GPU Acceleration

When `--use-gpu` is enabled (default), Stage H uses RAPIDS libraries for accelerated computation. The system falls back to CPU (scanpy) automatically on GPU errors.

**Note**: GPU Leiden is non-deterministic and may produce different cluster counts across runs even with fixed random_state. For reproducible results, use `--annotation-only` to reuse existing clusters.

---

### 3. Differential Expression

Differential expression analysis identifies marker genes that distinguish each cluster from the rest of the population.

#### Wilcoxon Rank-Sum Test

For each cluster, DE testing identifies markers significantly enriched compared to all other cells:

```python
sc.tl.rank_genes_groups(adata, groupby="cluster_lvl0", method="wilcoxon")
```

Results are stored in `adata.uns[f"de_wilcoxon_{layer}"]`.

#### Dynamic K Calculation

The number of top DE genes considered (K) adapts to panel size:

```text
K_target = ceil(panel_size x de_top_frac)
K = clamp(K_target, de_min_k, de_max_k)

Default parameters:
  de_top_frac = 0.2
  de_min_k = 3
  de_max_k = 12

Example (panel_size = 40):
  K_target = ceil(40 x 0.2) = 8
  K = clamp(8, 3, 12) = 8
```

#### Rank-Based DE Lookup

DE contribution uses rank-weighted scoring rather than binary hit/miss. Higher-ranked DE genes contribute more to the score, with a commonness penalty for markers that appear in many cell type definitions.

---

### 4. Marker Scoring

For detailed information about the scoring algorithm, see [Marker Scoring Algorithm](./marker-scoring-algorithm.md).

#### Summary

Each (cluster, cell_type) pair receives a composite score:

```text
score = mean_enrichment + mean_positive + de_component - anti_penalty
```

| Component | What It Measures |
|-----------|------------------|
| `mean_enrichment` | Z-score of cluster expression vs global median |
| `mean_positive` | Fraction of cells expressing markers above Q75 |
| `de_component` | Rank-weighted bonus for markers in top-K DE genes |
| `anti_penalty` | Penalty for anti-markers being expressed |

#### Score Interpretation

- **Typical winning scores**: 1.5 to 4.0 for well-matched cell types
- **Marginal matches**: 0.5 to 1.5, may indicate subtypes or transitional states
- **Poor matches**: < 0.5, usually indicates wrong cell type assignment
- **Negative scores**: Strong evidence against the cell type

---

### 5. Hierarchical Assignment

For detailed information about the gating algorithm, see [Hierarchical Gating Algorithm](./hierarchical-gating-algorithm.md).

#### Summary

The algorithm traverses the marker hierarchy top-down:

1. **Identify roots**: Extract labels that are never children (Epithelium, Immune, etc.)
2. **Evaluate roots**: Check hard requirements, veto markers, and standard gates
3. **Select best root**: Sort passing roots by score, check margin
4. **Descend hierarchy**: At each level, filter by gates and select best scoring child
5. **Record result**: Store assigned label, confidence, and stop reason

#### Stop Reasons

| Stop Reason | Meaning | Assigned Label |
|-------------|---------|----------------|
| `leaf_reached` | Best case - reached deepest level | Deepest node |
| `ambiguous_siblings` | Top 2 children too close | Parent node |
| `no_child_passed` | All children failed gates | Current node |
| `no_root_passed` | No root category fit | "Unassigned" |
| `ambiguous_root` | Top 2 roots too close | "Root1~Root2" |

#### Gating Parameters by Level

| Level | min_coverage | min_pos_frac | min_enrichment | min_gap |
|-------|--------------|--------------|----------------|---------|
| 0 (Root) | 0.50 | 0.30 | 0.00 | 0.50 |
| 1 | 0.40 | 0.20 | -0.50 | 0.30 |
| 2 | 0.30 | 0.15 | -1.00 | 0.20 |
| 3+ | 0.30 | 0.15 | -1.00 | 0.20 |

---

### 6. Per-Cell Voting

Per-cell voting provides additional evidence when cluster-level decisions are ambiguous.

#### Key Principle

```text
+=====================================================================+
|  CRITICAL: Voting is EVIDENCE ONLY                                  |
|                                                                     |
|  - assigned_label = PARENT (not any child)                          |
|  - stop_reason = "ambiguous_siblings" (not "voted")                 |
|  - Composition stored for audit and Stage I refinement hints        |
+=====================================================================+
```

Per-cell voting does NOT override cluster-level assignments. When siblings are too close, each cell is scored against the competing candidates, and the vote distribution is stored for downstream analysis.

#### Output Format

```python
{
    "cluster_id": 42,
    "assigned_label": "Lymphoid",
    "stop_reason": "ambiguous_siblings",
    "cell_votes": {
        "T Cell": 0.62,
        "B Cell": 0.38
    }
}
```

---

### 7. Outputs

Stage H produces multiple output files for analysis and downstream processing.

#### Primary Outputs

| File | Description |
|------|-------------|
| `coarse_clusters.h5ad` | AnnData with cluster assignments and annotation columns |
| `cluster_annotations.csv` | Per-cluster summary with assigned labels |
| `marker_scores.csv` | Full scoring matrix (clusters x cell types) |
| `decision_steps.csv` | Step-by-step hierarchical decision trace |

#### AnnData Contents

```text
adata.obs columns:
  cluster_lvl0          - Leiden cluster ID (string)
  cell_type_auto        - Assigned cell type label
  cell_type_auto_root   - Root category (Epithelium, Immune, etc.)
  cell_type_auto_conf   - Confidence score (min margin along path)

adata.uns:
  de_wilcoxon_{layer}   - Differential expression results
  stage_h_params        - Run parameters
  stage_h_gating_params - Gating thresholds used
```

---

## Stage I Refinement Flow

Stage I refines initial annotations from Stage H, addressing ambiguous assignments and discovering hidden subpopulations.

```text
                    Input: coarse_clusters.h5ad + refinement plan
                                       │
              ┌────────────────────────┴────────────────────────┐
              ▼                                                 ▼
┌──────────────────────────────────────┐    ┌──────────────────────────────────────┐
│  DIAGNOSTIC MODE (default)           │    │  EXECUTION MODE (--execute)          │
│  No modification to data             │    │  Applies refinement plan             │
├──────────────────────────────────────┤    ├──────────────────────────────────────┤
│                                      │    │                                      │
│  - Evaluate cluster quality          │    │  - Apply operations in order:        │
│  - Detect ambiguous assignments      │    │      override -> merge ->            │
│  - Analyze subtype signals           │    │      subcluster -> relabel           │
│  - Suggest operations                │    │  - Re-score affected clusters        │
│                                      │    │  - Track provenance                  │
│  Output:                             │    │                                      │
│    diagnostic_report.csv             │    │  Output:                             │
│    score_distribution.png            │    │    refined.h5ad                      │
│    decision_summary.json             │    │    refinement_log.csv                │
│                                      │    │                                      │
└──────────────────────────────────────┘    └──────────────────────────────────────┘
```

### Refinement Decision Logic

For detailed information, see [Refinement Decision Logic](./refinement-decision-logic.md).

#### Decision Criteria

| Criterion | Condition | Action |
|-----------|-----------|--------|
| Low Confidence | score < 1.0 AND n_cells >= 500 | SUBCLUSTER |
| Homogeneous Parent | clear child signal, large gap | RELABEL |
| Heterogeneous Parent | multiple competing children | SUBCLUSTER |
| Mixed Population | strong parent, weak children | SUBCLUSTER |

#### Operation Types

| Operation | Description | Source |
|-----------|-------------|--------|
| **OVERRIDE** | Direct label assignment | ManualPolicy (YAML) |
| **MERGE** | Combine multiple clusters | ManualPolicy |
| **SUBCLUSTER** | Re-cluster with finer resolution | AutoPolicy or ManualPolicy |
| **RELABEL** | Instant label change (no re-clustering) | AutoPolicy |

### Iterative Refinement

Stage I supports iterative chains until convergence:

```text
Stage H -> Stage I (iter 1) -> Stage I (iter 2) -> Stage I (iter 3) -> ...
```

Each iteration uses the previous output's clusters as input. Typically converges in 2-3 iterations.

#### Convergence Criteria

- No clusters meet refinement criteria (all SKIP)
- Maximum iteration count reached (default: 5)
- Score improvement falls below threshold (default: 0.05)

---

## End-to-End Example

```text
merged_data.h5ad (100,000 cells, 40 markers)
         │
         ▼
┌─────────────────────────────────────────────────────────────────┐
│  STAGE H                                                        │
│  - Preprocessing: 37 markers retained                           │
│  - Clustering: 42 clusters at resolution 0.6                    │
│  - Annotation: 38 assigned, 4 ambiguous                         │
└─────────────────────────────────────────────────────────────────┘
         │
         ▼
coarse_clusters.h5ad
         │
         ▼
┌─────────────────────────────────────────────────────────────────┐
│  STAGE I - Iteration 1                                          │
│  - 4 clusters flagged for refinement                            │
│  - 2 subclustered (mixed populations)                           │
│  - 2 relabeled (clear subtype signal)                           │
└─────────────────────────────────────────────────────────────────┘
         │
         ▼
refined_iter1.h5ad (48 clusters)
         │
         ▼
┌─────────────────────────────────────────────────────────────────┐
│  STAGE I - Iteration 2                                          │
│  - 1 cluster flagged for refinement                             │
│  - 1 subclustered                                               │
└─────────────────────────────────────────────────────────────────┘
         │
         ▼
refined_iter2.h5ad (50 clusters) - CONVERGED
```

---

## See Also

- [Methodology Overview](./index.md) - Conceptual background on the annotation approach
- [Marker Scoring Algorithm](./marker-scoring-algorithm.md) - Detailed scoring formula explanation
- [Hierarchical Gating Algorithm](./hierarchical-gating-algorithm.md) - Gate checks and assignment logic
- [Refinement Decision Logic](./refinement-decision-logic.md) - Stage I decision criteria
