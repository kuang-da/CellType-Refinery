# Stage H: Coarse Clustering and Cell-Type Annotation

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                    STAGE H: COARSE CLUSTERING & ANNOTATION                  │
│                                                                             │
│  Clusters cells from merged spatial data and assigns cell-type labels      │
│  using hierarchical marker-based scoring with differential expression.     │
└─────────────────────────────────────────────────────────────────────────────┘
```

## Table of Contents

1. [Overview](#overview)
2. [Pipeline Flow](#pipeline-flow)
3. [Input / Output](#input--output)
4. [Preprocessing](#1-preprocessing)
5. [Clustering](#2-clustering)
6. [Differential Expression](#3-differential-expression)
7. [Marker Scoring](#4-marker-scoring)
8. [Hierarchical Assignment](#5-hierarchical-assignment)
9. [Per-Cell Voting](#6-per-cell-voting)
10. [Gating Parameters](#7-gating-parameters)
11. [CLI Reference](#cli-reference)
12. [Output Files](#output-files)
13. [Execution Modes](#execution-modes)
14. [Configuration](#configuration)
15. [Design Principles](#design-principles)

---

## Overview

Stage H performs two major tasks:

1. **Coarse Clustering**: Leiden community detection on PCA-reduced expression data
2. **Cell-Type Annotation**: Hierarchical marker-based scoring with gating

The annotation algorithm traverses a marker hierarchy top-down, comparing only
sibling cell types at each level, using both expression enrichment and
differential expression evidence.

---

## Pipeline Flow

```
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
│  │ • Select expression layer (batchcorr / aligned / raw)                  │ │
│  │ • Filter low-variance markers (std < 1e-3)                             │ │
│  │ • Exclude technical markers (DAPI, Collagen IV, Beta-actin) → .obs     │ │
│  │ • Load and resolve marker sets from JSON hierarchy                     │ │
│  └────────────────────────────────────────────────────────────────────────┘ │
└─────────────────────────────────────────────────────────────────────────────┘
                                 │
                                 ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  2. CLUSTERING (GPU-accelerated if available)                               │
│  ┌────────────────────────────────────────────────────────────────────────┐ │
│  │ • Scale features (zero-center, clip at ±10 std)                        │ │
│  │ • PCA dimensionality reduction (default: 30 components)                │ │
│  │ • k-NN graph construction (default: k=15)                              │ │
│  │ • Leiden community detection (default: resolution=0.6)                 │ │
│  │ • Optional: UMAP embedding for visualization                           │ │
│  └────────────────────────────────────────────────────────────────────────┘ │
└─────────────────────────────────────────────────────────────────────────────┘
                                 │
                                 ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  3. DIFFERENTIAL EXPRESSION                                                 │
│  ┌────────────────────────────────────────────────────────────────────────┐ │
│  │ • Wilcoxon rank-sum test with tie correction (cluster vs rest)         │ │
│  │ • Dynamic K: ceil(panel_size × de_top_frac) clamped to [de_min_k, 12]  │ │
│  │ • Rank-based lookup for continuous DE bonus scoring                    │ │
│  │ • Commonness penalty for markers shared across many cell types         │ │
│  └────────────────────────────────────────────────────────────────────────┘ │
└─────────────────────────────────────────────────────────────────────────────┘
                                 │
                                 ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  4. MARKER SCORING (per cluster × cell type)                                │
│  ┌────────────────────────────────────────────────────────────────────────┐ │
│  │                                                                        │ │
│  │  score = mean_enrichment + mean_positive + de_component - anti_penalty │ │
│  │          ─────────────────────────────────────────────────────────     │ │
│  │               │              │           │            │                │ │
│  │               ▼              ▼           ▼            ▼                │ │
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
│  │                              │                                         │ │
│  │                         ┌────┴────┐                                    │ │
│  │                         ▼         ▼                                    │ │
│  │                    Neutrophils  Eosinophils                            │ │
│  │                                                                        │ │
│  │  Algorithm:                                                            │ │
│  │  1. Find best passing ROOT (score > threshold, coverage > min)         │ │
│  │  2. Descend: compare SIBLINGS only at each level                       │ │
│  │  3. Stop when: leaf reached, no child passes gate, or ambiguous        │ │
│  │                                                                        │ │
│  │  Stop Reasons:                                                         │ │
│  │  • leaf_reached          → best case, deepest annotation               │ │
│  │  • ambiguous_siblings    → close scores, needs review                  │ │
│  │  • no_child_passed       → stopped mid-hierarchy                       │ │
│  │  • no_root_passed        → completely unassigned                       │ │
│  │  • ambiguous_root        → multiple roots too close to distinguish     │ │
│  └────────────────────────────────────────────────────────────────────────┘ │
└─────────────────────────────────────────────────────────────────────────────┘
                                 │
                                 ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  6. PER-CELL VOTING (Evidence Collection)                                   │
│  ┌────────────────────────────────────────────────────────────────────────┐ │
│  │ • Triggered when sibling gap < min_gap threshold                       │ │
│  │ • Computes per-cell Z-score for each candidate's markers               │ │
│  │ • Records vote composition as EVIDENCE ONLY (does NOT override labels) │ │
│  │ • Stored for downstream audit and Stage I refinement                   │ │
│  └────────────────────────────────────────────────────────────────────────┘ │
└─────────────────────────────────────────────────────────────────────────────┘
                                 │
                                 ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  7. QC & DEBUG OUTPUTS                                                      │
│  ┌────────────────────────────────────────────────────────────────────────┐ │
│  │ • Cluster size distribution plot                                       │ │
│  │ • Confidence distribution plot                                         │ │
│  │ • Depth distribution plot                                              │ │
│  │ • Stop reason summary plot                                             │ │
│  │ • Near-miss report (close-call clusters)                               │ │
│  │ • Marker evidence table (--expand-markers)                             │ │
│  │ • Audit cards HTML (--generate-audit-cards)                            │ │
│  │ • Review checklist (--generate-review-checklist)                       │ │
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

## Input / Output

### Inputs

| Input | Description | Required |
|-------|-------------|----------|
| `merged_data.h5ad` | Stage F merged AnnData with expression layers | Yes |
| `cell_type_markers.json` | Hierarchical marker map with cell types and markers | Yes |
| CLI arguments | Resolution, n_pcs, layer selection, etc. | Optional |

### Outputs

| Output | Description |
|--------|-------------|
| `coarse_clusters.h5ad` | AnnData with `cluster_lvl0` and annotation columns |
| `cluster_annotations.csv` | Per-cluster summary with assigned labels |
| `marker_scores.csv` | Full scoring matrix (clusters × cell types) |
| `decision_steps.csv` | Step-by-step hierarchical decision trace |
| `marker_evidence.csv` | Per-marker contribution details (optional) |
| `figures/` | QC visualizations |

---

## 1. Preprocessing

### Layer Selection

Stage H operates on one expression layer from the merged AnnData:

```
Priority: batchcorr > aligned > X (raw)
```

The selected layer is copied to `adata.layers["stage_h_input"]` and pinned to
`adata.X` for downstream processing.

### Technical Marker Exclusion

Technical markers are excluded from clustering features but preserved for QC:

```python
TECHNICAL_MARKERS = ["DAPI", "Collagen IV", "Beta-actin"]
```

These are moved to `.obs` columns named `{marker}_intensity`.

### Low-Variance Filtering

Markers with standard deviation < 1e-3 are dropped to avoid numerical issues
in scaling and PCA.

### Marker Resolution

Marker names from the JSON hierarchy are resolved against the AnnData panel:

```
┌─────────────────────────────────────────────────────────────────────────┐
│  Marker Map JSON                    AnnData var_names                   │
│  ─────────────────                  ─────────────────                   │
│  "CD45" ──────────────────────────► "CD45"        ✓ Resolved            │
│  "Pan-Cytokeratin" ───────────────► "Pan-CK"      ✓ Alias matched       │
│  "CD3epsilon" ────────────────────► (not found)   ✗ Missing             │
└─────────────────────────────────────────────────────────────────────────┘
```

Missing markers are logged and excluded from scoring for that cell type.

---

## 2. Clustering

### Algorithm

```
┌─────────────────────────────────────────────────────────────────────────┐
│                         CLUSTERING PIPELINE                             │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                         │
│  Expression Matrix (n_cells × n_markers)                                │
│         │                                                               │
│         ▼                                                               │
│  ┌─────────────────────────────────────────────────────────────────┐   │
│  │ 1. SCALE                                                        │   │
│  │    • Zero-center each marker                                    │   │
│  │    • Clip values at ±10 standard deviations                     │   │
│  │    • Prevents outlier domination in PCA                         │   │
│  └─────────────────────────────────────────────────────────────────┘   │
│         │                                                               │
│         ▼                                                               │
│  ┌─────────────────────────────────────────────────────────────────┐   │
│  │ 2. PCA                                                          │   │
│  │    • n_components = 30 (default)                                │   │
│  │    • SVD solver: arpack                                         │   │
│  │    • Stores in adata.obsm["X_pca"]                              │   │
│  └─────────────────────────────────────────────────────────────────┘   │
│         │                                                               │
│         ▼                                                               │
│  ┌─────────────────────────────────────────────────────────────────┐   │
│  │ 3. NEIGHBORS                                                    │   │
│  │    • k = 15 (default)                                           │   │
│  │    • Uses PCA coordinates                                       │   │
│  │    • Stores connectivity in adata.obsp                          │   │
│  └─────────────────────────────────────────────────────────────────┘   │
│         │                                                               │
│         ▼                                                               │
│  ┌─────────────────────────────────────────────────────────────────┐   │
│  │ 4. LEIDEN                                                       │   │
│  │    • resolution = 0.6 (default)                                 │   │
│  │    • GPU: rapids_singlecell (RAPIDS/cuGraph)                    │   │
│  │    • CPU: igraph backend, n_iterations=2                        │   │
│  │    • random_state = 1337                                        │   │
│  │    • Stores in adata.obs["cluster_lvl0"]                        │   │
│  └─────────────────────────────────────────────────────────────────┘   │
│         │                                                               │
│         ▼                                                               │
│  ┌─────────────────────────────────────────────────────────────────┐   │
│  │ 5. UMAP (optional, --compute-umap)                              │   │
│  │    • 2D embedding for visualization                             │   │
│  │    • Stores in adata.obsm["X_umap"]                             │   │
│  └─────────────────────────────────────────────────────────────────┘   │
│                                                                         │
└─────────────────────────────────────────────────────────────────────────┘
```

### GPU Acceleration

When `--use-gpu` is enabled (default), Stage H uses RAPIDS libraries:

- `rapids_singlecell` for PCA, neighbors, Leiden, UMAP
- Falls back to CPU (scanpy) on GPU errors

### Non-Determinism Warning

```
┌─────────────────────────────────────────────────────────────────────────┐
│  ⚠️  GPU LEIDEN IS NON-DETERMINISTIC                                    │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                         │
│  GPU-accelerated Leiden clustering produces different cluster counts    │
│  across runs even with fixed random_state. This is due to:              │
│                                                                         │
│  • Parallel execution order differences                                 │
│  • Floating-point accumulation variations                               │
│                                                                         │
│  Typical variation: ±5 clusters (e.g., 34-42 for resolution=0.6)        │
│                                                                         │
│  RECOMMENDED WORKFLOW:                                                  │
│  1. Run Stage H once to establish baseline clusters                     │
│  2. Use --annotation-only for subsequent iterations                     │
│  3. Focus on composition similarity, not exact cluster counts           │
│                                                                         │
└─────────────────────────────────────────────────────────────────────────┘
```

---

## 3. Differential Expression

### Wilcoxon Rank-Sum Test

For each cluster, DE testing identifies markers significantly enriched
compared to all other cells:

```python
sc.tl.rank_genes_groups(adata, groupby="cluster_lvl0", method="wilcoxon")
```

Results are stored in `adata.uns[f"de_wilcoxon_{layer}"]`.

### Dynamic K Calculation

The number of top DE genes considered (K) is dynamically computed:

```
K_target = ceil(panel_size × de_top_frac)
K = clamp(K_target, de_min_k, de_max_k)

Default parameters:
  de_top_frac = 0.2
  de_min_k = 3
  de_max_k = 12

Example (panel_size = 40):
  K_target = ceil(40 × 0.2) = 8
  K = clamp(8, 3, 12) = 8
```

### Rank-Based DE Lookup

Instead of binary hit/miss, DE contribution uses rank-weighted scoring:

```
┌─────────────────────────────────────────────────────────────────────────┐
│  DE Rank Lookup Structure (per cluster)                                 │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                         │
│  {                                                                      │
│    "K": 8,                    # Effective K used                        │
│    "K_target": 8,             # Calculated K before clamping            │
│    "n_pos": 15,               # Number of upregulated markers           │
│    "rank_map": {              # Canonical marker → rank (1-indexed)     │
│      "CD45": 1,                                                         │
│      "CD3": 2,                                                          │
│      "CD4": 5,                                                          │
│      ...                                                                │
│    },                                                                   │
│    "topK_names": ["CD45", "CD3", ...],                                  │
│    "topK_scores": [8.5, 7.2, ...]                                       │
│  }                                                                      │
│                                                                         │
│  NOTE: Only UPREGULATED markers (Wilcoxon score > 0) are included       │
│                                                                         │
└─────────────────────────────────────────────────────────────────────────┘
```

---

## 4. Marker Scoring

### Scoring Formula

For each (cluster, cell_type) pair:

```
┌─────────────────────────────────────────────────────────────────────────┐
│                                                                         │
│  score = mean_enrichment + mean_positive + de_component - anti_penalty  │
│                                                                         │
└─────────────────────────────────────────────────────────────────────────┘
```

### Component Breakdown

#### Mean Enrichment

Z-score of cluster expression vs global expression:

```
For each marker m in cell_type:
    enrichment[m] = (cluster_median[m] - global_median[m]) / global_std[m]

mean_enrichment = average(enrichment[all markers])
```

#### Mean Positive Fraction

Fraction of cells expressing markers above threshold:

```
For each marker m in cell_type:
    threshold[m] = global_quantile[m, 0.75]  # Q75
    positive[m] = fraction of cells in cluster where expr >= threshold[m]

mean_positive = average(positive[all markers])
```

#### DE Component (Rank-Weighted)

```
┌─────────────────────────────────────────────────────────────────────────┐
│  DE BONUS CALCULATION                                                   │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                         │
│  For each marker m in cell_type's resolved markers:                     │
│                                                                         │
│    IF m is in cluster's top-K DE genes (upregulated):                   │
│      rank = de_rank_lookup[cluster][m]  # 1-indexed                     │
│      rank_weight = (K - rank + 1) / K   # Higher rank → higher weight   │
│                                                                         │
│      commonness = doc_freq[m] / n_cell_types                            │
│      commonness_penalty = commonness ^ de_commonness_alpha              │
│                                                                         │
│      marker_bonus = rank_weight × (1 - commonness_penalty)              │
│    ELSE:                                                                │
│      marker_bonus = 0                                                   │
│                                                                         │
│  de_component = de_bonus × mean(marker_bonus for all markers)           │
│                                                                         │
│  Default: de_bonus = 0.6, de_commonness_alpha = 0.5                     │
│                                                                         │
└─────────────────────────────────────────────────────────────────────────┘
```

**Example:**

```
Cell type: Immune Cells
Markers: [CD45, CD3, CD19]

Cluster 5 top-K DE: [CD45 (rank 1), VWF (rank 2), CD3 (rank 4), ...]
K = 8

CD45: rank=1, rank_weight=(8-1+1)/8=1.0, commonness=0.1, bonus=1.0×0.68=0.68
CD3:  rank=4, rank_weight=(8-4+1)/8=0.625, commonness=0.2, bonus=0.625×0.55=0.34
CD19: not in top-K, bonus=0

de_component = 0.6 × mean([0.68, 0.34, 0]) = 0.6 × 0.34 = 0.20
```

#### Anti-Marker Penalty

```
┌─────────────────────────────────────────────────────────────────────────┐
│  ANTI-MARKER PENALTY CALCULATION                                        │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                         │
│  For each anti-marker a in cell_type.anti_markers:                      │
│                                                                         │
│    anti_enrichment[a] = (cluster_median[a] - global_median[a]) / std[a] │
│    anti_enrichment[a] = max(anti_enrichment[a], 0)  # Clip negative     │
│                                                                         │
│    anti_positive[a] = fraction of cells >= Q75 threshold                │
│                                                                         │
│  Aggregation modes (--anti-agg):                                        │
│  ┌─────────────────┬─────────────────────────────────────────┐          │
│  │ "max"           │ Strictest: max(anti-markers)            │          │
│  ├─────────────────┼─────────────────────────────────────────┤          │
│  │ "top2mean"      │ DEFAULT: mean of top 2 anti-markers     │          │
│  ├─────────────────┼─────────────────────────────────────────┤          │
│  │ "mean"          │ Lenient: mean of all anti-markers       │          │
│  └─────────────────┴─────────────────────────────────────────┘          │
│                                                                         │
│  anti_penalty = anti_weight × (agg_enrichment + agg_positive)           │
│                                                                         │
│  Default: anti_weight = 0.8                                             │
│                                                                         │
│  HARD GATE: If anti_penalty > 1.0, the cell type is REJECTED            │
│                                                                         │
└─────────────────────────────────────────────────────────────────────────┘
```

**Example Scenarios:**

| Cluster Type | Anti-markers | anti_enrichment | anti_positive | anti_penalty |
|--------------|--------------|-----------------|---------------|--------------|
| True epithelial | CD45, CD3 | -0.5 → 0 (clipped) | 0.05 | 0.8 × (0 + 0.05) = 0.04 |
| Immune contaminated | CD45, CD3 | +2.0 | 0.60 | 0.8 × (2.0 + 0.6) = 2.08 **REJECTED** |
| Mixed cluster | CD45: +1.0, CD3: -0.2 | 0.5 (top2mean) | 0.25 | 0.8 × (0.5 + 0.25) = 0.60 |

### Additional Metrics

| Metric | Definition |
|--------|------------|
| `coverage` | Fraction of cell type's markers present in panel |
| `frac_markers_on` | Fraction of markers above expression threshold |
| `n_markers_on` | Count of markers above threshold |
| `n_cells` | Number of cells in cluster |

---

## 5. Hierarchical Assignment

### Algorithm Overview

```
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

### Gate Checking (`_passes_gate`)

For each candidate at level N:

```
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

### Stop Reasons

| Stop Reason | Condition | Assigned Label |
|-------------|-----------|----------------|
| `leaf_reached` | No children defined in hierarchy | Deepest node |
| `ambiguous_siblings` | Gap between top 2 children < min_gap | Parent node |
| `no_child_passed` | All children failed gates | Current node |
| `no_root_passed` | No root passed gates | "Unassigned" |
| `ambiguous_root` | Top 2 roots gap < root_gap_threshold | "Root1~Root2" |

### Confidence Calculation

```
confidence = min_margin_along_path

Where margin at each level = winner_score - runner_up_score

Special cases:
  • No competition (single candidate): margin = infinity → sentinel 1e6
  • no_root_passed: confidence = 0.0
  • margin_is_infinite flag distinguishes "no competition" from actual margins
```

---

## 6. Per-Cell Voting

### Purpose

When siblings are too close (ambiguous_siblings), per-cell voting provides
**evidence** for downstream refinement. It does NOT override the cluster-level
assignment.

### Algorithm

```
┌─────────────────────────────────────────────────────────────────────────┐
│  PER-CELL VOTING                                                        │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                         │
│  Data Source: adata.X (Z-scored matrix used for PCA/clustering)         │
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
│  Aggregate votes into composition:                                      │
│                                                                         │
│    composition = {                                                      │
│      "Ciliated Epithelium": 0.55,   # 55% voted this                    │
│      "Glandular Epithelium": 0.40,  # 40% voted this                    │
│      "Uncertain": 0.05,             # 5% uncertain                      │
│      "_evidence_only": True,        # Flag: no override                 │
│      "_top_by_score": "Ciliated...",                                    │
│      "_runner_up_by_score": "Gland..",                                  │
│    }                                                                    │
│                                                                         │
└─────────────────────────────────────────────────────────────────────────┘

╔═══════════════════════════════════════════════════════════════════════╗
║  CRITICAL: Voting is EVIDENCE ONLY                                    ║
║                                                                       ║
║  • assigned_label = PARENT (not any child)                            ║
║  • stop_reason = "ambiguous_siblings" (not "voted")                   ║
║  • Composition stored for audit and Stage I refinement hints          ║
╚═══════════════════════════════════════════════════════════════════════╝
```

---

## 7. Gating Parameters

### Default Parameters (Tissue-Agnostic)

```python
BASE_GATING_PARAMS = {
    # Level-dependent thresholds (0=root, 1=child, 2=grandchild, ...)
    "min_coverage":        {0: 0.50, 1: 0.40, 2: 0.30, 3: 0.30},
    "min_pos_frac":        {0: 0.30, 1: 0.20, 2: 0.15, 3: 0.15},
    "min_enrichment":      {0: 0.00, 1: -0.50, 2: -1.00, 3: -1.00},
    "min_gap":             {0: 0.50, 1: 0.30, 2: 0.20, 3: 0.20},
    "min_frac_markers_on": {0: 0.40, 1: 0.30, 2: 0.25, 3: 0.25},

    # Scalar thresholds
    "anti_penalty_hard_gate": 1.0,
    "per_cell_gap_threshold": 0.1,
    "root_gap_threshold": 0.25,
    "anti_agg": "top2mean",

    # Tissue-specific (empty by default)
    "root_hard_requirements": {},
    "root_veto_markers": {},
}
```

### Tissue-Specific Example (Fallopian Tube)

```python
FT_GATING_PARAMS = {
    "root_hard_requirements": {
        "Immune Cells": {
            "marker": "CD45",
            "min_pos_frac": 0.30,
            "min_enrichment": 0.0,
        }
    },
    "root_veto_markers": {
        "Immune Cells": {
            "markers": ["Pan-Cytokeratin", "E-cadherin"],
            "max_pos_frac": 0.20,
        }
    },
}
```

### Parameter Merging

Parameters are merged with priority: `base < tissue_params < override`

```python
from celltype_refinery.core.annotation.gating import merge_gating_params

# Load tissue params from marker map
tissue_params = marker_map.get("_gating_params", {})

# Merge with explicit overrides
params = merge_gating_params(
    base=BASE_GATING_PARAMS,
    tissue_params=tissue_params,
    override={"anti_penalty_hard_gate": 1.5},  # Relaxed threshold
)
```

---

## CLI Reference

### Basic Usage

```bash
python -m celltype_refinery.core.clustering \
    --input merged_data.h5ad \
    --marker-map cell_type_markers.json \
    --output output/stage_h
```

### Full Options

```bash
python -m celltype_refinery.core.clustering \
    --input merged_data.h5ad \
    --marker-map cell_type_markers.json \
    --output output/stage_h \
    --layer batchcorr \
    --n-pcs 30 \
    --neighbors-k 15 \
    --resolution 0.6 \
    --use-gpu \
    --compute-umap \
    --expand-markers \
    --generate-audit-cards \
    --generate-review-checklist \
    --log-dir logs/stage_h
```

### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | required | Path to Stage F merged AnnData |
| `--marker-map` | required | Path to marker map JSON |
| `--output` | required | Output directory |
| `--layer` | `batchcorr` | Expression layer to use |
| `--n-pcs` | `30` | PCA components |
| `--neighbors-k` | `15` | k-NN neighbors |
| `--resolution` | `0.6` | Leiden resolution |
| `--use-gpu` | `True` | Enable GPU acceleration |
| `--no-gpu` | - | Disable GPU, use CPU |
| `--compute-umap` | `False` | Compute UMAP embedding |
| `--annotation-only` | `False` | Skip clustering, use existing clusters |
| `--clustering-only` | `False` | Skip annotation, only cluster |
| `--expand-markers` | `False` | Output per-marker evidence table |
| `--de-bonus` | `0.6` | Maximum DE bonus |
| `--anti-weight` | `0.8` | Anti-marker penalty weight |
| `--anti-agg` | `top2mean` | Anti-marker aggregation mode |

---

## Output Files

### coarse_clusters.h5ad

AnnData with added columns:

```
adata.obs columns:
  cluster_lvl0          - Leiden cluster ID (string)
  cell_type_auto        - Assigned cell type label
  cell_type_auto_root   - Root category (Epithelium, Immune, etc.)
  cell_type_auto_conf   - Confidence score (min margin along path)
  DAPI_intensity        - Technical marker intensities
  ...

adata.uns:
  de_wilcoxon_{layer}   - Differential expression results
  stage_h_params        - Run parameters
  stage_h_gating_params - Gating thresholds used
```

### cluster_annotations.csv

| Column | Description |
|--------|-------------|
| `cluster_id` | Cluster identifier |
| `assigned_label` | Final cell type label |
| `assigned_path` | Full hierarchy path |
| `assigned_level` | Hierarchy depth reached |
| `assigned_score` | Score at assigned level |
| `root_label` | Root category |
| `confidence` | Min margin along path |
| `stop_reason` | Why descent stopped |
| `n_cells` | Cells in cluster |
| `coverage` | Marker coverage |
| `resolved_markers` | Markers used for scoring |
| `is_ambiguous_root` | True if root was ambiguous |
| `composition` | JSON: per-cell vote composition |

### marker_scores.csv

Full scoring matrix with columns:

```
cluster_id, label, path, level, score, coverage, mean_enrichment,
mean_positive_fraction, frac_markers_on, de_component, anti_penalty,
resolved_markers, missing_markers, n_cells, ...
```

### decision_steps.csv

Step-by-step trace:

```
cluster_id, step_idx, parent_label, child_label, child_passed_gate,
child_score, child_coverage, child_pos_frac, child_enrichment,
child_anti_penalty, selected, margin_to_runner_up, fail_reason
```

---

## Execution Modes

### Full Pipeline (Default)

Runs clustering + annotation:

```bash
python -m celltype_refinery.core.clustering \
    --input merged_data.h5ad \
    --marker-map markers.json \
    --output output/stage_h
```

### Annotation Only

Skip clustering, use existing `cluster_lvl0` from input:

```bash
python -m celltype_refinery.core.clustering \
    --input coarse_clusters.h5ad \
    --marker-map markers_v2.json \
    --annotation-only \
    --output output/stage_h_v2
```

**Use case:** Iterate on marker maps or gating parameters without re-clustering.

### Clustering Only

Run clustering without annotation:

```bash
python -m celltype_refinery.core.clustering \
    --input merged_data.h5ad \
    --clustering-only \
    --output output/stage_h
```

**Use case:** Establish clusters first, annotate separately.

---

## Configuration

### Marker Map JSON Structure

```json
{
    "_marker_map_metadata": {
        "version": "v5",
        "tissue": "uterus",
        "created": "2025-01-10"
    },
    "_gating_params": {
        "root_hard_requirements": {
            "Immune Cells": {
                "marker": "CD45",
                "min_pos_frac": 0.30
            }
        },
        "root_veto_markers": {
            "Immune Cells": {
                "markers": ["Pan-Cytokeratin"],
                "max_pos_frac": 0.20
            }
        }
    },
    "Epithelium": {
        "markers": ["EpCAM", "E-cadherin", "Pan-Cytokeratin"],
        "Anti_markers": ["CD45", "VWF"],
        "subtypes": {
            "Ciliated Epithelium": {
                "markers": ["FOXJ1", "Acetylated Tubulin"],
                "Anti_markers": []
            },
            "Glandular Epithelium": {
                "markers": ["PAX8", "MUC1"],
                "Anti_markers": []
            }
        }
    },
    "Immune Cells": {
        "markers": ["CD45"],
        "Anti_markers": ["EpCAM", "VWF"],
        "subtypes": {
            "Lymphoids": {
                "markers": ["CD3", "CD19", "CD56"],
                "subtypes": { ... }
            },
            "Myeloids": {
                "markers": ["CD14", "CD68"],
                "subtypes": { ... }
            }
        }
    }
}
```

### Per-Cell-Type Gating Overrides

Individual cell types can override default gating thresholds:

```json
{
    "Lymphatic Endothelium": {
        "markers": ["Podoplanin", "LYVE1"],
        "Anti_markers": ["CD146"],
        "gating_overrides": {
            "min_coverage": 0.3,
            "min_pos_frac": 0.15,
            "require_anti_markers_off": true
        }
    }
}
```

---

## Design Principles

### 1. Parents as Gates

A cluster must pass the parent's gate before comparing with siblings:

```
✗ WRONG: Compare Epithelium vs Immune vs Neutrophils directly
✓ RIGHT: Must pass Immune gate before comparing Neutrophils vs Eosinophils
```

### 2. Siblings-Only Competition

At each level, only compare nodes that share a parent:

```
Level 0: Epithelium vs Immune vs Mesenchymal vs Endothelium
Level 1 (under Immune): Lymphoids vs Myeloids
Level 2 (under Myeloids): Granulocytes vs Monocytes

Never compare: Neutrophils vs Ciliated Epithelium (different branches)
```

### 3. Traceable Decisions

Every assignment step is logged in `decision_steps.csv`:

- Which candidates were considered
- Which passed/failed gates (and why)
- Score margins between winners and runners-up
- Final selection rationale

### 4. Ambiguity Detection

Close calls are flagged for expert review rather than forced decisions:

- `ambiguous_root`: Multiple roots within 0.25 score margin
- `ambiguous_siblings`: Siblings within min_gap[level] margin

### 5. Evidence vs Override

Per-cell voting provides evidence but never overrides cluster-level assignments.
This maintains consistency while providing refinement hints for Stage I.

---

## Module Structure

```
celltype_refinery/core/
├── clustering/
│   ├── __init__.py          # Package exports
│   ├── __main__.py          # CLI entry point, Stage H orchestration
│   ├── config.py            # Configuration dataclasses
│   ├── engine.py            # ClusteringEngine (Leiden, PCA, etc.)
│   └── de.py                # DERunner (Wilcoxon testing)
│
└── annotation/
    ├── __init__.py          # Package exports
    ├── __main__.py          # Standalone annotation CLI
    ├── engine.py            # AnnotationEngine orchestrator
    ├── marker_loading.py    # MarkerSet, load_marker_sets()
    ├── scoring.py           # compute_marker_scores(), DE rank lookup
    ├── gating.py            # assign_labels_hierarchical(), gate checks
    ├── assignment.py        # annotate_obs() - cell-level mapping
    └── export.py            # CSV export functions
```

---

## See Also

- **Stage I Documentation**: Refinement and subclustering
- **Marker Map Specification**: JSON schema for marker hierarchies
- **Tissue Configuration Guide**: Setting up tissue-specific gating rules
