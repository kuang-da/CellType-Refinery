---
sidebar_position: 3
---

# Marker Scoring Algorithm

:::info Prerequisites
This page explains scoring details. For conceptual background, see [Methodology Overview](./index.md) first.
:::

## Overview

The CellType-Refinery scoring algorithm evaluates how well each cluster matches each candidate cell type in the marker map. The algorithm produces a single numeric score for each (cluster, cell_type) pair, enabling systematic ranking and selection of the best-fitting annotation.

### The Core Formula

The scoring formula for each (cluster, cell_type) pair:

```text
score = mean_enrichment + mean_positive + de_component - anti_penalty
```

Each component captures a different biological signal:

| Component | Range | What It Measures |
|-----------|-------|------------------|
| `mean_enrichment` | Unbounded (typically -2 to +4) | How much higher is expression vs. global average? |
| `mean_positive` | 0 to 1 | What fraction of cells clearly express markers? |
| `de_component` | 0 to ~0.6 | Are markers differentially expressed in this cluster? |
| `anti_penalty` | 0 to unbounded | Is conflicting marker expression present? |

### Score Interpretation

- **Typical winning scores**: 1.5 to 4.0 for well-matched cell types
- **Marginal matches**: 0.5 to 1.5, may indicate subtypes or transitional states
- **Poor matches**: < 0.5, usually indicates wrong cell type assignment
- **Negative scores**: Strong evidence against the cell type

The following sections explain each component in detail with worked examples.

---

## Component 1: Mean Enrichment

The mean enrichment component measures whether marker genes are expressed at higher levels in the cluster compared to the global population. This is the foundation of marker-based annotation.

### Calculation

```text
For each marker m in cell_type:
    enrichment[m] = (cluster_median[m] - global_median[m]) / global_std[m]

mean_enrichment = average(enrichment[all markers])
```

### Step-by-Step Breakdown

1. **Calculate cluster median**: For marker gene `m`, compute the median expression value across all cells in the cluster
2. **Calculate global median**: Compute the median expression across ALL cells in the dataset
3. **Calculate global standard deviation**: Compute std across all cells (measures natural variation)
4. **Compute Z-score**: Standardize the difference using the global std
5. **Average across markers**: Take the mean of all marker Z-scores

### Why Z-Scores?

Z-score normalization is essential because:

- **Different genes have different scales**: CD45 might range 0-10 while CD3 ranges 0-50
- **Comparability**: Z-scores put all markers on the same scale
- **Statistical interpretation**: A Z-score of 2 means "2 standard deviations above average"

### Interpretation Guide

| Z-Score Range | Interpretation | Typical Scenario |
|---------------|----------------|------------------|
| > 2.0 | Strong enrichment | Clear positive marker |
| 1.0 to 2.0 | Moderate enrichment | Typical marker expression |
| 0 to 1.0 | Weak enrichment | Subpopulation or dim expression |
| -1.0 to 0 | No enrichment | Marker not characteristic |
| < -1.0 | Depletion | Wrong cell type |

### Worked Example

```text
Cell type: T Cells
Markers: [CD3, CD4, CD8]

Cluster 3 statistics:
  CD3: cluster_median=4.2, global_median=1.8, global_std=1.5
  CD4: cluster_median=3.1, global_median=1.2, global_std=1.0
  CD8: cluster_median=0.8, global_median=0.9, global_std=0.8

Enrichment calculations:
  CD3: (4.2 - 1.8) / 1.5 = 1.60
  CD4: (3.1 - 1.2) / 1.0 = 1.90
  CD8: (0.8 - 0.9) / 0.8 = -0.13

mean_enrichment = (1.60 + 1.90 + (-0.13)) / 3 = 1.12
```

This cluster shows strong enrichment for CD3 and CD4, but not CD8, suggesting CD4+ T cells rather than CD8+ T cells.

### Edge Cases and Handling

- **Zero variance markers**: If `global_std = 0`, the marker is constant and skipped
- **Missing markers**: Markers not in the expression matrix are excluded from calculation
- **Sparse data**: Median is robust to zeros from dropout in single-cell data

---

## Component 2: Mean Positive Fraction

While mean enrichment measures *how much* higher expression is, mean positive fraction measures *how many* cells in the cluster clearly express the markers. This captures the consistency of marker expression.

### Calculation

```text
For each marker m in cell_type:
    threshold[m] = global_quantile[m, 0.75]  # Q75 by default
    positive[m] = fraction of cells in cluster where expr >= threshold[m]

mean_positive = average(positive[all markers])
```

### Step-by-Step Breakdown

1. **Determine threshold**: Calculate the 75th percentile of expression across ALL cells
2. **Count positive cells**: For each cell in the cluster, check if expression >= threshold
3. **Calculate fraction**: Divide positive count by total cells in cluster
4. **Average across markers**: Mean of all marker positive fractions

### Why Q75 (75th Percentile)?

The 75th percentile threshold is chosen because:

- **Noise filtering**: Single-cell data has substantial technical noise and dropout
- **Clear signal**: Q75 identifies cells with unambiguous expression
- **Conservative**: Avoids false positives from background noise
- **Biologically meaningful**: Captures the "clearly positive" population

### Tuning the Threshold

The `positive_quantile` parameter allows adjustment:

| Setting | Use Case | Effect |
|---------|----------|--------|
| 0.90 | Very stringent | Only the brightest cells count as positive |
| 0.75 | Default | Balanced, works for most datasets |
| 0.50 | Sparse data | More sensitive, useful for low-expression markers |
| 0.25 | Very sparse | Very sensitive, may include noise |

### Worked Example

```text
Cell type: B Cells
Markers: [CD19, CD20, CD79a]

Global Q75 thresholds:
  CD19: Q75 = 2.5
  CD20: Q75 = 1.8
  CD79a: Q75 = 3.0

Cluster 7 (100 cells):
  CD19: 85 cells >= 2.5 → positive = 0.85
  CD20: 78 cells >= 1.8 → positive = 0.78
  CD79a: 72 cells >= 3.0 → positive = 0.72

mean_positive = (0.85 + 0.78 + 0.72) / 3 = 0.78
```

A mean positive fraction of 0.78 indicates strong, consistent marker expression.

### Interpreting Mean Positive

| Range | Interpretation | Typical Scenario |
|-------|----------------|------------------|
| > 0.8 | Excellent | Clear cell type, homogeneous cluster |
| 0.6 - 0.8 | Good | Typical positive cluster |
| 0.4 - 0.6 | Moderate | Mixed population or dim expression |
| 0.2 - 0.4 | Weak | Subpopulation or wrong cell type |
| < 0.2 | Poor | Likely wrong cell type |

---

## Component 3: DE Component (Rank-Weighted)

The differential expression (DE) component rewards cell types whose markers appear in the cluster's top differentially expressed genes. This leverages unsupervised DE analysis to validate marker-based annotations.

### Conceptual Overview

```text
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
└─────────────────────────────────────────────────────────────────────────┘
```

### Step-by-Step Breakdown

#### Step 1: Determine K (Number of Top DE Genes)

```text
K = clip(n_markers × de_top_frac, de_min_k, de_max_k)
```

Where:
- `n_markers` = total unique markers in the panel
- `de_top_frac` = fraction of panel to consider (default 0.2)
- `de_min_k` = minimum K value (default 3)
- `de_max_k` = maximum K value (default 12)

Example: Panel with 40 markers → K = clip(40 × 0.2, 3, 12) = clip(8, 3, 12) = 8

#### Step 2: Check DE Membership

For each marker in the cell type, check if it appears in the cluster's top-K upregulated genes.

#### Step 3: Calculate Rank Weight

```text
rank_weight = (K - rank + 1) / K
```

This creates a linear decay from 1.0 (rank 1) to 1/K (rank K):

| Rank | K=8 Weight | K=12 Weight |
|------|------------|-------------|
| 1 | 1.000 | 1.000 |
| 2 | 0.875 | 0.917 |
| 3 | 0.750 | 0.833 |
| 4 | 0.625 | 0.750 |
| 5 | 0.500 | 0.667 |
| 6 | 0.375 | 0.583 |

#### Step 4: Calculate Commonness Penalty

Markers that appear in many cell type definitions are less informative:

```text
commonness = doc_freq[m] / n_cell_types
commonness_penalty = commonness ^ de_commonness_alpha
effective_weight = 1 - commonness_penalty
```

Example with `de_commonness_alpha = 0.5`:

| Marker Usage | Commonness | Penalty | Effective Weight |
|--------------|------------|---------|------------------|
| Unique (1/20) | 0.05 | 0.22 | 0.78 |
| Moderate (5/20) | 0.25 | 0.50 | 0.50 |
| Common (10/20) | 0.50 | 0.71 | 0.29 |
| Ubiquitous (20/20) | 1.00 | 1.00 | 0.00 |

#### Step 5: Combine and Average

```text
marker_bonus = rank_weight × (1 - commonness_penalty)
de_component = de_bonus × mean(marker_bonus for all markers)
```

### Worked Example

```text
Cell type: Immune Cells
Markers: [CD45, CD3, CD19]

Cluster 5 top-K DE: [CD45 (rank 1), VWF (rank 2), CD3 (rank 4), ...]
K = 8, de_bonus = 0.6, de_commonness_alpha = 0.5

Document frequencies (20 cell types total):
  CD45: appears in 2 cell types → commonness = 0.10
  CD3: appears in 4 cell types → commonness = 0.20
  CD19: appears in 2 cell types → commonness = 0.10

Calculations:
  CD45: rank=1
        rank_weight = (8-1+1)/8 = 1.0
        commonness_penalty = 0.10^0.5 = 0.316
        bonus = 1.0 × (1 - 0.316) = 0.68

  CD3:  rank=4
        rank_weight = (8-4+1)/8 = 0.625
        commonness_penalty = 0.20^0.5 = 0.447
        bonus = 0.625 × (1 - 0.447) = 0.34

  CD19: not in top-K
        bonus = 0

de_component = 0.6 × mean([0.68, 0.34, 0]) = 0.6 × 0.34 = 0.20
```

### Why Rank Weighting?

Rank weighting captures the intuition that:

1. **Top-ranked DE genes are most distinctive** for the cluster
2. **Lower-ranked genes add some evidence** but are less specific
3. **Non-DE genes contribute nothing** (they don't distinguish the cluster)

### Why Commonness Penalty?

The commonness penalty prevents ubiquitous markers from dominating:

- CD45 (pan-immune) shouldn't boost all immune types equally
- Cell-type-specific markers (e.g., CD19 for B cells) should matter more
- Prevents "generic" markers from creating false matches

---

## Component 4: Anti-Marker Penalty

Anti-markers define genes that should NOT be expressed in a cell type. The anti-marker penalty reduces scores when conflicting markers are detected, preventing misannotation.

### Conceptual Overview

```text
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
│  ┌─────────────────┬─────────────────────────────────────────────────┐  │
│  │ "max"           │ Strictest: max(anti-markers)                    │  │
│  ├─────────────────┼─────────────────────────────────────────────────┤  │
│  │ "top2mean"      │ DEFAULT: mean of top 2 anti-markers             │  │
│  ├─────────────────┼─────────────────────────────────────────────────┤  │
│  │ "mean"          │ Lenient: mean of all anti-markers               │  │
│  └─────────────────┴─────────────────────────────────────────────────┘  │
│                                                                         │
│  anti_penalty = anti_weight × (agg_enrichment + agg_positive)           │
│                                                                         │
│  HARD GATE: If anti_penalty > 1.0, the cell type is REJECTED            │
│                                                                         │
└─────────────────────────────────────────────────────────────────────────┘
```

### Step-by-Step Breakdown

#### Step 1: Calculate Anti-Enrichment

For each anti-marker, compute a Z-score similar to mean enrichment:

```text
anti_enrichment[a] = (cluster_median[a] - global_median[a]) / global_std[a]
anti_enrichment[a] = max(anti_enrichment[a], 0)  # Only positive values penalize
```

The clipping to zero is important: we only penalize when anti-markers are *enriched*, not when they're absent (which is expected).

#### Step 2: Calculate Anti-Positive Fraction

```text
anti_positive[a] = fraction of cells in cluster where expr >= Q75 threshold
```

This measures how many cells express the anti-marker.

#### Step 3: Aggregate Across Anti-Markers

The aggregation method determines how multiple anti-markers combine:

| Method | Formula | Use Case |
|--------|---------|----------|
| `max` | `max(anti_enrichment), max(anti_positive)` | Strict: any conflict rejects |
| `top2mean` | Mean of top 2 values | Balanced: robust to single outliers |
| `mean` | Mean of all values | Lenient: only strong conflicts penalize |

#### Step 4: Compute Final Penalty

```text
anti_penalty = anti_weight × (agg_enrichment + agg_positive)
```

### Hard Gate Behavior

When `anti_penalty > 1.0`, the cell type is completely rejected:

- The score is set to a large negative value (e.g., -999)
- This cell type cannot be selected regardless of other scores
- Prevents obvious biological conflicts

### Example Scenarios

| Scenario | Anti-Enrichment | Anti-Positive | Penalty (w=0.5) | Result |
|----------|-----------------|---------------|-----------------|--------|
| Clean | 0.0 | 0.05 | 0.025 | Minor |
| Mild conflict | 0.5 | 0.15 | 0.325 | Moderate |
| Moderate conflict | 1.2 | 0.35 | 0.775 | Significant |
| Strong conflict | 2.5 | 0.65 | 1.575 | **REJECTED** |

### Worked Example

```text
Cell type: CD4+ T Cells
Anti-markers: [CD8, CD19, CD14]

Cluster 8 statistics:
  CD8:  cluster_median=2.1, global_median=0.8, global_std=0.9
        enrichment = (2.1-0.8)/0.9 = 1.44, positive = 0.25
  CD19: cluster_median=0.2, global_median=0.5, global_std=0.6
        enrichment = max((0.2-0.5)/0.6, 0) = 0, positive = 0.02
  CD14: cluster_median=0.1, global_median=0.3, global_std=0.4
        enrichment = max((0.1-0.3)/0.4, 0) = 0, positive = 0.01

Using top2mean aggregation:
  Top 2 enrichments: [1.44, 0] → mean = 0.72
  Top 2 positives: [0.25, 0.02] → mean = 0.135

anti_penalty = 0.5 × (0.72 + 0.135) = 0.43
```

This penalty of 0.43 significantly reduces the score, flagging potential CD8 contamination.

### Choosing Anti-Aggregation Mode

| Mode | Best For | Behavior |
|------|----------|----------|
| `max` | Well-defined types with clear exclusions | Single strong conflict rejects |
| `top2mean` | Most datasets (default) | Robust to noise, catches real conflicts |
| `mean` | Complex panels, many anti-markers | Only consistent conflicts penalize |

---

## Putting It All Together

### Complete Scoring Example

```text
Cell type: NK Cells
Markers: [CD56, CD16, NKG2D]
Anti-markers: [CD3, CD19]

Cluster 12 (150 cells)

Step 1: Mean Enrichment
  CD56: Z = 2.1
  CD16: Z = 1.8
  NKG2D: Z = 1.5
  mean_enrichment = 1.80

Step 2: Mean Positive
  CD56: 0.82
  CD16: 0.75
  NKG2D: 0.68
  mean_positive = 0.75

Step 3: DE Component
  CD56: rank 2, bonus = 0.52
  CD16: rank 5, bonus = 0.28
  NKG2D: not in top-K, bonus = 0
  de_component = 0.6 × 0.27 = 0.16

Step 4: Anti-Penalty
  CD3: enrichment = 0.3, positive = 0.08
  CD19: enrichment = 0, positive = 0.01
  anti_penalty = 0.5 × (0.15 + 0.045) = 0.10

Final Score:
  score = 1.80 + 0.75 + 0.16 - 0.10 = 2.61
```

A score of 2.61 indicates a strong NK cell match.

---

## Parameter Reference

The following parameters control scoring behavior. They can be set via command-line arguments or configuration files.

### Positive Fraction Parameters

| Parameter | Default | Effect | When to Adjust |
|-----------|---------|--------|----------------|
| `positive_quantile` | 0.75 | Threshold for "positive" cells | Sparse data → lower to 0.5-0.6 |

### DE Component Parameters

| Parameter | Default | Effect | When to Adjust |
|-----------|---------|--------|----------------|
| `de_bonus` | 0.5-0.6 | Weight of DE component | Good DE results → increase to 0.7 |
| `de_top_frac` | 0.2 | Fraction of panel for K | Large panel (>50) → decrease to 0.15 |
| `de_min_k` | 3 | Minimum top-K | Small panels → keep at 3 |
| `de_max_k` | 12 | Maximum top-K | Large panels (>60) → increase to 15 |
| `de_commonness_alpha` | 0.5 | Commonness penalty strength | Many shared markers → increase to 0.7 |

### Anti-Marker Parameters

| Parameter | Default | Effect | When to Adjust |
|-----------|---------|--------|----------------|
| `anti_weight` | 0.5-0.8 | Penalty strength | Wrong assignments → increase to 0.8-1.0 |
| `anti_agg` | top2mean | Aggregation method | Strict → use max; lenient → use mean |

### Parameter Interaction Guidelines

1. **Sparse data**: Lower `positive_quantile`, may need higher `de_bonus`
2. **Dense data**: Default parameters usually work well
3. **Many similar cell types**: Increase `anti_weight`, consider `max` aggregation
4. **Large marker panel**: Decrease `de_top_frac`, increase `de_max_k`
5. **Poor DE results**: Decrease `de_bonus`, rely more on expression

---

## Debugging Scores

When scores seem incorrect, check these common issues:

### Low Scores for Expected Cell Types

1. **Check marker expression**: Are markers actually expressed in the cluster?
2. **Verify DE results**: Are markers appearing in DE gene lists?
3. **Inspect anti-markers**: Is an anti-marker unexpectedly enriched?
4. **Review thresholds**: Is `positive_quantile` too stringent?

### Multiple High-Scoring Cell Types

1. **Check for subtypes**: Parent and child types may both score well
2. **Review anti-markers**: Add more discriminating anti-markers
3. **Inspect markers**: Are marker sets too overlapping?
4. **Consider DE bonus**: Increase to differentiate similar types

### Unexpected Rejections (Hard Gate)

1. **Verify anti-marker enrichment**: Is rejection justified?
2. **Check for batch effects**: Technical artifacts can cause false conflicts
3. **Review aggregation mode**: Switch from `max` to `top2mean`
4. **Adjust anti_weight**: Lower if rejections seem too aggressive

---

## Mathematical Properties

### Score Range

- **Theoretical minimum**: Unbounded negative (due to enrichment Z-scores)
- **Practical minimum**: -2 to 0 for poor matches
- **Theoretical maximum**: Unbounded positive
- **Practical maximum**: 3 to 5 for excellent matches

### Component Independence

The four components are largely independent:

- A cluster can have high enrichment but low positive fraction (bimodal expression)
- DE bonus can be high even with moderate enrichment (distinct markers)
- Anti-penalty applies regardless of positive component values

### Normalization Considerations

The algorithm does NOT normalize scores across cell types. This is intentional:

- Absolute scores indicate match quality
- Enables setting quality thresholds (e.g., "minimum score of 1.0")
- Allows comparison across different runs

---

## See Also

- [Hierarchical Gating Algorithm](./hierarchical-gating-algorithm.md) - How scores determine labels at each hierarchy level
- [Tuning Guide](./tuning-guide.md) - Practical guidance for adjusting scoring parameters
- [Methodology Overview](./index.md) - Conceptual background on the annotation approach
