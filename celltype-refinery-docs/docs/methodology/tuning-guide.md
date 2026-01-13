---
sidebar_position: 6
---

# Tuning Guide

:::tip Expected Usage
Parameter tuning is a normal part of the annotation workflow. Almost every dataset benefits from some adjustment to match your tissue type and marker panel.
:::

## Overview

This guide helps you adjust parameters when results need improvement. Each section addresses a common scenario with specific recommendations.

Cell type annotation is not a one-size-fits-all process. Different tissues have different marker expression patterns, different sequencing depths produce different signal-to-noise ratios, and different marker panels have different coverage of the cell type hierarchy. This guide provides systematic approaches to diagnosing and fixing common issues.

### How to Use This Guide

1. **Identify your symptoms** - Look through the scenarios below to find matching issues
2. **Diagnose the cause** - Use the diagnostic checks to understand why the issue occurs
3. **Apply targeted fixes** - Start with the most likely fix and iterate
4. **Validate improvements** - Check that fixes improve results without creating new problems

### General Tuning Philosophy

- **Start with defaults** - The default parameters work well for many datasets
- **Change one thing at a time** - Makes it easier to understand what helps
- **Use biological knowledge** - Let your understanding of the tissue guide adjustments
- **Document changes** - Keep track of what you tried and why

---

## Scenario 1: Too Many "Unassigned" Clusters

**Symptoms**: Many clusters labeled "Unassigned" with stop_reason = "no_root_passed"

This is the most common issue users encounter. It indicates that clusters are failing to meet the minimum requirements for any root cell type in the hierarchy.

### Diagnostic Checklist

Before adjusting parameters, verify:

1. **Marker panel coverage**
   - Do you have markers for all expected cell types?
   - Are the markers in your panel present in the marker map?
   - Check for naming mismatches (e.g., "CD45" vs "PTPRC")

2. **Expression layer quality**
   - Is there sufficient signal in your expression data?
   - Are expression values in the expected range (log-normalized)?
   - Check for batch effects or technical artifacts

3. **Root requirements alignment**
   - Are `root_hard_requirements` markers present in your panel?
   - Are the thresholds appropriate for your data type?

### Primary Fixes

#### Fix 1.1: Lower Coverage Requirements

The `min_coverage` parameter controls what fraction of a cell type's markers must be expressed. Lower values allow clusters with partial marker expression to pass.

```bash
# Default at root level is 0.50
# Try lowering progressively
ctr annotate ... --gating-params '{"level_defaults": {"0": {"min_coverage": 0.40}}}'

# For sparse data (e.g., spatial transcriptomics)
ctr annotate ... --gating-params '{"level_defaults": {"0": {"min_coverage": 0.30}}}'

# For very sparse panels (< 50 markers)
ctr annotate ... --gating-params '{"level_defaults": {"0": {"min_coverage": 0.20}}}'
```

**When to use**: Your marker panel doesn't include all canonical markers for each cell type, or expression is generally sparse.

#### Fix 1.2: Lower Positive Fraction Requirements

The `min_pos_frac` parameter controls what fraction of cells must express the marker for it to count as "positive" for the cluster.

```bash
# Default is 0.30 (30% of cells must express marker)
# Lower for heterogeneous clusters
ctr annotate ... --gating-params '{"level_defaults": {"0": {"min_pos_frac": 0.20}}}'

# For tissues with high heterogeneity
ctr annotate ... --gating-params '{"level_defaults": {"0": {"min_pos_frac": 0.15}}}'
```

**When to use**: Clusters contain mixed populations or markers have bimodal expression.

#### Fix 1.3: Review Hard Requirements

If specific root types are consistently failing, check whether hard requirements are too strict.

```json
{
    "_gating_params": {
        "root_hard_requirements": {
            "Immune Cells": {
                "marker": "CD45",
                "min_pos_frac": 0.30,  // Lower from 0.50 if CD45 expression is low
                "min_enrichment": 0.0
            }
        }
    }
}
```

**When to use**: A specific cell type consistently fails despite clear marker expression.

#### Fix 1.4: Check Expression Layer

Verify your expression data has the expected characteristics:

```python
import scanpy as sc
import numpy as np

# Load data
adata = sc.read_h5ad("your_data.h5ad")

# Check expression range
print(f"Expression range: {adata.X.min():.2f} to {adata.X.max():.2f}")
print(f"Mean expression: {adata.X.mean():.2f}")

# Check for specific markers
markers_to_check = ["CD45", "CD3", "CD19", "CD14"]
for marker in markers_to_check:
    if marker in adata.var_names:
        expr = adata[:, marker].X.toarray().flatten()
        pos_frac = (expr > 0).mean()
        print(f"{marker}: {pos_frac:.1%} positive, mean={expr.mean():.2f}")
```

### Advanced Diagnostics

If the basic fixes don't help, examine the detailed outputs:

```python
import pandas as pd

# Check gate results for a specific cluster
results = pd.read_csv("output_dir/gating_details.csv")
cluster_results = results[results['cluster'] == 'problematic_cluster']

# Look at coverage scores for each root type
print(cluster_results[['root_type', 'coverage', 'pos_frac', 'passed']])
```

---

## Scenario 2: Annotation Stops Too Early (Parent-Level Only)

**Symptoms**: Clusters assigned at parent level (e.g., "Epithelium" instead of "Ciliated Epithelium")

This occurs when the algorithm successfully identifies a root or intermediate type but cannot confidently descend to more specific child types.

### Diagnostic Checklist

1. **Child marker definitions**
   - Are child types well-defined with unique markers in the marker map?
   - Do child types have sufficient marker coverage in your panel?

2. **Gap thresholds**
   - Is `min_gap` preventing descent when scores are close?
   - Are child types genuinely difficult to distinguish in your data?

3. **Expression patterns**
   - Do the expected child-specific markers show differential expression?

### Primary Fixes

#### Fix 2.1: Lower Gap Thresholds at Child Levels

The `min_gap` parameter controls the required score difference between the top candidate and alternatives. Lower values allow the algorithm to make decisions with smaller margins.

```bash
# Default is often 0.30-0.50 depending on level
# Lower for deeper annotation
ctr annotate ... --gating-params '{
    "level_defaults": {
        "1": {"min_gap": 0.20},
        "2": {"min_gap": 0.15}
    }
}'
```

**When to use**: Child types have similar marker profiles, or your panel doesn't strongly distinguish subtypes.

#### Fix 2.2: Adjust Coverage at Child Levels

Child types may need different coverage requirements than roots:

```bash
ctr annotate ... --gating-params '{
    "level_defaults": {
        "0": {"min_coverage": 0.50},
        "1": {"min_coverage": 0.35},
        "2": {"min_coverage": 0.25}
    }
}'
```

**When to use**: Child types have fewer markers defined, or markers are more specific/sparse.

#### Fix 2.3: Add Distinguishing Markers

Update your marker map with markers that distinguish child types:

```json
{
    "Epithelium": {
        "markers": ["EPCAM", "KRT8", "KRT18"],
        "children": {
            "Ciliated Epithelium": {
                "markers": ["FOXJ1", "TUBB4B", "DNAH5", "CCDC78"],
                "anti_markers": ["MUC5B", "SCGB1A1"]
            },
            "Secretory Epithelium": {
                "markers": ["MUC5B", "SCGB1A1", "MUC16"],
                "anti_markers": ["FOXJ1", "TUBB4B"]
            }
        }
    }
}
```

**When to use**: Child types are biologically distinct but markers don't adequately capture the differences.

#### Fix 2.4: Check Child Marker Expression

Verify that expected markers for child types are actually expressed:

```python
# For clusters stuck at "Epithelium"
epithelial_clusters = adata[adata.obs['annotation'] == 'Epithelium']

# Check ciliated markers
ciliated_markers = ['FOXJ1', 'TUBB4B', 'DNAH5']
for marker in ciliated_markers:
    if marker in adata.var_names:
        expr = epithelial_clusters[:, marker].X.toarray().flatten()
        print(f"{marker}: {(expr > 0).mean():.1%} positive")
```

---

## Scenario 3: Wrong Cell Types Assigned

**Symptoms**: High confidence but biologically incorrect assignments

This is the most concerning scenario because it can lead to incorrect biological conclusions. It typically indicates that the marker definitions or anti-marker constraints are insufficient.

### Diagnostic Checklist

1. **Marker specificity**
   - Are assigned markers actually specific to that cell type in your tissue?
   - Could markers be expressed by multiple cell types?

2. **Anti-marker definitions**
   - Are there anti-markers defined that should prevent the assignment?
   - Is `anti_weight` high enough to penalize conflicting expression?

3. **Root veto markers**
   - Should certain markers completely block a root assignment?

### Primary Fixes

#### Fix 3.1: Add Anti-Markers

Anti-markers penalize scores when markers that should NOT be expressed are present:

```json
{
    "T Cells": {
        "markers": ["CD3D", "CD3E", "CD3G"],
        "anti_markers": ["CD19", "CD20", "MS4A1", "CD14", "CD68"]
    },
    "B Cells": {
        "markers": ["CD19", "CD20", "MS4A1", "CD79A"],
        "anti_markers": ["CD3D", "CD3E", "CD14", "CD68"]
    }
}
```

#### Fix 3.2: Increase Anti-Weight

Make anti-markers more impactful on the score:

```bash
# Default is typically 0.5
# Increase for stricter conflict handling
ctr annotate ... --scoring-params '{"anti_weight": 0.8}'

# For very strict anti-marker enforcement
ctr annotate ... --scoring-params '{"anti_weight": 1.0}'
```

**When to use**: Anti-markers are defined but wrong assignments still occur.

#### Fix 3.3: Add Root Hard Requirements

Require specific markers for root assignment:

```json
{
    "_gating_params": {
        "root_hard_requirements": {
            "T Cells": {
                "marker": "CD3D",
                "min_pos_frac": 0.40,
                "min_enrichment": 0.0
            },
            "B Cells": {
                "marker": "CD19",
                "min_pos_frac": 0.30,
                "min_enrichment": 0.0
            }
        }
    }
}
```

**When to use**: Specific lineage markers are definitive for cell type identity.

#### Fix 3.4: Add Root Veto Markers

Completely block assignment when incompatible markers are expressed:

```json
{
    "_gating_params": {
        "root_veto_markers": {
            "T Cells": {
                "markers": ["CD19", "MS4A1"],
                "max_pos_frac": 0.15
            },
            "B Cells": {
                "markers": ["CD3D", "CD3E"],
                "max_pos_frac": 0.15
            },
            "Epithelium": {
                "markers": ["CD45"],
                "max_pos_frac": 0.20
            }
        }
    }
}
```

**When to use**: Certain marker combinations are biologically impossible.

#### Fix 3.5: Change Anti-Aggregation Method

For stricter anti-marker handling:

```bash
# Default is 'top2mean' (average of top 2 anti-markers)
# Use 'max' for strictest enforcement
ctr annotate ... --scoring-params '{"anti_agg": "max"}'
```

**When to use**: Even one strongly expressed anti-marker should disqualify the type.

---

## Scenario 4: Refinement Too Aggressive

**Symptoms**: Too many clusters being subclustered when they should be kept intact

The refinement step identifies clusters that may benefit from additional subclustering. When too aggressive, this creates unnecessary fragmentation.

### Diagnostic Checklist

1. **Score threshold**
   - Is the heterogeneity score threshold too low?
   - Are many clusters just barely triggering refinement?

2. **Minimum cell counts**
   - Are small clusters being subclustered unnecessarily?

3. **Heterogeneity interpretation**
   - Is biological heterogeneity being confused with technical noise?

### Primary Fixes

#### Fix 4.1: Raise Score Threshold

Higher thresholds mean only clearly heterogeneous clusters are flagged:

```bash
# Default is 1.0
# Raise for more conservative refinement
ctr annotate ... --refinement-params '{"score_threshold": 1.5}'

# For very conservative refinement
ctr annotate ... --refinement-params '{"score_threshold": 2.0}'
```

**When to use**: Many clusters are being flagged that appear homogeneous on inspection.

#### Fix 4.2: Increase Minimum Cell Requirement

Prevent subclustering of small clusters:

```bash
# Default is 500
# Raise for larger datasets
ctr annotate ... --refinement-params '{"min_cells": 1000}'

# For very large datasets (> 100k cells)
ctr annotate ... --refinement-params '{"min_cells": 2000}'
```

**When to use**: Small clusters are being unnecessarily fragmented.

#### Fix 4.3: Adjust Heterogeneity Gap

Control when relabeling is preferred over subclustering:

```bash
# Default is 0.3
# Higher values favor relabeling over subclustering
ctr annotate ... --refinement-params '{"heterogeneity_gap": 0.4}'
```

**When to use**: Clusters are being subclustered when a different label would be more appropriate.

#### Fix 4.4: Disable Refinement

For initial exploration or when refinement is not needed:

```bash
ctr annotate ... --skip-refinement
```

**When to use**: You want to focus on primary annotation without refinement complexity.

---

## Scenario 5: Refinement Too Conservative

**Symptoms**: Heterogeneous clusters not being split when they should be

The opposite of Scenario 4 - clusters that clearly contain multiple cell types are being kept together.

### Primary Fixes

#### Fix 5.1: Lower Score Threshold

```bash
# Default is 1.0
# Lower for more aggressive refinement
ctr annotate ... --refinement-params '{"score_threshold": 0.7}'
```

#### Fix 5.2: Lower Minimum Cell Requirement

```bash
# Allow smaller clusters to be flagged
ctr annotate ... --refinement-params '{"min_cells": 300}'
```

#### Fix 5.3: Lower Heterogeneity Gap

```bash
# Lower values favor subclustering over relabeling
ctr annotate ... --refinement-params '{"heterogeneity_gap": 0.2}'
```

---

## Scenario 6: Results Inconsistent Between Runs

**Cause**: GPU non-determinism in Leiden clustering

When using GPU-accelerated clustering (via rapids-singlecell), results may vary slightly between runs due to floating-point non-determinism.

### Primary Fixes

#### Fix 6.1: Disable GPU for Reproducibility

```bash
# Force CPU-based clustering
ctr annotate ... --no-gpu
```

**Trade-off**: Slower execution but perfectly reproducible results.

#### Fix 6.2: Fix Clustering, Iterate Annotation

Run clustering once, then iterate on annotation parameters:

```bash
# First run - includes clustering
ctr annotate input.h5ad --output-dir run1/

# Subsequent runs - use existing clustering
ctr annotate run1/output.h5ad --annotation-only --output-dir run2/
```

**When to use**: You want to compare different annotation parameters without re-clustering.

#### Fix 6.3: Set Random Seeds

While not guaranteed for GPU operations, setting seeds improves reproducibility:

```bash
ctr annotate ... --random-seed 42
```

---

## Scenario 7: High Confidence but Biologically Wrong

**Issue**: Markers may not be specific for your tissue

Different tissues have different expression patterns. A marker highly specific in one tissue may be expressed by multiple types in another.

### Primary Fixes

#### Fix 7.1: Review Marker Map for Tissue Appropriateness

Audit your marker map against known tissue biology:

```python
# Load marker map
import json
with open("marker_map.json") as f:
    markers = json.load(f)

# List all markers for review
def list_all_markers(node, path=""):
    all_markers = []
    if "markers" in node:
        for m in node["markers"]:
            all_markers.append((path, m, "positive"))
    if "anti_markers" in node:
        for m in node["anti_markers"]:
            all_markers.append((path, m, "anti"))
    if "children" in node:
        for child_name, child_node in node["children"].items():
            all_markers.extend(list_all_markers(child_node, f"{path}/{child_name}"))
    return all_markers

# Review marker list
for path, marker, marker_type in list_all_markers(markers):
    print(f"{path}: {marker} ({marker_type})")
```

#### Fix 7.2: Create Tissue-Specific Marker Map

Start with the default and customize for your tissue:

```bash
# Export default marker map
ctr export-markers --output my_tissue_markers.json

# Edit for your tissue
# Then use your customized map
ctr annotate ... --marker-map my_tissue_markers.json
```

#### Fix 7.3: Add Tissue-Specific Gating Parameters

Include custom gating in your marker map:

```json
{
    "_gating_params": {
        "root_hard_requirements": { ... },
        "root_veto_markers": { ... }
    },
    "Immune Cells": { ... },
    "Epithelium": { ... }
}
```

---

## Parameter Reference

### Scoring Parameters

These parameters control how marker scores are calculated.

| Parameter | Default | Range | Effect | When to Adjust |
|-----------|---------|-------|--------|----------------|
| `positive_quantile` | 0.75 | 0.5-0.95 | Percentile for positive expression threshold | Sparse data: lower to 0.5-0.6 |
| `de_bonus` | 0.5-0.6 | 0.0-1.0 | Weight for differential expression | Good DE results: increase to 0.7 |
| `anti_weight` | 0.5-0.8 | 0.0-1.0 | Penalty weight for anti-markers | Wrong assignments: increase |
| `anti_agg` | top2mean | top2mean/max/mean | How to aggregate anti-markers | Strict: use max |
| `coverage_weight` | 0.3 | 0.0-1.0 | Weight for marker coverage in score | Low panel coverage: decrease |

#### Detailed Parameter Explanations

**`positive_quantile`**: Determines the expression threshold for considering a marker "positive." A value of 0.75 means the threshold is set at the 75th percentile of non-zero expression. Lower values make it easier to be positive, higher values are more stringent.

```bash
# For sparse spatial data
ctr annotate ... --scoring-params '{"positive_quantile": 0.50}'

# For deep scRNA-seq
ctr annotate ... --scoring-params '{"positive_quantile": 0.85}'
```

**`de_bonus`**: Controls how much differential expression contributes to the score. Higher values emphasize markers that are specifically enriched in the cluster compared to background.

```bash
# Emphasize DE
ctr annotate ... --scoring-params '{"de_bonus": 0.7}'

# De-emphasize DE (rely more on absolute expression)
ctr annotate ... --scoring-params '{"de_bonus": 0.3}'
```

**`anti_weight`**: Controls the penalty applied when anti-markers are expressed. Higher values more strongly penalize "forbidden" marker expression.

```bash
# Moderate anti-marker penalty
ctr annotate ... --scoring-params '{"anti_weight": 0.5}'

# Strong anti-marker penalty
ctr annotate ... --scoring-params '{"anti_weight": 0.9}'
```

### Gating Parameters

These parameters control the hierarchical decision-making process.

| Parameter | Level | Default | Effect |
|-----------|-------|---------|--------|
| `min_coverage` | 0 | 0.50 | Fraction of markers that must be positive |
| `min_coverage` | 1+ | 0.30-0.40 | Generally lower at deeper levels |
| `min_pos_frac` | 0 | 0.30 | Fraction of cells positive for counted markers |
| `min_gap` | all | 0.20-0.50 | Required score gap between top candidates |
| `max_ambiguity` | all | 3 | Maximum types within gap to still proceed |

#### Level-Specific Configuration

```bash
# Configure each level separately
ctr annotate ... --gating-params '{
    "level_defaults": {
        "0": {
            "min_coverage": 0.50,
            "min_pos_frac": 0.30,
            "min_gap": 0.40
        },
        "1": {
            "min_coverage": 0.40,
            "min_pos_frac": 0.25,
            "min_gap": 0.30
        },
        "2": {
            "min_coverage": 0.30,
            "min_pos_frac": 0.20,
            "min_gap": 0.20
        }
    }
}'
```

### Refinement Parameters

These parameters control the refinement step that identifies heterogeneous clusters.

| Parameter | Default | Range | Effect |
|-----------|---------|-------|--------|
| `score_threshold` | 1.0 | 0.5-3.0 | Heterogeneity score to trigger refinement |
| `min_cells` | 500 | 100-5000 | Minimum cluster size for refinement |
| `heterogeneity_gap` | 0.3 | 0.1-0.5 | Gap to prefer relabel over subcluster |
| `max_iterations` | 3 | 1-5 | Maximum refinement iterations |

#### Combined Refinement Configuration

```bash
# Conservative refinement
ctr annotate ... --refinement-params '{
    "score_threshold": 1.5,
    "min_cells": 1000,
    "heterogeneity_gap": 0.4,
    "max_iterations": 2
}'

# Aggressive refinement
ctr annotate ... --refinement-params '{
    "score_threshold": 0.7,
    "min_cells": 300,
    "heterogeneity_gap": 0.2,
    "max_iterations": 5
}'
```

---

## Tissue-Specific Configuration

### Adding Gating Parameters to Marker Map

The marker map can include a `_gating_params` section for tissue-specific configuration:

```json
{
    "_gating_params": {
        "root_hard_requirements": {
            "Immune Cells": {
                "marker": "CD45",
                "min_pos_frac": 0.30,
                "min_enrichment": 0.0
            },
            "Epithelium": {
                "marker": "EPCAM",
                "min_pos_frac": 0.40,
                "min_enrichment": 0.0
            }
        },
        "root_veto_markers": {
            "Immune Cells": {
                "markers": ["Pan-Cytokeratin", "E-cadherin", "EPCAM"],
                "max_pos_frac": 0.20
            },
            "Epithelium": {
                "markers": ["CD45", "CD3D", "CD19"],
                "max_pos_frac": 0.15
            }
        },
        "level_defaults": {
            "0": {"min_coverage": 0.45},
            "1": {"min_coverage": 0.35}
        }
    },
    "Immune Cells": { ... },
    "Epithelium": { ... }
}
```

### Example: Fallopian Tube Configuration

The fallopian tube presents unique challenges for cell type annotation:

1. **Mixed epithelial populations**: Ciliated and secretory cells intermingle
2. **Stromal complexity**: Multiple fibroblast subtypes with overlapping markers
3. **Immune infiltration**: Variable immune cell populations
4. **Hormonal influence**: Expression patterns vary with menstrual cycle

#### Recommended Configuration

```json
{
    "_gating_params": {
        "root_hard_requirements": {
            "Epithelium": {
                "marker": "EPCAM",
                "min_pos_frac": 0.50,
                "min_enrichment": 0.0
            },
            "Immune Cells": {
                "marker": "PTPRC",
                "min_pos_frac": 0.35,
                "min_enrichment": 0.0
            }
        },
        "root_veto_markers": {
            "Epithelium": {
                "markers": ["PTPRC", "VIM"],
                "max_pos_frac": 0.25
            },
            "Fibroblasts": {
                "markers": ["EPCAM", "PTPRC"],
                "max_pos_frac": 0.15
            }
        },
        "level_defaults": {
            "0": {
                "min_coverage": 0.45,
                "min_pos_frac": 0.30,
                "min_gap": 0.35
            },
            "1": {
                "min_coverage": 0.35,
                "min_pos_frac": 0.25,
                "min_gap": 0.25
            }
        }
    },
    "Epithelium": {
        "markers": ["EPCAM", "KRT8", "KRT18", "CDH1"],
        "children": {
            "Ciliated Epithelium": {
                "markers": ["FOXJ1", "TUBB4B", "DNAH5", "CAPS", "PIFO"],
                "anti_markers": ["MUC5B", "SCGB2A1", "OVGP1"]
            },
            "Secretory Epithelium": {
                "markers": ["PAX8", "MUC16", "OVGP1", "SCGB2A1"],
                "anti_markers": ["FOXJ1", "TUBB4B", "CAPS"]
            }
        }
    }
}
```

#### Why These Settings Work

- **EPCAM hard requirement for Epithelium**: Ensures epithelial clusters truly express this definitive marker
- **PTPRC (CD45) for Immune**: Standard pan-immune marker
- **VIM veto for Epithelium**: Prevents stromal contamination (though some epithelial-mesenchymal transition cells express VIM)
- **Lower coverage at level 1**: Allows ciliated/secretory distinction even with partial marker panels
- **Anti-markers for ciliated/secretory**: Provides mutual exclusion between these related types

### Example: Lung Configuration

Lung tissue has distinct challenges:

1. **Diverse epithelium**: AT1, AT2, club, ciliated, basal cells
2. **Immune complexity**: Alveolar macrophages, interstitial macrophages
3. **Endothelial diversity**: Arterial, venous, capillary subtypes

```json
{
    "_gating_params": {
        "root_hard_requirements": {
            "Epithelium": {
                "marker": "EPCAM",
                "min_pos_frac": 0.40,
                "min_enrichment": 0.0
            },
            "Endothelium": {
                "marker": "PECAM1",
                "min_pos_frac": 0.50,
                "min_enrichment": 0.0
            }
        },
        "root_veto_markers": {
            "Epithelium": {
                "markers": ["PTPRC", "PECAM1"],
                "max_pos_frac": 0.15
            },
            "Macrophages": {
                "markers": ["EPCAM", "CD3D"],
                "max_pos_frac": 0.10
            }
        }
    }
}
```

### Example: Tumor Microenvironment

Tumor samples present additional challenges:

1. **Aberrant expression**: Cancer cells may express unexpected markers
2. **Mixed populations**: Tumor heterogeneity creates complex clusters
3. **Immune exhaustion**: Immune markers may be downregulated

```json
{
    "_gating_params": {
        "root_hard_requirements": {
            "T Cells": {
                "marker": "CD3D",
                "min_pos_frac": 0.25,
                "min_enrichment": 0.0
            }
        },
        "level_defaults": {
            "0": {
                "min_coverage": 0.35,
                "min_pos_frac": 0.20
            }
        }
    }
}
```

---

## Iterative Tuning Workflow

### Step 1: Initial Run with Defaults

```bash
ctr annotate input.h5ad --output-dir run_v1/
```

### Step 2: Evaluate Results

```python
import pandas as pd

# Load results
results = pd.read_csv("run_v1/cluster_annotations.csv")

# Check assignment distribution
print(results['final_annotation'].value_counts())

# Check stop reasons
print(results['stop_reason'].value_counts())

# Identify problematic clusters
unassigned = results[results['final_annotation'] == 'Unassigned']
print(f"Unassigned clusters: {len(unassigned)}")
```

### Step 3: Diagnose Issues

```python
# For unassigned clusters, check coverage scores
details = pd.read_csv("run_v1/gating_details.csv")
for cluster in unassigned['cluster']:
    cluster_details = details[details['cluster'] == cluster]
    print(f"\n{cluster}:")
    print(cluster_details[['root_type', 'coverage', 'pos_frac', 'score']])
```

### Step 4: Apply Targeted Fixes

Based on diagnosis, apply appropriate fixes (see scenarios above).

### Step 5: Re-run and Compare

```bash
# Use annotation-only to skip re-clustering
ctr annotate run_v1/output.h5ad \
    --annotation-only \
    --gating-params '{"level_defaults": {"0": {"min_coverage": 0.35}}}' \
    --output-dir run_v2/
```

### Step 6: Validate Improvements

```python
# Compare versions
v1 = pd.read_csv("run_v1/cluster_annotations.csv")
v2 = pd.read_csv("run_v2/cluster_annotations.csv")

# Check if unassigned reduced
print(f"V1 unassigned: {(v1['final_annotation'] == 'Unassigned').sum()}")
print(f"V2 unassigned: {(v2['final_annotation'] == 'Unassigned').sum()}")

# Check for new issues
changed = v1.merge(v2, on='cluster', suffixes=('_v1', '_v2'))
changed = changed[changed['final_annotation_v1'] != changed['final_annotation_v2']]
print(f"Changed assignments: {len(changed)}")
```

---

## Common Parameter Combinations

### Sparse Data (Spatial Transcriptomics)

```bash
ctr annotate ... \
    --scoring-params '{"positive_quantile": 0.50}' \
    --gating-params '{
        "level_defaults": {
            "0": {"min_coverage": 0.30, "min_pos_frac": 0.15},
            "1": {"min_coverage": 0.20, "min_pos_frac": 0.10}
        }
    }'
```

### High-Quality Deep scRNA-seq

```bash
ctr annotate ... \
    --scoring-params '{"positive_quantile": 0.80, "de_bonus": 0.7}' \
    --gating-params '{
        "level_defaults": {
            "0": {"min_coverage": 0.55, "min_gap": 0.45},
            "1": {"min_coverage": 0.45, "min_gap": 0.35}
        }
    }'
```

### Heterogeneous Tumor Samples

```bash
ctr annotate ... \
    --scoring-params '{"anti_weight": 0.8}' \
    --gating-params '{
        "level_defaults": {
            "0": {"min_coverage": 0.35, "min_pos_frac": 0.20}
        }
    }' \
    --refinement-params '{"score_threshold": 0.8, "min_cells": 300}'
```

### Conservative Annotation (High Precision)

```bash
ctr annotate ... \
    --scoring-params '{"anti_weight": 0.9, "anti_agg": "max"}' \
    --gating-params '{
        "level_defaults": {
            "0": {"min_coverage": 0.55, "min_gap": 0.50}
        }
    }' \
    --refinement-params '{"score_threshold": 1.5}'
```

---

## See Also

- [Marker Scoring Algorithm](./marker-scoring-algorithm.md) - Understand scoring components
- [Hierarchical Gating Algorithm](./hierarchical-gating-algorithm.md) - Understand gating thresholds
- [Refinement Decision Logic](./refinement-decision-logic.md) - Understand heterogeneity detection
- [Annotation Pipeline](./annotation-pipeline.md) - End-to-end workflow understanding
