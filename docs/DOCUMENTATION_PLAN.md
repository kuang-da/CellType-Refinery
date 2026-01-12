# Documentation Improvement Plan (Revised)

## Executive Summary

This plan improves documentation for Stages H (Clustering & Initial Annotation) and I (Refinement), which together implement the core cell type annotation methodology - the project's key contribution.

**Key Revisions Based on Review:**
- Restructured phases to follow user journey (not algorithm-first)
- Added getting-started enhancements
- Added API reference for developers
- Added explicit cross-references between sections
- Reduced redundancy through link-heavy approach

---

## Goals

1. **Make the annotation methodology clear and accessible** - Users understand HOW cell types are assigned
2. **Enable parameter tuning** - Document scoring formulas and gating thresholds with tuning guidance
3. **Support the user journey** - From "what is this?" to "how do I tune this?"
4. **Serve developers** - API reference, source code links, extension guides
5. **Integrate rich technical content** - Bring stage_h.md ASCII charts and Stage I docs into user docs

---

## Proposed Sidebar Structure

```javascript
const sidebars = {
  docsSidebar: [
    'intro',
    {
      type: 'category',
      label: 'Getting Started',
      items: [
        'getting-started/installation',
        'getting-started/quickstart',
        'getting-started/first-annotation',
        'getting-started/key-concepts',        // NEW: Terminology glossary
        'getting-started/configuration',
      ],
    },
    {
      type: 'category',
      label: 'Methodology',                    // NEW SECTION
      items: [
        'methodology/index',                   // Philosophy & design principles
        'methodology/annotation-pipeline',     // Complete H→I pipeline
        'methodology/marker-scoring-algorithm',// Full formula breakdown
        'methodology/hierarchical-gating-algorithm', // Gating decision tree
        'methodology/refinement-decision-logic',// When to subcluster vs relabel
        'methodology/tuning-guide',            // Parameter adjustment guide
      ],
    },
    {
      type: 'category',
      label: 'Core Workflows',
      items: [/* keep as-is */],
    },
    {
      type: 'category',
      label: 'Modules',
      items: [
        /* Enhanced with cross-references */
      ],
    },
    {
      type: 'category',
      label: 'Reference',                      // NEW SECTION
      items: [
        'reference/glossary',                  // Term definitions
        'reference/faq',                       // Common questions
        'reference/parameter-reference',       // All parameters in one place
        'reference/output-files',              // CSV/H5AD schema reference
      ],
    },
    {
      type: 'category',
      label: 'API Reference',                  // NEW SECTION (for developers)
      items: [
        'api/overview',
        'api/annotation-engine',
        'api/refinement-engine',
        'api/data-structures',
      ],
    },
    {
      type: 'category',
      label: 'CLI Reference',
      items: [/* keep as-is */],
    },
    {
      type: 'category',
      label: 'Configuration',
      items: [/* keep as-is */],
    },
    {
      type: 'category',
      label: 'Examples',
      items: [/* keep as-is */],
    },
  ],
};
```

---

## User Paths After Restructuring

### Path 1: New User (Wants to Run)
```
intro → getting-started/quickstart → core-workflows/automated-workflow
→ examples/fallopian-tube
```

### Path 2: Researcher (Wants to Understand the Algorithm)
```
intro → methodology/index → methodology/annotation-pipeline
→ methodology/marker-scoring-algorithm → methodology/hierarchical-gating-algorithm
→ methodology/tuning-guide (for customization)
```

### Path 3: Developer (Extending the System)
```
intro → methodology/ (understand algorithms) → api/overview
→ api/annotation-engine → modules/ (implementation details)
```

---

## Implementation Phases (User-Journey First)

### Phase 1: Foundation & Quick Wins (Priority: Critical)
Enable users to get started and understand core concepts.

| File | Lines | Description |
|------|-------|-------------|
| `getting-started/key-concepts.md` | 200 | Terminology glossary (gating, enrichment, stop_reason, etc.) |
| `methodology/index.md` | 200 | Design philosophy, why this approach works |
| `reference/glossary.md` | 250 | Comprehensive term definitions |
| `reference/output-files.md` | 200 | Schema for all CSV/H5AD outputs |

### Phase 2: Core Methodology (Priority: High)
Document the cell type annotation algorithm - the key contribution.

| File | Lines | Description |
|------|-------|-------------|
| `methodology/annotation-pipeline.md` | 300 | Complete H→I flow with ASCII diagrams |
| `methodology/marker-scoring-algorithm.md` | 400 | Full formula with all 4 components |
| `methodology/hierarchical-gating-algorithm.md` | 350 | Gate checks, descent, stop reasons |
| `methodology/refinement-decision-logic.md` | 300 | 4 criteria for subcluster vs relabel |

### Phase 3: Practical Support (Priority: High)
Help users tune parameters and understand common scenarios.

| File | Lines | Description |
|------|-------|-------------|
| `methodology/tuning-guide.md` | 500 | Parameter adjustment by scenario (expanded with practical tips) |
| `reference/faq.md` | 250 | Common questions answered |

### Phase 4: Enhanced Module Docs (Priority: Medium)
Improve existing pages with cross-references and detail.

| File | Action | Lines |
|------|--------|-------|
| `modules/annotation/overview.md` | ENHANCE: Add ASCII diagram, links to methodology | +150 |
| `modules/annotation/marker-scoring.md` | ENHANCE: Link to methodology, add quick reference | +100 |
| `modules/annotation/hierarchical-gating.md` | ENHANCE: Link to methodology, add stop reasons table | +100 |
| `modules/annotation/per-cell-voting.md` | NEW: Evidence collection explained | 150 |
| `modules/annotation/design-principles.md` | NEW: Parents-as-gates, sibling comparison | 150 |
| `modules/refinement/overview.md` | ENHANCE: Full architecture diagram | +200 |
| `modules/refinement/autopolicy-criteria.md` | NEW: 4 selection criteria details | 200 |
| `modules/refinement/focus-controls.md` | NEW: --focus-labels/clusters usage | 150 |
| `modules/refinement/smart-rescoring.md` | NEW: Optimization explanation | 150 |
| `modules/refinement/iterative-refinement.md` | NEW: Multi-iteration workflows | 200 |
| `modules/clustering/overview.md` | ENHANCE: Add preprocessing, GPU notes | +150 |

### Phase 5: Developer Documentation (Priority: Medium)
API reference and extension guides.

| File | Lines | Description |
|------|-------|-------------|
| `api/overview.md` | 150 | API introduction, when to use |
| `api/annotation-engine.md` | 300 | AnnotationEngine class, params, methods |
| `api/refinement-engine.md` | 300 | RefinementEngine class, operations |
| `api/data-structures.md` | 250 | Key dataclasses (ScoringContext, MarkerSet, etc.) |

---

## Content Details for Key Files

### `methodology/index.md` - The Core Philosophy

**Purpose**: Explain WHY this approach works for cell type annotation.

**Content Outline**:
1. **Problem Statement**
   - Why traditional clustering + manual annotation fails at scale
   - The challenge of marker ambiguity and overlapping cell types

2. **Key Insight: Hierarchical Marker-Based Gating**
   - Not just clustering - structured traversal through marker hierarchy
   - Sibling-only comparison ensures fair evaluation

3. **Design Principles** (integrate from stage_h.md lines 1004-1046)
   - Parents as gates (not just labels)
   - Sibling-only competition
   - Evidence collection without override
   - Traceable decisions at every step

4. **Visual: High-Level Pipeline**
   - ASCII diagram showing Stage H → Stage I flow

**Cross-references**:
- → `methodology/marker-scoring-algorithm.md` for formula details
- → `methodology/hierarchical-gating-algorithm.md` for decision logic
- → `reference/glossary.md` for terminology

---

### `methodology/marker-scoring-algorithm.md` - The Formula

**Purpose**: Complete technical explanation so users can tune parameters.

**Content Outline**:

```
score = mean_enrichment + mean_positive + de_component - anti_penalty
```

1. **Mean Enrichment** (Z-score calculation)
   ```
   enrichment[m] = (cluster_median[m] - global_median[m]) / global_std[m]
   ```
   - Why Z-scores? Normalizes across markers with different scales
   - Interpretation: >0 = enriched, <0 = depleted

2. **Mean Positive Fraction** (Q75 threshold)
   ```
   positive[m] = fraction of cells >= global_quantile[m, 0.75]
   ```
   - Why Q75? Filters noise, focuses on clearly expressing cells
   - Tunable via `positive_quantile` parameter

3. **DE Component** (rank-weighted bonus)
   ```
   K = clamp(ceil(panel_size × de_top_frac), de_min_k, de_max_k)
   rank_weight = (K - rank + 1) / K
   commonness_penalty = (doc_freq[m] / n_cell_types) ^ de_commonness_alpha
   marker_bonus = rank_weight × (1 - commonness_penalty)
   de_component = de_bonus × mean(marker_bonus)
   ```
   - Why rank-weighted? Continuous bonus, not binary
   - Why commonness penalty? Prevents ubiquitous markers from dominating

4. **Anti-Penalty** (conflict detection)
   ```
   anti_enrichment = max(cluster_median[m] - global_median[m], 0) / std[m]
   anti_positive = fraction of cells >= Q75
   anti_penalty = anti_weight × aggregated(anti_enrichment + anti_positive)
   ```
   - Aggregation modes: `max`, `top2mean` (default), `mean`
   - Hard gate: `anti_penalty > 1.0 → REJECT`

5. **Parameter Reference Table**
   | Parameter | Default | Effect | When to Adjust |
   |-----------|---------|--------|----------------|
   | `positive_quantile` | 0.75 | Higher = stricter positivity | Sparse data → lower to 0.5 |
   | `de_bonus` | 0.5-0.6 | Higher = DE matters more | Good DE results → increase |
   | `de_top_frac` | 0.2 | Fraction of panel for K | Large panel → decrease |
   | `anti_weight` | 0.5-0.8 | Higher = stricter conflict | Wrong assignments → increase |
   | `anti_agg` | top2mean | Aggregation method | Strict → max, lenient → mean |

6. **ASCII Formula Visualization** (from stage_h.md lines 98-109)

**Cross-references**:
- → `methodology/hierarchical-gating-algorithm.md` for how scores are used
- → `methodology/tuning-guide.md` for adjustment scenarios

---

### `methodology/hierarchical-gating-algorithm.md` - The Decision Logic

**Purpose**: Explain exactly how labels are assigned.

**Content Outline**:

1. **Overview: Top-Down Traversal**
   - Start at roots, descend through hierarchy
   - Only siblings compete at each level

2. **Gate Check Sequence** (per level)
   ```
   Level 0 (Roots only):
   ├─ CHECK 1: root_hard_requirements
   │  - marker pos_frac < required → REJECT
   └─ CHECK 2: root_veto_markers
      - veto marker pos_frac > max_allowed → REJECT

   All levels:
   ├─ CHECK 3: coverage >= min_coverage[level]
   ├─ CHECK 4: pos_frac >= min_pos_frac[level] OR frac_markers_on >= threshold
   ├─ CHECK 5: enrichment >= min_enrichment[level]
   └─ CHECK 6: anti_penalty <= 1.0
   ```

3. **Base Gating Parameters Table**
   | Level | min_coverage | min_pos_frac | min_enrichment | min_gap |
   |-------|--------------|--------------|----------------|---------|
   | 0 (Root) | 0.50 | 0.30 | 0.00 | 0.50 |
   | 1 | 0.40 | 0.20 | -0.50 | 0.30 |
   | 2 | 0.30 | 0.15 | -1.00 | 0.20 |
   | 3+ | 0.30 | 0.15 | -1.00 | 0.20 |

4. **Descent Algorithm** (pseudocode with ASCII flowchart)
   ```
   function assign_label(cluster, hierarchy):
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

5. **Stop Reasons Explained**
   | Stop Reason | Meaning | What to Do |
   |-------------|---------|------------|
   | `leaf_reached` | Best case - reached deepest level | None needed |
   | `ambiguous_siblings` | Top 2 children too close | Consider refinement |
   | `no_child_passed` | All children failed gates | Check marker coverage |
   | `no_root_passed` | No root category fit | Check panel, lower thresholds |
   | `ambiguous_root` | Top 2 roots too close | Check root markers |

6. **ASCII Decision Tree** (from stage_h.md lines 550-585)

**Cross-references**:
- → `methodology/marker-scoring-algorithm.md` for score calculation
- → `methodology/refinement-decision-logic.md` for what to do with ambiguous clusters

---

### `methodology/refinement-decision-logic.md` - When to Refine

**Purpose**: Explain AutoPolicy criteria for subclustering vs relabeling.

**Content Outline**:

1. **Overview: The Refinement Question**
   - Stage H produces initial annotations
   - Some clusters are confident, others need refinement
   - Question: SUBCLUSTER (re-cluster) or RELABEL (instant reassignment)?

2. **Criterion 1: Low Confidence**
   ```
   IF score < score_threshold (default: 1.0)
   AND n_cells >= min_cells (default: 500)
   → SUBCLUSTER
   ```
   - Rationale: Large, poorly-scoring clusters likely contain mixed populations
   - Subclustering may reveal hidden subpopulations

3. **Criterion 2a: Homogeneous Parent**
   ```
   IF assigned at parent level (e.g., "Epithelium" not "Ciliated Epithelium")
   AND best_child_score > subtype_signal_threshold (default: 0.5)
   AND (only_one_passing_child OR gap_to_runner_up >= heterogeneity_gap)
   → RELABEL to best child
   ```
   - Rationale: Clear subtype signal, no need to re-cluster
   - Fast: instant relabel without clustering

4. **Criterion 2b: Heterogeneous Parent**
   ```
   IF assigned at parent level
   AND best_child_score > subtype_signal_threshold
   AND multiple_passing_children
   AND gap < heterogeneity_gap (default: 0.3)
   → SUBCLUSTER
   ```
   - Rationale: Competing child signals suggest mixed population
   - Subclustering should separate the populations

5. **Criterion 3: Mixed Population (Extension)**
   ```
   IF parent has high score
   AND all children below subtype_signal_threshold
   → SUBCLUSTER
   ```
   - Rationale: Strong parent signal but weak children = novel subtype or noise

6. **Criterion 4: Weak Leaf (Extension)**
   ```
   IF at leaf level
   AND score low
   AND marker heterogeneity high
   → SUBCLUSTER
   ```
   - Rationale: Leaf with high variance may contain distinct subpopulations

7. **Decision Flowchart** (ASCII)
   ```
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

**Cross-references**:
- → `modules/refinement/autopolicy-criteria.md` for implementation details
- → `modules/refinement/operations.md` for operation types

---

### `methodology/tuning-guide.md` - Practical Parameter Tuning

**Purpose**: Help users adjust parameters for their specific data and tissue types. This is expected - almost every dataset requires some tuning.

**Content Outline**:

1. **Too Many "Unassigned" Clusters**
   - Check: Marker panel coverage (do you have markers for all cell types?)
   - Check: Root requirements (are they too strict?)
   - Fix: Lower `min_coverage` (0.5 → 0.3)
   - Fix: Lower `min_pos_frac` (0.3 → 0.15)
   - Fix: Review `root_hard_requirements` in marker map

2. **Annotation Stops Too Early (Parent-Level Only)**
   - Check: Child marker definitions (are children well-defined?)
   - Check: min_gap thresholds (are they too strict?)
   - Fix: Lower `min_gap` at level 1+ (0.3 → 0.2)
   - Fix: Add more specific markers for child types

3. **Wrong Cell Types Assigned**
   - Check: Anti-markers (are they defined for conflicting types?)
   - Check: Root veto markers (should block incorrect roots)
   - Fix: Add anti-markers to marker map
   - Fix: Increase `anti_weight` (0.5 → 0.8)
   - Fix: Add `root_hard_requirements` for key markers

4. **Refinement Too Aggressive (Too Much Subclustering)**
   - Check: `score_threshold` (is it too high?)
   - Check: `min_cells` (should filter small clusters)
   - Fix: Raise `score_threshold` (1.0 → 1.5)
   - Fix: Raise `min_cells` (500 → 1000)

5. **Results Inconsistent Between Runs**
   - Cause: GPU non-determinism in Leiden clustering
   - Fix: Use `--no-gpu` for reproducibility
   - Alternative: Run clustering once, then use annotation-only mode for iterations

6. **High Confidence but Biologically Wrong**
   - Issue: Markers may not be specific for your tissue
   - Fix: Review marker map for tissue appropriateness
   - Fix: Add tissue-specific `gating_params` section to marker map

7. **Parameter Reference Quick Table**
   - All tunable parameters organized by category
   - Default values, valid ranges, and when to adjust

8. **Tissue-Specific Configuration**
   - How to create custom `_gating_params` in marker map JSON
   - Example: Fallopian tube defaults explained
   - Template for new tissues

**Cross-references**:
- → `methodology/marker-scoring-algorithm.md` for formula details
- → `methodology/hierarchical-gating-algorithm.md` for gating thresholds
- → `reference/parameter-reference.md` for complete parameter list

---

### `reference/glossary.md` - Terminology

Key terms to define:
- **Enrichment**: Z-score of cluster expression vs global mean
- **Positive Fraction**: % of cells expressing above Q75 threshold
- **Coverage**: Fraction of expected markers present in panel
- **Gating**: Sequential filters where cluster must meet threshold to proceed
- **Hierarchical Descent**: Start at root, drill down to specific cell type
- **Sibling Comparison**: Only compare cell types that share a parent
- **Stop Reason**: Why the algorithm stopped descending
- **Confidence**: Minimum margin (score gap) along decision path
- **Anti-Marker**: Marker that should NOT be expressed in a cell type
- **DE Bonus**: Extra score for markers in top differential expression results
- **Commonness Penalty**: Reduction for markers expressed in many cell types

---

## ASCII Diagrams to Integrate

From `docs/stage_h.md`, integrate:

1. **Main Pipeline Flow** (lines 47-192) → `methodology/annotation-pipeline.md`
2. **Marker Hierarchy Example** (lines 116-148) → `configuration/marker-maps.md`
3. **Clustering Pipeline** (lines 269-318) → `modules/clustering/overview.md`
4. **Scoring Formula Breakdown** (lines 98-109) → `methodology/marker-scoring-algorithm.md`
5. **Hierarchical Assignment Algorithm** (lines 550-585) → `methodology/hierarchical-gating-algorithm.md`

From Stage I documentation (user-provided), integrate:
1. **Policy Architecture** → `modules/refinement/overview.md`
2. **Execution Modes** → `modules/refinement/execution-mode.md`
3. **Candidate Selection Criteria** → `methodology/refinement-decision-logic.md`
4. **Operation Types** → `modules/refinement/operations.md`
5. **Focus Controls** → `modules/refinement/focus-controls.md`

---

## Cross-Reference Strategy

Every page should include:

1. **"See Also" Section** at bottom with related pages
2. **Inline Links** where concepts are mentioned
3. **"Prerequisites" Note** at top if page depends on understanding another

Example for `methodology/marker-scoring-algorithm.md`:
```markdown
---
sidebar_position: 3
---

# Marker Scoring Algorithm

:::info Prerequisites
This page explains scoring details. For conceptual background, see
[Methodology Overview](./index.md) first.
:::

[... content ...]

## See Also

- [Hierarchical Gating Algorithm](./hierarchical-gating-algorithm.md) - How scores determine labels
- [Tuning Guide](./tuning-guide.md) - Adjusting scoring parameters
- [Parameter Reference](../reference/parameter-reference.md) - All parameters in one place
- [Glossary: Enrichment](../reference/glossary.md#enrichment) - Term definition
```

---

## Success Criteria

After implementation, users should be able to:

1. **New User**: Install, run quickstart, get first results in <30 minutes
2. **Understand Methodology**: Read methodology/ section and understand HOW annotation works
3. **Tune Parameters**: Know which parameters affect what behavior and adjust confidently
4. **Interpret Results**: Understand stop_reason, confidence, and decision_steps
5. **Developer**: Find API reference, understand data structures, extend the system

---

## Estimated Effort

| Phase | New Content | Enhanced Content | Total Lines |
|-------|-------------|------------------|-------------|
| Phase 1 | 650 lines | 0 | ~650 |
| Phase 2 | 1,350 lines | 0 | ~1,350 |
| Phase 3 | 750 lines | 0 | ~750 |
| Phase 4 | 850 lines | 600 | ~1,450 |
| Phase 5 | 1,000 lines | 0 | ~1,000 |
| **Total** | **4,600 lines** | **600 lines** | **~5,200 lines** |

---

## Notes

### Format Decision: ASCII vs Mermaid
- **ASCII**: Used for complex algorithm flowcharts, formula breakdowns, multi-level decisions
- **Mermaid**: Used for simple high-level overviews (4-5 nodes max)
- Rationale: ASCII renders consistently everywhere, better for technical detail

### Link-Heavy Approach
- Avoid duplicating content across pages
- Each concept documented in ONE place, linked from others
- Example: Scoring formula in `methodology/marker-scoring-algorithm.md`, referenced from `modules/annotation/marker-scoring.md`

### Consolidation of Existing Redundancies
- `configuration/marker-maps.md` and `modules/annotation/marker-maps.md` → Consolidate to annotation/, link from configuration/
- Multiple mentions of scoring without detail → All point to methodology/marker-scoring-algorithm.md
