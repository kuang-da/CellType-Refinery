---
sidebar_position: 1
---

# Refinement Overview

The refinement module improves cell type annotations through automatic and manual policies. It provides a flexible framework for correcting, merging, splitting, and relabeling clusters based on statistical evidence and expert knowledge.

:::tip Deep Dive
For the decision logic behind refinement, see [Refinement Decision Logic](/docs/methodology/refinement-decision-logic.md).
:::

## Architecture

The refinement system uses a two-policy architecture that allows automatic recommendations to be combined with manual overrides:

```text
┌─────────────────────────────────────────────────────────────────────────────┐
│                        REFINEMENT ARCHITECTURE                               │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│    ┌─────────────────┐         ┌─────────────────┐                          │
│    │   AutoPolicy    │         │  ManualPolicy   │                          │
│    │  (--auto flag)  │         │ (--config YAML) │                          │
│    └────────┬────────┘         └────────┬────────┘                          │
│             │                           │                                   │
│             ▼                           ▼                                   │
│        base_plan                   overlay_plan                             │
│             │                           │                                   │
│             └───────────┬───────────────┘                                   │
│                         ▼                                                   │
│              ┌─────────────────────┐                                        │
│              │    merge_plans()    │                                        │
│              │  (overlay wins on   │                                        │
│              │   conflicts)        │                                        │
│              └──────────┬──────────┘                                        │
│                         ▼                                                   │
│                    merged_plan                                              │
│                         │                                                   │
│                         ▼                                                   │
│              ┌─────────────────────┐                                        │
│              │  RefinementEngine   │                                        │
│              │  (execute plan)     │                                        │
│              └──────────┬──────────┘                                        │
│                         ▼                                                   │
│                   refined.h5ad                                              │
└─────────────────────────────────────────────────────────────────────────────┘
```

### Key Components

- **AutoPolicy**: Analyzes scoring metrics and automatically generates refinement recommendations based on configurable thresholds
- **ManualPolicy**: Reads user-defined YAML configuration files specifying exact operations to perform
- **merge_plans()**: Combines base and overlay plans, with manual overrides taking precedence on conflicts
- **RefinementEngine**: Executes the merged plan and produces the refined output

## Execution Modes

The refinement module supports two primary execution modes that control whether changes are applied to the data:

### Diagnostic Mode (Default)

When run without the `--execute` flag, the refinement module operates in diagnostic mode:

- Analyzes the input data and generates recommendations
- Produces a `diagnostic_report.csv` with suggested refinements
- **Does not modify the input data**
- Useful for reviewing proposed changes before committing

```bash
# Diagnostic mode - review recommendations first
celltype-refinery refine input.h5ad --auto
```

### Execution Mode

When the `--execute` flag is provided, refinements are applied:

- Executes all planned operations on the data
- Outputs a new `refined.h5ad` file with updated annotations
- Original input file remains unchanged
- Generates an execution log for audit purposes

```bash
# Execution mode - apply refinements
celltype-refinery refine input.h5ad --auto --execute
```

### Policy Combinations

The refinement module supports three policy configurations:

#### Auto-only Mode

Uses only the AutoPolicy to generate and apply refinements based on statistical analysis:

```bash
celltype-refinery refine input.h5ad --auto --execute
```

Best for:
- Initial annotation passes
- Large datasets requiring automated curation
- When manual review is not feasible

#### Manual-only Mode

Uses only the ManualPolicy with a user-defined YAML configuration:

```bash
celltype-refinery refine input.h5ad --config curation.yaml --execute
```

Best for:
- Expert-driven curation workflows
- Applying known corrections from literature
- Reproducible, version-controlled refinements

#### Hybrid Mode

Combines AutoPolicy recommendations with manual overrides. Manual configurations take precedence when conflicts occur:

```bash
celltype-refinery refine input.h5ad --auto --config overrides.yaml --execute
```

Best for:
- Leveraging automation while maintaining expert control
- Correcting specific auto-generated recommendations
- Iterative refinement workflows

## Operation Types

The refinement module supports five operation types, each serving a distinct purpose:

| Operation | Description | Source |
|-----------|-------------|--------|
| Override | Direct label assignment | ManualPolicy |
| Merge | Combine similar clusters | ManualPolicy |
| Subcluster | Re-cluster at finer resolution | AutoPolicy or ManualPolicy |
| Relabel | Update label without re-clustering | AutoPolicy |
| Rescore | Recompute scores | AutoPolicy (automatic) |

### Operation Details

**Override**: Directly assigns a new cell type label to a cluster, bypassing all scoring logic. Used when expert knowledge contradicts automated predictions.

**Merge**: Combines two or more clusters into a single cluster. Typically used when clusters represent the same cell type but were over-split during initial clustering.

**Subcluster**: Re-clusters cells within a cluster at a finer resolution. Used when a cluster contains heterogeneous cell populations that should be separated.

**Relabel**: Updates the assigned label based on re-evaluation of scores without modifying cluster membership. Used when the top-scoring label changes after threshold adjustments.

**Rescore**: Recomputes confidence scores for all clusters. Automatically triggered after structural changes (merge, subcluster) to ensure score consistency.

### Execution Order

Operations are executed in a specific order to ensure consistency and avoid conflicts:

```
override → merge → subcluster → relabel → rescore
```

This ordering ensures that:
1. Direct overrides are applied first, preventing unnecessary computation
2. Merges consolidate clusters before subclustering decisions
3. Subclustering creates new clusters that can then be relabeled
4. Relabeling uses the final cluster structure
5. Rescoring reflects all structural changes

## AutoPolicy Selection Criteria

The AutoPolicy uses a rules-based system to determine which operations to recommend for each cluster. Key decision factors include:

- **Confidence Score**: Primary metric for label reliability
- **Delta Score**: Difference between top two candidate labels
- **Cluster Size**: Number of cells in the cluster
- **Marker Expression**: Presence of canonical markers for predicted types
- **Entropy**: Distribution of cells across candidate labels

### Decision Thresholds

| Metric | Low | Medium | High |
|--------|-----|--------|------|
| Confidence | < 0.3 | 0.3 - 0.7 | > 0.7 |
| Delta | < 0.1 | 0.1 - 0.3 | > 0.3 |

For detailed decision trees and threshold tuning, see the [Tuning Guide](/docs/methodology/tuning-guide.md) and [Refinement Decision Logic](/docs/methodology/refinement-decision-logic.md) documentation.

## CLI Examples

### Basic Usage

```bash
# View help for refine command
celltype-refinery refine --help

# Diagnostic run with auto policy
celltype-refinery refine input.h5ad --auto

# Execute auto refinements
celltype-refinery refine input.h5ad --auto --execute
```

### Manual Configuration

```bash
# Diagnostic run with manual config
celltype-refinery refine input.h5ad --config curation.yaml

# Execute manual refinements
celltype-refinery refine input.h5ad --config curation.yaml --execute

# Hybrid mode with overrides
celltype-refinery refine input.h5ad --auto --config overrides.yaml --execute
```

### Output Control

```bash
# Specify output path
celltype-refinery refine input.h5ad --auto --execute --output refined_output.h5ad

# Generate detailed diagnostic report
celltype-refinery refine input.h5ad --auto --report detailed_report.csv

# Verbose logging
celltype-refinery refine input.h5ad --auto --execute --verbose
```

### Advanced Options

```bash
# Custom confidence threshold
celltype-refinery refine input.h5ad --auto --execute --min-confidence 0.5

# Limit operations to specific clusters
celltype-refinery refine input.h5ad --auto --execute --clusters 0,1,5,12

# Dry run to preview execution plan
celltype-refinery refine input.h5ad --auto --dry-run
```

## See Also

- [Refinement Decision Logic](/docs/methodology/refinement-decision-logic.md) - Detailed explanation of the decision framework
- [Policies](./policies.md) - Policy types and configuration options
- [Operations](./operations.md) - Available refinement operations
- [Tuning Guide](/docs/methodology/tuning-guide.md) - Practical tips for optimizing refinement parameters
- [Marker Scoring](/docs/modules/annotation/marker-scoring.md) - Understanding the scores used in refinement decisions
