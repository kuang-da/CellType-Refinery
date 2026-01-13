---
sidebar_position: 1
---

# Glossary

Quick reference for terminology used in CellType-Refinery documentation.

---

## A

### Anti-Marker

A marker that should NOT be expressed in a cell type. Used to penalize incorrect assignments. For example, CD45 is an anti-marker for epithelial cells.

### Anti-Penalty

The score deduction applied when a cluster expresses anti-markers for a cell type. Higher anti-marker expression results in larger penalties.

### Ambiguous Root

A stop reason indicating the top two root categories scored too close together (gap < root_gap_threshold). The cluster cannot be confidently assigned to a single root.

### Ambiguous Siblings

A stop reason indicating children at the current hierarchy level scored too close together (gap < min_gap[level]). Assignment stops at the parent level.

### Assignment Path

The sequence of nodes traversed during hierarchical descent, from root to the final assigned cell type.

---

## B

### Base Score

The fundamental score component derived from marker expression, before applying DE bonuses or penalties.

### Batch Effect

Technical variation between experimental batches that can affect marker expression measurements and clustering results.

---

## C

### Cell State

A transient or activation-dependent phenotype of a cell, such as "activated," "exhausted," or "proliferating." Distinguished from cell type identity.

### Cell Type

A distinct biological category of cells defined by characteristic marker expression patterns and functional properties.

### Cluster

A group of cells identified by unsupervised clustering algorithms based on expression similarity.

### Commonness Penalty

A score reduction applied to markers that are frequently expressed across many cell types, reducing their discriminative power.

### Confidence

The minimum margin (score gap) along the entire hierarchical assignment path. Higher confidence indicates more decisive assignments. Range: 0 to ~5+.

### Coverage

The fraction of a cell type's expected markers that are present in the marker panel. Coverage of 0.8 means 80% of markers are available for scoring.

---

## D

### DE Bonus

Extra score contribution for markers that appear in the cluster's top differentially expressed genes. Rank-weighted to give higher weight to top-ranked DE genes.

### DE Component

The differential expression portion of the marker score. Calculated from rank-weighted DE hits with commonness penalty.

### DE Hits

The count of a cell type's markers that appear in the cluster's differentially expressed gene list.

### Differential Expression (DE)

Statistical comparison identifying genes with significantly different expression between a cluster and other cells.

---

## E

### Enrichment

Z-score of a cluster's marker expression compared to the global mean. Positive enrichment indicates the cluster expresses the marker more than average.

### Enrichment Threshold

The minimum Z-score required for a marker to be considered enriched in a cluster.

### Expression Matrix

The data structure containing gene expression values for each cell, with genes as rows and cells as columns.

---

## F

### Final Score

The composite score combining base score, DE component, and anti-penalty used for ranking cell type candidates.

### Frac Markers On

The fraction of a cell type's markers that are expressed above threshold in the cluster.

---

## G

### Gap Threshold

The minimum score difference required between candidates for confident assignment. Different thresholds may apply at different hierarchy levels.

### Gating

The process of checking whether a cluster meets threshold requirements to be considered for a cell type. Gates include coverage, positive fraction, enrichment, and anti-penalty checks.

### Gap

The score difference between the best candidate and the runner-up at any hierarchy level. Larger gaps indicate more confident assignments.

### Global Mean

The average expression of a marker across all cells in the dataset, used as reference for enrichment calculations.

---

## H

### Hard Requirements

Mandatory marker expression criteria that must be met for a cell type to be considered. Failure to meet hard requirements results in immediate exclusion.

### Hierarchical Descent

The top-down traversal through the marker hierarchy, starting from root categories and descending to specific cell types.

### Hierarchy Level

The depth in the marker map tree structure. Level 0 is root, with increasing levels representing more specific cell type categories.

---

## I

### Intermediate Node

A non-leaf node in the marker hierarchy that has children. Assignment may stop at intermediate nodes when children are ambiguous.

---

## L

### Leaf Node

A terminal node in the marker hierarchy with no children. Represents the most specific cell type category available.

### Leaf Reached

A stop reason indicating the algorithm successfully descended to the deepest level of the hierarchy. This is the ideal outcome.

---

## M

### Marker

A gene whose expression is characteristic of a particular cell type or category.

### Marker Map

The hierarchical data structure defining cell types, their relationships, and associated markers.

### Marker Panel

The set of genes measured in the experiment that can be used for cell type annotation.

### Mean Enrichment

Average Z-score of a cluster's expression for a cell type's markers compared to global expression.

### Mean Positive Fraction

Average fraction of cells in a cluster expressing each of a cell type's markers above the Q75 threshold.

### Min Gap

The minimum score difference required between the best and second-best candidates at each hierarchy level. Below this threshold, assignment stops as "ambiguous."

---

## N

### No Child Passed

A stop reason indicating none of the children at the current hierarchy level passed the gating requirements. Assignment stops at the current (parent) level.

### No Root Passed

A stop reason indicating no root category passed the gating requirements. The cluster is labeled "Unassigned."

### Node

An entry in the marker hierarchy representing a cell type category at any level.

---

## O

### Overclustering

Excessive subdivision of cell populations resulting in biologically redundant clusters.

---

## P

### Parent Node

The immediate ancestor of a node in the marker hierarchy. All siblings share the same parent.

### Per-Cell Voting

A process that calculates cell-level scores when siblings are ambiguous. Provides evidence for refinement but does NOT override cluster-level assignments.

### Positive Fraction

The fraction of cells in a cluster expressing a marker above the Q75 (75th percentile) threshold.

### Positive Fraction Threshold

The minimum fraction of cells required to express a marker for the cluster to be considered positive for that marker.

---

## Q

### Q75 Threshold

The 75th percentile expression value across all cells, used as the cutoff for determining positive marker expression.

---

## R

### Rank Weight

The weight applied to DE genes based on their rank position, giving higher influence to top-ranked genes.

### Refinement

The process of improving cell type assignments through subclustering, merging, or manual adjustment.

### Root Category

A top-level category in the marker hierarchy (e.g., Immune, Epithelial, Stromal). All hierarchical descent begins at root level.

### Root Gap Threshold

The minimum score difference required between the top two root categories for confident root assignment.

### Root Hard Requirements

Mandatory marker expression requirements for root categories. For example, Immune Cells may require CD45 expression above a threshold.

### Root Veto Markers

Markers that, if expressed above threshold, block assignment to a root category. For example, Pan-Cytokeratin expression may veto Immune Cells assignment.

### Runner-Up

The second-highest scoring candidate at any hierarchy level. Used to calculate the gap for confidence assessment.

---

## S

### Score

The numerical value representing how well a cluster matches a cell type based on marker expression patterns.

### Sibling Comparison

The principle that only cell types sharing the same parent are compared at each hierarchy level. This ensures fair comparison within biological categories.

### Sibling Node

A node sharing the same parent as another node in the marker hierarchy.

### Stop Reason

The explanation for why hierarchical descent stopped at a particular level. Values: leaf_reached, ambiguous_siblings, no_child_passed, no_root_passed, ambiguous_root.

### Subcluster

A refinement operation that re-clusters a cluster at finer resolution to separate mixed populations.

---

## T

### Technical Marker

A marker used for quality control or batch correction rather than biological cell type identification.

### Top-Down Traversal

The annotation strategy of starting at broad root categories and descending through the hierarchy to specific cell types, rather than comparing all cell types directly.

### Total Score

Synonym for final score; the aggregate score used for ranking candidates.

---

## U

### Unassigned

The label given to clusters that fail to pass gating requirements for any root category.

### Underclustering

Insufficient subdivision of cell populations resulting in mixed clusters containing multiple cell types.

---

## V

### Veto

The blocking of a cell type assignment due to expression of exclusionary markers.

### Voting Fraction

In per-cell voting, the proportion of cells within a cluster that vote for a particular cell type.

---

## W

### Weight

A multiplier applied to marker scores to adjust their influence on final scoring. Some markers may be weighted higher due to specificity.

---

## Z

### Z-Score

A statistical measure indicating how many standard deviations a value is from the mean. Used for enrichment calculations.

---

## See Also

- [Marker Scoring Algorithm](../methodology/marker-scoring-algorithm.md) - Detailed explanation of how scores are calculated
- [Hierarchical Gating Algorithm](../methodology/hierarchical-gating-algorithm.md) - How the algorithm traverses the marker map
- [Output Files](./output-files.md) - Output file format specifications
