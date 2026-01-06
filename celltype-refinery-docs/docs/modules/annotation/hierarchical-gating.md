---
sidebar_position: 3
---

# Hierarchical Gating

Top-down traversal through the marker hierarchy.

```mermaid
flowchart TD
    START[Start: Cluster] --> ROOT{Any ROOT<br/>passes?}
    ROOT -->|No| UNAS[Unassigned]
    ROOT -->|Yes| BEST[Select Best<br/>ROOT]

    BEST --> CHILD{Has children?}
    CHILD -->|No| LEAF[Leaf Reached]
    CHILD -->|Yes| CHECK{Any child<br/>passes?}

    CHECK -->|No| STOP1[No Child Passed]
    CHECK -->|Yes| SIBS{Ambiguous<br/>siblings?}

    SIBS -->|Yes| STOP2[Ambiguous<br/>Siblings]
    SIBS -->|No| NEXT[Select Best<br/>Child]
    NEXT --> CHILD

    style UNAS fill:#FFB6C1
    style LEAF fill:#90EE90
    style STOP1 fill:#FFFACD
    style STOP2 fill:#FFFACD
```

## Algorithm

1. Find best passing ROOT (score > threshold)
2. Descend: compare SIBLINGS at each level
3. Stop when: leaf reached, no child passes, or ambiguous

## Stop Reasons

| Reason | Description |
|--------|-------------|
| `leaf_reached` | Deepest annotation |
| `ambiguous_siblings` | Close scores |
| `no_child_passed` | Stopped mid-hierarchy |
| `no_root_passed` | Unassigned |

## Scoring Formula

```
score = mean_enrichment + mean_positive + de_bonus - anti_penalty
```
