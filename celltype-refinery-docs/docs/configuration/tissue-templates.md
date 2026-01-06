---
sidebar_position: 2
---

# Tissue Templates

YAML templates for tissue-specific configuration.

```mermaid
flowchart TD
    subgraph Template["Tissue Template Structure"]
        subgraph Gating["Gating Parameters"]
            G1[score_threshold]
            G2[ambiguity_gap]
            G3[min_coverage]
        end

        subgraph Patterns["Pattern Matching"]
            P1[Expected proportions]
            P2[Regional patterns]
            P3[Known artifacts]
        end

        subgraph Biology["Biology Metrics"]
            B1[E:S ratio range]
            B2[Regional gradients]
            B3[Cell type lists]
        end
    end

    style Gating fill:#FFD700
    style Patterns fill:#ADD8E6
    style Biology fill:#90EE90
```

## Location

`configs/tissues/template.yaml`

## Sections

- Gating parameters
- Pattern matching
- Regional configuration
- Biology metrics
