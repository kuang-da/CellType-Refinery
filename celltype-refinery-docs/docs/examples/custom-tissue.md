---
sidebar_position: 2
---

# Custom Tissue Example

How to configure CellType-Refinery for a new tissue.

```mermaid
flowchart TD
    subgraph Setup["Initial Setup"]
        S1["1. Create marker<br/>map JSON"] --> S2["2. Create tissue<br/>template YAML"]
    end

    subgraph Iterate["Iterative Development"]
        S3["3. Run initial<br/>clustering"] --> S4["4. Review<br/>results"]
        S4 --> S5{"Satisfied?"}
        S5 -->|No| S6["Refine marker<br/>definitions"]
        S6 --> S3
        S5 -->|Yes| S7["5. Validate with<br/>domain experts"]
    end

    S2 --> S3
    S7 --> DONE[Production<br/>Ready]

    style DONE fill:#90EE90
```

## Steps

1. Create marker map JSON
2. Create tissue template YAML
3. Run initial clustering
4. Iterate on marker definitions
5. Validate with domain experts

## Tips

- Start with major lineages
- Add subtypes incrementally
- Use anti-markers to distinguish similar types
- Validate with known marker expression
