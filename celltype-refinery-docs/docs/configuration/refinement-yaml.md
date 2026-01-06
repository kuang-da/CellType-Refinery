---
sidebar_position: 4
---

# Refinement YAML

Configuration for manual refinement.

```mermaid
flowchart TD
    subgraph Config["refinement.yaml"]
        V[version: 1.0]
        I[iteration: 1]
    end

    subgraph Operations["Operations"]
        subgraph OV["overrides"]
            O1[cluster_id: 3]
            O2[cell_type: X]
        end

        subgraph MR["merge"]
            M1[source: 5, 6]
            M2[target: Combined]
        end

        subgraph SC["subcluster"]
            S1[cluster_id: 12]
            S2[resolution: 0.3]
        end
    end

    V --> Operations
    I --> Operations

    style OV fill:#FFD700
    style MR fill:#ADD8E6
    style SC fill:#90EE90
```

## Example

```yaml
version: "1.0"
iteration: 1

overrides:
  - cluster_id: "3"
    cell_type: "My_Cell_Type"
    reason: "Explanation"

merge:
  - source_clusters: ["5", "6"]
    target_label: "Combined"
    reason: "Explanation"

subcluster:
  - cluster_id: "12"
    resolution: 0.3
    reason: "Explanation"
```
