---
sidebar_position: 2
---

# Orphan Rescue

Rescue cells with strong subtype markers but weak root markers.

```mermaid
flowchart TD
    A[Unassigned Cell] --> B{Strong subtype<br/>markers?}
    B -->|No| STAY[Remains<br/>Unassigned]
    B -->|Yes| C{Weak root<br/>markers?}

    C -->|No| STAY
    C -->|Yes| D{Clear lineage<br/>evidence?}

    D -->|No| STAY
    D -->|Yes| RESCUE[Rescue as<br/>"Subtype (orphan)"]

    style STAY fill:#FFB6C1
    style RESCUE fill:#90EE90
```

**Example**: Cell with strong Podoplanin (lymphatic marker) but weak CD31 (endothelium root) → Rescue as "Lymphatic Endothelium (orphan)"

## Detection

Identify cells that:
- Have strong subtype scores
- Failed parent gating
- Show clear lineage evidence
