---
sidebar_position: 1
---

# Marker Maps

JSON files defining cell-type hierarchy and markers.

```mermaid
flowchart TD
    subgraph JSON["markers.json Structure"]
        META["_marker_map_metadata"]
        ROOT1["Epithelium"]
        ROOT2["Immune Cells"]
        ROOT3["Mesenchymal"]

        SUB1["Ciliated"]
        SUB2["Secretory"]
        SUB3["T Cells"]
        SUB4["Macrophages"]
    end

    META -.-> ROOT1
    META -.-> ROOT2
    META -.-> ROOT3

    ROOT1 --> SUB1
    ROOT1 --> SUB2
    ROOT2 --> SUB3
    ROOT2 --> SUB4

    style META fill:#E6E6FA
    style ROOT1 fill:#F0FFF0
    style ROOT2 fill:#F0FFF0
    style ROOT3 fill:#F0FFF0
```

## Format

```json
{
  "_marker_map_metadata": {
    "version": "1.0",
    "tissue": "my_tissue",
    "description": "Description"
  },
  "Cell_Type": {
    "markers": ["Marker1", "Marker2"],
    "anti_markers": ["AntiMarker"],
    "subtypes": {}
  }
}
```

## Fields

See [Annotation: Marker Maps](../modules/annotation/marker-maps) for details.
