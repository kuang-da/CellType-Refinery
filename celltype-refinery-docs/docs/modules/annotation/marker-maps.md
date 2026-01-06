---
sidebar_position: 2
---

# Marker Maps

JSON files defining the cell-type hierarchy and markers.

```mermaid
flowchart TD
    subgraph MarkerMap["Marker Map Structure"]
        META["_marker_map_metadata<br/>version, tissue"]

        subgraph ROOT["ROOT Cell Types"]
            EPI["Epithelium<br/>markers: [EPCAM, KRT]"]
            IMM["Immune Cells<br/>markers: [CD45]"]
            MES["Mesenchymal<br/>markers: [VIM]"]
        end

        subgraph Subtypes["Nested Subtypes"]
            CIL["Ciliated<br/>markers: [FOXJ1]"]
            SEC["Secretory<br/>markers: [MUC5B]"]
        end
    end

    META -.-> ROOT
    EPI --> CIL
    EPI --> SEC

    style META fill:#E6E6FA
    style ROOT fill:#F0FFF0
    style Subtypes fill:#FFFACD
```

## Format

```json
{
  "_marker_map_metadata": {
    "version": "1.0",
    "tissue": "my_tissue"
  },
  "Cell_Type": {
    "markers": ["Marker1", "Marker2"],
    "anti_markers": ["AntiMarker"],
    "subtypes": {
      "Subtype_A": {
        "markers": ["SubMarker1"]
      }
    }
  }
}
```

## Fields

| Field | Required | Description |
|-------|----------|-------------|
| `markers` | Yes | Positive markers |
| `anti_markers` | No | Negative markers |
| `subtypes` | No | Child cell types |
