# Feature Requests

## FEAT-001: Stage F Marker Map Validation

**Status**: Planned
**Priority**: Medium
**Created**: 2025-12-29

### Problem

Stage F sanitizes column names to be AnnData-compatible, replacing characters like `/` and `\` with `_`. For example:
- `Keratin 8/18` becomes `Keratin 8_18`

However, downstream marker maps (used in Stage H annotation) may still reference the original unsanitized names. This causes silent marker resolution failures during annotation, resulting in reduced scoring accuracy.

### Example

```
# Stage F sanitization (stage_f_data_merge.py:38)
_COLUMN_SANITIZE_MAP: Dict[str, str] = {"/": "_", "\\": "_"}

# Marker map (FT_cell_type_markers_v9.json) - BEFORE fix
"Glandular Epithelium": {
    "markers": ["HLA-DR", "EPCAM", "CD105", "PAX8", "RUNX3", "PGR", "Keratin 8/18"]
}

# AnnData var_names after Stage F
# ["HLA-DR", "EPCAM", "CD105", "PAX8", "RUNX3", "PGR", "Keratin 8_18"]

# Result: "Keratin 8/18" not found in AnnData, marker silently ignored
```

### Proposed Solution

Add a marker map validation step in Stage F that:

1. **Input**: Optional path to marker map JSON file (`--marker-map` argument)
2. **When**: After column sanitization, before writing AnnData
3. **Action**:
   - Extract all marker names from the JSON hierarchy
   - Compare against sanitized `var_names`
   - Report any mismatches with suggested fixes

### Implementation Sketch

```python
def validate_marker_map_names(
    adata: anndata.AnnData,
    marker_map_path: Path,
    logger,
) -> List[str]:
    """Validate marker map names match AnnData var_names after sanitization.

    Returns list of warnings for mismatched marker names.
    """
    import json

    with open(marker_map_path) as f:
        marker_map = json.load(f)

    # Extract all marker names from hierarchy
    all_markers = set()
    def extract_markers(node):
        if isinstance(node, dict):
            if "markers" in node:
                all_markers.update(node["markers"])
            if "Anti_markers" in node:
                all_markers.update(node["Anti_markers"])
            if "subtypes" in node:
                for subtype in node["subtypes"].values():
                    extract_markers(subtype)

    for key, value in marker_map.items():
        if not key.startswith("_"):  # Skip metadata
            extract_markers(value)

    var_names = set(adata.var_names)
    warnings = []

    for marker in all_markers:
        if marker not in var_names:
            # Check if sanitized version exists
            sanitized = marker
            for bad, replacement in _COLUMN_SANITIZE_MAP.items():
                sanitized = sanitized.replace(bad, replacement)

            if sanitized in var_names:
                warnings.append(
                    f"Marker '{marker}' in marker map should be '{sanitized}' "
                    f"(sanitized during Stage F)"
                )
            else:
                warnings.append(
                    f"Marker '{marker}' not found in AnnData (neither original nor sanitized)"
                )

    return warnings
```

### CLI Addition

```bash
python -m celltype_refinery.core.preprocessing --stage F \
    --marker-map configs/tissues/fallopian_tube/markers_v9.json \
    ...

# Output:
# [WARNING] Marker map validation found 1 issue:
#   - Marker 'Keratin 8/18' in marker map should be 'Keratin 8_18' (sanitized during Stage F)
```

### Acceptance Criteria

1. Stage F accepts optional `--marker-map` argument
2. When provided, validates all marker names after sanitization
3. Logs clear warnings with exact original → sanitized mapping
4. Does not fail the pipeline (warning only, unless `--strict` is set)
5. Works with any marker map JSON format (handles nested subtypes)

### Related

- Fixed in marker maps: `Keratin 8/18` → `Keratin 8_18` (2025-12-29)
- Verification copies preserved: `*_verification.json` files retain original names for testing deterministic stages
