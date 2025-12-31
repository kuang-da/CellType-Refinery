---
sidebar_position: 1
---

# Refinement Overview

The refinement module improves annotations through automatic and manual policies.

## Architecture

```mermaid
flowchart TB
    AUTO[AutoPolicy] --> BP[base_plan]
    MANUAL[ManualPolicy] --> OP[overlay_plan]
    BP --> MERGE[merge_plans]
    OP --> MERGE
    MERGE --> ENGINE[RefinementEngine]
    ENGINE --> OUTPUT[refined.h5ad]
```

## Execution Modes

- **Diagnostic**: Generate recommendations only
- **Auto-only**: Automatic refinement
- **Manual-only**: YAML configuration
- **Hybrid**: Auto + manual overrides
