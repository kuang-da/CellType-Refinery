# Claude Code Guidelines for celltype-refinery-docs

## Documentation Standards

### Link Validation

**CRITICAL: Only link to pages that exist.**

Before adding any link in documentation:
1. Verify the target file exists in the `docs/` directory
2. Use the correct relative path from the source file
3. Do NOT create placeholder links to non-existent pages

### Valid Link Paths

The documentation uses these established paths:

```
docs/
├── methodology/
│   ├── index.md
│   ├── annotation-pipeline.md
│   ├── marker-scoring-algorithm.md
│   ├── hierarchical-gating-algorithm.md
│   ├── refinement-decision-logic.md
│   └── tuning-guide.md
├── reference/
│   ├── glossary.md
│   └── output-files.md
├── modules/
│   ├── annotation/
│   │   ├── overview.md
│   │   ├── marker-maps.md
│   │   ├── hierarchical-gating.md
│   │   ├── marker-scoring.md
│   │   ├── per-cell-voting.md
│   │   └── exports.md
│   ├── refinement/
│   │   ├── overview.md
│   │   ├── policies.md
│   │   ├── operations.md
│   │   ├── diagnostic-mode.md
│   │   ├── execution-mode.md
│   │   ├── focus-controls.md
│   │   └── iterative-refinement.md
│   └── clustering/
│       ├── overview.md
│       ├── leiden-clustering.md
│       └── differential-expression.md
├── configuration/
│   ├── marker-maps.md
│   ├── tissue-templates.md
│   ├── pipeline-yaml.md
│   └── refinement-yaml.md
└── cli/
    ├── overview.md
    ├── cluster.md
    ├── annotate.md
    └── refine.md
```

### Link Format

Use relative paths from the current file:
- Same directory: `./other-page.md`
- Parent directory: `../other-section/page.md`
- Absolute from docs root: `/docs/section/page.md`

### Common Mistakes to Avoid

1. **Don't link to planned but non-existent pages**
   - Bad: `[API Reference](/docs/api/overview.md)` (if api/ doesn't exist)
   - Good: Remove the link or create the page first

2. **Don't use wrong paths**
   - Bad: `/docs/guides/tuning-guide.md`
   - Good: `/docs/methodology/tuning-guide.md`

3. **Don't assume page names**
   - Bad: `./refinement-algorithm.md`
   - Good: `./refinement-decision-logic.md`

### CI Build Validation

The Docusaurus build will fail if broken links are detected. Always run:

```bash
cd celltype-refinery-docs
npm run build
```

This validates all links before pushing.

## ASCII Diagrams

Use ASCII art for complex diagrams (renders consistently everywhere). Keep ASCII in code blocks:

```text
┌─────────────┐
│  Example    │
└─────────────┘
```

Use Mermaid only for simple flowcharts (4-5 nodes max).
