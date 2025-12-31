// @ts-check

/**
 * @type {import('@docusaurus/plugin-content-docs').SidebarsConfig}
 */
const sidebars = {
  docsSidebar: [
    'intro',
    {
      type: 'category',
      label: 'Getting Started',
      items: [
        'getting-started/installation',
        'getting-started/quickstart',
        'getting-started/configuration',
        'getting-started/first-annotation',
      ],
    },
    {
      type: 'category',
      label: 'Core Workflows',
      items: [
        'core-workflows/workflow-overview',
        'core-workflows/automated-workflow',
        'core-workflows/expert-curation',
        'core-workflows/hybrid-workflow',
        'core-workflows/iterative-refinement',
      ],
    },
    {
      type: 'category',
      label: 'Modules',
      items: [
        {
          type: 'category',
          label: 'Preprocessing',
          items: [
            'modules/preprocessing/overview',
            'modules/preprocessing/data-loading',
            'modules/preprocessing/cell-qc',
            'modules/preprocessing/normalization',
            'modules/preprocessing/alignment',
            'modules/preprocessing/batch-correction',
            'modules/preprocessing/data-merging',
          ],
        },
        {
          type: 'category',
          label: 'Clustering',
          items: [
            'modules/clustering/overview',
            'modules/clustering/leiden-clustering',
            'modules/clustering/differential-expression',
          ],
        },
        {
          type: 'category',
          label: 'Annotation',
          items: [
            'modules/annotation/overview',
            'modules/annotation/marker-maps',
            'modules/annotation/hierarchical-gating',
            'modules/annotation/marker-scoring',
            'modules/annotation/exports',
          ],
        },
        {
          type: 'category',
          label: 'Refinement',
          items: [
            'modules/refinement/overview',
            'modules/refinement/policies',
            'modules/refinement/operations',
            'modules/refinement/diagnostic-mode',
            'modules/refinement/execution-mode',
          ],
        },
        {
          type: 'category',
          label: 'Consolidation',
          items: [
            'modules/consolidation/overview',
            'modules/consolidation/orphan-rescue',
            'modules/consolidation/iel-rescue',
            'modules/consolidation/harmonization',
          ],
        },
        {
          type: 'category',
          label: 'Composition',
          items: [
            'modules/composition/overview',
            'modules/composition/diversity-metrics',
            'modules/composition/biology-metrics',
            'modules/composition/enrichment',
          ],
        },
        {
          type: 'category',
          label: 'Spatial',
          items: [
            'modules/spatial/overview',
            'modules/spatial/neighborhood-enrichment',
            'modules/spatial/morans-i',
            'modules/spatial/local-diversity',
          ],
        },
        {
          type: 'category',
          label: 'Review',
          items: [
            'modules/review/overview',
            'modules/review/flagging-rules',
            'modules/review/quality-metrics',
          ],
        },
      ],
    },
    {
      type: 'category',
      label: 'CLI Reference',
      items: [
        'cli/overview',
        'cli/preprocess',
        'cli/cluster',
        'cli/annotate',
        'cli/refine',
        'cli/consolidate',
        'cli/analyze',
        'cli/pipeline',
      ],
    },
    {
      type: 'category',
      label: 'Configuration',
      items: [
        'configuration/marker-maps',
        'configuration/tissue-templates',
        'configuration/pipeline-yaml',
        'configuration/refinement-yaml',
      ],
    },
    {
      type: 'category',
      label: 'Examples',
      items: [
        'examples/fallopian-tube',
        'examples/custom-tissue',
      ],
    },
  ],
};

export default sidebars;
