import clsx from 'clsx';
import Heading from '@theme/Heading';
import styles from './styles.module.css';

const FeatureList = [
  {
    title: 'Hierarchical Annotation',
    emoji: '🧬',
    description: (
      <>
        Marker-based hierarchical gating with multi-level scoring.
        Configure for any tissue via JSON marker maps.
      </>
    ),
  },
  {
    title: 'Iterative Refinement',
    emoji: '🔄',
    description: (
      <>
        Policy-based refinement with automatic and manual modes.
        Subcluster, merge, and override operations with full provenance.
      </>
    ),
  },
  {
    title: 'Spatial Analysis',
    emoji: '🗺️',
    description: (
      <>
        Neighborhood enrichment, Moran's I, and cell-type interactions.
        Permutation-based statistical tests with FDR correction.
      </>
    ),
  },
  {
    title: 'Composition Metrics',
    emoji: '📊',
    description: (
      <>
        Shannon diversity, regional enrichment, and tissue-specific
        biology metrics with comprehensive export formats.
      </>
    ),
  },
];

function Feature({emoji, title, description}) {
  return (
    <div className={clsx('col col--3')}>
      <div className="text--center">
        <span className={styles.featureEmoji} role="img" aria-label={title}>
          {emoji}
        </span>
      </div>
      <div className="text--center padding-horiz--md">
        <Heading as="h3">{title}</Heading>
        <p>{description}</p>
      </div>
    </div>
  );
}

export default function HomepageFeatures() {
  return (
    <section className={styles.features}>
      <div className="container">
        <div className="row">
          {FeatureList.map((props, idx) => (
            <Feature key={idx} {...props} />
          ))}
        </div>
      </div>
    </section>
  );
}
