"""CellType-Refinery: Cell-type annotation and refinement for spatial proteomics data.

This package provides tools for:
- Cell-type annotation using hierarchical marker-based gating
- Iterative refinement with automatic and manual policies
- Composition analysis with diversity metrics
- Spatial analysis including neighborhood enrichment
- Quality review with configurable flagging rules

The package is designed to be tissue-agnostic, with tissue-specific
parameters loaded from YAML configuration files.

Example usage:
    >>> from celltype_refinery.core.annotation import AnnotationEngine
    >>> from celltype_refinery.core.refinement import RefinementEngine
    >>>
    >>> # Load and annotate
    >>> engine = AnnotationEngine(marker_map_path)
    >>> result = engine.run(adata, output_dir)
    >>>
    >>> # Refine annotations
    >>> refiner = RefinementEngine()
    >>> refiner.execute(adata, plan)
"""

__version__ = "0.1.0"

# Package-level exports will be added as modules are implemented
