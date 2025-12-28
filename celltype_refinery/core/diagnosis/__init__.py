"""Diagnostic module for cluster refinement decisions.

This module provides shared diagnostic functionality:

- DiagnosticEngine: Generates diagnostic reports with SUBCLUSTER/RELABEL/SKIP recommendations
- CriteriaConfig: Configuration for diagnostic thresholds
- CriteriaEvaluator: Evaluates refinement criteria for individual clusters
- StoppingCriteria: Controls iteration stopping
- StoppingStatus: Result of stopping criteria evaluation

Example Usage
-------------
Basic diagnostic report:

    >>> from celltype_refinery.core.diagnosis import DiagnosticEngine, CriteriaConfig
    >>> engine = DiagnosticEngine(config=CriteriaConfig(score_threshold=1.0))
    >>> report = engine.diagnose(cluster_annotations, marker_scores)
    >>> report.to_csv("diagnostic_report.csv", index=False)

Iterative refinement with stopping criteria:

    >>> from celltype_refinery.core.diagnosis import DiagnosticEngine, StoppingCriteria
    >>> engine = DiagnosticEngine()
    >>> stopping = StoppingCriteria(max_iterations=10)
    >>>
    >>> for iteration in range(1, 11):
    ...     # Execute refinement
    ...     execute_refinement(...)
    ...
    ...     # Generate diagnostic report for next iteration
    ...     report = engine.diagnose(cluster_annotations, marker_scores)
    ...
    ...     # Check if we should stop
    ...     status = stopping.evaluate(report, iteration)
    ...     if status.should_stop:
    ...         print(f"Stopping: {status.reason}")
    ...         break
"""

from .criteria import (
    CriteriaConfig,
    CriteriaEvaluator,
    build_hierarchy_from_scores,
)
from .engine import DiagnosticEngine
from .stopping import StoppingCriteria, StoppingStatus, check_convergence

__all__ = [
    # Core classes
    "DiagnosticEngine",
    "CriteriaConfig",
    "CriteriaEvaluator",
    # Stopping criteria
    "StoppingCriteria",
    "StoppingStatus",
    # Utilities
    "build_hierarchy_from_scores",
    "check_convergence",
]
