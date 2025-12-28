"""ReviewEngine - main orchestrator for cell-type annotation review.

Coordinates rule execution and result aggregation.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional
import logging
import time

from .config import ReviewConfig, TissueTemplate
from .aggregator import MetricsAggregator
from .rules import RuleRegistry
from .rules.base import FlaggedIssue

logger = logging.getLogger(__name__)


@dataclass
class ReviewResult:
    """Result from single-column review.

    Attributes
    ----------
    cell_type_col : str
        Cell type column reviewed
    issues : List[FlaggedIssue]
        All flagged issues
    n_issues : int
        Total issue count
    n_critical : int
        Critical issue count
    n_warning : int
        Warning issue count
    n_note : int
        Note count
    summary : Dict
        Summary statistics
    rules_run : List[str]
        Rules that were executed
    rules_skipped : List[str]
        Rules that were skipped
    execution_time_seconds : float
        Execution time
    """

    cell_type_col: str = ""
    issues: List[FlaggedIssue] = field(default_factory=list)
    n_issues: int = 0
    n_critical: int = 0
    n_warning: int = 0
    n_note: int = 0
    summary: Dict[str, Any] = field(default_factory=dict)
    rules_run: List[str] = field(default_factory=list)
    rules_skipped: List[str] = field(default_factory=list)
    execution_time_seconds: float = 0.0

    def __post_init__(self):
        self._update_counts()

    def _update_counts(self):
        """Update issue counts from issues list."""
        self.n_issues = len(self.issues)
        self.n_critical = sum(1 for i in self.issues if i.severity == "critical")
        self.n_warning = sum(1 for i in self.issues if i.severity == "warning")
        self.n_note = sum(1 for i in self.issues if i.severity == "note")

    def add_issue(self, issue: FlaggedIssue) -> None:
        """Add an issue and update counts."""
        self.issues.append(issue)
        self._update_counts()

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            "cell_type_col": self.cell_type_col,
            "n_issues": self.n_issues,
            "n_critical": self.n_critical,
            "n_warning": self.n_warning,
            "n_note": self.n_note,
            "summary": self.summary,
            "rules_run": self.rules_run,
            "rules_skipped": self.rules_skipped,
            "execution_time_seconds": round(self.execution_time_seconds, 2),
            "issues": [i.to_dict() for i in self.issues],
        }


@dataclass
class MultiColumnReviewResult:
    """Result from multi-column review.

    Attributes
    ----------
    results : Dict[str, ReviewResult]
        Per-column results
    columns_processed : List[str]
        Columns that were processed
    columns_skipped : List[str]
        Columns that were skipped
    combined_issues : List[FlaggedIssue]
        All issues with column attribution
    total_execution_time_seconds : float
        Total execution time
    """

    results: Dict[str, ReviewResult] = field(default_factory=dict)
    columns_processed: List[str] = field(default_factory=list)
    columns_skipped: List[str] = field(default_factory=list)
    combined_issues: List[FlaggedIssue] = field(default_factory=list)
    total_execution_time_seconds: float = 0.0

    @property
    def n_total_issues(self) -> int:
        return len(self.combined_issues)

    @property
    def n_total_critical(self) -> int:
        return sum(1 for i in self.combined_issues if i.severity == "critical")

    def get_result(self, col: str) -> Optional[ReviewResult]:
        """Get result for a specific column."""
        return self.results.get(col)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            "columns_processed": self.columns_processed,
            "columns_skipped": self.columns_skipped,
            "n_total_issues": self.n_total_issues,
            "n_total_critical": self.n_total_critical,
            "total_execution_time_seconds": round(self.total_execution_time_seconds, 2),
            "per_column": {
                col: result.to_dict() for col, result in self.results.items()
            },
        }


class ReviewEngine:
    """Engine for cell-type annotation review.

    Parameters
    ----------
    config : ReviewConfig, optional
        Review configuration
    template : TissueTemplate, optional
        Tissue template with expectations

    Example
    -------
    >>> from celltype_refinery.core.review import ReviewEngine, ReviewConfig, TissueTemplate
    >>> config = ReviewConfig.from_yaml("review.yaml")
    >>> template = TissueTemplate(Path("template.yaml"))
    >>> engine = ReviewEngine(config, template)
    >>> result = engine.execute(
    ...     composition_dir=Path("composition/"),
    ...     spatial_dir=Path("spatial/"),
    ... )
    """

    def __init__(
        self,
        config: Optional[ReviewConfig] = None,
        template: Optional[TissueTemplate] = None,
    ):
        self.config = config or ReviewConfig.default()
        self.template = template or TissueTemplate.empty()

        # Load template from config if specified
        if self.config.template_path and template is None:
            self.template = TissueTemplate(self.config.template_path)

    def execute(
        self,
        composition_dir: Optional[Path] = None,
        spatial_dir: Optional[Path] = None,
        consolidation_dir: Optional[Path] = None,
        cell_type_col: Optional[str] = None,
    ) -> ReviewResult:
        """Execute review for a single cell type column.

        Parameters
        ----------
        composition_dir : Path, optional
            Directory with composition outputs
        spatial_dir : Path, optional
            Directory with spatial outputs
        consolidation_dir : Path, optional
            Directory with consolidation outputs
        cell_type_col : str, optional
            Cell type column name (for labeling)

        Returns
        -------
        ReviewResult
            Review results
        """
        start_time = time.time()
        result = ReviewResult(cell_type_col=cell_type_col or "")

        # Load data
        aggregator = MetricsAggregator(
            composition_dir=composition_dir,
            spatial_dir=spatial_dir,
            consolidation_dir=consolidation_dir,
        )

        logger.info(f"Loaded data: {aggregator.summary()}")

        # Build dataset summary
        result.summary["dataset"] = self._build_dataset_summary(aggregator)

        # Instantiate and run rules
        rules = RuleRegistry.instantiate_all(
            template=self.template,
            configs=self.config.rules,
            skip_rules=self.config.skip_rules,
        )

        for rule in rules:
            if not self.config.is_rule_enabled(rule.rule_id):
                result.rules_skipped.append(rule.rule_id)
                continue

            if not rule.is_applicable(aggregator):
                result.rules_skipped.append(rule.rule_id)
                continue

            try:
                issues = rule.check(aggregator)
                for issue in issues:
                    result.add_issue(issue)
                result.rules_run.append(rule.rule_id)
            except Exception as e:
                logger.warning(f"Rule {rule.rule_id} failed: {e}")
                result.rules_skipped.append(rule.rule_id)

        result.execution_time_seconds = time.time() - start_time
        logger.info(
            f"Review complete: {result.n_issues} issues "
            f"({result.n_critical} critical, {result.n_warning} warning, {result.n_note} note) "
            f"in {result.execution_time_seconds:.2f}s"
        )

        return result

    def execute_multi(
        self,
        composition_dirs: Optional[Dict[str, Path]] = None,
        spatial_dirs: Optional[Dict[str, Path]] = None,
        consolidation_dirs: Optional[Dict[str, Path]] = None,
        cell_type_cols: Optional[List[str]] = None,
    ) -> MultiColumnReviewResult:
        """Execute review for multiple cell type columns.

        Parameters
        ----------
        composition_dirs : Dict[str, Path], optional
            Per-column composition directories
        spatial_dirs : Dict[str, Path], optional
            Per-column spatial directories
        consolidation_dirs : Dict[str, Path], optional
            Per-column consolidation directories
        cell_type_cols : List[str], optional
            Cell type column names

        Returns
        -------
        MultiColumnReviewResult
            Multi-column review results
        """
        start_time = time.time()
        result = MultiColumnReviewResult()

        cols = cell_type_cols or self.config.cell_type_cols
        if not cols:
            logger.warning("No cell type columns specified")
            return result

        composition_dirs = composition_dirs or {}
        spatial_dirs = spatial_dirs or {}
        consolidation_dirs = consolidation_dirs or {}

        for col in cols:
            logger.info(f"\n{'='*60}\nReviewing column: {col}\n{'='*60}")

            col_result = self.execute(
                composition_dir=composition_dirs.get(col),
                spatial_dir=spatial_dirs.get(col),
                consolidation_dir=consolidation_dirs.get(col),
                cell_type_col=col,
            )

            result.results[col] = col_result
            result.columns_processed.append(col)

            # Add issues with column attribution
            for issue in col_result.issues:
                issue.column = col
                result.combined_issues.append(issue)

        result.total_execution_time_seconds = time.time() - start_time

        logger.info(
            f"\nMulti-column review complete: "
            f"{len(result.columns_processed)} columns, "
            f"{result.n_total_issues} total issues "
            f"in {result.total_execution_time_seconds:.2f}s"
        )

        return result

    def _build_dataset_summary(self, aggregator: MetricsAggregator) -> Dict[str, Any]:
        """Build dataset summary from aggregator."""
        summary = {
            "n_cell_types": len(aggregator.get_cell_types()),
            "n_regions": len(aggregator.get_regions()),
        }

        # Try to get total cells
        if aggregator.composition_global is not None:
            if "count" in aggregator.composition_global.columns:
                summary["total_cells"] = int(
                    aggregator.composition_global["count"].sum()
                )
            elif "n_cells" in aggregator.composition_global.columns:
                summary["total_cells"] = int(
                    aggregator.composition_global["n_cells"].sum()
                )

        # Calculate assignment rate
        unassigned_pct = aggregator.get_global_pct("Unassigned")
        if unassigned_pct is not None:
            summary["unassigned_pct"] = round(unassigned_pct, 2)
            summary["assignment_rate"] = round(100 - unassigned_pct, 2)

        # Calculate hybrid percentage
        hybrid_pct = 0.0
        for ct in aggregator.get_cell_types():
            if "~" in ct:
                pct = aggregator.get_global_pct(ct)
                if pct:
                    hybrid_pct += pct
        summary["hybrid_pct"] = round(hybrid_pct, 2)

        # Sample count
        if aggregator.composition_by_sample is not None:
            if "sample_id" in aggregator.composition_by_sample.columns:
                summary["n_samples"] = aggregator.composition_by_sample[
                    "sample_id"
                ].nunique()

        return summary

    def list_available_rules(self) -> Dict[str, List[str]]:
        """List all available rules by category.

        Returns
        -------
        Dict[str, List[str]]
            Map of category to list of rule IDs
        """
        return RuleRegistry.summary()
