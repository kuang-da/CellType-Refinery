"""Export functions for the Review module.

Provides:
- export_json: Full structured output
- export_csv: Flattened issues table
- export_markdown: Human-readable report
- export_provenance: Audit trail
"""

import json
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import pandas as pd

from .engine import ReviewResult, MultiColumnReviewResult
from .rules.base import FlaggedIssue


def export_json(
    result: Union[ReviewResult, MultiColumnReviewResult],
    output_path: Path,
) -> Path:
    """Export review result to JSON file.

    Parameters
    ----------
    result : ReviewResult or MultiColumnReviewResult
        Review result
    output_path : Path
        Output file path

    Returns
    -------
    Path
        Path to created file
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    data = result.to_dict()
    data["export_timestamp"] = datetime.now().isoformat()
    data["export_version"] = "2.0"

    with open(output_path, "w") as f:
        json.dump(data, f, indent=2, default=str)

    return output_path


def export_csv(
    issues: List[FlaggedIssue],
    output_path: Path,
) -> Path:
    """Export flagged issues to CSV file.

    Parameters
    ----------
    issues : List[FlaggedIssue]
        List of flagged issues
    output_path : Path
        Output file path

    Returns
    -------
    Path
        Path to created file
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    rows = [issue.to_dict() for issue in issues]
    df = pd.DataFrame(rows)

    # Sort by severity then rule
    severity_order = {"critical": 0, "warning": 1, "note": 2}
    if len(df) > 0:
        df["severity_order"] = df["severity"].map(severity_order)
        df = df.sort_values(["severity_order", "rule_id", "cell_type"])
        df = df.drop(columns=["severity_order"])

    df.to_csv(output_path, index=False)
    return output_path


def export_markdown(
    result: Union[ReviewResult, MultiColumnReviewResult],
    output_path: Path,
    include_recommendations: bool = True,
) -> Path:
    """Export review result to Markdown report.

    Parameters
    ----------
    result : ReviewResult or MultiColumnReviewResult
        Review result
    output_path : Path
        Output file path
    include_recommendations : bool
        Whether to include recommendations

    Returns
    -------
    Path
        Path to created file
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    lines = []

    if isinstance(result, MultiColumnReviewResult):
        lines.extend(_format_multi_column_report(result, include_recommendations))
    else:
        lines.extend(_format_single_column_report(result, include_recommendations))

    with open(output_path, "w") as f:
        f.write("\n".join(lines))

    return output_path


def _format_single_column_report(
    result: ReviewResult,
    include_recommendations: bool,
) -> List[str]:
    """Format single-column report."""
    lines = [
        "# Cell-Type Annotation Review Report",
        "",
        f"**Column**: `{result.cell_type_col}`",
        f"**Generated**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"**Execution Time**: {result.execution_time_seconds:.2f}s",
        "",
    ]

    # Summary section
    lines.extend([
        "## Summary",
        "",
        "| Metric | Value |",
        "|--------|-------|",
        f"| Total Issues | {result.n_issues} |",
        f"| Critical | {result.n_critical} |",
        f"| Warning | {result.n_warning} |",
        f"| Note | {result.n_note} |",
        "",
    ])

    # Dataset summary
    if "dataset" in result.summary:
        ds = result.summary["dataset"]
        lines.extend([
            "### Dataset",
            "",
            f"- Total cells: {ds.get('total_cells', 0):,}",
            f"- Samples: {ds.get('n_samples', 0)}",
            f"- Cell types: {ds.get('n_cell_types', 0)}",
            f"- Assignment rate: {ds.get('assignment_rate', 0):.1f}%",
            f"- Unassigned: {ds.get('unassigned_pct', 0):.1f}%",
            f"- Hybrids: {ds.get('hybrid_pct', 0):.1f}%",
            "",
        ])

    # Issues by severity
    if result.issues:
        lines.extend([
            "## Issues",
            "",
        ])

        # Group by severity
        grouped = {"critical": [], "warning": [], "note": []}
        for issue in result.issues:
            grouped[issue.severity].append(issue)

        for severity in ["critical", "warning", "note"]:
            issues = grouped[severity]
            if issues:
                icon = {"critical": "!!!", "warning": "!!", "note": "!"}[severity]
                lines.append(f"### {severity.title()} ({len(issues)})")
                lines.append("")

                for issue in issues:
                    lines.append(f"**[{icon}] {issue.rule_id}**: {issue.cell_type}")
                    lines.append(f"- {issue.description}")
                    lines.append(f"- Evidence: {issue.evidence}")
                    if include_recommendations and issue.recommendation:
                        lines.append(f"- Recommendation: {issue.recommendation}")
                    lines.append("")

    else:
        lines.extend([
            "## Issues",
            "",
            "No issues found.",
            "",
        ])

    return lines


def _format_multi_column_report(
    result: MultiColumnReviewResult,
    include_recommendations: bool,
) -> List[str]:
    """Format multi-column comparison report."""
    lines = [
        "# Cell-Type Annotation Review Report (Multi-Column)",
        "",
        f"**Generated**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"**Execution Time**: {result.total_execution_time_seconds:.2f}s",
        "",
    ]

    # Comparison table
    lines.extend([
        "## Column Comparison",
        "",
        "| Column | Issues | Critical | Warning | Note |",
        "|--------|--------|----------|---------|------|",
    ])

    for col in result.columns_processed:
        r = result.results[col]
        lines.append(
            f"| `{col}` | {r.n_issues} | {r.n_critical} | {r.n_warning} | {r.n_note} |"
        )

    lines.append("")

    # Per-column details
    for col in result.columns_processed:
        r = result.results[col]
        lines.extend([
            "---",
            "",
            f"## {col}",
            "",
        ])

        if r.issues:
            # Just show counts by rule
            rule_counts: Dict[str, int] = {}
            for issue in r.issues:
                rule_counts[issue.rule_id] = rule_counts.get(issue.rule_id, 0) + 1

            lines.append("| Rule | Count |")
            lines.append("|------|-------|")
            for rule, count in sorted(rule_counts.items(), key=lambda x: -x[1]):
                lines.append(f"| {rule} | {count} |")
            lines.append("")
        else:
            lines.append("No issues found.")
            lines.append("")

    # Combined issues (top 10)
    if result.combined_issues:
        lines.extend([
            "---",
            "",
            "## Top Issues (All Columns)",
            "",
        ])

        critical = [i for i in result.combined_issues if i.severity == "critical"]
        for issue in critical[:10]:
            lines.append(
                f"- **[!!!] {issue.rule_id}** ({issue.column}): "
                f"{issue.cell_type} - {issue.description}"
            )

        lines.append("")

    return lines


def export_provenance(
    result: Union[ReviewResult, MultiColumnReviewResult],
    output_path: Path,
    config_path: Optional[Path] = None,
    template_path: Optional[Path] = None,
) -> Path:
    """Export provenance information for audit trail.

    Parameters
    ----------
    result : ReviewResult or MultiColumnReviewResult
        Review result
    output_path : Path
        Output file path
    config_path : Path, optional
        Config file used
    template_path : Path, optional
        Template file used

    Returns
    -------
    Path
        Path to created file
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    provenance = {
        "export_timestamp": datetime.now().isoformat(),
        "export_version": "2.0",
        "config_path": str(config_path) if config_path else None,
        "template_path": str(template_path) if template_path else None,
    }

    if isinstance(result, MultiColumnReviewResult):
        provenance["type"] = "multi_column"
        provenance["columns_processed"] = result.columns_processed
        provenance["total_issues"] = len(result.combined_issues)
        provenance["execution_time_seconds"] = result.total_execution_time_seconds
    else:
        provenance["type"] = "single_column"
        provenance["cell_type_col"] = result.cell_type_col
        provenance["total_issues"] = result.n_issues
        provenance["execution_time_seconds"] = result.execution_time_seconds

    with open(output_path, "w") as f:
        json.dump(provenance, f, indent=2)

    return output_path


def export_all(
    result: Union[ReviewResult, MultiColumnReviewResult],
    output_dir: Path,
    config_path: Optional[Path] = None,
    template_path: Optional[Path] = None,
) -> Dict[str, Path]:
    """Export all review outputs.

    Parameters
    ----------
    result : ReviewResult or MultiColumnReviewResult
        Review result
    output_dir : Path
        Output directory
    config_path : Path, optional
        Config file used
    template_path : Path, optional
        Template file used

    Returns
    -------
    Dict[str, Path]
        Map of output type to file path
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    outputs = {}

    # JSON summary
    if isinstance(result, MultiColumnReviewResult):
        outputs["json"] = export_json(result, output_dir / "review_summary_multi.json")
        outputs["csv"] = export_csv(
            result.combined_issues, output_dir / "flagged_issues.csv"
        )
    else:
        outputs["json"] = export_json(result, output_dir / "review_summary.json")
        outputs["csv"] = export_csv(result.issues, output_dir / "flagged_issues.csv")

    # Markdown report
    outputs["markdown"] = export_markdown(result, output_dir / "review_report.md")

    # Provenance
    outputs["provenance"] = export_provenance(
        result, output_dir / "provenance.json", config_path, template_path
    )

    return outputs
