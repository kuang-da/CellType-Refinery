"""
Audit Card Generation - Per-cluster visual evidence cards.

This module generates HTML audit cards for each cluster showing:
- Decision path trace through the hierarchy
- Top competing cell types and their scores
- Marker evidence (if marker_evidence table is available)
- QC flags and warnings

Ported from ft/src/stage_h_audit_cards.py
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd


# CSS styles for audit cards
AUDIT_CARD_CSS = """
<style>
    body {
        font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
        background-color: #f5f5f5;
        margin: 0;
        padding: 20px;
    }
    .header {
        text-align: center;
        padding: 20px;
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        border-radius: 10px;
        margin-bottom: 20px;
    }
    .header h1 {
        margin: 0;
        font-size: 2em;
    }
    .header p {
        margin: 10px 0 0 0;
        opacity: 0.9;
    }
    .summary-stats {
        display: flex;
        justify-content: center;
        gap: 20px;
        flex-wrap: wrap;
        margin-top: 15px;
    }
    .stat-box {
        background: rgba(255,255,255,0.2);
        padding: 10px 20px;
        border-radius: 5px;
        text-align: center;
    }
    .stat-box .value {
        font-size: 1.5em;
        font-weight: bold;
    }
    .stat-box .label {
        font-size: 0.9em;
        opacity: 0.8;
    }
    .card-grid {
        display: grid;
        grid-template-columns: repeat(auto-fill, minmax(450px, 1fr));
        gap: 20px;
    }
    .audit-card {
        background: white;
        border-radius: 10px;
        box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        overflow: hidden;
    }
    .card-header {
        background: linear-gradient(135deg, #4A90D9 0%, #357ABD 100%);
        color: white;
        padding: 15px;
        display: flex;
        justify-content: space-between;
        align-items: center;
        flex-wrap: wrap;
    }
    .card-header h3 {
        margin: 0;
        font-size: 1.1em;
    }
    .card-header .badges {
        display: flex;
        gap: 8px;
        flex-wrap: wrap;
    }
    .badge {
        padding: 4px 10px;
        border-radius: 15px;
        font-size: 0.85em;
        font-weight: 500;
    }
    .badge.score { background: rgba(255,255,255,0.3); }
    .badge.cells { background: rgba(255,255,255,0.2); }
    .badge.high { background: #4CAF50; }
    .badge.medium { background: #FFC107; color: #333; }
    .badge.low { background: #FF9800; }
    .badge.very_low { background: #F44336; }
    .card-body {
        padding: 15px;
    }
    .section {
        margin-bottom: 15px;
    }
    .section h4 {
        margin: 0 0 8px 0;
        font-size: 0.95em;
        color: #555;
        border-bottom: 1px solid #eee;
        padding-bottom: 5px;
    }
    .decision-trace {
        display: flex;
        align-items: center;
        flex-wrap: wrap;
        gap: 5px;
        font-size: 0.9em;
    }
    .trace-step {
        background: #e3f2fd;
        padding: 4px 10px;
        border-radius: 15px;
        color: #1565c0;
    }
    .trace-step.current {
        background: #1565c0;
        color: white;
    }
    .trace-arrow {
        color: #999;
    }
    .competing-table {
        width: 100%;
        border-collapse: collapse;
        font-size: 0.85em;
    }
    .competing-table th, .competing-table td {
        padding: 8px;
        text-align: left;
        border-bottom: 1px solid #eee;
    }
    .competing-table th {
        background: #f9f9f9;
        font-weight: 600;
    }
    .competing-table tr:hover {
        background: #f5f5f5;
    }
    .competing-table .winner-row {
        background: #e8f5e9;
        font-weight: 600;
    }
    .competing-table .winner-row:hover {
        background: #c8e6c9;
    }
    .margin-badge {
        display: inline-block;
        padding: 2px 8px;
        border-radius: 4px;
        font-size: 0.8em;
        margin-left: 8px;
    }
    .margin-badge.high {
        background: #c8e6c9;
        color: #2e7d32;
    }
    .margin-badge.medium {
        background: #fff9c4;
        color: #f57f17;
    }
    .margin-badge.low {
        background: #ffcdd2;
        color: #c62828;
    }
    .note {
        color: #666;
        font-size: 0.85em;
        margin-bottom: 8px;
    }
    .marker-evidence {
        font-size: 0.85em;
    }
    .marker-row {
        display: flex;
        justify-content: space-between;
        padding: 4px 0;
        border-bottom: 1px solid #f0f0f0;
    }
    .marker-name {
        font-weight: 500;
    }
    .marker-stats {
        color: #666;
    }
    .qc-flags {
        display: flex;
        flex-wrap: wrap;
        gap: 8px;
    }
    .flag {
        padding: 4px 10px;
        border-radius: 5px;
        font-size: 0.8em;
    }
    .flag.warning { background: #fff3e0; color: #e65100; }
    .flag.critical { background: #ffebee; color: #c62828; }
    .flag.info { background: #e3f2fd; color: #1565c0; }
    .flag.success { background: #e8f5e9; color: #2e7d32; }
    .no-issues {
        color: #4CAF50;
        font-style: italic;
    }
    /* Score breakdown styles */
    .score-breakdown {
        display: grid;
        grid-template-columns: repeat(2, 1fr);
        gap: 8px;
        font-size: 0.85em;
    }
    .score-component {
        display: flex;
        justify-content: space-between;
        padding: 6px 10px;
        background: #f8f9fa;
        border-radius: 5px;
    }
    .score-component.positive { background: #e8f5e9; }
    .score-component.negative { background: #ffebee; }
    .score-component .label { color: #666; }
    .score-component .value { font-weight: 600; }
    .score-formula {
        margin-top: 8px;
        padding: 8px;
        background: #f5f5f5;
        border-radius: 5px;
        font-family: monospace;
        font-size: 0.8em;
        text-align: center;
    }
    /* Vote composition styles */
    .vote-composition {
        margin-top: 8px;
    }
    .vote-bar-container {
        display: flex;
        height: 24px;
        border-radius: 4px;
        overflow: hidden;
        margin-bottom: 8px;
    }
    .vote-bar {
        display: flex;
        align-items: center;
        justify-content: center;
        color: white;
        font-size: 0.75em;
        font-weight: 500;
        min-width: 30px;
    }
    .vote-legend {
        display: flex;
        flex-wrap: wrap;
        gap: 8px;
        font-size: 0.8em;
    }
    .vote-legend-item {
        display: flex;
        align-items: center;
        gap: 4px;
    }
    .vote-color-box {
        width: 12px;
        height: 12px;
        border-radius: 2px;
    }
    /* Ambiguous root styles */
    .ambiguous-analysis {
        background: #fff3e0;
        padding: 10px;
        border-radius: 5px;
        border-left: 3px solid #ff9800;
    }
    .ambiguous-candidates {
        display: flex;
        gap: 10px;
        margin-top: 8px;
    }
    .candidate-box {
        flex: 1;
        padding: 8px;
        background: white;
        border-radius: 5px;
        text-align: center;
    }
    .candidate-box .name { font-weight: 600; font-size: 0.9em; }
    .candidate-box .score { color: #666; font-size: 0.85em; }
    .gap-info {
        margin-top: 8px;
        font-size: 0.85em;
        color: #666;
    }
    /* Confidence explanation styles */
    .confidence-explain {
        font-size: 0.85em;
        padding: 8px;
        background: #f5f5f5;
        border-radius: 5px;
        margin-top: 8px;
    }
    .confidence-explain .metric {
        display: flex;
        justify-content: space-between;
        margin-bottom: 4px;
    }
    /* Gate summary styles */
    .gate-summary {
        font-size: 0.85em;
    }
    .gate-row {
        display: flex;
        justify-content: space-between;
        padding: 4px 0;
        border-bottom: 1px solid #f0f0f0;
    }
    .gate-row .gate-name { color: #555; }
    .gate-row .gate-status { font-weight: 500; }
    .gate-row .gate-status.pass { color: #4CAF50; }
    .gate-row .gate-status.fail { color: #F44336; }
    /* Trace step with score */
    .trace-step .trace-score {
        font-size: 0.8em;
        opacity: 0.8;
        margin-left: 4px;
    }
    .trace-step.ambiguous {
        background: #fff3e0;
        color: #e65100;
    }
    .trace-step.current.ambiguous {
        background: #ff9800;
        color: white;
    }
</style>
"""


def _render_score_breakdown(
    marker_scores: pd.DataFrame,
    cluster_id: str,
    assigned_label: str,
) -> str:
    """Render score breakdown showing all components of the final score."""
    row = marker_scores[
        (marker_scores["cluster_id"].astype(str) == str(cluster_id)) &
        (marker_scores["label"] == assigned_label)
    ]

    if row.empty:
        return "<p>Score breakdown not available</p>"

    row = row.iloc[0]

    mean_enrichment = row.get("mean_enrichment", 0.0)
    mean_positive = row.get("mean_positive_fraction", 0.0)
    de_component = row.get("de_component", 0.0)
    anti_penalty = row.get("anti_penalty", 0.0)
    frac_markers_on = row.get("frac_markers_on", 0.0)
    final_score = row.get("score", 0.0)

    def fmt(val, is_pct=False):
        if pd.isna(val):
            return "N/A"
        if is_pct:
            return f"{val:.0%}"
        return f"{val:.2f}"

    return f"""
    <div class="score-breakdown">
        <div class="score-component positive">
            <span class="label">Mean Enrichment</span>
            <span class="value">+{fmt(mean_enrichment)}</span>
        </div>
        <div class="score-component positive">
            <span class="label">Mean Positive</span>
            <span class="value">+{fmt(mean_positive)}</span>
        </div>
        <div class="score-component positive">
            <span class="label">DE Bonus</span>
            <span class="value">+{fmt(de_component)}</span>
        </div>
        <div class="score-component negative">
            <span class="label">Anti-Penalty</span>
            <span class="value">-{fmt(anti_penalty)}</span>
        </div>
    </div>
    <div class="score-formula">
        {fmt(mean_enrichment)} + {fmt(mean_positive)} + {fmt(de_component)} - {fmt(anti_penalty)} = <strong>{fmt(final_score)}</strong>
    </div>
    <div style="margin-top: 6px; font-size: 0.8em; color: #666;">
        Frac markers on: {fmt(frac_markers_on, is_pct=True)}
    </div>
    """


def _render_vote_composition(cluster_annotation: pd.Series) -> str:
    """Render vote composition for ambiguous cases."""
    stop_reason = cluster_annotation.get("stop_reason", "")
    composition_json = cluster_annotation.get("composition", "")

    if stop_reason not in ("ambiguous_siblings",) or not composition_json:
        return ""

    try:
        if isinstance(composition_json, str):
            composition = json.loads(composition_json)
        else:
            composition = composition_json

        if not composition or not isinstance(composition, dict):
            return ""
    except (json.JSONDecodeError, TypeError):
        return ""

    votes = {k: v for k, v in composition.items() if not k.startswith("_")}
    if not votes:
        return ""

    sorted_votes = sorted(votes.items(), key=lambda x: -x[1])
    colors = ["#4A90D9", "#67B26F", "#FFC107", "#FF9800", "#9C27B0", "#607D8B"]

    bar_segments = []
    legend_items = []
    for i, (label, pct) in enumerate(sorted_votes):
        color = colors[i % len(colors)]
        width_pct = pct * 100
        display_pct = "" if width_pct < 5 else f"{pct:.0%}"
        bar_segments.append(
            f'<div class="vote-bar" style="width: {width_pct}%; background: {color};">{display_pct}</div>'
        )
        legend_items.append(
            f'<div class="vote-legend-item">'
            f'<div class="vote-color-box" style="background: {color};"></div>'
            f'<span>{label}: {pct:.0%}</span>'
            f'</div>'
        )

    evidence_note = ""
    if composition.get("_evidence_only"):
        evidence_note = '<div style="font-size: 0.8em; color: #666; margin-top: 4px;">Info: Evidence only (label stays at parent)</div>'

    return f"""
    <div class="vote-composition">
        <p class="note">Per-cell voting results:</p>
        <div class="vote-bar-container">
            {"".join(bar_segments)}
        </div>
        <div class="vote-legend">
            {"".join(legend_items)}
        </div>
        {evidence_note}
    </div>
    """


def _render_ambiguous_root_analysis(cluster_annotation: pd.Series) -> str:
    """Render detailed analysis for ambiguous root cases."""
    if not cluster_annotation.get("is_ambiguous_root", False):
        return ""

    trace_json = cluster_annotation.get("decision_trace", "")
    if not trace_json:
        return ""

    try:
        if isinstance(trace_json, str):
            trace = json.loads(trace_json)
        else:
            trace = trace_json

        if not isinstance(trace, dict) or trace.get("type") != "ambiguous_root":
            return ""
    except (json.JSONDecodeError, TypeError):
        return ""

    candidates = trace.get("candidates", [])
    gap = trace.get("gap", 0)
    threshold = trace.get("threshold", 0.25)

    if len(candidates) < 2:
        return ""

    candidate_html = []
    for i, cand in enumerate(candidates[:3]):
        label = cand.get("label", "?")
        score = cand.get("score", 0)
        rank = "1st" if i == 0 else ("2nd" if i == 1 else "3rd")
        candidate_html.append(f"""
            <div class="candidate-box">
                <div class="name">{label}</div>
                <div class="score">{rank}: {score:.2f}</div>
            </div>
        """)

    return f"""
    <div class="ambiguous-analysis">
        <strong>Warning: Ambiguous Root Detection</strong>
        <div class="ambiguous-candidates">
            {"".join(candidate_html)}
        </div>
        <div class="gap-info">
            Gap between top 2: <strong>{gap:.3f}</strong> (threshold: {threshold})
            <br>Since gap &lt; threshold, both roots are considered viable candidates.
        </div>
    </div>
    """


def _render_confidence_explanation(cluster_annotation: pd.Series) -> str:
    """Render confidence/margin explanation."""
    margin = cluster_annotation.get("min_margin_along_path", 0)
    margin_is_infinite = cluster_annotation.get("margin_is_infinite", False)
    stop_reason = cluster_annotation.get("stop_reason", "")
    assigned_level = cluster_annotation.get("assigned_level", 0)

    if stop_reason == "no_root_passed":
        conf_text = "No assignment (no root passed gates)"
        conf_class = "critical"
    elif margin_is_infinite:
        conf_text = "No competition at any level (highest confidence)"
        conf_class = "success"
    elif margin >= 0.5:
        conf_text = "High margin to runner-up"
        conf_class = "success"
    elif margin >= 0.2:
        conf_text = "Moderate margin"
        conf_class = "info"
    else:
        conf_text = "Low margin (close competition)"
        conf_class = "warning"

    margin_display = "Inf (no competitors)" if margin_is_infinite else f"{margin:.3f}"

    return f"""
    <div class="confidence-explain">
        <div class="metric">
            <span>Min margin along path:</span>
            <span><strong>{margin_display}</strong></span>
        </div>
        <div class="metric">
            <span>Assignment depth:</span>
            <span>Level {assigned_level}</span>
        </div>
        <div class="metric">
            <span>Interpretation:</span>
            <span class="flag {conf_class}" style="padding: 2px 6px; font-size: 0.85em;">{conf_text}</span>
        </div>
    </div>
    """


def _render_gate_summary(
    cluster_id: str,
    decision_steps: Optional[pd.DataFrame],
) -> str:
    """Render gate check summary showing pass/fail at each level."""
    if decision_steps is None or decision_steps.empty:
        return "<p>Gate details not available</p>"

    steps = decision_steps[decision_steps["cluster_id"].astype(str) == str(cluster_id)]
    if steps.empty:
        return "<p>No gate checks recorded</p>"

    gate_rows = []
    parents = steps["parent_label"].unique()
    for parent in parents:
        level_steps = steps[steps["parent_label"] == parent]
        n_passed = level_steps["child_passed_gate"].sum()
        n_total = len(level_steps)
        n_failed = n_total - n_passed

        failed_steps = level_steps[~level_steps["child_passed_gate"]]
        fail_reasons = failed_steps["fail_reason"].dropna().unique()
        fail_reason_str = ", ".join(fail_reasons[:3]) if len(fail_reasons) > 0 else ""

        detail = ""
        if n_failed > 0 and fail_reason_str:
            detail = f'<span style="font-size: 0.8em; color: #666;">({fail_reason_str})</span>'

        gate_rows.append(f"""
            <div class="gate-row">
                <span class="gate-name">{parent} -></span>
                <span>
                    <span class="gate-status pass">Pass {n_passed}</span> /
                    <span class="gate-status fail">Fail {n_failed}</span>
                    {detail}
                </span>
            </div>
        """)

    return f'<div class="gate-summary">{"".join(gate_rows)}</div>'


def _render_decision_breadcrumb(
    cluster_annotation: pd.Series,
    decision_steps: Optional[pd.DataFrame],
    cluster_id: str,
) -> str:
    """Render decision path as HTML breadcrumb trail."""
    trace_json = cluster_annotation.get("decision_trace", "")
    if trace_json and isinstance(trace_json, str):
        try:
            trace = json.loads(trace_json)
            if trace:
                if isinstance(trace, dict) and trace.get("type") == "ambiguous_root":
                    candidates = trace.get("candidates", [])
                    if candidates:
                        labels_with_scores = [
                            f"{c.get('label', '?')} ({c.get('score', 0):.2f})"
                            for c in candidates[:2]
                        ]
                        combined = " ~ ".join(labels_with_scores)
                        return f'<span class="trace-step current ambiguous">{combined}</span>'
                elif isinstance(trace, list):
                    steps_html = []
                    for i, (label, score) in enumerate(trace):
                        is_current = i == len(trace) - 1
                        css_class = "trace-step current" if is_current else "trace-step"
                        score_html = f'<span class="trace-score">({score:.2f})</span>'
                        steps_html.append(f'<span class="{css_class}">{label}{score_html}</span>')
                        if not is_current:
                            steps_html.append('<span class="trace-arrow">-></span>')
                    return " ".join(steps_html)
        except (json.JSONDecodeError, TypeError, ValueError):
            pass

    path = cluster_annotation.get("assigned_path", "")
    if path and isinstance(path, str):
        parts = [p.strip() for p in path.split("/")]
        steps_html = []
        for i, part in enumerate(parts):
            is_current = i == len(parts) - 1
            css_class = "trace-step current" if is_current else "trace-step"
            steps_html.append(f'<span class="{css_class}">{part}</span>')
            if not is_current:
                steps_html.append('<span class="trace-arrow">-></span>')
        return " ".join(steps_html)

    label = cluster_annotation.get("assigned_label", "Unknown")
    return f'<span class="trace-step current">{label}</span>'


def _render_competing_types_table(
    marker_scores: pd.DataFrame,
    cluster_id: str,
    assigned_path: str = "",
    winner_label: str = "",
    decision_steps: Optional[pd.DataFrame] = None,
    top_n: int = 5,
) -> str:
    """Render table of top competing cell types."""
    if decision_steps is not None and not decision_steps.empty and assigned_path:
        sibling_html = _render_sibling_competition(
            cluster_id, assigned_path, winner_label, decision_steps, top_n
        )
        if sibling_html:
            return sibling_html

    cluster_scores = marker_scores[marker_scores["cluster_id"].astype(str) == str(cluster_id)].copy()
    if cluster_scores.empty:
        return "<p>No marker scores available</p>"

    top_scores = cluster_scores.nlargest(top_n, "score")

    rows = []
    for _, row in top_scores.iterrows():
        label = row.get("label", "Unknown")
        score = row.get("score", 0)
        coverage = row.get("coverage", 0)
        pos_frac = row.get("mean_positive_fraction", 0)

        rows.append(f"""
            <tr>
                <td>{label}</td>
                <td>{score:.2f}</td>
                <td>{coverage:.0%}</td>
                <td>{pos_frac:.0%}</td>
            </tr>
        """)

    return f"""
        <p class="note">Global ranking (non-hierarchical mode)</p>
        <table class="competing-table">
            <tr>
                <th>Cell Type</th>
                <th>Score</th>
                <th>Coverage</th>
                <th>Pos Frac</th>
            </tr>
            {"".join(rows)}
        </table>
    """


def _render_sibling_competition(
    cluster_id: str,
    assigned_path: str,
    winner_label: str,
    decision_steps: pd.DataFrame,
    top_n: int = 5,
) -> Optional[str]:
    """Render sibling competition table from decision_steps."""
    if " / " in str(assigned_path):
        parts = str(assigned_path).split(" / ")
        parent = parts[-2] if len(parts) >= 2 else "ROOT"
    else:
        parent = "ROOT"

    siblings = decision_steps[
        (decision_steps["cluster_id"].astype(str) == str(cluster_id)) &
        (decision_steps["parent_label"] == parent) &
        (decision_steps["child_passed_gate"] == True)
    ].sort_values("child_score", ascending=False)

    if siblings.empty:
        return None

    if len(siblings) >= 2:
        margin = siblings.iloc[0]["child_score"] - siblings.iloc[1]["child_score"]
        if margin >= 0.5:
            margin_class = "high"
            margin_label = "High confidence"
        elif margin >= 0.2:
            margin_class = "medium"
            margin_label = "Moderate confidence"
        else:
            margin_class = "low"
            margin_label = "Low confidence (review)"
        margin_html = f'<span class="margin-badge {margin_class}">{margin_label} (margin: {margin:.2f})</span>'
    else:
        margin_html = '<span class="margin-badge high">No competitors</span>'

    rows = []
    for _, row in siblings.head(top_n).iterrows():
        label = row.get("child_label", "Unknown")
        score = row.get("child_score", 0)
        coverage = row.get("child_coverage", 0)
        pos_frac = row.get("child_pos_frac", 0)

        if label == winner_label:
            row_class = "winner-row"
            selected = "Check"
        else:
            row_class = ""
            selected = ""

        rows.append(f"""
            <tr class="{row_class}">
                <td>{label} {selected}</td>
                <td>{score:.2f}</td>
                <td>{coverage:.0%}</td>
                <td>{pos_frac:.0%}</td>
            </tr>
        """)

    return f"""
        <p class="note">Sibling competition at level: <strong>{parent}</strong> {margin_html}</p>
        <table class="competing-table">
            <tr>
                <th>Cell Type</th>
                <th>Score</th>
                <th>Coverage</th>
                <th>Pos Frac</th>
            </tr>
            {"".join(rows)}
        </table>
    """


def _render_marker_evidence(
    marker_evidence: Optional[pd.DataFrame],
    cluster_id: str,
    assigned_label: str,
    max_markers: int = 8,
) -> str:
    """Render marker evidence mini-table."""
    if marker_evidence is None or marker_evidence.empty:
        return "<p>Marker evidence not available (use --expand-markers to generate)</p>"

    evidence = marker_evidence[
        (marker_evidence["cluster_id"].astype(str) == str(cluster_id)) &
        (marker_evidence["cell_type"] == assigned_label)
    ].copy()

    if evidence.empty:
        return "<p>No marker evidence for assigned cell type</p>"

    positive = evidence[evidence["role"] == "positive"].nlargest(max_markers, "enrichment")
    anti = evidence[evidence["role"] == "anti"].head(3)

    rows = []
    for _, row in positive.iterrows():
        marker = row["marker"]
        enrichment = row.get("enrichment", 0)
        pos_frac = row.get("pos_frac", 0)
        is_on = row.get("is_on", False)
        status = "Check" if is_on else "o"

        rows.append(f"""
            <div class="marker-row">
                <span class="marker-name">{status} {marker}</span>
                <span class="marker-stats">enrich={enrichment:.2f}, pos={pos_frac:.0%}</span>
            </div>
        """)

    if not anti.empty:
        rows.append('<div style="margin-top: 8px; font-weight: 500; color: #c62828;">Anti-markers:</div>')
        for _, row in anti.iterrows():
            marker = row["marker"]
            enrichment = row.get("enrichment", 0)
            pos_frac = row.get("pos_frac", 0)
            rows.append(f"""
                <div class="marker-row">
                    <span class="marker-name">X {marker}</span>
                    <span class="marker-stats">enrich={enrichment:.2f}, pos={pos_frac:.0%}</span>
                </div>
            """)

    return f'<div class="marker-evidence">{"".join(rows)}</div>'


def _render_qc_flags(cluster_annotation: pd.Series) -> str:
    """Render QC flag badges."""
    flags = []

    stop_reason = cluster_annotation.get("stop_reason", "")
    if stop_reason == "no_root_passed":
        flags.append('<span class="flag critical">Unassigned</span>')
    elif stop_reason == "ambiguous_root":
        candidates = cluster_annotation.get("ambiguous_root_candidates", "")
        if candidates:
            flags.append(f'<span class="flag warning">Ambiguous root: {candidates.replace(";", " vs ")}</span>')
        else:
            flags.append('<span class="flag warning">Ambiguous root</span>')
    elif stop_reason == "ambiguous_siblings":
        flags.append('<span class="flag warning">Ambiguous siblings</span>')
    elif stop_reason == "no_child_passed":
        flags.append('<span class="flag warning">No child passed gate</span>')
    elif stop_reason == "leaf_reached":
        flags.append('<span class="flag success">Leaf reached</span>')

    margin_is_infinite = cluster_annotation.get("margin_is_infinite", False)
    if margin_is_infinite:
        flags.append('<span class="flag success">No competition</span>')

    if cluster_annotation.get("stopped_before_leaf", False):
        flags.append('<span class="flag warning">Needs refinement</span>')

    coverage = cluster_annotation.get("coverage", 1.0)
    if coverage < 0.5:
        flags.append('<span class="flag critical">Very low coverage</span>')
    elif coverage < 0.75:
        flags.append('<span class="flag warning">Low coverage</span>')

    score = cluster_annotation.get("assigned_score", 0)
    if score < 0.5:
        flags.append('<span class="flag critical">Very low score</span>')
    elif score < 1.0:
        flags.append('<span class="flag warning">Low score</span>')

    if not flags:
        return '<span class="no-issues">No issues detected</span>'

    return " ".join(flags)


def generate_audit_card_html(
    cluster_id: str,
    cluster_annotation: pd.Series,
    decision_steps: Optional[pd.DataFrame],
    marker_scores: pd.DataFrame,
    marker_evidence: Optional[pd.DataFrame] = None,
) -> str:
    """Generate HTML for a single cluster audit card."""
    assigned_label = cluster_annotation.get("assigned_label", "Unknown")
    assigned_score = cluster_annotation.get("assigned_score", 0)
    assigned_path = cluster_annotation.get("assigned_path", "")
    n_cells = cluster_annotation.get("n_cells", 0)

    if assigned_score >= 2.0:
        conf_class = "high"
        conf_label = "High"
    elif assigned_score >= 1.0:
        conf_class = "medium"
        conf_label = "Medium"
    elif assigned_score >= 0.5:
        conf_class = "low"
        conf_label = "Low"
    else:
        conf_class = "very_low"
        conf_label = "Very Low"

    breadcrumb = _render_decision_breadcrumb(cluster_annotation, decision_steps, cluster_id)
    competing_table = _render_competing_types_table(
        marker_scores,
        cluster_id,
        assigned_path=assigned_path,
        winner_label=assigned_label,
        decision_steps=decision_steps,
    )
    marker_evidence_html = _render_marker_evidence(marker_evidence, cluster_id, assigned_label)
    qc_flags = _render_qc_flags(cluster_annotation)
    score_breakdown = _render_score_breakdown(marker_scores, cluster_id, assigned_label)
    vote_composition = _render_vote_composition(cluster_annotation)
    ambiguous_root_analysis = _render_ambiguous_root_analysis(cluster_annotation)
    confidence_explanation = _render_confidence_explanation(cluster_annotation)
    gate_summary = _render_gate_summary(cluster_id, decision_steps)

    vote_section = ""
    if vote_composition:
        vote_section = f"""
            <div class="section">
                <h4>Vote Composition (Ambiguous)</h4>
                {vote_composition}
            </div>
        """

    ambiguous_root_section = ""
    if ambiguous_root_analysis:
        ambiguous_root_section = f"""
            <div class="section">
                <h4>Ambiguous Root Analysis</h4>
                {ambiguous_root_analysis}
            </div>
        """

    return f"""
    <div class="audit-card" id="cluster-{cluster_id}">
        <div class="card-header">
            <h3>Cluster {cluster_id}: {assigned_label}</h3>
            <div class="badges">
                <span class="badge score">{assigned_score:.2f}</span>
                <span class="badge cells">{n_cells:,} cells</span>
                <span class="badge {conf_class}">{conf_label}</span>
            </div>
        </div>
        <div class="card-body">
            <div class="section">
                <h4>Decision Path</h4>
                <div class="decision-trace">
                    {breadcrumb}
                </div>
            </div>

            {ambiguous_root_section}

            <div class="section">
                <h4>Score Breakdown</h4>
                {score_breakdown}
            </div>

            <div class="section">
                <h4>Confidence Analysis</h4>
                {confidence_explanation}
            </div>

            {vote_section}

            <div class="section">
                <h4>Top Competing Cell Types</h4>
                {competing_table}
            </div>

            <div class="section">
                <h4>Gate Checks</h4>
                {gate_summary}
            </div>

            <div class="section">
                <h4>Marker Evidence</h4>
                {marker_evidence_html}
            </div>

            <div class="section">
                <h4>QC Flags</h4>
                <div class="qc-flags">
                    {qc_flags}
                </div>
            </div>
        </div>
    </div>
    """


def generate_all_audit_cards(
    cluster_annotations: pd.DataFrame,
    decision_steps: Optional[pd.DataFrame],
    marker_scores: pd.DataFrame,
    marker_evidence: Optional[pd.DataFrame],
    output_dir: Path,
    logger: logging.Logger,
    stage_name: str = "Stage H",
) -> Path:
    """Generate standalone HTML file with all audit cards."""
    from datetime import datetime

    if cluster_annotations.empty:
        logger.warning("Empty cluster annotations; generating empty audit cards page.")
        output_path = output_dir / "audit_cards.html"
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(f"""
        <!DOCTYPE html>
        <html><head><title>{stage_name} Audit Cards</title></head>
        <body><h1>No clusters to display</h1></body></html>
        """)
        return output_path

    n_clusters = len(cluster_annotations)

    mean_score = cluster_annotations["assigned_score"].mean() if "assigned_score" in cluster_annotations.columns else 0
    n_high_conf = len(cluster_annotations[cluster_annotations.get("assigned_score", 0) >= 2.0]) if "assigned_score" in cluster_annotations.columns else 0
    n_needs_review = len(cluster_annotations[cluster_annotations.get("stopped_before_leaf", False) == True]) if "stopped_before_leaf" in cluster_annotations.columns else 0

    cards_html = []
    for _, row in cluster_annotations.iterrows():
        cluster_id = str(row.get("cluster_id", "?"))
        card_html = generate_audit_card_html(
            cluster_id=cluster_id,
            cluster_annotation=row,
            decision_steps=decision_steps,
            marker_scores=marker_scores,
            marker_evidence=marker_evidence,
        )
        cards_html.append(card_html)

    html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>{stage_name} Audit Cards</title>
        {AUDIT_CARD_CSS}
    </head>
    <body>
        <div class="header">
            <h1>{stage_name} Audit Cards</h1>
            <p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            <div class="summary-stats">
                <div class="stat-box">
                    <div class="value">{n_clusters}</div>
                    <div class="label">Clusters</div>
                </div>
                <div class="stat-box">
                    <div class="value">{mean_score:.2f}</div>
                    <div class="label">Avg Score</div>
                </div>
                <div class="stat-box">
                    <div class="value">{n_high_conf}</div>
                    <div class="label">High Confidence</div>
                </div>
                <div class="stat-box">
                    <div class="value">{n_needs_review}</div>
                    <div class="label">Needs Review</div>
                </div>
            </div>
        </div>

        <div class="card-grid">
            {"".join(cards_html)}
        </div>
    </body>
    </html>
    """

    output_path = output_dir / "audit_cards.html"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(html)

    logger.info("Generated audit cards for %d clusters: %s", n_clusters, output_path)
    return output_path
