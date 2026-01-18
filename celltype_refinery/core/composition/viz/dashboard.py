"""Interactive HTML dashboard for composition analysis results.

Generates a comprehensive single-page HTML dashboard with:
- Summary metrics cards
- Global composition bar chart
- Diversity scatter plot (Shannon vs Simpson by region)
- Composition by region stacked bars
- Biology metrics by region (configurable)
- Composition heatmap (samples x cell types)
- Regional enrichment table
- Full cell type statistics table

Usage:
    from celltype_refinery.core.composition.viz.dashboard import generate_composition_dashboard

    generate_composition_dashboard(
        output_dir=Path("output/composition"),
        output_path=Path("output/composition/dashboard.html"),
    )
"""

from pathlib import Path
from typing import Dict, List, Optional, Union
import json
import logging

import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


# Default colors for regions (fallback if not specified in config)
DEFAULT_REGION_COLORS = [
    "#3498db",  # blue
    "#e74c3c",  # red
    "#2ecc71",  # green
    "#9b59b6",  # purple
    "#f39c12",  # orange
    "#1abc9c",  # teal
    "#e67e22",  # dark orange
    "#34495e",  # dark blue-gray
    "#95a5a6",  # gray
    "#d35400",  # burnt orange
]


def _sort_regions(regions, region_order: Optional[List[str]] = None) -> List[str]:
    """Sort regions by specified order or alphabetically."""
    if not region_order:
        return sorted(regions)

    order_map = {r.lower(): i for i, r in enumerate(region_order)}

    def get_order(r):
        r_lower = r.lower() if isinstance(r, str) else str(r).lower()
        return order_map.get(r_lower, len(region_order))

    return sorted(regions, key=get_order)


def _get_region_colors(regions: List[str], config_colors: Optional[Dict[str, str]] = None) -> Dict[str, str]:
    """Get color mapping for regions."""
    colors = {}

    # Use config colors if provided
    if config_colors:
        for region in regions:
            r_lower = region.lower()
            if r_lower in config_colors:
                colors[region] = config_colors[r_lower]
            elif region in config_colors:
                colors[region] = config_colors[region]

    # Fill in missing regions with default colors
    color_idx = 0
    for region in regions:
        if region not in colors:
            colors[region] = DEFAULT_REGION_COLORS[color_idx % len(DEFAULT_REGION_COLORS)]
            color_idx += 1

    return colors


def generate_composition_dashboard(
    output_dir: Union[str, Path],
    output_path: Optional[Union[str, Path]] = None,
    html_dir: Optional[Union[str, Path]] = None,
    title: str = "Cell-Type Composition Analysis",
    region_order: Optional[List[str]] = None,
    region_colors: Optional[Dict[str, str]] = None,
) -> Optional[Path]:
    """Generate interactive HTML dashboard from composition CSVs.

    Supports both single-column and multi-column output structures.
    In multi-column mode, generates an index dashboard at root linking
    to individual column dashboards.

    Parameters
    ----------
    output_dir : Path
        Directory containing composition output CSVs
    output_path : Path, optional
        Path to save HTML. If None, saves to output_dir/dashboard.html
        Ignored if html_dir is specified.
    html_dir : Path, optional
        Output directory for HTML dashboards (flat structure).
        If specified, all HTML files are written to this directory:
        - Multi-column: index.html + <column>.html files
        - Single-column: dashboard.html
        This makes the HTML output easy to share.
    title : str
        Dashboard title
    region_order : List[str], optional
        Canonical region order for consistent display
    region_colors : Dict[str, str], optional
        Color mapping for regions (region name -> hex color)

    Returns
    -------
    Optional[Path]
        Path to saved HTML file, or None if plotly not available
    """
    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
    except ImportError:
        logger.warning("plotly not available, skipping dashboard generation")
        return None

    output_dir = Path(output_dir)

    # Handle html_dir for flat output structure
    if html_dir is not None:
        html_dir = Path(html_dir)
        html_dir.mkdir(parents=True, exist_ok=True)

    # Determine output path
    if html_dir is not None:
        output_path = html_dir / "dashboard.html"
    elif output_path is None:
        output_path = output_dir / "dashboard.html"
    else:
        output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Check for multi-column mode
    multi_summary_path = output_dir / "multi_column_summary.json"
    if multi_summary_path.exists():
        logger.info("Detected multi-column mode, generating dashboards for each column")
        return _generate_multi_column_dashboard(
            output_dir, output_path, html_dir, title, region_order, region_colors
        )

    # Single-column mode
    # Load data files
    data = _load_composition_data(output_dir)
    if data is None:
        logger.error("Failed to load composition data")
        return None

    # Store region config in data for use by build functions
    data["_region_order"] = region_order
    data["_region_colors"] = region_colors

    # Build dashboard HTML
    html_parts = []

    # Header and styles
    html_parts.append(_build_header(title))

    # Summary metrics row
    html_parts.append(_build_summary_metrics(data))

    # Row 1: Global composition + Diversity scatter
    html_parts.append(_build_row_composition_diversity(data))

    # Row 2: Composition by region + biology
    html_parts.append(_build_row_region_biology(data))

    # Row 3: Composition heatmap (full width)
    html_parts.append(_build_row_heatmap(data))

    # Row 4: Regional enrichment table
    html_parts.append(_build_row_enrichment(data))

    # Row 5: Cell type statistics table
    html_parts.append(_build_row_cell_type_table(data))

    # Footer
    html_parts.append(_build_footer(data))

    # Write HTML
    html_content = "\n".join(html_parts)
    output_path.write_text(html_content)

    logger.info(f"Saved composition dashboard to {output_path}")
    return output_path


def _generate_multi_column_dashboard(
    output_dir: Path,
    output_path: Path,
    html_dir: Optional[Path],
    title: str,
    region_order: Optional[List[str]] = None,
    region_colors: Optional[Dict[str, str]] = None,
) -> Optional[Path]:
    """Generate dashboards for multi-column mode.

    Creates:
    - Individual dashboard for each column in subdirectories (or flat in html_dir)
    - Index dashboard at root with comparison and links

    Parameters
    ----------
    output_dir : Path
        Root output directory containing column subdirectories
    output_path : Path
        Path for index dashboard
    html_dir : Path, optional
        If provided, write all HTML files to this directory (flat structure):
        - index.html for the comparison dashboard
        - <column>.html for each individual dashboard
    title : str
        Dashboard title
    region_order : List[str], optional
        Canonical region order
    region_colors : Dict[str, str], optional
        Color mapping for regions

    Returns
    -------
    Optional[Path]
        Path to index dashboard
    """
    import plotly.graph_objects as go

    # Load multi-column summary
    with open(output_dir / "multi_column_summary.json") as f:
        multi_summary = json.load(f)

    columns_processed = multi_summary.get("columns_processed", [])
    comparison = multi_summary.get("comparison", {})

    if not columns_processed:
        logger.error("No columns processed in multi-column summary")
        return None

    # Determine if using flat output structure
    flat_output = html_dir is not None

    # Generate individual dashboards for each column
    column_dashboards = {}
    for col in columns_processed:
        col_dir = output_dir / col
        if col_dir.exists():
            # Determine output path for this column's dashboard
            if flat_output:
                col_dashboard_path = html_dir / f"{col}.html"
            else:
                col_dashboard_path = col_dir / "dashboard.html"

            col_title = f"{title} - {col}"

            # Load column data and generate dashboard
            data = _load_composition_data(col_dir)
            if data is not None:
                # Store region config
                data["_region_order"] = region_order
                data["_region_colors"] = region_colors

                # Generate single-column dashboard for this column
                html_parts = []
                html_parts.append(_build_header(col_title))
                html_parts.append(_build_summary_metrics(data))
                html_parts.append(_build_row_composition_diversity(data))
                html_parts.append(_build_row_region_biology(data))
                html_parts.append(_build_row_heatmap(data))
                html_parts.append(_build_row_enrichment(data))
                html_parts.append(_build_row_cell_type_table(data))
                html_parts.append(_build_footer(data))

                html_content = "\n".join(html_parts)
                col_dashboard_path.write_text(html_content)
                column_dashboards[col] = col_dashboard_path
                logger.info(f"Generated dashboard for {col}")

    # Determine index dashboard output path
    if flat_output:
        index_path = html_dir / "index.html"
    else:
        index_path = output_path

    # Generate index dashboard with comparison
    html_parts = []
    html_parts.append(_build_multi_column_header(title, multi_summary))
    html_parts.append(_build_multi_column_comparison(comparison, columns_processed))
    html_parts.append(_build_multi_column_links(columns_processed, column_dashboards, flat_output))
    html_parts.append(_build_multi_column_footer(multi_summary))

    html_content = "\n".join(html_parts)
    index_path.write_text(html_content)

    logger.info(f"Saved multi-column index dashboard to {index_path}")
    return index_path


def _build_multi_column_header(title: str, multi_summary: Dict) -> str:
    """Build header for multi-column index dashboard."""
    n_columns = len(multi_summary.get("columns_processed", []))
    n_skipped = len(multi_summary.get("columns_skipped", []))
    total_time = multi_summary.get("total_execution_time_seconds", 0)

    return f"""
<!DOCTYPE html>
<html>
<head>
    <title>{title}</title>
    <style>
        * {{
            box-sizing: border-box;
        }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
            margin: 0;
            padding: 20px;
            background: #f8f9fa;
            color: #333;
        }}
        .dashboard {{
            max-width: 1400px;
            margin: 0 auto;
        }}
        .header {{
            background: linear-gradient(135deg, #8e44ad 0%, #3498db 100%);
            color: white;
            padding: 25px 35px;
            border-radius: 12px;
            margin-bottom: 25px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }}
        .header h1 {{
            margin: 0;
            font-size: 28px;
            font-weight: 600;
        }}
        .header .subtitle {{
            opacity: 0.9;
            font-size: 14px;
            margin-top: 8px;
        }}
        .header .stats {{
            display: flex;
            gap: 30px;
            margin-top: 15px;
            font-size: 13px;
            opacity: 0.9;
        }}
        .row {{
            display: flex;
            gap: 20px;
            margin-bottom: 20px;
            flex-wrap: wrap;
        }}
        .card {{
            background: white;
            border-radius: 12px;
            padding: 20px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.08);
            flex: 1;
            min-width: 300px;
        }}
        .card.full-width {{
            flex: none;
            width: 100%;
        }}
        .card h3 {{
            margin: 0 0 15px 0;
            color: #2c3e50;
            font-size: 16px;
            font-weight: 600;
            border-bottom: 2px solid #8e44ad;
            padding-bottom: 10px;
        }}
        .chart-container {{
            width: 100%;
            min-height: 350px;
            overflow: hidden;
        }}
        .column-links {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
            gap: 15px;
        }}
        .column-link {{
            display: block;
            padding: 20px;
            background: linear-gradient(135deg, #3498db 0%, #2980b9 100%);
            color: white;
            text-decoration: none;
            border-radius: 10px;
            transition: transform 0.2s, box-shadow 0.2s;
        }}
        .column-link:hover {{
            transform: translateY(-2px);
            box-shadow: 0 4px 12px rgba(52, 152, 219, 0.4);
        }}
        .column-link .col-name {{
            font-size: 16px;
            font-weight: 600;
            margin-bottom: 8px;
        }}
        .column-link .col-stats {{
            font-size: 12px;
            opacity: 0.9;
        }}
        .footer {{
            text-align: center;
            color: #95a5a6;
            font-size: 12px;
            margin-top: 30px;
            padding: 20px;
        }}
    </style>
</head>
<body>
    <div class="dashboard">
        <div class="header">
            <h1>{title}</h1>
            <div class="subtitle">Multi-column comparison across {n_columns} cell-type annotation schemes</div>
            <div class="stats">
                <span>Columns processed: {n_columns}</span>
                <span>Columns skipped: {n_skipped}</span>
                <span>Total time: {total_time:.1f}s</span>
            </div>
        </div>
"""


def _build_multi_column_comparison(comparison: Dict, columns: list) -> str:
    """Build comparison charts for multi-column dashboard."""
    import plotly.graph_objects as go

    # Extract comparison metrics
    n_cell_types = comparison.get("n_cell_types", {})
    mean_shannon = comparison.get("mean_shannon_entropy", {})
    mean_simpson = comparison.get("mean_simpson_index", {})
    n_enrichments = comparison.get("n_significant_enrichments", {})

    # Build bar chart for cell type counts
    bar_fig = go.Figure(data=[
        go.Bar(
            x=[col.replace("cell_type_", "") for col in columns],
            y=[n_cell_types.get(col, 0) for col in columns],
            marker=dict(color=["#3498db", "#e74c3c", "#2ecc71", "#9b59b6"][:len(columns)]),
            text=[n_cell_types.get(col, 0) for col in columns],
            textposition="outside",
        )
    ])
    bar_fig.update_layout(
        margin=dict(l=50, r=20, t=50, b=50),
        height=350,
        autosize=True,
        xaxis=dict(title="Annotation Column", automargin=True),
        yaxis=dict(title="Number of Cell Types", automargin=True),
        plot_bgcolor="white",
    )
    bar_html = bar_fig.to_html(full_html=False, include_plotlyjs=True, config={'responsive': True})

    # Build diversity comparison
    diversity_fig = go.Figure()
    diversity_fig.add_trace(go.Bar(
        name="Shannon Entropy",
        x=[col.replace("cell_type_", "") for col in columns],
        y=[mean_shannon.get(col, 0) for col in columns],
        marker=dict(color="#3498db"),
    ))
    diversity_fig.add_trace(go.Bar(
        name="Simpson Index",
        x=[col.replace("cell_type_", "") for col in columns],
        y=[mean_simpson.get(col, 0) for col in columns],
        marker=dict(color="#2ecc71"),
    ))
    diversity_fig.update_layout(
        barmode="group",
        margin=dict(l=50, r=20, t=30, b=50),
        height=350,
        autosize=True,
        xaxis=dict(title="Annotation Column", automargin=True),
        yaxis=dict(title="Diversity Metric", automargin=True),
        plot_bgcolor="white",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )
    diversity_html = diversity_fig.to_html(full_html=False, include_plotlyjs=False, config={'responsive': True})

    return f"""
        <div class="row">
            <div class="card">
                <h3>Cell Type Counts by Column</h3>
                <div class="chart-container">
                    {bar_html}
                </div>
            </div>
            <div class="card">
                <h3>Diversity Metrics Comparison</h3>
                <div class="chart-container">
                    {diversity_html}
                </div>
            </div>
        </div>
"""


def _build_multi_column_links(columns: list, dashboards: Dict, flat_output: bool = False) -> str:
    """Build links to individual column dashboards.

    Parameters
    ----------
    columns : list
        List of column names
    dashboards : Dict
        Mapping of column names to dashboard paths
    flat_output : bool
        If True, use flat paths (e.g., cell_type_broad.html)
        If False, use nested paths (e.g., cell_type_broad/dashboard.html)
    """
    links_html = []

    for col in columns:
        dashboard_path = dashboards.get(col)
        if dashboard_path:
            # Use relative path based on output structure
            if flat_output:
                rel_path = f"{col}.html"
            else:
                rel_path = f"{col}/dashboard.html"
            col_display = col.replace("cell_type_", "").replace("_", " ").title()

            links_html.append(f"""
                <a href="{rel_path}" class="column-link">
                    <div class="col-name">{col_display}</div>
                    <div class="col-stats">Click to view detailed dashboard</div>
                </a>
            """)

    return f"""
        <div class="row">
            <div class="card full-width">
                <h3>Individual Column Dashboards</h3>
                <div class="column-links">
                    {''.join(links_html)}
                </div>
            </div>
        </div>
"""


def _build_multi_column_footer(multi_summary: Dict) -> str:
    """Build footer for multi-column dashboard."""
    from datetime import datetime

    timestamp = multi_summary.get("timestamp", datetime.now().isoformat())
    version = multi_summary.get("version", "1.1.0")

    try:
        dt = datetime.fromisoformat(timestamp)
        formatted_time = dt.strftime("%Y-%m-%d %H:%M:%S")
    except (ValueError, TypeError):
        formatted_time = timestamp

    return f"""
        <div class="footer">
            <p>Generated: {formatted_time} | Module version: {version}</p>
            <p>Multi-Column Composition Analysis | celltype_refinery.core.composition</p>
        </div>
    </div>
    <script>
        // Force Plotly charts to resize after page fully loads
        window.addEventListener('load', function() {{
            setTimeout(function() {{
                window.dispatchEvent(new Event('resize'));
            }}, 100);
        }});
    </script>
</body>
</html>
"""


def _load_composition_data(output_dir: Path) -> Optional[Dict]:
    """Load all composition CSV files."""
    data = {}

    # Required files
    required_files = {
        "summary": "composition_summary.json",
        "global": "composition_global.csv",
        "by_sample": "composition_by_sample.csv",
        "by_region": "composition_by_region.csv",
        "diversity": "diversity_by_sample.csv",
    }

    # Optional files (biology files are tissue-agnostic now)
    optional_files = {
        "by_donor": "composition_by_donor.csv",
        "wide": "composition_wide.csv",
        "diversity_summary": "diversity_summary.csv",
        "biology_sample": "biology_by_sample.csv",
        "biology_region": "biology_by_region.csv",
        "enrichment": "regional_enrichment.csv",
    }

    # Load required files
    for key, filename in required_files.items():
        filepath = output_dir / filename
        if not filepath.exists():
            logger.error(f"Required file not found: {filepath}")
            return None

        if filename.endswith(".json"):
            with open(filepath) as f:
                data[key] = json.load(f)
        else:
            data[key] = pd.read_csv(filepath)

    # Load optional files
    for key, filename in optional_files.items():
        filepath = output_dir / filename
        if filepath.exists():
            data[key] = pd.read_csv(filepath)
        else:
            data[key] = None

    return data


def _build_header(title: str) -> str:
    """Build HTML header with styles."""
    return f"""
<!DOCTYPE html>
<html>
<head>
    <title>{title}</title>
    <style>
        * {{
            box-sizing: border-box;
        }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
            margin: 0;
            padding: 20px;
            background: #f8f9fa;
            color: #333;
        }}
        .dashboard {{
            max-width: 1600px;
            margin: 0 auto;
        }}
        .header {{
            background: linear-gradient(135deg, #2c3e50 0%, #3498db 100%);
            color: white;
            padding: 25px 35px;
            border-radius: 12px;
            margin-bottom: 25px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }}
        .header h1 {{
            margin: 0;
            font-size: 28px;
            font-weight: 600;
        }}
        .header .subtitle {{
            opacity: 0.9;
            font-size: 14px;
            margin-top: 8px;
        }}
        .row {{
            display: flex;
            gap: 20px;
            margin-bottom: 20px;
            flex-wrap: wrap;
        }}
        .card {{
            background: white;
            border-radius: 12px;
            padding: 20px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.08);
            flex: 1;
            min-width: 300px;
        }}
        .card.full-width {{
            flex: none;
            width: 100%;
        }}
        .card h3 {{
            margin: 0 0 15px 0;
            color: #2c3e50;
            font-size: 16px;
            font-weight: 600;
            border-bottom: 2px solid #3498db;
            padding-bottom: 10px;
        }}
        .metrics-grid {{
            display: grid;
            grid-template-columns: repeat(6, 1fr);
            gap: 15px;
        }}
        @media (max-width: 1200px) {{
            .metrics-grid {{
                grid-template-columns: repeat(3, 1fr);
            }}
        }}
        @media (max-width: 768px) {{
            .metrics-grid {{
                grid-template-columns: repeat(2, 1fr);
            }}
        }}
        .metric {{
            text-align: center;
            padding: 15px 10px;
            background: #f8f9fa;
            border-radius: 8px;
        }}
        .metric .value {{
            font-size: 28px;
            font-weight: 700;
            color: #2c3e50;
        }}
        .metric .label {{
            font-size: 11px;
            color: #7f8c8d;
            margin-top: 5px;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }}
        .chart-container {{
            width: 100%;
            min-height: 350px;
            overflow: hidden;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            font-size: 13px;
        }}
        th, td {{
            padding: 10px 12px;
            text-align: left;
            border-bottom: 1px solid #ecf0f1;
        }}
        th {{
            background: #f8f9fa;
            font-weight: 600;
            color: #2c3e50;
            position: sticky;
            top: 0;
        }}
        tr:hover {{
            background: #f8f9fa;
        }}
        .sig-yes {{
            color: #27ae60;
            font-weight: 600;
        }}
        .sig-no {{
            color: #95a5a6;
        }}
        .table-scroll {{
            max-height: 400px;
            overflow-y: auto;
        }}
        .footer {{
            text-align: center;
            color: #95a5a6;
            font-size: 12px;
            margin-top: 30px;
            padding: 20px;
        }}
    </style>
</head>
<body>
    <div class="dashboard">
        <div class="header">
            <h1>{title}</h1>
            <div class="subtitle">Comprehensive analysis of cell-type composition, diversity, and regional patterns</div>
        </div>
"""


def _build_summary_metrics(data: Dict) -> str:
    """Build summary metrics cards."""
    summary = data["summary"]
    meta = summary.get("metadata", {})
    diversity = summary.get("diversity", {})

    n_cells = meta.get("n_cells_total", 0)
    n_types = meta.get("n_cell_types", 0)
    n_samples = meta.get("n_samples", 0)
    n_donors = meta.get("n_donors", 0)
    n_regions = meta.get("n_regions", 0)
    mean_shannon = diversity.get("mean_shannon", 0)

    return f"""
        <div class="row">
            <div class="card full-width">
                <h3>Summary Metrics</h3>
                <div class="metrics-grid">
                    <div class="metric">
                        <div class="value">{n_cells:,}</div>
                        <div class="label">Total Cells</div>
                    </div>
                    <div class="metric">
                        <div class="value">{n_types}</div>
                        <div class="label">Cell Types</div>
                    </div>
                    <div class="metric">
                        <div class="value">{n_samples}</div>
                        <div class="label">Samples</div>
                    </div>
                    <div class="metric">
                        <div class="value">{n_donors}</div>
                        <div class="label">Donors</div>
                    </div>
                    <div class="metric">
                        <div class="value">{n_regions}</div>
                        <div class="label">Regions</div>
                    </div>
                    <div class="metric">
                        <div class="value">{mean_shannon:.2f}</div>
                        <div class="label">Mean Shannon H'</div>
                    </div>
                </div>
            </div>
        </div>
"""


def _build_row_composition_diversity(data: Dict) -> str:
    """Build row with global composition and diversity scatter."""
    import plotly.graph_objects as go

    # Global composition bar chart
    df_global = data["global"].head(15).copy()

    bar_fig = go.Figure(data=[go.Bar(
        x=[float(v) for v in df_global["global_pct"].values],
        y=[str(v) for v in df_global["cell_type"].values],
        orientation="h",
        marker=dict(color="#3498db"),
        text=[f"{v:.1f}%" for v in df_global["global_pct"].values],
        textposition="outside",
        hovertemplate="%{y}<br>%{x:.2f}% (%{customdata:,} cells)<extra></extra>",
        customdata=[int(v) for v in df_global["total_count"].values],
    )])
    bar_fig.update_layout(
        margin=dict(l=20, r=60, t=30, b=40),
        height=400,
        yaxis=dict(categoryorder="total ascending", title=""),
        xaxis=dict(title="Proportion (%)", tickformat=".1f"),
        plot_bgcolor="white",
    )
    bar_html = bar_fig.to_html(full_html=False, include_plotlyjs=True)

    # Diversity scatter plot
    df_div = data["diversity"]

    # Get region order and colors
    region_order = data.get("_region_order")
    region_colors_config = data.get("_region_colors")

    scatter_traces = []

    # Check if region column exists
    if "region" in df_div.columns:
        regions = _sort_regions(df_div["region"].unique(), region_order)
        region_colors = _get_region_colors(regions, region_colors_config)

        for region in regions:
            mask = df_div["region"] == region
            df_r = df_div[mask]
            scatter_traces.append(go.Scatter(
                x=[float(v) for v in df_r["shannon_entropy"].values],
                y=[float(v) for v in df_r["simpson_index"].values],
                mode="markers",
                name=region.title(),
                marker=dict(
                    size=10,
                    color=region_colors.get(region, "#95a5a6"),
                    opacity=0.8,
                ),
                text=[str(v) for v in df_r["sample_id"].values],
                hovertemplate="%{text}<br>Shannon: %{x:.3f}<br>Simpson: %{y:.3f}<extra></extra>",
            ))
    else:
        # No region column - show all samples with single color
        scatter_traces.append(go.Scatter(
            x=[float(v) for v in df_div["shannon_entropy"].values],
            y=[float(v) for v in df_div["simpson_index"].values],
            mode="markers",
            name="All Samples",
            marker=dict(
                size=10,
                color="#3498db",
                opacity=0.8,
            ),
            text=[str(v) for v in df_div["sample_id"].values],
            hovertemplate="%{text}<br>Shannon: %{x:.3f}<br>Simpson: %{y:.3f}<extra></extra>",
        ))

    scatter_fig = go.Figure(data=scatter_traces)
    scatter_fig.update_layout(
        margin=dict(l=50, r=20, t=30, b=50),
        height=400,
        xaxis=dict(title="Shannon Entropy (H')", gridcolor="#ecf0f1"),
        yaxis=dict(title="Simpson Index (1-D)", gridcolor="#ecf0f1"),
        plot_bgcolor="white",
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1,
        ),
    )
    scatter_html = scatter_fig.to_html(full_html=False, include_plotlyjs=False)

    return f"""
        <div class="row">
            <div class="card">
                <h3>Global Composition (Top 15)</h3>
                <div class="chart-container">
                    {bar_html}
                </div>
            </div>
            <div class="card">
                <h3>Diversity by Sample</h3>
                <div class="chart-container">
                    {scatter_html}
                </div>
            </div>
        </div>
"""


def _build_row_region_biology(data: Dict) -> str:
    """Build row with composition by region and biology metrics."""
    import plotly.graph_objects as go

    # Get region order
    region_order = data.get("_region_order")

    # Composition by region - stacked bar
    df_region = data["by_region"]

    # Check if region data is available
    if df_region is None or len(df_region) == 0 or "region" not in df_region.columns:
        return """
        <div class="row">
            <div class="card">
                <h3>Composition by Region</h3>
                <p style='color: #95a5a6; text-align: center; padding: 50px;'>Regional data not available</p>
            </div>
            <div class="card">
                <h3>Biology Metrics by Region</h3>
                <p style='color: #95a5a6; text-align: center; padding: 50px;'>Biology metrics not available</p>
            </div>
        </div>
"""

    # Get top 10 cell types by global count
    top_types = data["global"].head(10)["cell_type"].tolist()

    # Filter and pivot
    df_filtered = df_region[df_region["cell_type"].isin(top_types)]

    regions = _sort_regions(df_region["region"].unique(), region_order)

    # Build stacked bar traces
    colors = [
        "#3498db", "#e74c3c", "#2ecc71", "#9b59b6", "#f39c12",
        "#1abc9c", "#e67e22", "#34495e", "#95a5a6", "#d35400"
    ]

    stacked_traces = []
    for i, ct in enumerate(top_types):
        ct_data = df_filtered[df_filtered["cell_type"] == ct]
        # Align with regions
        values = []
        for r in regions:
            row = ct_data[ct_data["region"] == r]
            if len(row) > 0:
                values.append(float(row["mean_pct"].values[0]))
            else:
                values.append(0)

        stacked_traces.append(go.Bar(
            name=ct,
            x=[r.title() for r in regions],
            y=values,
            marker=dict(color=colors[i % len(colors)]),
            hovertemplate=f"{ct}<br>%{{x}}: %{{y:.1f}}%<extra></extra>",
        ))

    stacked_fig = go.Figure(data=stacked_traces)
    stacked_fig.update_layout(
        barmode="stack",
        margin=dict(l=50, r=20, t=30, b=50),
        height=400,
        xaxis=dict(title="Region"),
        yaxis=dict(title="Mean Proportion (%)", tickformat=".0f"),
        plot_bgcolor="white",
        legend=dict(
            orientation="v",
            yanchor="top",
            y=1,
            xanchor="left",
            x=1.02,
            font=dict(size=10),
        ),
    )
    stacked_html = stacked_fig.to_html(full_html=False, include_plotlyjs=False)

    # Biology metrics (tissue-agnostic)
    biology_html = ""
    if data.get("biology_region") is not None:
        df_bio = data["biology_region"]

        # Detect the region/group column name (support both 'region' and 'group_id')
        bio_region_col = None
        for col_name in ["region", "group_id"]:
            if col_name in df_bio.columns:
                bio_region_col = col_name
                break

        if bio_region_col is None:
            # No region column found, skip biology metrics
            pass
        else:
            # Auto-detect available metrics from columns
            available_metrics = []
            for col in df_bio.columns:
                if col.endswith("_pct") and not col.endswith("_pct_valid"):
                    label = col.replace("_pct", "").replace("_", " ").title() + " %"
                    available_metrics.append((col, label))
                elif col.endswith("_ratio") and not col.endswith("_ratio_valid"):
                    label = col.replace("_ratio", "").replace("_", " ").title() + " Ratio"
                    available_metrics.append((col, label))

            if available_metrics:
                bio_traces = []
                bar_colors = ["#3498db", "#2ecc71", "#e74c3c", "#9b59b6", "#f39c12", "#1abc9c"]

                for i, (col, label) in enumerate(available_metrics[:6]):
                    if col in df_bio.columns:
                        bio_traces.append(go.Bar(
                            name=label,
                            x=[str(r).title() for r in df_bio[bio_region_col].values],
                            y=[float(v) for v in df_bio[col].values],
                            marker=dict(color=bar_colors[i % len(bar_colors)]),
                            hovertemplate=f"{label}<br>%{{x}}: %{{y:.2f}}<extra></extra>",
                        ))

            if bio_traces:
                bio_fig = go.Figure(data=bio_traces)
                bio_fig.update_layout(
                    barmode="group",
                    margin=dict(l=50, r=20, t=30, b=50),
                    height=400,
                    xaxis=dict(title="Region"),
                    yaxis=dict(title="Value"),
                    plot_bgcolor="white",
                    legend=dict(
                        orientation="h",
                        yanchor="bottom",
                        y=1.02,
                        xanchor="right",
                        x=1,
                    ),
                )
                biology_html = bio_fig.to_html(full_html=False, include_plotlyjs=False)

    if not biology_html:
        biology_html = "<p style='color: #95a5a6; text-align: center; padding: 50px;'>Biology metrics not available</p>"

    return f"""
        <div class="row">
            <div class="card">
                <h3>Composition by Region (Top 10 Types)</h3>
                <div class="chart-container">
                    {stacked_html}
                </div>
            </div>
            <div class="card">
                <h3>Biology Metrics by Region</h3>
                <div class="chart-container">
                    {biology_html}
                </div>
            </div>
        </div>
"""


def _build_row_heatmap(data: Dict) -> str:
    """Build composition heatmap row with hierarchical clustering and region colors."""
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    if data.get("wide") is None:
        return ""

    df_wide = data["wide"].copy()

    # Get sample_id column
    if "sample_id" not in df_wide.columns:
        return ""

    # Limit to top 20 cell types by global count
    top_types = data["global"].head(20)["cell_type"].tolist()
    cell_types = [c for c in top_types if c in df_wide.columns]

    if len(cell_types) == 0:
        return ""

    # Get region for each sample (if available)
    df_div = data["diversity"]
    has_region = "region" in df_div.columns
    if has_region:
        sample_to_region = dict(zip(df_div["sample_id"], df_div["region"]))
        df_wide["region"] = df_wide["sample_id"].map(sample_to_region)
    else:
        df_wide["region"] = "Unknown"  # Default region when not available

    # Perform hierarchical clustering on rows
    try:
        from scipy.cluster.hierarchy import linkage, leaves_list
        from scipy.spatial.distance import pdist

        # Get the data matrix for clustering
        data_matrix = df_wide[cell_types].values
        # Handle any NaN values
        data_matrix = np.nan_to_num(data_matrix, nan=0.0)

        if len(data_matrix) > 1:
            # Compute linkage
            distances = pdist(data_matrix, metric='euclidean')
            Z = linkage(distances, method='ward')
            # Get row order from clustering
            row_order = leaves_list(Z)
            df_wide = df_wide.iloc[row_order].reset_index(drop=True)
    except ImportError:
        # If scipy not available, skip clustering
        pass

    sample_ids = df_wide["sample_id"].tolist()
    regions = df_wide["region"].tolist()

    # Build z matrix for heatmap
    z_values = df_wide[cell_types].values.tolist()

    # Get region colors
    region_colors_config = data.get("_region_colors")
    unique_regions = list(set(regions))
    region_color_map = _get_region_colors(unique_regions, region_colors_config)

    region_colors = [region_color_map.get(r, "#95a5a6") for r in regions]

    # Create subplots: region annotation column + main heatmap
    fig = make_subplots(
        rows=1, cols=2,
        column_widths=[0.03, 0.97],
        horizontal_spacing=0.005,
        shared_yaxes=True,
    )

    # Region color annotation (left column)
    # Create a numeric mapping for regions to use as z values
    region_to_num = {r: i for i, r in enumerate(unique_regions)}
    region_z = [[region_to_num.get(r, len(unique_regions))] for r in regions]

    # Custom colorscale for regions
    region_colorscale = []
    for i, r in enumerate(unique_regions):
        pos = i / max(len(unique_regions) - 1, 1)
        region_colorscale.append([pos, region_color_map[r]])
    if len(region_colorscale) == 1:
        region_colorscale = [[0, region_colorscale[0][1]], [1, region_colorscale[0][1]]]

    fig.add_trace(
        go.Heatmap(
            z=region_z,
            x=["Region"],
            y=sample_ids,
            colorscale=region_colorscale,
            showscale=False,
            hovertemplate="%{y}<br>Region: %{customdata}<extra></extra>",
            customdata=[[r] for r in regions],
        ),
        row=1, col=1
    )

    # Main heatmap (right column)
    fig.add_trace(
        go.Heatmap(
            z=z_values,
            x=cell_types,
            y=sample_ids,
            colorscale="YlOrRd",
            colorbar=dict(
                title="Proportion (%)",
                x=1.02,
                len=0.5,
                y=0.75,
            ),
            hovertemplate="%{y}<br>%{x}: %{z:.2f}%<extra></extra>",
        ),
        row=1, col=2
    )

    # Update layout
    n_samples = len(sample_ids)
    fig_height = max(500, min(800, n_samples * 18))

    fig.update_layout(
        height=fig_height,
        margin=dict(l=100, r=80, t=30, b=150),
        plot_bgcolor="white",
        yaxis=dict(
            title="",
            tickfont=dict(size=9),
            automargin=True,
        ),
        yaxis2=dict(
            showticklabels=False,
        ),
        xaxis=dict(
            title="",
            showticklabels=False,
        ),
        xaxis2=dict(
            title="",
            tickangle=45,
            tickfont=dict(size=9),
            automargin=True,
        ),
    )

    # Add region legend as annotations at the bottom
    region_order = data.get("_region_order")
    sorted_regions = _sort_regions(unique_regions, region_order)

    legend_y = -0.22
    legend_x_start = 0.3
    legend_spacing = 0.15

    annotations = []
    for i, region in enumerate(sorted_regions):
        x_pos = legend_x_start + i * legend_spacing
        # Color box (using a shape would be better, but annotations work)
        annotations.append(dict(
            x=x_pos - 0.02,
            y=legend_y,
            xref="paper",
            yref="paper",
            text=f"<b>■</b>",
            showarrow=False,
            font=dict(size=16, color=region_color_map[region]),
        ))
        annotations.append(dict(
            x=x_pos + 0.02,
            y=legend_y,
            xref="paper",
            yref="paper",
            text=region,
            showarrow=False,
            font=dict(size=11),
            xanchor="left",
        ))

    # Add "Region" title for legend
    annotations.append(dict(
        x=legend_x_start - 0.08,
        y=legend_y,
        xref="paper",
        yref="paper",
        text="<b>Region:</b>",
        showarrow=False,
        font=dict(size=11),
        xanchor="right",
    ))

    fig.update_layout(annotations=annotations)

    heatmap_html = fig.to_html(full_html=False, include_plotlyjs=False)

    return f"""
        <div class="row">
            <div class="card full-width">
                <h3>Composition Heatmap</h3>
                <div class="chart-container" style="min-height: {fig_height}px;">
                    {heatmap_html}
                </div>
            </div>
        </div>
"""


def _build_row_enrichment(data: Dict) -> str:
    """Build regional enrichment table row."""
    if data.get("enrichment") is None:
        return ""

    df_enrich = data["enrichment"]

    if len(df_enrich) == 0:
        return ""

    # Sort by fold_change descending, show top 30
    df_sorted = df_enrich.sort_values("fold_change", ascending=False).head(30)

    # Build table HTML
    table_rows = []
    for _, row in df_sorted.iterrows():
        sig_class = "sig-yes" if row.get("significant", False) else "sig-no"
        sig_text = "Yes" if row.get("significant", False) else "No"

        pval = row.get("pvalue", 1)
        pval_adj = row.get("pvalue_adj", 1)

        table_rows.append(f"""
            <tr>
                <td>{row.get('cell_type', '')}</td>
                <td>{str(row.get('region', '')).title()}</td>
                <td>{row.get('mean_pct_in_region', 0):.2f}%</td>
                <td>{row.get('mean_pct_outside', 0):.2f}%</td>
                <td><strong>{row.get('fold_change', 0):.2f}x</strong></td>
                <td>{pval:.2e}</td>
                <td>{pval_adj:.2e}</td>
                <td class="{sig_class}">{sig_text}</td>
            </tr>
        """)

    table_html = f"""
        <div class="table-scroll">
            <table>
                <thead>
                    <tr>
                        <th>Cell Type</th>
                        <th>Region</th>
                        <th>Mean % (in)</th>
                        <th>Mean % (out)</th>
                        <th>Fold Change</th>
                        <th>P-value</th>
                        <th>Adj. P-value</th>
                        <th>Significant</th>
                    </tr>
                </thead>
                <tbody>
                    {''.join(table_rows)}
                </tbody>
            </table>
        </div>
    """

    return f"""
        <div class="row">
            <div class="card full-width">
                <h3>Regional Enrichment Analysis (Top 30 by Fold Change)</h3>
                {table_html}
            </div>
        </div>
"""


def _build_row_cell_type_table(data: Dict) -> str:
    """Build cell type statistics table row."""
    df_global = data["global"]

    # Build table HTML
    table_rows = []
    for _, row in df_global.iterrows():
        table_rows.append(f"""
            <tr>
                <td><strong>{row.get('cell_type', '')}</strong></td>
                <td>{int(row.get('total_count', 0)):,}</td>
                <td>{row.get('global_pct', 0):.2f}%</td>
                <td>{row.get('mean_pct', 0):.2f}%</td>
                <td>{row.get('std_pct', 0):.2f}%</td>
                <td>{row.get('cv', 0):.2f}</td>
                <td>{int(row.get('n_samples_present', 0))}</td>
                <td>{row.get('min_pct', 0):.2f}% - {row.get('max_pct', 0):.2f}%</td>
            </tr>
        """)

    table_html = f"""
        <div class="table-scroll">
            <table>
                <thead>
                    <tr>
                        <th>Cell Type</th>
                        <th>Total Count</th>
                        <th>Global %</th>
                        <th>Mean %</th>
                        <th>Std %</th>
                        <th>CV</th>
                        <th>Samples Present</th>
                        <th>Range</th>
                    </tr>
                </thead>
                <tbody>
                    {''.join(table_rows)}
                </tbody>
            </table>
        </div>
    """

    return f"""
        <div class="row">
            <div class="card full-width">
                <h3>Cell Type Statistics</h3>
                {table_html}
            </div>
        </div>
"""


def _build_footer(data: Dict) -> str:
    """Build footer with metadata."""
    from datetime import datetime

    summary = data["summary"]
    config = summary.get("config", {})
    execution = summary.get("execution", {})

    tissue = config.get("tissue", "unknown").replace("_", " ").title()
    time_sec = execution.get("time_seconds", 0)
    timestamp = summary.get("timestamp", datetime.now().isoformat())

    # Parse timestamp
    try:
        dt = datetime.fromisoformat(timestamp)
        formatted_time = dt.strftime("%Y-%m-%d %H:%M:%S")
    except (ValueError, TypeError):
        formatted_time = timestamp

    return f"""
        <div class="footer">
            <p>Tissue: {tissue} | Analysis time: {time_sec:.1f}s | Generated: {formatted_time}</p>
            <p>Composition Analysis Module v1.0 | celltype_refinery.core.composition</p>
        </div>
    </div>
    <script>
        // Force Plotly charts to resize after page fully loads
        window.addEventListener('load', function() {{
            setTimeout(function() {{
                window.dispatchEvent(new Event('resize'));
            }}, 100);
        }});
    </script>
</body>
</html>
"""
