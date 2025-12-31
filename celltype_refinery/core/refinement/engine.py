"""
RefinementEngine: Executes RefinePlan operations on AnnData objects.

This engine provides the execution layer for cell-type annotation refinement.
It takes a RefinePlan and applies operations in deterministic order:
    overrides → merges → subclusters → relabels → rescore

Key design decision: DE test functions are passed via dependency injection
rather than imported directly, allowing decoupling from specific clustering
implementations.
"""

from __future__ import annotations

import logging
import time
import warnings
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd

try:
    import scanpy as sc
except ImportError:
    sc = None

from .plan import (
    RefineOp,
    RefinePlan,
    SubclusterOp,
    MergeOp,
    OverrideOp,
    RelabelOp,
    RescoreOp,
)

# Type alias for DE runner functions
DERunner = Callable[..., Dict[str, List[str]]]


def _sanitize_df_for_h5ad(df: pd.DataFrame) -> pd.DataFrame:
    """Sanitize DataFrame for H5AD serialization.

    H5AD/anndata has strict requirements for DataFrame columns:
    - String columns must not contain mixed types or NaN floats
    - All string columns should have NaN replaced with empty strings

    This function fixes common issues:
    - Converts cluster_id to string
    - Replaces NaN in object/string columns with empty strings
    - Converts float columns that should be strings (e.g., marker_weights with NaN)

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to sanitize

    Returns
    -------
    pd.DataFrame
        Sanitized DataFrame safe for H5AD storage
    """
    df = df.copy()

    # Known string columns that may have NaN (from marker_scores)
    string_columns = [
        "cluster_id",
        "label",
        "path",
        "marker_weights",
        "missing_markers",
        "resolved_markers",
        "missing_anti_markers",
        "resolved_anti_markers",
    ]

    for col in string_columns:
        if col in df.columns:
            # Convert to string and replace NaN/None with empty string
            df[col] = df[col].fillna("").astype(str)

    # For any remaining object columns, also sanitize
    for col in df.columns:
        if df[col].dtype == "object" and col not in string_columns:
            df[col] = df[col].fillna("").astype(str)

    return df


def default_scanpy_de_runner(
    adata: "sc.AnnData",
    cluster_key: str = "cluster_lvl0",
    method: str = "wilcoxon",
    n_genes: int = 20,
    layer: Optional[str] = None,
    **kwargs,
) -> Dict[str, List[str]]:
    """Default DE runner using scanpy.

    This is the default implementation used when no custom DE runner is provided.

    Parameters
    ----------
    adata : sc.AnnData
        AnnData object with cluster annotations
    cluster_key : str
        Column in adata.obs with cluster assignments
    method : str
        DE method (default: "wilcoxon")
    n_genes : int
        Number of top genes per cluster (default: 20)
    layer : str, optional
        Layer to use for DE (None = adata.X)
    **kwargs
        Additional arguments passed to sc.tl.rank_genes_groups

    Returns
    -------
    Dict[str, List[str]]
        Mapping from cluster ID to list of top DE genes
    """
    import scanpy as sc

    # Run DE
    use_layer = layer if layer and layer in adata.layers else None
    sc.tl.rank_genes_groups(
        adata,
        groupby=cluster_key,
        method=method,
        use_raw=False,
        layer=use_layer,
        **kwargs,
    )

    # Extract top genes per cluster
    de_lookup = {}
    groups = adata.uns["rank_genes_groups"]["names"].dtype.names

    for group in groups:
        genes = adata.uns["rank_genes_groups"]["names"][group][:n_genes].tolist()
        de_lookup[str(group)] = genes

    return de_lookup


def default_scanpy_de_within_parent_runner(
    adata: "sc.AnnData",
    subcluster_key: str = "cluster_lvl1",
    parent_key: str = "cluster_lvl0",
    method: str = "wilcoxon",
    n_genes: int = 20,
    n_workers: int = 8,
    layer: Optional[str] = None,
    tie_correct: bool = True,
    logger: Optional[logging.Logger] = None,
    subclustered_parents: Optional[List[str]] = None,
) -> Dict[str, List[str]]:
    """Within-parent DE runner: compare subclusters within each parent cluster.

    This is biologically appropriate for subclustering:
    - "3:0" is compared against "3:1" + "3:2" (not entire dataset)
    - Much faster: ~100K cells vs 3.4M cells
    - Parallelizable: each parent is independent

    Parameters
    ----------
    adata : sc.AnnData
        AnnData with cluster annotations
    subcluster_key : str
        Column with subcluster IDs (e.g., "3:0", "3:1")
    parent_key : str
        Column with parent cluster IDs (unused, kept for API compat)
    method : str
        DE method (default: "wilcoxon")
    n_genes : int
        Number of top DE genes per cluster
    n_workers : int
        Number of parallel workers
    layer : str, optional
        Data layer to use
    tie_correct : bool
        Apply tie correction for Wilcoxon
    logger : logging.Logger
        Logger instance
    subclustered_parents : List[str], optional
        List of parent cluster IDs that were subclustered in this iteration.
        If provided, only these parents will be processed for DE.

    Returns
    -------
    Dict[str, List[str]]
        Mapping from subcluster ID to list of top DE genes
    """
    from concurrent.futures import ThreadPoolExecutor, as_completed
    import threading

    if logger is None:
        logger = logging.getLogger(__name__)

    # Identify parent clusters that have subclusters
    subcluster_ids = adata.obs[subcluster_key].astype(str).unique()
    all_parents_with_subclusters = set()
    for sc_id in subcluster_ids:
        if ":" in sc_id:
            parts = sc_id.split(":")
            parent = ":".join(parts[:-1])
            all_parents_with_subclusters.add(parent)

    # Filter to only subclustered parents if provided
    if subclustered_parents is not None and len(subclustered_parents) > 0:
        subclustered_set = set(str(p) for p in subclustered_parents)
        parents_with_subclusters = all_parents_with_subclusters & subclustered_set
        n_skipped = len(all_parents_with_subclusters) - len(parents_with_subclusters)
        if n_skipped > 0:
            logger.info(
                "Filtered to %d newly-subclustered parents (skipping %d unchanged)",
                len(parents_with_subclusters), n_skipped
            )
    else:
        parents_with_subclusters = all_parents_with_subclusters

    n_parents = len(parents_with_subclusters)
    if n_parents == 0:
        logger.info("No subclusters found (no ':' in cluster IDs), falling back to global DE")
        return default_scanpy_de_runner(
            adata, subcluster_key, method, n_genes, layer
        )

    logger.info(
        "Within-parent DE: %d parent clusters with subclusters, %d workers",
        n_parents, n_workers
    )

    all_results: Dict[str, List[str]] = {}
    results_lock = threading.Lock()
    completed_count = [0]

    def process_parent(parent_id: str) -> Tuple[Dict[str, List[str]], str, int, float]:
        """Run DE for one parent cluster's subclusters."""
        try:
            # Extract cells for this parent
            parent_mask = adata.obs[subcluster_key].astype(str).str.startswith(f"{parent_id}:")
            n_cells = parent_mask.sum()

            if n_cells < 100:
                return {}, parent_id, n_cells, 0.0

            # Create subset AnnData
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=pd.errors.SettingWithCopyWarning)
                warnings.filterwarnings("ignore", message=".*SettingWithCopyWarning.*")
                parent_adata = adata[parent_mask].copy()

            # Check we have at least 2 subclusters
            subclusters_in_parent = parent_adata.obs[subcluster_key].astype(str).unique()
            if len(subclusters_in_parent) < 2:
                return {}, parent_id, n_cells, 0.0

            # Run DE within this parent
            start = time.time()
            sc.tl.rank_genes_groups(
                parent_adata,
                groupby=subcluster_key,
                method=method,
                n_genes=n_genes,
                layer=layer if layer and layer in parent_adata.layers else None,
                use_raw=False,
                tie_correct=tie_correct,
                pts=True,
            )
            elapsed = time.time() - start

            # Extract results
            result = {}
            for sc_id in subclusters_in_parent:
                try:
                    df = sc.get.rank_genes_groups_df(parent_adata, group=str(sc_id))
                    df = df.dropna(subset=["names"])
                    result[str(sc_id)] = df.head(n_genes)["names"].astype(str).tolist()
                except Exception:
                    result[str(sc_id)] = []

            return result, parent_id, n_cells, elapsed

        except Exception as e:
            if logger:
                logger.warning("DE failed for parent %s: %s", parent_id, e)
            return {}, parent_id, 0, 0.0

    # Process parents in parallel
    start_total = time.time()

    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        futures = {executor.submit(process_parent, p): p for p in parents_with_subclusters}

        for future in as_completed(futures):
            result, parent_id, n_cells, elapsed = future.result()
            with results_lock:
                all_results.update(result)
                completed_count[0] += 1
                if elapsed > 0:
                    logger.info(
                        "DE parent %s: %d cells, %.1fs (%d/%d)",
                        parent_id, n_cells, elapsed, completed_count[0], n_parents
                    )

    # Handle non-subclustered clusters
    non_subclustered = [c for c in subcluster_ids if ":" not in c]
    if non_subclustered:
        logger.info("Adding %d non-subclustered clusters (no DE needed)", len(non_subclustered))
        for cluster_id in non_subclustered:
            all_results[str(cluster_id)] = []

    total_elapsed = time.time() - start_total
    logger.info("Within-parent DE complete: %d clusters, %.1f sec total", len(all_results), total_elapsed)

    return all_results


@dataclass
class EngineResult:
    """Result of RefinementEngine execution."""
    success: bool
    n_cells_modified: int
    operations_executed: Dict[str, int]
    subclusters_created: List[str]
    errors: List[str]
    warnings: List[str]
    execution_time_seconds: float


class RefinementEngine:
    """Executes RefinePlan operations on AnnData objects.

    This engine applies refinement operations in a deterministic order:
    1. Overrides: Direct cell-type label assignments
    2. Merges: Combine multiple clusters into one label
    3. Subclusters: Re-cluster specific clusters
    4. Relabels: Instant relabeling without clustering
    5. Rescore: Recompute marker scores for modified clusters

    Key design: DE test functions are provided via dependency injection,
    allowing this engine to work independently of any specific clustering
    implementation.

    Example:
        >>> from celltype_refinery.core.refinement import RefinePlan, RefinementEngine
        >>> plan = RefinePlan()
        >>> plan.add_override("3", "Ciliated_Epithelium", "High FOXJ1")
        >>> plan.add_subcluster("5", resolution=0.3, reason="Heterogeneous")
        >>> engine = RefinementEngine(label_key_out="cell_type_curated")
        >>> adata = engine.execute(adata, plan)
    """

    def __init__(
        self,
        marker_map: Optional[Dict] = None,
        label_key_out: str = "cell_type_curated",
        reason_key_out: str = "curation_reason",
        logger: Optional[logging.Logger] = None,
        # Input adapter for column name normalization
        adapter: Optional["InputAdapter"] = None,
        # DE function injection (key abstraction for tissue-agnostic operation)
        de_runner: Optional[DERunner] = None,
        de_within_parent_runner: Optional[DERunner] = None,
        # Rescoring parameters
        stage_h_scores_path: Optional[Path] = None,
        stage_h_annotations_path: Optional[Path] = None,
        positive_quantile: float = 0.75,
        de_bonus: float = 0.5,
        anti_weight: float = 0.5,
        expand_markers: bool = False,
        skip_clustering: bool = False,
        # DE parameters
        de_layer: str = "X",
        de_tie_correct: bool = True,
        de_workers: int = 8,
        de_within_parent: bool = True,
        # Scoring layer (for marker score computation)
        scoring_layer: str = "batchcorr",
        # Rank-weighted DE bonus parameters
        de_commonness_alpha: float = 0.5,
        de_top_frac: float = 0.2,
        de_min_k: int = 3,
        de_max_k: int = 12,
        # Debug output
        output_dir: Optional[Path] = None,
        # Parallel execution
        n_workers: int = 1,
        scoring_workers: int = 1,
        scoring_batch_size: int = 50,
        gating_workers: int = 1,
    ):
        """Initialize RefinementEngine.

        Parameters
        ----------
        marker_map : Dict, optional
            Marker map for rescoring operations (JSON dict or Path to JSON file)
        label_key_out : str
            Output column name for cell type labels (default: "cell_type_curated")
        reason_key_out : str
            Output column name for assignment reasons (default: "curation_reason")
        logger : logging.Logger, optional
            Logger for tracking operations
        adapter : InputAdapter, optional
            Adapter for input column name normalization (auto-created if not provided).
        de_runner : Callable, optional
            Function to run DE tests. Signature:
            de_runner(adata, cluster_key, method, n_genes, layer, **kwargs) -> Dict[str, List[str]]
            If None, uses default_scanpy_de_runner.
        de_within_parent_runner : Callable, optional
            Function to run within-parent DE tests. If None, uses de_runner.
        stage_h_scores_path : Path, optional
            Path to Stage H marker_scores.csv for smart rescoring
        stage_h_annotations_path : Path, optional
            Path to Stage H cluster_annotations.csv
        positive_quantile : float
            Quantile threshold for marker positivity (default: 0.75)
        de_bonus : float
            Score bonus for differential expression overlap (default: 0.5)
        anti_weight : float
            Weight for anti-marker penalty (default: 0.5)
        expand_markers : bool
            If True, generate per-marker evidence table (default: False)
        skip_clustering : bool
            If True, skip SubclusterOp operations (default: False)
        de_layer : str
            Layer to use for DE analysis (default: 'X')
        de_tie_correct : bool
            Apply tie correction for Wilcoxon test (default: True)
        de_workers : int
            Number of parallel workers for DE computation (default: 8)
        de_within_parent : bool
            If True, run DE within each parent cluster (default: True)
        de_commonness_alpha : float
            Exponent for commonness penalty in DE scoring (default: 0.5)
        de_top_frac : float
            Fraction of panel size for target K in rank-weighted DE (default: 0.2)
        de_min_k : int
            Minimum K for rank-weighted DE top-K selection (default: 3)
        de_max_k : int
            Maximum K for rank-weighted DE top-K selection (default: 12)
        output_dir : Path, optional
            Output directory for debug files
        n_workers : int
            Number of parallel workers for subclustering (default: 1)
        scoring_workers : int
            Number of parallel workers for marker scoring (default: 1)
        scoring_batch_size : int
            Number of clusters per worker batch for scoring (default: 50)
        gating_workers : int
            Number of parallel workers for hierarchical gating (default: 1)
        """
        self.marker_map = marker_map
        self.label_key_out = label_key_out
        self.reason_key_out = reason_key_out
        self.logger = logger or logging.getLogger(__name__)

        # Store adapter (lazy-create if not provided)
        self._adapter = adapter

        # DE function injection - key abstraction for decoupling
        self._de_runner = de_runner or default_scanpy_de_runner
        self._de_within_parent_runner = de_within_parent_runner

        # Rescoring parameters
        self.stage_h_scores_path = Path(stage_h_scores_path) if stage_h_scores_path else None
        self.stage_h_annotations_path = Path(stage_h_annotations_path) if stage_h_annotations_path else None
        self.positive_quantile = positive_quantile
        self.de_bonus = de_bonus
        self.anti_weight = anti_weight
        self.expand_markers = expand_markers
        self.skip_clustering = skip_clustering

        # DE parameters
        self.de_layer = de_layer
        self.de_tie_correct = de_tie_correct
        self.de_workers = de_workers
        self.de_within_parent = de_within_parent
        self.scoring_layer = scoring_layer

        # Rank-weighted DE bonus parameters
        self.de_commonness_alpha = de_commonness_alpha
        self.de_top_frac = de_top_frac
        self.de_min_k = de_min_k
        self.de_max_k = de_max_k

        # Debug output
        self.output_dir = Path(output_dir) if output_dir else None

        # Parallel execution
        self.n_workers = n_workers
        self.scoring_workers = scoring_workers
        self.scoring_batch_size = scoring_batch_size
        self.gating_workers = gating_workers

        # Track state during execution
        self._subclustered_clusters: List[str] = []
        self._result: Optional[EngineResult] = None
        self._marker_sets: Optional[Sequence] = None  # Cached marker sets

        # Iteration tracking for additive levels design
        self._iteration: int = 0
        self._cluster_col: str = "cluster_lvl1"
        self._label_col: str = "cell_type_lvl1"

        # Lazy-loaded annotation functions
        self._scoring_loaded = False
        self._compute_marker_scores = None
        self._compute_marker_scores_parallel = None
        self._load_marker_sets = None
        self._assign_labels_hierarchical = None
        self._DEFAULT_GATING_PARAMS = None
        self._build_de_rank_lookup = None
        self._compute_marker_idf = None

    def _ensure_scoring_functions(self):
        """Lazy-load scoring functions from annotation module."""
        if self._scoring_loaded:
            return

        from ..annotation import (
            compute_marker_scores,
            compute_marker_scores_parallel,
            load_marker_sets,
            assign_labels_hierarchical,
            DEFAULT_GATING_PARAMS,
            build_de_rank_lookup,
            compute_marker_idf,
            merge_gating_params,
            extract_gating_params_from_marker_map,
            log_gating_params,
        )

        self._compute_marker_scores = compute_marker_scores
        self._compute_marker_scores_parallel = compute_marker_scores_parallel
        self._load_marker_sets = load_marker_sets
        self._assign_labels_hierarchical = assign_labels_hierarchical
        self._build_de_rank_lookup = build_de_rank_lookup
        self._compute_marker_idf = compute_marker_idf
        self._merge_gating_params = merge_gating_params
        self._extract_gating_params_from_marker_map = extract_gating_params_from_marker_map
        self._log_gating_params = log_gating_params

        # Build gating params: base defaults + marker map tissue-specific rules
        self._DEFAULT_GATING_PARAMS = self._build_gating_params(DEFAULT_GATING_PARAMS)
        self._scoring_loaded = True

    def _build_gating_params(self, default_params: Dict[str, Any]) -> Dict[str, Any]:
        """Build gating params by merging defaults with marker map tissue-specific rules.

        Priority (lowest to highest):
        1. default_params (from annotation module, includes FT defaults for backward compat)
        2. marker_map._gating_params (tissue-specific from marker map)

        Parameters
        ----------
        default_params : dict
            Default gating parameters from annotation module

        Returns
        -------
        dict
            Merged gating parameters
        """
        # If no marker map, use defaults
        if self.marker_map is None:
            return default_params

        # Load marker map if it's a path
        marker_map_dict = None
        if isinstance(self.marker_map, (str, Path)):
            try:
                import json
                with open(self.marker_map, "r") as f:
                    marker_map_dict = json.load(f)
            except Exception as e:
                self.logger.warning("Failed to load marker map for gating params: %s", e)
                return default_params
        elif isinstance(self.marker_map, dict):
            marker_map_dict = self.marker_map

        if marker_map_dict is None:
            return default_params

        # Extract tissue-specific gating params from marker map
        tissue_params = self._extract_gating_params_from_marker_map(marker_map_dict)

        if tissue_params:
            self.logger.info(
                "[Stage I] Found _gating_params in marker map: %s",
                list(tissue_params.keys())
            )
            # Log details of tissue-specific rules
            if "root_hard_requirements" in tissue_params:
                for root, req in tissue_params["root_hard_requirements"].items():
                    self.logger.info(
                        "  root_hard_requirements[%s]: marker=%s, min_pos_frac=%.2f",
                        root, req.get("marker", "?"), req.get("min_pos_frac", 0)
                    )
            if "root_veto_markers" in tissue_params:
                for root, veto in tissue_params["root_veto_markers"].items():
                    self.logger.info(
                        "  root_veto_markers[%s]: markers=%s, max_pos_frac=%.2f",
                        root, veto.get("markers", []), veto.get("max_pos_frac", 0)
                    )
            # Merge: defaults < tissue_params
            self.logger.info("  Merging with DEFAULT_GATING_PARAMS")
            merged_params = self._merge_gating_params(
                base=default_params,
                tissue_params=tissue_params
            )
            # Log full merged gating parameters
            self._log_gating_params(merged_params, self.logger, "[Stage I]")
            return merged_params

        self.logger.info("[Stage I] No _gating_params in marker map, using DEFAULT_GATING_PARAMS")
        # Log full default gating parameters
        self._log_gating_params(default_params, self.logger, "[Stage I]")
        return default_params

    def execute(self, adata: "sc.AnnData", plan: RefinePlan) -> "sc.AnnData":
        """Execute RefinePlan operations on AnnData.

        Operations are executed in deterministic order:
        overrides → merges → subclusters → relabels → rescore

        Parameters
        ----------
        adata : sc.AnnData
            AnnData object to modify (modified in place)
        plan : RefinePlan
            Plan with operations to execute

        Returns
        -------
        sc.AnnData
            Modified AnnData object
        """
        start_time = time.time()

        self.logger.info("=" * 75)
        self.logger.info("RefinementEngine executing plan")
        self.logger.info("=" * 75)
        self.logger.info("Plan summary: %s", plan.summary())
        self.logger.info("Output label column: %s", self.label_key_out)

        # Reset tracking state
        self._subclustered_clusters = []
        operations_executed = {"override": 0, "merge": 0, "subcluster": 0, "relabel": 0, "rescore": 0}
        errors = []
        warnings_list = []
        total_cells_modified = 0

        # Ensure plan is sorted
        plan.sort_operations()

        # Initialize output columns
        self._init_output_columns(adata)

        # Separate operations by type
        override_ops = [op for op in plan.operations if isinstance(op, OverrideOp)]
        merge_ops = [op for op in plan.operations if isinstance(op, MergeOp)]
        subcluster_ops = [op for op in plan.operations if isinstance(op, SubclusterOp)]
        relabel_ops = [op for op in plan.operations if isinstance(op, RelabelOp)]
        rescore_ops = [op for op in plan.operations if isinstance(op, RescoreOp)]

        # Phase 1: Execute overrides
        for op in override_ops:
            try:
                n_cells = self._execute_override(adata, op)
                total_cells_modified += n_cells
                operations_executed["override"] += 1
            except Exception as e:
                error_msg = f"Error executing override for cluster {op.cluster_id}: {e}"
                self.logger.error(error_msg)
                errors.append(error_msg)

        # Phase 2: Execute merges
        for op in merge_ops:
            try:
                n_cells = self._execute_merge(adata, op)
                total_cells_modified += n_cells
                operations_executed["merge"] += 1
            except Exception as e:
                error_msg = f"Error executing merge for cluster {op.cluster_id}: {e}"
                self.logger.error(error_msg)
                errors.append(error_msg)

        # Phase 3: Execute subclustering
        if subcluster_ops:
            if self.skip_clustering:
                for op in subcluster_ops:
                    self.logger.info(
                        "SKIP: Subclustering for cluster %s (skip_clustering enabled)",
                        op.cluster_id
                    )
                operations_executed["subcluster_skipped"] = len(subcluster_ops)
            elif self.n_workers > 1:
                # Parallel path
                self.logger.info(
                    "Running %d subclustering operations in parallel (%d workers)",
                    len(subcluster_ops), self.n_workers
                )
                try:
                    from .parallel import run_subclusters_parallel

                    source_cluster_col = "cluster_lvl0" if self._iteration == 0 else f"cluster_lvl{self._iteration - 1}"

                    subclustered_ids, parallel_errors = run_subclusters_parallel(
                        adata=adata,
                        subcluster_ops=subcluster_ops,
                        n_workers=self.n_workers,
                        cluster_col=source_cluster_col,
                        output_cluster_col=self._cluster_col,
                        logger=self.logger,
                    )

                    for cid in subclustered_ids:
                        mask = adata.obs[self._cluster_col].astype(str).str.startswith(f"{cid}:")
                        total_cells_modified += mask.sum()

                    operations_executed["subcluster"] = len(subclustered_ids)
                    self._subclustered_clusters.extend(subclustered_ids)

                    for err in parallel_errors:
                        warnings_list.append(err)
                except Exception as e:
                    error_msg = f"Parallel subclustering failed: {e}"
                    self.logger.error(error_msg)
                    errors.append(error_msg)
            else:
                # Sequential path
                for op in subcluster_ops:
                    try:
                        n_cells = self._execute_subcluster(adata, op)
                        if n_cells > 0:
                            total_cells_modified += n_cells
                            operations_executed["subcluster"] += 1
                            self._subclustered_clusters.append(op.cluster_id)
                        else:
                            warnings_list.append(f"Subcluster failed for cluster {op.cluster_id}")
                    except Exception as e:
                        error_msg = f"Error executing subcluster for cluster {op.cluster_id}: {e}"
                        self.logger.error(error_msg)
                        errors.append(error_msg)

        # Phase 4: Execute relabels
        for op in relabel_ops:
            try:
                n_cells = self._execute_relabel(adata, op)
                total_cells_modified += n_cells
                operations_executed["relabel"] += 1
            except Exception as e:
                error_msg = f"Error executing relabel for cluster {op.cluster_id}: {e}"
                self.logger.error(error_msg)
                errors.append(error_msg)

        # Phase 5: Execute rescore
        for op in rescore_ops:
            try:
                self._execute_rescore(adata, op)
                operations_executed["rescore"] += 1
            except Exception as e:
                error_msg = f"Error executing rescore: {e}"
                self.logger.error(error_msg)
                errors.append(error_msg)

        # Calculate execution time
        execution_time = time.time() - start_time

        # Store result
        self._result = EngineResult(
            success=len(errors) == 0,
            n_cells_modified=total_cells_modified,
            operations_executed=operations_executed,
            subclusters_created=self._subclustered_clusters,
            errors=errors,
            warnings=warnings_list,
            execution_time_seconds=execution_time,
        )

        self.logger.info("-" * 75)
        self.logger.info("Execution complete")
        self.logger.info("  Cells modified: %d", total_cells_modified)
        self.logger.info("  Operations: %s", operations_executed)
        self.logger.info("  Subclusters created: %s", self._subclustered_clusters)
        self.logger.info("  Execution time: %.2f seconds", execution_time)
        if warnings_list:
            self.logger.warning("  Warnings: %s", warnings_list)
        if errors:
            self.logger.error("  Errors: %s", errors)

        return adata

    def get_result(self) -> Optional[EngineResult]:
        """Get result of last execution."""
        return self._result

    def _get_adapter(self) -> "InputAdapter":
        """Get or create InputAdapter for column name normalization."""
        if self._adapter is None:
            from .schema import InputAdapter
            self._adapter = InputAdapter(logger=self.logger)
        return self._adapter

    def _init_output_columns(self, adata: "sc.AnnData") -> None:
        """Initialize output columns if not present.

        Creates cell type columns based on iteration level:
        - cell_type: Original labels (preserved, never modified)
        - cell_type_lvl{N}: Refined labels for iteration N

        Detects current iteration from existing cluster_lvl columns.
        """
        from .operations import get_current_iteration

        # Get adapter and validate adata
        adapter = self._get_adapter()
        adata_result = adapter.validate_adata(adata)

        # Log warnings but don't fail on warnings-only
        adata_result.log_warnings(self.logger)

        # Raise if validation failed (missing required columns)
        adata_result.raise_if_invalid("input AnnData")

        # Detect current iteration from existing columns
        self._iteration = get_current_iteration(adata)

        # Set dynamic column names based on iteration
        if self._iteration == 0:
            self._cluster_col = "cluster_lvl1"
            self._label_col = "cell_type_lvl1"
        else:
            self._cluster_col = f"cluster_lvl{self._iteration}"
            self._label_col = f"cell_type_lvl{self._iteration}"

        self.logger.info(
            "Detected iteration %d: using %s and %s",
            self._iteration, self._cluster_col, self._label_col
        )

        # Get actual column names from adapter
        cluster_col = adata_result.adapted_columns.get(
            adapter.schema.cluster_key, "cluster_lvl0"
        )
        label_col = adata_result.adapted_columns.get(adapter.schema.label_key_in)

        # Initialize cell_type column from input (preserved original labels)
        if "cell_type" not in adata.obs:
            if label_col and label_col in adata.obs:
                adata.obs["cell_type"] = adata.obs[label_col].astype(str).copy()
                self.logger.info("Initialized cell_type from %s", label_col)
            else:
                adata.obs["cell_type"] = "Unknown"
                self.logger.warning(
                    "Label column not found; initialized cell_type to 'Unknown'"
                )
        else:
            if adata.obs["cell_type"].dtype.name == "category":
                adata.obs["cell_type"] = adata.obs["cell_type"].astype(str)

        # Initialize the target label column for this iteration
        if self._label_col not in adata.obs:
            if self._iteration == 0:
                source_col = "cell_type"
            else:
                prev_label_col = f"cell_type_lvl{self._iteration - 1}" if self._iteration > 1 else "cell_type_lvl1"
                source_col = prev_label_col if prev_label_col in adata.obs else "cell_type"

            adata.obs[self._label_col] = adata.obs[source_col].astype(str).copy()
            self.logger.info("Initialized %s from %s", self._label_col, source_col)
        else:
            if adata.obs[self._label_col].dtype.name == "category":
                adata.obs[self._label_col] = adata.obs[self._label_col].astype(str)

        # Also initialize legacy label_key_out if different
        if self.label_key_out not in ["cell_type", self._label_col]:
            if self.label_key_out not in adata.obs:
                adata.obs[self.label_key_out] = adata.obs[self._label_col].copy()
            elif adata.obs[self.label_key_out].dtype.name == "category":
                adata.obs[self.label_key_out] = adata.obs[self.label_key_out].astype(str)

        # Initialize reason column
        if self.reason_key_out not in adata.obs:
            adata.obs[self.reason_key_out] = ""
        else:
            if adata.obs[self.reason_key_out].dtype.name == "category":
                adata.obs[self.reason_key_out] = adata.obs[self.reason_key_out].astype(str)

        # Initialize cluster column for this iteration if not present
        if self._cluster_col not in adata.obs:
            if self._iteration == 0:
                if cluster_col in adata.obs:
                    adata.obs[self._cluster_col] = adata.obs[cluster_col].astype(str)
                    self.logger.info(
                        "Initialized %s from '%s'",
                        self._cluster_col, cluster_col,
                    )
                else:
                    raise ValueError(
                        f"Cluster column '{cluster_col}' not found in adata.obs. "
                        f"Available columns: {list(adata.obs.columns)[:20]}..."
                    )
            else:
                raise ValueError(
                    f"Cluster column '{self._cluster_col}' not found. "
                    f"Did you call initialize_iteration_level() before execute()?"
                )
        else:
            if adata.obs[self._cluster_col].dtype.name == "category":
                adata.obs[self._cluster_col] = adata.obs[self._cluster_col].astype(str)
                self.logger.info("Converted %s from Categorical to string dtype", self._cluster_col)

    def _execute_override(self, adata: "sc.AnnData", op: OverrideOp) -> int:
        """Execute override operation.

        Returns number of cells affected.
        """
        cluster_id = str(op.cluster_id)
        new_label = op.new_label
        reason = op.reason or f"Override to {new_label}"

        mask = adata.obs[self._cluster_col].astype(str) == cluster_id

        n_cells = mask.sum()
        if n_cells == 0:
            self.logger.warning("Override: Cluster %s has 0 cells; skipping", cluster_id)
            return 0

        adata.obs.loc[mask, self._label_col] = new_label
        adata.obs.loc[mask, self.reason_key_out] = reason

        if self.label_key_out != self._label_col:
            adata.obs.loc[mask, self.label_key_out] = new_label

        self.logger.info(
            "Override applied: Cluster %s → '%s' (%d cells) [%s]",
            cluster_id, new_label, n_cells, reason[:50],
        )
        return n_cells

    def _execute_merge(self, adata: "sc.AnnData", op: MergeOp) -> int:
        """Execute merge operation.

        Returns number of cells affected.
        """
        source_clusters = [str(c) for c in op.source_clusters]
        target_label = op.target_label
        reason = op.reason or f"Merged from {source_clusters}"

        mask = pd.Series(False, index=adata.obs.index)
        for cluster_id in source_clusters:
            mask = mask | (adata.obs[self._cluster_col].astype(str) == cluster_id)

        n_cells = mask.sum()
        if n_cells == 0:
            self.logger.warning("Merge: Source clusters %s have 0 cells; skipping", source_clusters)
            return 0

        adata.obs.loc[mask, self._label_col] = target_label
        adata.obs.loc[mask, self.reason_key_out] = reason

        if self.label_key_out != self._label_col:
            adata.obs.loc[mask, self.label_key_out] = target_label

        self.logger.info(
            "Merge applied: %s → '%s' (%d cells) [%s]",
            source_clusters, target_label, n_cells, reason[:50],
        )
        return n_cells

    def _execute_subcluster(self, adata: "sc.AnnData", op: SubclusterOp) -> int:
        """Execute subcluster operation.

        Returns number of cells modified (0 if failed).
        """
        from .operations import run_subcluster

        cluster_id = str(op.cluster_id)
        self.logger.info(
            "Subclustering cluster %s (resolution=%.2f) → %s",
            cluster_id, op.resolution, self._cluster_col
        )

        try:
            n_cells = run_subcluster(
                adata=adata,
                parent_cluster=cluster_id,
                cluster_col=self._cluster_col,
                resolution=op.resolution,
                n_pcs=op.n_pcs,
                neighbors_k=op.neighbors_k,
                focus_markers=op.focus_markers,
                min_cells=op.min_cells,
                logger=self.logger,
            )
            return n_cells
        except Exception as e:
            self.logger.error("Subcluster failed for cluster %s: %s", cluster_id, e)
            return 0

    def _execute_relabel(self, adata: "sc.AnnData", op: RelabelOp) -> int:
        """Execute relabel operation (instant, no clustering).

        NOTE: RELABEL is not supported in iterative mode.
        """
        if self._iteration > 0:
            raise NotImplementedError(
                f"RELABEL operations are not supported in iterative mode (iteration={self._iteration}). "
                "Use SUBCLUSTER to create new cluster levels, or modify the marker map."
            )

        cluster_id = str(op.cluster_id)
        new_label = op.new_label
        old_label = op.old_label
        reason = op.reason or f"Relabeled from {old_label} to {new_label}"

        mask = adata.obs[self._cluster_col].astype(str) == cluster_id

        n_cells = mask.sum()
        if n_cells == 0:
            self.logger.warning("Relabel: Cluster %s has 0 cells; skipping", cluster_id)
            return 0

        adata.obs.loc[mask, self._label_col] = new_label
        adata.obs.loc[mask, self.reason_key_out] = reason

        if self.label_key_out != self._label_col:
            adata.obs.loc[mask, self.label_key_out] = new_label

        self.logger.info(
            "Relabel applied: Cluster %s: '%s' → '%s' (%d cells, score=%.3f)",
            cluster_id, old_label, new_label, n_cells, op.confidence_score,
        )
        return n_cells

    def _execute_rescore(self, adata: "sc.AnnData", op: RescoreOp) -> pd.DataFrame:
        """Execute rescore operation.

        Recomputes marker scores for modified clusters.
        """
        self.logger.info("Rescore requested (mode=%s)", op.mode)

        if op.mode == "none":
            self.logger.info("Rescore mode is 'none'; skipping")
            return pd.DataFrame()

        if self.marker_map is None:
            self.logger.warning("Rescore requested but no marker_map provided; skipping")
            return pd.DataFrame()

        # Ensure scoring functions are loaded
        self._ensure_scoring_functions()

        # Load marker sets if not cached
        if self._marker_sets is None:
            self._marker_sets = self._load_marker_sets_from_map(adata.var_names)
            if self._marker_sets is None:
                self.logger.error("Failed to load marker sets; cannot rescore")
                return pd.DataFrame()

        # Determine which clusters need rescoring
        if op.target_clusters:
            changed_clusters = [str(c) for c in op.target_clusters]
        else:
            changed_clusters = [str(c) for c in self._subclustered_clusters]

        self.logger.info(
            "Rescore: mode=%s, changed_clusters=%s, recompute_de=%s",
            op.mode,
            changed_clusters or "all",
            op.recompute_de,
        )

        # Get DE lookup
        de_lookup, de_rank_lookup, marker_doc_freq = self._get_de_lookup(
            adata, op.recompute_de, subclustered_parents=self._subclustered_clusters
        )

        # Execute rescoring based on mode
        if op.mode == "smart" and self.stage_h_scores_path and self.stage_h_scores_path.exists():
            scores = self._smart_rescore(
                adata, changed_clusters, de_lookup, de_rank_lookup, marker_doc_freq
            )
        else:
            scores = self._full_rescore(adata, de_lookup, de_rank_lookup, marker_doc_freq)

        # Store in AnnData
        if not scores.empty:
            scores = _sanitize_df_for_h5ad(scores)
            adata.uns["marker_scores_refined"] = scores
            self.logger.info("Stored %d score records in adata.uns['marker_scores_refined']", len(scores))

            # Assign labels to subclustered clusters
            self._assign_labels_to_subclusters(adata, scores, changed_clusters)

        return scores

    def _load_marker_sets_from_map(self, var_names: Sequence[str]) -> Optional[Sequence]:
        """Load marker sets from marker_map (dict or path)."""
        if self.marker_map is None:
            return None

        try:
            if isinstance(self.marker_map, (str, Path)):
                marker_path = Path(self.marker_map)
                if not marker_path.exists():
                    self.logger.error("Marker map file not found: %s", marker_path)
                    return None
                import json
                with open(marker_path, "r") as f:
                    marker_map_dict = json.load(f)
                return self._load_marker_sets(marker_map_dict, var_names, self.logger)
            elif isinstance(self.marker_map, dict):
                return self._load_marker_sets(self.marker_map, var_names, self.logger)
            else:
                self.logger.error("Invalid marker_map type: %s", type(self.marker_map))
                return None
        except Exception as e:
            self.logger.error("Failed to load marker sets: %s", e)
            return None

    def _get_de_lookup(
        self, adata: "sc.AnnData", recompute: bool, subclustered_parents: Optional[List[str]] = None
    ) -> Tuple[Dict[str, List[str]], Optional[Dict[str, Dict]], Optional[Dict[str, int]]]:
        """Get differential expression lookup table."""
        self.logger.info("Getting DE lookup (recompute=%s)", recompute)
        de_lookup_start = time.time()

        cluster_key = self._cluster_col
        de_key = f"de_wilcoxon_{self.de_layer}"

        de_lookup = None
        if not recompute:
            if de_key in adata.uns:
                self.logger.info("Reusing cached DE results from adata.uns['%s']", de_key)
                de_results = adata.uns[de_key]
                if isinstance(de_results, dict):
                    de_lookup = de_results
            elif "de_results" in adata.uns:
                de_results = adata.uns["de_results"]
                if isinstance(de_results, dict):
                    de_lookup = de_results

        # Recompute DE tests if needed
        if de_lookup is None:
            has_subclusters = any(
                ":" in str(c) for c in adata.obs[cluster_key].astype(str).unique()
            )
            use_within_parent = self.de_within_parent and has_subclusters

            self.logger.info(
                "Computing differential expression tests (layer=%s, within_parent=%s)...",
                self.de_layer,
                use_within_parent,
            )
            try:
                if use_within_parent and self._de_within_parent_runner is not None:
                    de_lookup = self._de_within_parent_runner(
                        adata,
                        subcluster_key=cluster_key,
                        method="wilcoxon",
                        n_genes=20,
                        n_workers=self.de_workers,
                        layer=self.de_layer,
                        tie_correct=self.de_tie_correct,
                        logger=self.logger,
                        subclustered_parents=subclustered_parents,
                    )
                else:
                    de_lookup = self._de_runner(
                        adata,
                        cluster_key=cluster_key,
                        method="wilcoxon",
                        n_genes=20,
                        layer=self.de_layer,
                    )
                adata.uns["de_results"] = de_lookup
            except Exception as e:
                self.logger.warning("DE computation failed: %s; proceeding without DE bonus", e)
                de_lookup = {str(c): [] for c in adata.obs[cluster_key].unique()}

        # Build rank-weighted DE lookup
        de_rank_lookup = None
        if de_key in adata.uns or "rank_genes_groups" in adata.uns:
            try:
                de_rank_lookup = self._build_de_rank_lookup(
                    adata,
                    cluster_key=cluster_key,
                    de_key=de_key,
                    de_top_frac=self.de_top_frac,
                    de_min_k=self.de_min_k,
                    de_max_k=self.de_max_k,
                    logger=self.logger,
                )
            except Exception as e:
                self.logger.warning("Failed to build DE rank lookup: %s", e)

        # Compute marker document frequency
        marker_doc_freq = None
        if self._marker_sets is not None:
            try:
                marker_doc_freq = self._compute_marker_idf(self._marker_sets, logger=self.logger)
            except Exception as e:
                self.logger.warning("Failed to compute marker IDF: %s", e)

        de_lookup_elapsed = time.time() - de_lookup_start
        self.logger.info("DE lookup complete in %.1f sec", de_lookup_elapsed)
        return de_lookup, de_rank_lookup, marker_doc_freq

    def _smart_rescore(
        self,
        adata: "sc.AnnData",
        changed_clusters: List[str],
        de_lookup: Dict[str, List[str]],
        de_rank_lookup: Optional[Dict[str, Dict]] = None,
        marker_doc_freq: Optional[Dict[str, int]] = None,
    ) -> pd.DataFrame:
        """Smart rescoring: reuse unchanged cluster scores, recompute only modified."""
        self.logger.info("Smart rescoring: reusing unchanged clusters, recomputing subclusters")

        cluster_key = self._cluster_col
        h_scores = pd.read_csv(self.stage_h_scores_path)
        self.logger.info("Loaded %d Stage H score records", len(h_scores))

        all_cluster_ids = set(adata.obs[cluster_key].astype(str).unique())
        changed_set = set(str(c) for c in changed_clusters)

        unchanged = []
        new_cluster_ids = []

        for cid in all_cluster_ids:
            if ":" in str(cid):
                new_cluster_ids.append(str(cid))
            elif str(cid) not in changed_set:
                unchanged.append(str(cid))

        self.logger.info("Smart rescoring breakdown:")
        self.logger.info("  Total clusters: %d", len(all_cluster_ids))
        self.logger.info("  Unchanged (reuse): %d", len(unchanged))
        self.logger.info("  New subclusters (compute): %d", len(new_cluster_ids))

        reused_scores = h_scores[h_scores["cluster_id"].astype(str).isin(unchanged)].copy()
        if not reused_scores.empty:
            reused_scores["cluster_id"] = reused_scores["cluster_id"].astype(str)
        self.logger.info("Reused %d score records from Stage H", len(reused_scores))

        new_marker_evidence = None
        if new_cluster_ids:
            self.logger.info("Computing scores for %d new subclusters", len(new_cluster_ids))
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=pd.errors.SettingWithCopyWarning)
                new_subset = adata[adata.obs[cluster_key].astype(str).isin(new_cluster_ids)].copy()

            new_subset.obs["cluster_lvl0"] = new_subset.obs[cluster_key]

            # Determine scoring layer (use batchcorr/stage_h_input if available)
            use_layer = self.scoring_layer if self.scoring_layer in new_subset.layers else None
            if use_layer:
                self.logger.info("Using layer '%s' for subcluster scoring", use_layer)

            if self.expand_markers:
                new_scores, new_marker_evidence = self._compute_marker_scores_parallel(
                    new_subset,
                    marker_sets=self._marker_sets,
                    cluster_key="cluster_lvl0",
                    positive_quantile=self.positive_quantile,
                    de_lookup=de_lookup,
                    de_bonus=self.de_bonus,
                    anti_weight=self.anti_weight,
                    logger=self.logger,
                    expand_markers=True,
                    de_rank_lookup=de_rank_lookup,
                    marker_doc_freq=marker_doc_freq,
                    de_commonness_alpha=self.de_commonness_alpha,
                    n_workers=self.scoring_workers,
                    batch_size=self.scoring_batch_size,
                    layer=use_layer,
                )
            else:
                new_scores = self._compute_marker_scores_parallel(
                    new_subset,
                    marker_sets=self._marker_sets,
                    cluster_key="cluster_lvl0",
                    positive_quantile=self.positive_quantile,
                    de_lookup=de_lookup,
                    de_bonus=self.de_bonus,
                    anti_weight=self.anti_weight,
                    logger=self.logger,
                    expand_markers=False,
                    de_rank_lookup=de_rank_lookup,
                    marker_doc_freq=marker_doc_freq,
                    de_commonness_alpha=self.de_commonness_alpha,
                    n_workers=self.scoring_workers,
                    batch_size=self.scoring_batch_size,
                    layer=use_layer,
                )
            self.logger.info("Computed %d new score records", len(new_scores))
        else:
            new_scores = pd.DataFrame()

        combined = pd.concat([reused_scores, new_scores], ignore_index=True)
        if "cluster_id" in combined.columns:
            combined["cluster_id"] = combined["cluster_id"].astype(str)
        self.logger.info("Combined scores: %d total records", len(combined))

        if self.expand_markers and new_marker_evidence is not None and not new_marker_evidence.empty:
            new_marker_evidence = _sanitize_df_for_h5ad(new_marker_evidence)
            adata.uns["marker_evidence_refined"] = new_marker_evidence

        return combined

    def _full_rescore(
        self,
        adata: "sc.AnnData",
        de_lookup: Dict[str, List[str]],
        de_rank_lookup: Optional[Dict[str, Dict]] = None,
        marker_doc_freq: Optional[Dict[str, int]] = None,
    ) -> pd.DataFrame:
        """Full rescoring: recompute all scores from scratch."""
        self.logger.info("Full rescoring: computing all scores from scratch")

        cluster_key = self._cluster_col

        # Determine scoring layer (use batchcorr/stage_h_input if available)
        use_layer = self.scoring_layer if self.scoring_layer in adata.layers else None
        if use_layer:
            self.logger.info("Using layer '%s' for marker scoring", use_layer)

        if self.expand_markers:
            scores, marker_evidence = self._compute_marker_scores_parallel(
                adata,
                marker_sets=self._marker_sets,
                cluster_key=cluster_key,
                positive_quantile=self.positive_quantile,
                de_lookup=de_lookup,
                de_bonus=self.de_bonus,
                anti_weight=self.anti_weight,
                logger=self.logger,
                expand_markers=True,
                de_rank_lookup=de_rank_lookup,
                marker_doc_freq=marker_doc_freq,
                de_commonness_alpha=self.de_commonness_alpha,
                n_workers=self.scoring_workers,
                batch_size=self.scoring_batch_size,
                layer=use_layer,
            )
            if marker_evidence is not None and not marker_evidence.empty:
                marker_evidence = _sanitize_df_for_h5ad(marker_evidence)
                adata.uns["marker_evidence_refined"] = marker_evidence
        else:
            scores = self._compute_marker_scores_parallel(
                adata,
                marker_sets=self._marker_sets,
                cluster_key=cluster_key,
                positive_quantile=self.positive_quantile,
                de_lookup=de_lookup,
                de_bonus=self.de_bonus,
                anti_weight=self.anti_weight,
                logger=self.logger,
                expand_markers=False,
                de_rank_lookup=de_rank_lookup,
                marker_doc_freq=marker_doc_freq,
                de_commonness_alpha=self.de_commonness_alpha,
                n_workers=self.scoring_workers,
                batch_size=self.scoring_batch_size,
                layer=use_layer,
            )

        self.logger.info("Computed %d total score records", len(scores))
        return scores

    def _assign_labels_to_subclusters(
        self,
        adata: "sc.AnnData",
        scores: pd.DataFrame,
        changed_clusters: List[str],
    ) -> None:
        """Assign cell type labels to subclustered clusters based on marker scores."""
        if scores.empty:
            return

        cluster_col = self._cluster_col

        cluster_assignments, decision_steps = self._assign_labels_hierarchical(
            marker_scores=scores,
            adata=adata,
            cluster_key=cluster_col,
            layer=None,
            params=self._DEFAULT_GATING_PARAMS,
            logger=self.logger,
            n_workers=self.gating_workers,
        )

        if cluster_assignments.empty:
            self.logger.warning("No cluster assignments from hierarchical scoring; skipping label assignment")
            return

        if not decision_steps.empty:
            decision_steps = _sanitize_df_for_h5ad(decision_steps)
            adata.uns["decision_steps_subcluster"] = decision_steps

        if not cluster_assignments.empty:
            cluster_assignments = _sanitize_df_for_h5ad(cluster_assignments)
            adata.uns["cluster_annotations_subcluster"] = cluster_assignments

        # Build mappings: cluster_id -> (assigned_label, reason_string)
        label_mapping = {}
        reason_mapping = {}

        for _, row in cluster_assignments.iterrows():
            cid = str(row["cluster_id"])
            label_mapping[cid] = row["assigned_label"]

            # Build meaningful reason from scoring metadata (reference format)
            score = row.get("assigned_score", 0)
            confidence = row.get("confidence", 0)
            root = row.get("root_label", "")
            stop_reason = row.get("stop_reason", "")
            assigned_label = row["assigned_label"]

            # Format reason based on stop_reason and hierarchy
            if stop_reason == "ambiguous_root":
                reason_mapping[cid] = f"Ambiguous roots (gap={confidence:.2f})"
            elif stop_reason == "ambiguous_siblings":
                reason_mapping[cid] = f"Mixed population (score={score:.2f})"
            elif stop_reason == "no_root_passed":
                # Include detailed fail reasons if available (e.g., veto info)
                root_fail_reasons = row.get("root_fail_reasons", "")
                if root_fail_reasons:
                    reason_mapping[cid] = f"No root gate passed: {root_fail_reasons}"
                else:
                    reason_mapping[cid] = "No root gate passed"
            elif stop_reason == "no_child_passed":
                # Stopped at parent level
                reason_mapping[cid] = f"score={score:.2f}, margin={confidence:.2f}"
            elif root and root != assigned_label and "~" not in root:
                # Hierarchical descent: root -> subtype
                reason_mapping[cid] = f"{root} -> {assigned_label} (score={score:.2f}, margin={confidence:.2f})"
            else:
                # Default: just score and margin
                reason_mapping[cid] = f"score={score:.2f}, margin={confidence:.2f}"

        # Convert categorical columns to string
        if hasattr(adata.obs[self._label_col], "cat"):
            adata.obs[self._label_col] = adata.obs[self._label_col].astype(str)
        if self.reason_key_out in adata.obs and hasattr(adata.obs[self.reason_key_out], "cat"):
            adata.obs[self.reason_key_out] = adata.obs[self.reason_key_out].astype(str)

        # Filter to only subcluster mappings (contain ":")
        subcluster_label_map = {
            cid: label for cid, label in label_mapping.items() if ":" in str(cid)
        }
        subcluster_reason_map = {
            cid: reason_mapping.get(cid, "score unknown")
            for cid in subcluster_label_map.keys()
        }

        if not subcluster_label_map:
            self.logger.info("No subclusters to update")
            return

        cluster_ids = adata.obs[cluster_col].astype(str)

        # Apply new labels using vectorized map
        new_labels = cluster_ids.map(subcluster_label_map)
        update_mask = new_labels.notna()
        adata.obs.loc[update_mask, self._label_col] = new_labels[update_mask]

        # Apply reasons
        new_reasons = cluster_ids.map(subcluster_reason_map)
        adata.obs.loc[update_mask, self.reason_key_out] = new_reasons[update_mask]

        n_cells_updated = update_mask.sum()
        n_subclusters_updated = len([c for c in subcluster_label_map.keys()])

        self.logger.info("Updated %s: %d cells (%d subclusters)",
                        self._label_col, n_cells_updated, n_subclusters_updated)
