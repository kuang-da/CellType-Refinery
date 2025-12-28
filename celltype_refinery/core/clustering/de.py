"""Differential expression testing for clustering module.

Provides Wilcoxon rank-sum and other DE methods for identifying
cluster-specific markers.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple
import logging
import time
import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading

import numpy as np
import pandas as pd

from .config import DEConfig, StageHConfig


@dataclass
class DEResult:
    """Result from differential expression testing.

    Attributes
    ----------
    cluster_de_genes : Dict[str, List[str]]
        Map of cluster ID to list of top DE gene names
    key_added : str
        Key in adata.uns containing full DE results
    elapsed_seconds : float
        Time taken for DE computation
    """

    cluster_de_genes: Dict[str, List[str]] = field(default_factory=dict)
    key_added: str = ""
    elapsed_seconds: float = 0.0


class DERunner:
    """Differential expression test runner.

    Provides methods for computing differential expression between
    clusters, both globally and within parent clusters.

    Parameters
    ----------
    config : StageHConfig, optional
        Stage H configuration. If None, uses defaults.
    logger : logging.Logger, optional
        Logger instance. If None, creates default logger.

    Example
    -------
    >>> from celltype_refinery.core.clustering import DERunner, StageHConfig
    >>> config = StageHConfig()
    >>> runner = DERunner(config)
    >>> result = runner.run_de_tests(adata, cluster_key="leiden")
    """

    def __init__(
        self,
        config: Optional[StageHConfig] = None,
        logger: Optional[logging.Logger] = None,
    ):
        self.config = config or StageHConfig()
        self.logger = logger or logging.getLogger(__name__)
        self._check_dependencies()

    def _check_dependencies(self) -> None:
        """Check for required dependencies."""
        try:
            import scanpy
        except ImportError:
            raise RuntimeError(
                "Differential expression requires scanpy. "
                "Install with: pip install scanpy"
            )

    def run_de_tests(
        self,
        adata: Any,  # AnnData
        cluster_key: str,
        method: Optional[str] = None,
        n_genes: Optional[int] = None,
        layer: Optional[str] = None,
        tie_correct: Optional[bool] = None,
        output_dir: Optional[Path] = None,
        key_added: Optional[str] = None,
    ) -> DEResult:
        """Run differential expression tests between clusters.

        Parameters
        ----------
        adata : AnnData
            AnnData object with expression data and cluster assignments
        cluster_key : str
            Column name in adata.obs with cluster labels
        method : str, optional
            DE method ('wilcoxon', 't-test', etc.). Uses config default if None.
        n_genes : int, optional
            Number of top DE genes to return per cluster. Uses config default if None.
        layer : str, optional
            Layer to use for DE analysis. Uses config default if None.
        tie_correct : bool, optional
            Apply tie correction for Wilcoxon test. Uses config default if None.
        output_dir : Path, optional
            Directory to save debug output
        key_added : str, optional
            Key to store results in adata.uns

        Returns
        -------
        DEResult
            DE result with per-cluster gene lists
        """
        import scanpy as sc

        # Use config defaults where not specified
        cfg = self.config.de
        method = method if method is not None else cfg.method
        n_genes = n_genes if n_genes is not None else cfg.n_genes
        layer = layer if layer is not None else cfg.layer
        tie_correct = tie_correct if tie_correct is not None else cfg.tie_correct

        # Determine the key for storing results
        if key_added is None:
            layer_suffix = layer if layer else "X"
            key_added = f"de_{method}_{layer_suffix}"

        # Validate layer exists
        if layer and layer not in adata.layers:
            available = list(adata.layers.keys())
            self.logger.warning(
                "Layer '%s' not found in adata.layers (available: %s). "
                "Falling back to adata.X",
                layer,
                available,
            )
            layer = None

        n_cells = adata.n_obs
        n_clusters = adata.obs[cluster_key].nunique()
        self.logger.info(
            "Computing differential expression (method=%s, layer=%s, "
            "tie_correct=%s, top_n=%d, key=%s)",
            method,
            layer if layer else "X",
            tie_correct,
            n_genes,
            key_added,
        )
        self.logger.info(
            "Starting %s DE: %d cells, %d clusters. This may take several minutes...",
            method,
            n_cells,
            n_clusters,
        )

        de_start_time = time.time()

        # Run rank_genes_groups
        sc.tl.rank_genes_groups(
            adata,
            groupby=cluster_key,
            method=method,
            n_genes=n_genes,
            layer=layer,
            use_raw=False,
            tie_correct=tie_correct,
            key_added=key_added,
            pts=True,
        )

        de_elapsed = time.time() - de_start_time
        self.logger.info(
            "%s DE completed in %.1f seconds (%.1f min)",
            method,
            de_elapsed,
            de_elapsed / 60,
        )

        # Extract results
        clusters = sorted(adata.obs[cluster_key].astype(str).unique())
        result = DEResult(key_added=key_added, elapsed_seconds=de_elapsed)

        for cluster in clusters:
            df = sc.get.rank_genes_groups_df(adata, group=cluster, key=key_added)
            df = df.dropna(subset=["names"])
            top = df.head(n_genes)["names"].astype(str).tolist()
            result.cluster_de_genes[cluster] = top

        # Save debug output if requested
        if output_dir is not None:
            debug_dir = Path(output_dir) / "tmp" / "de_debug"
            debug_dir.mkdir(parents=True, exist_ok=True)
            for cluster in clusters:
                df = sc.get.rank_genes_groups_df(adata, group=cluster, key=key_added)
                debug_path = debug_dir / f"de_results_cluster_{cluster}.csv"
                df.to_csv(debug_path, index=False)
                self.logger.debug(
                    "Saved DE debug for cluster %s: %d rows, columns=%s",
                    cluster,
                    len(df),
                    list(df.columns),
                )
            self.logger.info("Saved DE debug files to %s", debug_dir)

        self.logger.info(
            "Computed differential expression for %d clusters", len(clusters)
        )
        return result

    def run_de_within_parent(
        self,
        adata: Any,  # AnnData
        subcluster_key: str = "cluster_lvl1",
        parent_key: str = "cluster_lvl0",
        method: Optional[str] = None,
        n_genes: Optional[int] = None,
        layer: Optional[str] = None,
        tie_correct: Optional[bool] = None,
        n_workers: int = 8,
        subclustered_parents: Optional[List[str]] = None,
    ) -> DEResult:
        """Run within-parent DE for subclusters.

        For each parent cluster, compares subclusters against each other
        (not entire dataset). This is biologically appropriate and faster.

        Parameters
        ----------
        adata : AnnData
            AnnData with cluster annotations
        subcluster_key : str
            Column with subcluster IDs (e.g., "3:0", "3:1")
        parent_key : str
            Column with parent cluster IDs (e.g., "3")
        method : str, optional
            DE method. Uses config default if None.
        n_genes : int, optional
            Number of top DE genes per cluster. Uses config default if None.
        layer : str, optional
            Data layer to use. Uses config default if None.
        tie_correct : bool, optional
            Apply tie correction for Wilcoxon. Uses config default if None.
        n_workers : int
            Number of parallel workers
        subclustered_parents : List[str], optional
            List of parent cluster IDs that were subclustered. If provided,
            only these parents will be processed for DE.

        Returns
        -------
        DEResult
            DE result with per-subcluster gene lists
        """
        import scanpy as sc

        # Use config defaults
        cfg = self.config.de
        method = method if method is not None else cfg.method
        n_genes = n_genes if n_genes is not None else cfg.n_genes
        layer = layer if layer is not None else cfg.layer
        tie_correct = tie_correct if tie_correct is not None else cfg.tie_correct

        # Identify parent clusters with subclusters
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
                self.logger.info(
                    "Filtered to %d newly-subclustered parents (skipping %d unchanged)",
                    len(parents_with_subclusters),
                    n_skipped,
                )
        else:
            parents_with_subclusters = all_parents_with_subclusters

        n_parents = len(parents_with_subclusters)
        if n_parents == 0:
            self.logger.info(
                "No subclusters found (no ':' in cluster IDs), falling back to global DE"
            )
            return self.run_de_tests(
                adata, subcluster_key, method, n_genes, layer, tie_correct
            )

        self.logger.info(
            "Within-parent DE: %d parent clusters with subclusters, %d workers",
            n_parents,
            n_workers,
        )

        start_time = time.time()
        all_results: Dict[str, List[str]] = {}
        results_lock = threading.Lock()

        def process_parent(parent_id: str) -> Tuple[Dict[str, List[str]], str, int, float]:
            """Run DE for one parent cluster's subclusters."""
            try:
                parent_mask = (
                    adata.obs[subcluster_key]
                    .astype(str)
                    .str.startswith(f"{parent_id}:")
                )
                n_cells = parent_mask.sum()

                if n_cells < 100:
                    return {}, parent_id, n_cells, 0.0

                # Create subset for this parent
                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore", category=pd.errors.SettingWithCopyWarning)
                    parent_adata = adata[parent_mask].copy()

                n_subclusters = parent_adata.obs[subcluster_key].nunique()
                if n_subclusters < 2:
                    return {}, parent_id, n_cells, 0.0

                parent_start = time.time()

                # Run DE within this parent's subclusters
                sc.tl.rank_genes_groups(
                    parent_adata,
                    groupby=subcluster_key,
                    method=method,
                    n_genes=n_genes,
                    layer=layer,
                    use_raw=False,
                    tie_correct=tie_correct,
                    pts=True,
                )

                # Extract results
                parent_results = {}
                for sc_id in parent_adata.obs[subcluster_key].unique():
                    df = sc.get.rank_genes_groups_df(parent_adata, group=sc_id)
                    df = df.dropna(subset=["names"])
                    parent_results[str(sc_id)] = df.head(n_genes)["names"].astype(str).tolist()

                elapsed = time.time() - parent_start
                return parent_results, parent_id, n_cells, elapsed

            except Exception as e:
                self.logger.warning(
                    "DE failed for parent %s: %s", parent_id, str(e)
                )
                return {}, parent_id, 0, 0.0

        # Run in parallel
        with ThreadPoolExecutor(max_workers=n_workers) as executor:
            futures = {
                executor.submit(process_parent, p): p
                for p in sorted(parents_with_subclusters)
            }

            for future in as_completed(futures):
                parent_results, parent_id, n_cells, elapsed = future.result()
                if parent_results:
                    with results_lock:
                        all_results.update(parent_results)
                    self.logger.debug(
                        "Completed DE for parent %s: %d subclusters in %.1fs",
                        parent_id,
                        len(parent_results),
                        elapsed,
                    )

        total_elapsed = time.time() - start_time
        self.logger.info(
            "Within-parent DE completed: %d subclusters in %.1f seconds",
            len(all_results),
            total_elapsed,
        )

        return DEResult(
            cluster_de_genes=all_results,
            key_added=f"de_{method}_{layer or 'X'}_within_parent",
            elapsed_seconds=total_elapsed,
        )

    @staticmethod
    def extract_de_results(
        adata: Any,  # AnnData
        key: str,
        n_genes: int = 20,
    ) -> Dict[str, pd.DataFrame]:
        """Extract DE results from adata.uns as DataFrames.

        Parameters
        ----------
        adata : AnnData
            AnnData with DE results stored
        key : str
            Key in adata.uns containing DE results
        n_genes : int
            Maximum number of genes to extract per cluster

        Returns
        -------
        Dict[str, pd.DataFrame]
            Map of cluster ID to DataFrame with DE statistics
        """
        import scanpy as sc

        if key not in adata.uns:
            raise KeyError(f"DE results not found at key '{key}' in adata.uns")

        results = {}
        groups = adata.uns[key]["names"].dtype.names

        for group in groups:
            df = sc.get.rank_genes_groups_df(adata, group=group, key=key)
            results[str(group)] = df.head(n_genes)

        return results

    @staticmethod
    def build_de_rank_lookup(
        de_result: DEResult,
        adata: Any = None,
        key: Optional[str] = None,
    ) -> Dict[str, Dict[str, int]]:
        """Build lookup table for DE gene ranks.

        Parameters
        ----------
        de_result : DEResult
            DE result object
        adata : AnnData, optional
            AnnData with full DE results (for more detailed lookup)
        key : str, optional
            Key in adata.uns for full results

        Returns
        -------
        Dict[str, Dict[str, int]]
            Map of cluster -> gene -> rank (0-indexed)
        """
        rank_lookup: Dict[str, Dict[str, int]] = {}

        for cluster, genes in de_result.cluster_de_genes.items():
            rank_lookup[cluster] = {gene: i for i, gene in enumerate(genes)}

        return rank_lookup


# Convenience function for dependency injection pattern
def default_de_runner(
    adata: Any,
    cluster_key: str,
    method: str = "wilcoxon",
    n_genes: int = 12,
    layer: Optional[str] = "batchcorr",
    tie_correct: bool = True,
    logger: Optional[logging.Logger] = None,
) -> Dict[str, List[str]]:
    """Default DE runner function for dependency injection.

    This function provides a simple interface for running DE tests,
    suitable for use in dependency injection patterns.

    Parameters
    ----------
    adata : AnnData
        AnnData object
    cluster_key : str
        Cluster column name
    method : str
        DE method
    n_genes : int
        Number of top genes
    layer : str, optional
        Data layer
    tie_correct : bool
        Apply tie correction
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    Dict[str, List[str]]
        Map of cluster to top DE genes
    """
    runner = DERunner(logger=logger)
    result = runner.run_de_tests(
        adata,
        cluster_key=cluster_key,
        method=method,
        n_genes=n_genes,
        layer=layer,
        tie_correct=tie_correct,
    )
    return result.cluster_de_genes
