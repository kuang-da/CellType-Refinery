"""Annotation engine for cell-type classification.

This module provides the main AnnotationEngine class that orchestrates
the annotation pipeline: marker loading, scoring, hierarchical gating,
and cell assignment.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import pandas as pd

try:
    import scanpy as sc
except ImportError:
    sc = None

from .assignment import annotate_obs
from .gating import (
    DEFAULT_GATING_PARAMS,
    assign_labels_hierarchical,
    merge_gating_params,
    extract_gating_params_from_marker_map,
    log_gating_params,
)
from .marker_loading import (
    MarkerSet,
    compute_marker_document_frequency,
    load_marker_sets,
)
from .scoring import (
    build_de_rank_lookup,
    compute_marker_idf,
    compute_marker_scores,
    extract_existing_de_results,
)


@dataclass
class AnnotationResult:
    """Result from annotation pipeline.

    Attributes:
        adata: AnnData with cell-type columns added to obs
        marker_scores: Full scoring matrix (clusters × cell types)
        marker_evidence: Per-marker evidence table (if requested)
        cluster_annotations: Cluster → cell type mapping with metadata
        decision_steps: Detailed decision trace for each cluster
    """

    adata: "sc.AnnData"
    marker_scores: pd.DataFrame
    marker_evidence: Optional[pd.DataFrame]
    cluster_annotations: pd.DataFrame
    decision_steps: pd.DataFrame


@dataclass
class AnnotationParams:
    """Parameters for annotation pipeline.

    Attributes:
        cluster_key: Column in adata.obs with cluster assignments
        label_col: Output column name for cell type labels
        layer: Layer to use for expression data
        positive_quantile: Quantile threshold for positive cells
        de_key: Key in adata.uns for DE results
        de_bonus: Maximum DE bonus
        de_top_frac: Fraction of panel for DE top-K
        de_min_k: Minimum K for DE
        de_max_k: Maximum K for DE
        de_commonness_alpha: Exponent for commonness penalty
        anti_weight: Weight for anti-marker penalty
        anti_agg: Aggregation mode for anti-markers
        use_idf: Whether to use IDF weighting
        expand_markers: Whether to generate per-marker evidence table
        gating_params: Override gating parameters
    """

    cluster_key: str = "cluster_lvl0"
    label_col: str = "cell_type_auto"
    layer: str = "batchcorr"
    positive_quantile: float = 0.75
    de_key: str = "de_wilcoxon"
    de_bonus: float = 0.6
    de_top_frac: float = 0.2
    de_min_k: int = 3
    de_max_k: int = 12
    de_commonness_alpha: float = 0.5
    anti_weight: float = 0.8
    anti_agg: str = "top2mean"
    use_idf: bool = True
    expand_markers: bool = True
    gating_params: Optional[Dict[str, Any]] = None


class AnnotationEngine:
    """Standalone annotation engine for pre-clustered AnnData.

    This engine performs cell-type annotation on AnnData objects that
    already have cluster assignments and DE results. It:
    1. Loads and resolves marker sets against the panel
    2. Scores clusters against marker sets
    3. Assigns labels hierarchically with gating
    4. Maps labels to individual cells

    Example:
        >>> engine = AnnotationEngine(
        ...     marker_map=Path("markers.json"),
        ... )
        >>> result = engine.run(
        ...     adata=adata,
        ...     cluster_key="cluster_lvl0",
        ...     output_dir=Path("output/"),
        ... )
        >>> # result.adata now has cell_type columns
        >>> # result.cluster_annotations has the mapping
    """

    def __init__(
        self,
        marker_map: Union[Dict, Path, str],
        params: Optional[AnnotationParams] = None,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize annotation engine.

        Args:
            marker_map: Marker map dict, Path, or string path to JSON
            params: Annotation parameters (uses defaults if None)
            logger: Logger instance
        """
        self.marker_map_source = marker_map
        self.params = params or AnnotationParams()
        self.logger = logger or logging.getLogger(__name__)

        # Load marker map if path
        if isinstance(marker_map, (str, Path)):
            path = Path(marker_map)
            if not path.exists():
                raise FileNotFoundError(f"Marker map not found: {path}")
            with open(path, "r") as f:
                self._marker_map = json.load(f)
            self._marker_map_path = path
        else:
            self._marker_map = marker_map
            self._marker_map_path = None

        # Extract gating params from marker map if present
        self._marker_map_gating_params = extract_gating_params_from_marker_map(self._marker_map)
        if self._marker_map_gating_params:
            self.logger.info(
                "[Stage H] Found _gating_params in marker map: %s",
                list(self._marker_map_gating_params.keys())
            )
            # Log details of tissue-specific rules
            if "root_hard_requirements" in self._marker_map_gating_params:
                for root, req in self._marker_map_gating_params["root_hard_requirements"].items():
                    self.logger.info(
                        "  root_hard_requirements[%s]: marker=%s, min_pos_frac=%.2f",
                        root, req.get("marker", "?"), req.get("min_pos_frac", 0)
                    )
            if "root_veto_markers" in self._marker_map_gating_params:
                for root, veto in self._marker_map_gating_params["root_veto_markers"].items():
                    self.logger.info(
                        "  root_veto_markers[%s]: markers=%s, max_pos_frac=%.2f",
                        root, veto.get("markers", []), veto.get("max_pos_frac", 0)
                    )

    def run(
        self,
        adata: "sc.AnnData",
        cluster_key: Optional[str] = None,
        output_dir: Optional[Path] = None,
    ) -> AnnotationResult:
        """Run full annotation pipeline.

        Args:
            adata: AnnData with cluster_key in obs + DE results in uns
            cluster_key: Column name for cluster IDs (overrides params)
            output_dir: Where to write CSVs (None = don't write)

        Returns:
            AnnotationResult with marker_scores, cluster_annotations,
            and modified adata (with cell_type columns)
        """
        cluster_key = cluster_key or self.params.cluster_key

        self.logger.info("=" * 70)
        self.logger.info("ANNOTATION ENGINE")
        self.logger.info("=" * 70)
        if self._marker_map_path:
            self.logger.info("Marker map: %s", self._marker_map_path)
        self.logger.info("Cluster key: %s", cluster_key)
        self.logger.info("Layer: %s", self.params.layer)
        self.logger.info("")

        # 1. Load and resolve marker sets
        self.logger.info("Phase 1: Loading marker sets...")
        marker_sets = load_marker_sets(
            self._marker_map,
            var_names=list(adata.var_names),
            logger=self.logger,
        )

        if not marker_sets:
            raise ValueError("No marker sets loaded from marker map")

        # 2. Compute IDF weights if enabled
        idf_weights = None
        if self.params.use_idf:
            self.logger.info("Phase 2: Computing IDF weights...")
            idf_weights = compute_marker_idf(marker_sets, logger=self.logger)

        # 3. Extract DE results
        self.logger.info("Phase 3: Extracting DE results...")
        de_lookup = extract_existing_de_results(
            adata, self.params.de_key, self.logger
        )

        de_rank_lookup = None
        if de_lookup:
            de_rank_lookup = build_de_rank_lookup(
                adata,
                cluster_key=cluster_key,
                de_key=self.params.de_key,
                de_top_frac=self.params.de_top_frac,
                de_min_k=self.params.de_min_k,
                de_max_k=self.params.de_max_k,
                logger=self.logger,
            )

        # 4. Compute marker document frequency for commonness penalty
        marker_doc_freq = compute_marker_document_frequency(marker_sets)

        # 5. Score clusters
        self.logger.info("Phase 4: Scoring clusters against marker sets...")
        score_result = compute_marker_scores(
            adata,
            marker_sets,
            cluster_key=cluster_key,
            positive_quantile=self.params.positive_quantile,
            de_lookup=de_lookup,
            de_bonus=self.params.de_bonus,
            anti_weight=self.params.anti_weight,
            logger=self.logger,
            idf_weights=idf_weights,
            expand_markers=self.params.expand_markers,
            anti_agg=self.params.anti_agg,
            de_rank_lookup=de_rank_lookup,
            marker_doc_freq=marker_doc_freq,
            de_commonness_alpha=self.params.de_commonness_alpha,
            layer=self.params.layer,
        )

        if self.params.expand_markers:
            marker_scores, marker_evidence = score_result
        else:
            marker_scores = score_result
            marker_evidence = None

        if marker_scores.empty:
            raise ValueError("No marker scores computed - check marker map vs panel")

        # 6. Hierarchical assignment
        self.logger.info("Phase 5: Hierarchical label assignment...")

        # Build gating params: explicit params > marker map params > defaults
        if self.params.gating_params:
            gating_params = self.params.gating_params
            self.logger.info("  Using explicit gating_params from AnnotationParams")
        elif self._marker_map_gating_params:
            gating_params = merge_gating_params(
                base=DEFAULT_GATING_PARAMS,
                tissue_params=self._marker_map_gating_params
            )
            self.logger.info("  Using gating_params merged from marker map _gating_params")
        else:
            gating_params = DEFAULT_GATING_PARAMS
            self.logger.info("  Using DEFAULT_GATING_PARAMS (no _gating_params in marker map)")

        # Log full gating parameters for debugging
        log_gating_params(gating_params, self.logger, "[Stage H]")

        cluster_annotations, decision_steps = assign_labels_hierarchical(
            marker_scores,
            adata=adata,
            cluster_key=cluster_key,
            layer=self.params.layer,
            params=gating_params,
            logger=self.logger,
            marker_sets=marker_sets,
        )

        # 7. Map to cells
        self.logger.info("Phase 6: Mapping annotations to cells...")
        annotate_obs(
            adata,
            cluster_annotations,
            cluster_key=cluster_key,
            label_col=self.params.label_col,
            logger=self.logger,
        )

        # 8. Export CSVs if output_dir provided
        if output_dir:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)

            marker_scores.to_csv(output_dir / "marker_scores.csv", index=False)
            self.logger.info("Wrote marker_scores.csv")

            cluster_annotations.to_csv(
                output_dir / "cluster_annotations.csv", index=False
            )
            self.logger.info("Wrote cluster_annotations.csv")

            decision_steps.to_csv(output_dir / "decision_steps.csv", index=False)
            self.logger.info("Wrote decision_steps.csv")

            if marker_evidence is not None and not marker_evidence.empty:
                marker_evidence.to_csv(
                    output_dir / "marker_evidence.csv", index=False
                )
                self.logger.info("Wrote marker_evidence.csv")

        self.logger.info("")
        self.logger.info("Annotation complete!")
        self.logger.info(
            "  Clusters: %d, Cell types: %d",
            len(cluster_annotations),
            len(cluster_annotations["assigned_label"].unique()),
        )

        return AnnotationResult(
            adata=adata,
            marker_scores=marker_scores,
            marker_evidence=marker_evidence,
            cluster_annotations=cluster_annotations,
            decision_steps=decision_steps,
        )

    def validate_input(self, adata: "sc.AnnData", cluster_key: str) -> List[str]:
        """Validate input AnnData has required data.

        Args:
            adata: AnnData to validate
            cluster_key: Expected cluster column

        Returns:
            List of validation errors (empty if valid)
        """
        errors = []

        # Check cluster column
        if cluster_key not in adata.obs.columns:
            errors.append(f"Missing cluster column: {cluster_key}")

        # Check layer
        if self.params.layer != "X" and self.params.layer not in adata.layers:
            errors.append(f"Missing layer: {self.params.layer}")

        # Check DE results (warning, not error)
        if self.params.de_key not in adata.uns:
            self.logger.warning(
                "DE results not found at uns['%s'] - DE bonus disabled",
                self.params.de_key,
            )

        return errors
