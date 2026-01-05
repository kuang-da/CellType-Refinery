"""ConsolidationEngine - main orchestrator for final cell-type consolidation.

This engine takes refinement output (refined.h5ad + diagnostic_report.csv)
and produces a single cell_type_phenocycler column with clean, production-ready labels.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple
import logging
import time

import pandas as pd

try:
    import scanpy as sc
except ImportError:
    sc = None

from .config import ConsolidationConfig
from .rules import (
    classify_label,
    select_final_label,
    apply_override,
    apply_relabel_rules,
    build_reason_chain,
    simplify_label,
    LabelCategory,
    ROOT_TYPES,
)
from .orphan_detection import (
    detect_orphaned_subtypes,
    apply_orphan_rescue,
    get_orphan_rescue_map,
    OrphanCandidate,
)
from .iel_rescue import (
    detect_iel_candidates,
    apply_iel_rescue,
    get_iel_rescue_map,
    IELCandidate,
    IELRescueConfig,
)
from .harmonize import (
    HarmonizeConfig,
    harmonize_adata,
)

logger = logging.getLogger(__name__)


@dataclass
class ConsolidationResult:
    """Result from consolidation execution.

    Attributes
    ----------
    success : bool
        Whether consolidation succeeded
    n_cells_total : int
        Total number of cells processed
    n_clusters_total : int
        Total number of clusters processed
    n_overrides_applied : int
        Number of manual overrides applied
    n_relabels_applied : int
        Number of relabel rules applied
    n_orphans_rescued : int
        Number of orphan clusters rescued
    n_orphans_flagged : int
        Number of orphan clusters flagged for review
    n_iel_rescued : int
        Number of IEL clusters rescued
    label_distribution : Dict[str, int]
        Simplified label distribution (label -> count)
    label_distribution_detailed : Dict[str, int]
        Detailed label distribution with orphan suffix and unsorted hybrids
    confidence_breakdown : Dict[str, int]
        Cell counts by confidence band
    category_breakdown : Dict[str, int]
        Cell counts by label category
    orphan_candidates : List[OrphanCandidate]
        All detected orphan candidates
    iel_candidates : List[IELCandidate]
        All detected IEL candidates
    errors : List[str]
        Any errors encountered
    warnings : List[str]
        Any warnings encountered
    execution_time_seconds : float
        Total execution time
    """

    success: bool = True
    n_cells_total: int = 0
    n_clusters_total: int = 0
    n_overrides_applied: int = 0
    n_relabels_applied: int = 0
    n_orphans_rescued: int = 0
    n_orphans_flagged: int = 0
    n_iel_rescued: int = 0
    label_distribution: Dict[str, int] = field(default_factory=dict)
    label_distribution_detailed: Dict[str, int] = field(default_factory=dict)
    confidence_breakdown: Dict[str, int] = field(default_factory=dict)
    category_breakdown: Dict[str, int] = field(default_factory=dict)
    orphan_candidates: List[OrphanCandidate] = field(default_factory=list)
    iel_candidates: List[IELCandidate] = field(default_factory=list)
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    execution_time_seconds: float = 0.0


class ConsolidationEngine:
    """Engine for final cell-type label consolidation.

    Takes refinement output (refined.h5ad + diagnostic_report.csv) and
    produces a single cell_type_phenocycler column with clean labels.

    Decision Logic:
    1. Start with assigned_label from diagnostic_report
    2. Apply orphan rescue for biologically plausible cases
    3. Keep root labels as honest uncertainty (don't force to subtypes)
    4. Keep hybrid labels (X~Y) as biological ambiguity
    5. Apply manual overrides from config (if provided)
    6. Apply global relabel rules from config (if provided)

    Parameters
    ----------
    config : ConsolidationConfig, optional
        Configuration with manual overrides and relabel rules
    output_col : str
        Output column name for final labels (default: cell_type_phenocycler)
    reason_col : str
        Output column name for decision reasons (default: consolidation_reason)
    enable_orphan_rescue : bool
        Whether to detect and rescue orphaned subtypes (default: True)
    orphan_suffix : str
        Suffix for rescued orphan labels (default: "(orphan)")
    enable_iel_rescue : bool
        Whether to detect and rescue IELs (default: False)
    harmonize_config : HarmonizeConfig, optional
        Configuration for two-level label harmonization (fine + broad vocab)
        If provided, adds harmonized columns to adata.obs after consolidation
    logger : logging.Logger, optional
        Logger instance
    """

    def __init__(
        self,
        config: Optional[ConsolidationConfig] = None,
        output_col: str = "cell_type_phenocycler",
        reason_col: str = "consolidation_reason",
        enable_orphan_rescue: bool = True,
        orphan_suffix: str = "(orphan)",
        enable_iel_rescue: bool = False,
        harmonize_config: Optional[HarmonizeConfig] = None,
        custom_logger: Optional[logging.Logger] = None,
    ):
        self.config = config or ConsolidationConfig()
        self.output_col = output_col
        self.reason_col = reason_col
        self.enable_orphan_rescue = enable_orphan_rescue
        self.orphan_suffix = orphan_suffix
        self.enable_iel_rescue = enable_iel_rescue or self.config.iel_rescue.enabled
        self.harmonize_config = harmonize_config
        self.logger = custom_logger or logger

        # Precompute lookup maps
        self._override_map = self.config.get_override_map()
        self._relabel_map = self.config.get_relabel_map()

    def execute(
        self,
        adata: "sc.AnnData",
        diagnostic_report: pd.DataFrame,
        marker_scores: Optional[pd.DataFrame] = None,
        cluster_col: str = "cluster_lvl1",
    ) -> ConsolidationResult:
        """Execute consolidation on AnnData.

        Modifies adata in place, adding output_col and reason_col columns.

        Parameters
        ----------
        adata : sc.AnnData
            AnnData with refined annotations
        diagnostic_report : pd.DataFrame
            Diagnostic report with columns: cluster_id, assigned_label,
            assigned_score, n_cells, confidence_band, recommendation
        marker_scores : pd.DataFrame, optional
            Marker scores for orphan detection. Required if enable_orphan_rescue=True.
        cluster_col : str
            Column in adata.obs containing cluster IDs (default: cluster_lvl1)

        Returns
        -------
        ConsolidationResult
            Execution result with statistics and any errors/warnings
        """
        start_time = time.time()
        result = ConsolidationResult()

        try:
            self.logger.info("Starting consolidation...")

            # Validate inputs
            self._validate_inputs(adata, diagnostic_report, cluster_col, result)
            if not result.success:
                return result

            # Build cluster to final mapping
            cluster_mapping = self._build_cluster_mapping(
                diagnostic_report,
                marker_scores,
                result,
                adata=adata,
                cluster_col=cluster_col,
            )

            # Apply mapping to AnnData
            self._apply_mapping_to_adata(
                adata,
                cluster_mapping,
                cluster_col,
                result,
            )

            # Apply harmonization if configured
            if self.harmonize_config is not None:
                self.logger.info("Applying two-level label harmonization...")
                harmonize_adata(
                    adata,
                    diagnostic_report,
                    self.harmonize_config,
                    label_col=self.output_col,
                    cluster_col=cluster_col,
                )

            # Compute statistics
            self._compute_statistics(adata, result)

            result.success = True
            self.logger.info(
                f"Consolidation complete: {result.n_cells_total:,} cells, "
                f"{result.n_clusters_total} clusters"
            )

        except Exception as e:
            result.success = False
            result.errors.append(f"Consolidation failed: {str(e)}")
            self.logger.error(f"Consolidation failed: {e}", exc_info=True)

        result.execution_time_seconds = time.time() - start_time
        return result

    def _validate_inputs(
        self,
        adata: "sc.AnnData",
        diagnostic_report: pd.DataFrame,
        cluster_col: str,
        result: ConsolidationResult,
    ) -> None:
        """Validate input data."""
        # Check required columns in diagnostic report
        required_cols = {"cluster_id", "assigned_label"}
        missing = required_cols - set(diagnostic_report.columns)
        if missing:
            result.success = False
            result.errors.append(f"Missing columns in diagnostic_report: {missing}")
            return

        # Check cluster column in adata
        if cluster_col not in adata.obs.columns:
            # Try to find alternative
            alternatives = ["cluster_lvl0", "cluster_id", "leiden"]
            found = None
            for alt in alternatives:
                if alt in adata.obs.columns:
                    found = alt
                    break
            if found:
                result.warnings.append(
                    f"Column '{cluster_col}' not found, using '{found}' instead"
                )
            else:
                result.success = False
                result.errors.append(
                    f"Cluster column '{cluster_col}' not found in adata.obs"
                )
                return

    def _build_cluster_mapping(
        self,
        diagnostic_report: pd.DataFrame,
        marker_scores: Optional[pd.DataFrame],
        result: ConsolidationResult,
        adata: Optional["sc.AnnData"] = None,
        cluster_col: str = "cluster_lvl1",
    ) -> Dict[str, Tuple[str, str, str]]:
        """Build mapping from cluster_id to (final_label, reason, confidence_band).

        Parameters
        ----------
        diagnostic_report : pd.DataFrame
            Diagnostic report
        marker_scores : pd.DataFrame, optional
            Marker scores for orphan detection
        result : ConsolidationResult
            Result object for tracking statistics
        adata : sc.AnnData, optional
            AnnData for IEL rescue (requires decision_steps in uns)
        cluster_col : str
            Column name for cluster IDs in adata.obs

        Returns
        -------
        Dict[str, Tuple[str, str, str]]
            Mapping: cluster_id -> (final_label, reason, confidence_band)
        """
        cluster_mapping: Dict[str, Tuple[str, str, str]] = {}

        # Step 1: Detect orphans if enabled
        orphan_rescue_map: Dict[str, Tuple[str, str]] = {}
        if self.enable_orphan_rescue and marker_scores is not None:
            orphan_rescue_map = self._detect_and_rescue_orphans(
                diagnostic_report,
                marker_scores,
                result,
            )

        # Step 1b: Detect IELs if enabled
        iel_rescue_map: Dict[str, Tuple[str, str]] = {}
        if self.enable_iel_rescue and adata is not None and marker_scores is not None:
            iel_rescue_map = self._detect_and_rescue_iels(
                adata,
                diagnostic_report,
                marker_scores,
                cluster_col,
                result,
            )

        # Step 2: Build base mapping from diagnostic report
        for _, row in diagnostic_report.iterrows():
            cluster_id = str(row["cluster_id"])
            assigned_label = str(row.get("assigned_label", "Unassigned"))
            confidence_band = str(row.get("confidence_band", "unknown"))
            recommendation = str(row.get("recommendation", "SKIP"))

            # Check if this cluster was rescued (IEL takes priority over orphan)
            if cluster_id in iel_rescue_map:
                final_label, iel_reason = iel_rescue_map[cluster_id]
                base_reason = iel_reason
            elif cluster_id in orphan_rescue_map:
                final_label, orphan_reason = orphan_rescue_map[cluster_id]
                base_reason = orphan_reason
            else:
                # Apply base decision logic
                final_label, base_reason = select_final_label(
                    cluster_id,
                    assigned_label,
                    confidence_band,
                    recommendation,
                )

            # Step 3: Apply manual override
            final_label, override_reason = apply_override(
                cluster_id,
                final_label,
                self._override_map,
            )
            if override_reason:
                result.n_overrides_applied += 1

            # Step 4: Apply relabel rules
            final_label, relabel_reason = apply_relabel_rules(
                final_label,
                self._relabel_map,
            )
            if relabel_reason:
                result.n_relabels_applied += 1

            # Build combined reason
            reason = build_reason_chain(
                base_reason,
                override_reason=override_reason,
                relabel_reason=relabel_reason,
            )

            cluster_mapping[cluster_id] = (final_label, reason, confidence_band)

        result.n_clusters_total = len(cluster_mapping)
        return cluster_mapping

    def _detect_and_rescue_orphans(
        self,
        diagnostic_report: pd.DataFrame,
        marker_scores: pd.DataFrame,
        result: ConsolidationResult,
    ) -> Dict[str, Tuple[str, str]]:
        """Detect and rescue orphaned subtypes.

        Parameters
        ----------
        diagnostic_report : pd.DataFrame
            Diagnostic report
        marker_scores : pd.DataFrame
            Marker scores
        result : ConsolidationResult
            Result object for tracking statistics

        Returns
        -------
        Dict[str, Tuple[str, str]]
            Mapping: cluster_id -> (final_label, reason) for rescued orphans
        """
        # Find Unassigned clusters
        unassigned_mask = diagnostic_report["assigned_label"].str.lower() == "unassigned"
        unassigned_ids = diagnostic_report.loc[unassigned_mask, "cluster_id"].astype(str).tolist()

        if not unassigned_ids:
            self.logger.info("No Unassigned clusters found for orphan detection")
            return {}

        self.logger.info(f"Detecting orphans in {len(unassigned_ids)} Unassigned clusters...")

        # Detect orphans
        candidates = detect_orphaned_subtypes(
            marker_scores,
            unassigned_ids,
            config=self.config.orphan_rescue,
        )

        if not candidates:
            self.logger.info("No orphan candidates detected")
            return {}

        self.logger.info(f"Found {len(candidates)} orphan candidates")

        # Apply rescue logic
        candidates = apply_orphan_rescue(
            candidates,
            suffix=self.orphan_suffix,
            config=self.config.orphan_rescue,
        )

        # Track results
        result.orphan_candidates = candidates
        result.n_orphans_rescued = sum(
            1 for c in candidates if c.final_label and "Unassigned" not in c.final_label
        )
        result.n_orphans_flagged = sum(
            1 for c in candidates if c.flag
        )

        self.logger.info(
            f"Orphan rescue: {result.n_orphans_rescued} rescued, "
            f"{result.n_orphans_flagged} flagged for review"
        )

        # Return rescue map
        return get_orphan_rescue_map(candidates)

    def _detect_and_rescue_iels(
        self,
        adata: "sc.AnnData",
        diagnostic_report: pd.DataFrame,
        marker_scores: pd.DataFrame,
        cluster_col: str,
        result: ConsolidationResult,
    ) -> Dict[str, Tuple[str, str]]:
        """Detect and rescue intraepithelial immune cells (IELs).

        IELs are immune cells that were vetoed by epithelial markers
        but are likely legitimate intraepithelial lymphocytes or
        tissue-resident macrophages.

        Parameters
        ----------
        adata : sc.AnnData
            AnnData with decision_steps_subcluster in uns
        diagnostic_report : pd.DataFrame
            Diagnostic report
        marker_scores : pd.DataFrame
            Marker scores
        cluster_col : str
            Column name for cluster IDs in adata.obs
        result : ConsolidationResult
            Result object for tracking statistics

        Returns
        -------
        Dict[str, Tuple[str, str]]
            Mapping: cluster_id -> (final_label, reason) for rescued IELs
        """
        self.logger.info("Detecting IEL candidates...")

        # Use config from consolidation config
        iel_config = IELRescueConfig(
            enabled=True,
            cd45_min_pos_frac=self.config.iel_rescue.cd45_min_pos_frac,
            lymphoid_score_threshold=self.config.iel_rescue.lymphoid_score_threshold,
            myeloid_score_threshold=self.config.iel_rescue.myeloid_score_threshold,
            suffix=self.config.iel_rescue.suffix,
        )

        # Detect IEL candidates
        candidates = detect_iel_candidates(
            adata=adata,
            diagnostic_report=diagnostic_report,
            marker_scores=marker_scores,
            cluster_col=cluster_col,
            layer=None,
            config=iel_config,
        )

        if not candidates:
            self.logger.info("No IEL candidates detected")
            return {}

        self.logger.info(f"Found {len(candidates)} IEL candidates")

        # Apply rescue logic
        candidates = apply_iel_rescue(candidates, config=iel_config)

        # Track results
        result.iel_candidates = candidates
        result.n_iel_rescued = len(candidates)

        # Log summary by type
        for c in candidates:
            self.logger.info(
                f"  IEL: {c.cluster_id} ({c.n_cells} cells) -> {c.final_label}"
            )

        self.logger.info(f"IEL rescue: {result.n_iel_rescued} clusters rescued")

        # Return rescue map
        return get_iel_rescue_map(candidates)

    def _apply_mapping_to_adata(
        self,
        adata: "sc.AnnData",
        cluster_mapping: Dict[str, Tuple[str, str, str]],
        cluster_col: str,
        result: ConsolidationResult,
    ) -> None:
        """Apply cluster mapping to AnnData.

        Parameters
        ----------
        adata : sc.AnnData
            AnnData object to modify
        cluster_mapping : Dict[str, Tuple[str, str, str]]
            Mapping: cluster_id -> (final_label, reason, confidence_band)
        cluster_col : str
            Column containing cluster IDs
        result : ConsolidationResult
            Result object for tracking
        """
        # Get cluster IDs from adata
        cluster_ids = adata.obs[cluster_col].astype(str)

        # Create label, reason, and confidence arrays
        final_labels = []
        reasons = []
        confidence_bands = []
        missing_clusters = set()

        for cid in cluster_ids:
            if cid in cluster_mapping:
                label, reason, conf_band = cluster_mapping[cid]
            else:
                label = "Unassigned"
                reason = "Cluster not in diagnostic report"
                conf_band = "unknown"
                missing_clusters.add(cid)
            final_labels.append(label)
            reasons.append(reason)
            confidence_bands.append(conf_band)

        # Create simplified labels (strip orphan suffix, sort hybrids)
        simplified_labels = [simplify_label(lbl) for lbl in final_labels]

        # Assign to adata - DETAILED column first (original labels with orphan suffix)
        adata.obs[f"{self.output_col}_detailed"] = pd.Categorical(final_labels)
        # Then simplified column (atlas-ready labels)
        adata.obs[self.output_col] = pd.Categorical(simplified_labels)
        adata.obs[self.reason_col] = reasons
        adata.obs["confidence_band"] = pd.Categorical(confidence_bands)

        result.n_cells_total = len(adata)

        if missing_clusters:
            result.warnings.append(
                f"{len(missing_clusters)} clusters in AnnData not found in diagnostic_report"
            )
            self.logger.warning(
                f"{len(missing_clusters)} clusters not in diagnostic_report: "
                f"{list(missing_clusters)[:5]}..."
            )

    def _compute_statistics(
        self,
        adata: "sc.AnnData",
        result: ConsolidationResult,
    ) -> None:
        """Compute summary statistics.

        Parameters
        ----------
        adata : sc.AnnData
            AnnData with final labels
        result : ConsolidationResult
            Result object to populate
        """
        # Simplified label distribution (for cell_type_phenocycler)
        label_counts = adata.obs[self.output_col].value_counts()
        result.label_distribution = label_counts.to_dict()

        # Detailed label distribution (for cell_type_phenocycler_detailed)
        detailed_col = f"{self.output_col}_detailed"
        if detailed_col in adata.obs.columns:
            detailed_counts = adata.obs[detailed_col].value_counts()
            result.label_distribution_detailed = detailed_counts.to_dict()

        # Category breakdown
        category_counts: Dict[str, int] = {
            "subtype": 0,
            "root": 0,
            "hybrid": 0,
            "unassigned": 0,
            "mixed": 0,
        }
        for label, count in label_counts.items():
            cat = classify_label(label)
            category_counts[cat.value] += count
        result.category_breakdown = category_counts

        # Confidence breakdown (if available in adata)
        if "confidence_band" in adata.obs.columns:
            conf_counts = adata.obs["confidence_band"].value_counts()
            result.confidence_breakdown = conf_counts.to_dict()

    def get_mapping_table(
        self,
        diagnostic_report: pd.DataFrame,
        marker_scores: Optional[pd.DataFrame] = None,
    ) -> pd.DataFrame:
        """Get the cluster to final label mapping table without modifying AnnData.

        Useful for previewing changes before applying.

        Parameters
        ----------
        diagnostic_report : pd.DataFrame
            Diagnostic report
        marker_scores : pd.DataFrame, optional
            Marker scores for orphan detection

        Returns
        -------
        pd.DataFrame
            Mapping table with columns: cluster_id, assigned_label, final_label,
            reason, category, confidence_band, n_cells
        """
        result = ConsolidationResult()
        cluster_mapping = self._build_cluster_mapping(
            diagnostic_report,
            marker_scores,
            result,
        )

        rows = []
        for _, row in diagnostic_report.iterrows():
            cluster_id = str(row["cluster_id"])
            assigned_label = str(row.get("assigned_label", "Unassigned"))
            n_cells = int(row.get("n_cells", 0))

            if cluster_id in cluster_mapping:
                final_label, reason, conf_band = cluster_mapping[cluster_id]
            else:
                final_label = "Unassigned"
                reason = "Not mapped"
                conf_band = "unknown"

            category = classify_label(final_label)

            rows.append({
                "cluster_id": cluster_id,
                "assigned_label": assigned_label,
                "final_label": final_label,
                "reason": reason,
                "category": category.value,
                "confidence_band": conf_band,
                "n_cells": n_cells,
            })

        return pd.DataFrame(rows)
