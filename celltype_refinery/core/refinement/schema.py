"""
CanonicalSchema and InputAdapter for refinement input normalization.

This module provides:
- CanonicalSchema: Defines expected column names and keys for inputs
- InputAdapter: Normalizes clustering outputs to match canonical schema

The adapter supports:
- Auto-detection of column names via alias matching
- Explicit CLI overrides for edge cases
- Path separator normalization
- Provenance type detection
- Actionable validation error messages
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Union

import pandas as pd

try:
    import scanpy as sc
except ImportError:
    sc = None

from .validation import (
    ValidationResult,
    ValidationError,
    MissingColumnError,
    InvalidClusterIdError,
    MissingProvenanceError,
    PathSeparatorError,
)


@dataclass
class CanonicalSchema:
    """Canonical schema for refinement inputs.

    Defines the expected column names for AnnData.obs, adata.uns, and CSV files.
    All downstream refinement code should use these canonical names after adaptation.

    Attributes
    ----------
    cluster_key : str
        Canonical name for cluster column in adata.obs (default: "cluster_lvl0")
    cluster_key_refined : str
        Canonical name for refined clusters (default: "cluster_lvl1")
    label_key_in : str
        Input label column (default: "cell_type_auto")
    label_key_out : str
        Output label column for refinement (default: "cell_type_lvl1")
    score_key : str
        Confidence score column (default: "cell_type_auto_score")
    root_key : str
        Root category column (default: "cell_type_auto_root")
    provenance_key_cluster : str
        Clustering stage provenance key in adata.uns
    provenance_key_refine_prefix : str
        Prefix for refinement iteration provenance keys
    csv_* : str
        Canonical column names for cluster_annotations.csv and marker_scores.csv
    path_separator : str
        Canonical separator for hierarchical paths (default: " / ")
    """

    # AnnData.obs columns
    cluster_key: str = "cluster_lvl0"
    cluster_key_refined: str = "cluster_lvl1"
    label_key_in: str = "cell_type_auto"
    label_key_out: str = "cell_type_lvl1"
    score_key: str = "cell_type_auto_score"
    root_key: str = "cell_type_auto_root"

    # Provenance keys in adata.uns
    provenance_key_cluster: str = "clustering_provenance"
    provenance_key_refine_prefix: str = "refinement_iteration_"

    # cluster_annotations.csv columns
    csv_cluster_id: str = "cluster_id"
    csv_assigned_label: str = "assigned_label"
    csv_assigned_score: str = "assigned_score"
    csv_n_cells: str = "n_cells"

    # marker_scores.csv columns
    csv_ms_cluster_id: str = "cluster_id"
    csv_ms_label: str = "label"
    csv_ms_path: str = "path"
    csv_ms_score: str = "score"

    # Path hierarchy separator
    path_separator: str = " / "

    @classmethod
    def with_overrides(cls, **kwargs) -> "CanonicalSchema":
        """Create schema with custom overrides.

        Parameters
        ----------
        **kwargs
            Attribute overrides

        Returns
        -------
        CanonicalSchema
            Schema with overridden values
        """
        return cls(**kwargs)


class InputAdapter:
    """Adapter to normalize clustering outputs to canonical schema.

    The adapter auto-detects column names via alias matching and normalizes
    them to the canonical schema. It supports explicit CLI overrides for
    edge cases where auto-detection fails.

    Features:
    - Auto-detect column names via alias matching (first match wins)
    - Normalize path separators ("|", "::", ">" → " / ")
    - Validate required columns with actionable errors
    - Support explicit CLI overrides
    - Detect provenance type (clustering vs refinement iteration)

    Example:
        >>> adapter = InputAdapter(cluster_key="leiden")  # explicit override
        >>> result = adapter.validate_adata(adata)
        >>> if result.is_valid:
        ...     cluster_col = result.adapted_columns["cluster_lvl0"]
        ... else:
        ...     result.raise_if_invalid("input AnnData")
        >>>
        >>> # Adapt CSV files
        >>> annotations_df = adapter.adapt_cluster_annotations(raw_df)
        >>> scores_df = adapter.adapt_marker_scores(raw_df)
    """

    # Alias maps for auto-detection (order matters - first match wins)
    CLUSTER_KEY_ALIASES: List[str] = [
        "cluster_lvl0",
        "cluster",
        "leiden",
        "leiden_h",
        "cluster_h",
        "louvain",
        "cluster_id",
    ]

    CLUSTER_KEY_REFINED_ALIASES: List[str] = [
        "cluster_lvl1",
        "cluster_refined",
        "subcluster",
        "leiden_refined",
    ]

    LABEL_KEY_ALIASES: List[str] = [
        "cell_type_auto",
        "cell_type_lvl0_auto",
        "cell_type",
        "assigned_label",
        "cell_type_lvl0",
        "annotation",
    ]

    SCORE_KEY_ALIASES: List[str] = [
        "cell_type_auto_score",
        "cell_type_lvl0_auto_score",
        "assigned_score",
        "score",
        "confidence",
        "confidence_score",
        "annotation_score",
    ]

    # CSV column aliases for cluster_annotations.csv
    CSV_CLUSTER_ID_ALIASES: List[str] = [
        "cluster_id",
        "cluster",
        "Cluster",
        "ClusterID",
        "cluster_idx",
    ]

    CSV_ASSIGNED_LABEL_ALIASES: List[str] = [
        "assigned_label",
        "label",
        "cell_type",
        "annotation",
        "AssignedLabel",
        "assigned_type",
    ]

    CSV_ASSIGNED_SCORE_ALIASES: List[str] = [
        "assigned_score",
        "score",
        "confidence",
        "AssignedScore",
        "annotation_score",
    ]

    CSV_N_CELLS_ALIASES: List[str] = [
        "n_cells",
        "ncells",
        "num_cells",
        "cell_count",
        "size",
        "count",
    ]

    # CSV column aliases for marker_scores.csv
    CSV_MS_LABEL_ALIASES: List[str] = [
        "label",
        "cell_type",
        "Label",
        "CellType",
    ]

    CSV_MS_PATH_ALIASES: List[str] = [
        "path",
        "hierarchy",
        "hierarchy_path",
        "full_path",
        "Path",
    ]

    CSV_MS_SCORE_ALIASES: List[str] = [
        "score",
        "Score",
        "marker_score",
        "final_score",
    ]

    # Supported path separators for normalization
    PATH_SEPARATORS: List[str] = [" / ", " | ", "|", "::", " > ", ">", " → ", "→"]
    CANONICAL_PATH_SEPARATOR: str = " / "

    def __init__(
        self,
        schema: Optional[CanonicalSchema] = None,
        cluster_key: Optional[str] = None,
        label_key_in: Optional[str] = None,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize adapter with optional explicit overrides.

        Parameters
        ----------
        schema : CanonicalSchema, optional
            Custom schema (uses default if not provided)
        cluster_key : str, optional
            Explicit cluster column name (overrides auto-detection).
            Use this when auto-detection fails or you want to force a specific column.
        label_key_in : str, optional
            Explicit label column name (overrides auto-detection).
        logger : logging.Logger, optional
            Logger for warnings (uses module logger if not provided)
        """
        self.schema = schema or CanonicalSchema()
        self.explicit_cluster_key = cluster_key
        self.explicit_label_key_in = label_key_in
        self.logger = logger or logging.getLogger(__name__)
        self._warnings: List[str] = []

    @property
    def warnings(self) -> List[str]:
        """Get accumulated warnings from adaptation."""
        return self._warnings.copy()

    def clear_warnings(self) -> None:
        """Clear accumulated warnings."""
        self._warnings = []

    # -------------------------------------------------------------------------
    # Core validation methods
    # -------------------------------------------------------------------------

    def validate_adata(self, adata: "sc.AnnData") -> ValidationResult:
        """Validate AnnData.obs has required columns.

        Checks for required columns (cluster_key) and optional columns
        (label_key_in, cluster_key_refined). Returns adapted column mapping.

        Parameters
        ----------
        adata : sc.AnnData
            AnnData object to validate

        Returns
        -------
        ValidationResult
            Result with adapted_columns mapping canonical → actual names,
            plus any errors or warnings

        Example
        -------
        >>> result = adapter.validate_adata(adata)
        >>> if result.is_valid:
        ...     cluster_col = result.adapted_columns["cluster_lvl0"]
        ... else:
        ...     for error in result.errors:
        ...         print(error)
        """
        errors: List[ValidationError] = []
        warnings: List[str] = []
        adapted_columns: Dict[str, str] = {}
        available = list(adata.obs.columns)

        # Required: cluster key
        cluster_col = self._find_column_by_alias(
            adata.obs,
            self.schema.cluster_key,
            self.CLUSTER_KEY_ALIASES,
            explicit_override=self.explicit_cluster_key,
        )
        if cluster_col is None:
            errors.append(
                MissingColumnError(
                    message="Required cluster column not found in adata.obs",
                    expected=f"One of: {self.CLUSTER_KEY_ALIASES}",
                    found=f"Available columns: {available[:20]}{'...' if len(available) > 20 else ''}",
                    available_columns=available,
                    suggestion="Use --cluster-key to specify the cluster column explicitly.",
                )
            )
        else:
            adapted_columns[self.schema.cluster_key] = cluster_col
            if cluster_col != self.schema.cluster_key:
                warnings.append(
                    f"Auto-mapped cluster column: '{cluster_col}' → '{self.schema.cluster_key}'"
                )

        # Optional: label key (may not exist for first-time input)
        label_col = self._find_column_by_alias(
            adata.obs,
            self.schema.label_key_in,
            self.LABEL_KEY_ALIASES,
            explicit_override=self.explicit_label_key_in,
        )
        if label_col is not None:
            adapted_columns[self.schema.label_key_in] = label_col
            if label_col != self.schema.label_key_in:
                warnings.append(
                    f"Auto-mapped label column: '{label_col}' → '{self.schema.label_key_in}'"
                )

        # Optional: score key
        score_col = self._find_column_by_alias(
            adata.obs,
            self.schema.score_key,
            self.SCORE_KEY_ALIASES,
        )
        if score_col is not None:
            adapted_columns[self.schema.score_key] = score_col

        # Optional: refined cluster key
        refined_col = self._find_column_by_alias(
            adata.obs,
            self.schema.cluster_key_refined,
            self.CLUSTER_KEY_REFINED_ALIASES,
        )
        if refined_col is not None:
            adapted_columns[self.schema.cluster_key_refined] = refined_col

        # Check provenance
        provenance_type = self.detect_provenance_type(adata)
        if provenance_type == "unknown":
            warnings.append(
                "No clustering or refinement provenance found in adata.uns. "
                "This may be acceptable for fresh inputs."
            )

        self._warnings.extend(warnings)

        return ValidationResult(
            is_valid=len(errors) == 0,
            errors=errors,
            warnings=warnings,
            adapted_columns=adapted_columns,
            missing_columns=[self.schema.cluster_key] if cluster_col is None else [],
            suggested_fixes=[e.suggestion for e in errors if e.suggestion],
        )

    def validate_cluster_annotations(self, df: pd.DataFrame) -> ValidationResult:
        """Validate cluster_annotations.csv has required columns.

        Required columns: cluster_id, assigned_label, assigned_score, n_cells

        Parameters
        ----------
        df : pd.DataFrame
            Raw cluster_annotations DataFrame to validate

        Returns
        -------
        ValidationResult
            Result with adapted_columns mapping and any errors
        """
        errors: List[ValidationError] = []
        warnings: List[str] = []
        adapted_columns: Dict[str, str] = {}
        available = list(df.columns)

        # Required columns mapping: (canonical_name, aliases, required)
        required_mappings = [
            (self.schema.csv_cluster_id, self.CSV_CLUSTER_ID_ALIASES, True),
            (self.schema.csv_assigned_label, self.CSV_ASSIGNED_LABEL_ALIASES, True),
            (self.schema.csv_assigned_score, self.CSV_ASSIGNED_SCORE_ALIASES, True),
            (self.schema.csv_n_cells, self.CSV_N_CELLS_ALIASES, True),
        ]

        for canonical, aliases, required in required_mappings:
            actual = self._find_column_by_alias(df, canonical, aliases)
            if actual is None:
                if required:
                    errors.append(
                        MissingColumnError(
                            message=f"Required column '{canonical}' not found in cluster_annotations.csv",
                            expected=f"One of: {aliases}",
                            found=f"Available: {available}",
                            available_columns=available,
                        )
                    )
            else:
                adapted_columns[canonical] = actual
                if actual != canonical:
                    warnings.append(
                        f"Auto-mapped CSV column: '{actual}' → '{canonical}'"
                    )

        self._warnings.extend(warnings)

        return ValidationResult(
            is_valid=len(errors) == 0,
            errors=errors,
            warnings=warnings,
            adapted_columns=adapted_columns,
        )

    def validate_marker_scores(self, df: pd.DataFrame) -> ValidationResult:
        """Validate marker_scores.csv has required columns.

        Required columns: cluster_id, label, path, score

        Parameters
        ----------
        df : pd.DataFrame
            Raw marker_scores DataFrame to validate

        Returns
        -------
        ValidationResult
            Result with adapted_columns mapping and any errors
        """
        errors: List[ValidationError] = []
        warnings: List[str] = []
        adapted_columns: Dict[str, str] = {}
        available = list(df.columns)

        # Required columns mapping
        required_mappings = [
            (self.schema.csv_ms_cluster_id, self.CSV_CLUSTER_ID_ALIASES, True),
            (self.schema.csv_ms_label, self.CSV_MS_LABEL_ALIASES, True),
            (self.schema.csv_ms_path, self.CSV_MS_PATH_ALIASES, True),
            (self.schema.csv_ms_score, self.CSV_MS_SCORE_ALIASES, True),
        ]

        for canonical, aliases, required in required_mappings:
            actual = self._find_column_by_alias(df, canonical, aliases)
            if actual is None:
                if required:
                    errors.append(
                        MissingColumnError(
                            message=f"Required column '{canonical}' not found in marker_scores.csv",
                            expected=f"One of: {aliases}",
                            found=f"Available: {available}",
                            available_columns=available,
                        )
                    )
            else:
                adapted_columns[canonical] = actual
                if actual != canonical:
                    warnings.append(
                        f"Auto-mapped CSV column: '{actual}' → '{canonical}'"
                    )

        # Check path separator
        if self.schema.csv_ms_path in adapted_columns:
            path_col = adapted_columns[self.schema.csv_ms_path]
            sample_paths = df[path_col].dropna().head(10).tolist()
            detected_sep = self._detect_path_separator(sample_paths)
            if detected_sep and detected_sep != self.CANONICAL_PATH_SEPARATOR:
                warnings.append(
                    f"Non-canonical path separator detected: '{detected_sep}' "
                    f"(will normalize to '{self.CANONICAL_PATH_SEPARATOR}')"
                )

        self._warnings.extend(warnings)

        return ValidationResult(
            is_valid=len(errors) == 0,
            errors=errors,
            warnings=warnings,
            adapted_columns=adapted_columns,
        )

    def validate_cluster_ids(
        self,
        cluster_ids: List[str],
        adata: "sc.AnnData",
        context: str = "",
    ) -> ValidationResult:
        """Validate cluster IDs exist in AnnData.

        Checks if all provided cluster IDs are present in either cluster_lvl0
        or cluster_lvl1 columns.

        Parameters
        ----------
        cluster_ids : List[str]
            Cluster IDs to validate
        adata : sc.AnnData
            AnnData object to check against
        context : str
            Context for error messages (e.g., "subcluster rule")

        Returns
        -------
        ValidationResult
            Result with any invalid cluster ID errors
        """
        errors: List[ValidationError] = []

        # Get valid cluster IDs from both levels
        adata_result = self.validate_adata(adata)
        if not adata_result.is_valid:
            return adata_result

        cluster_col = adata_result.adapted_columns.get(
            self.schema.cluster_key, self.schema.cluster_key
        )
        refined_col = adata_result.adapted_columns.get(self.schema.cluster_key_refined)

        valid_ids: Set[str] = set(adata.obs[cluster_col].astype(str).unique())
        if refined_col and refined_col in adata.obs:
            valid_ids |= set(adata.obs[refined_col].astype(str).unique())

        invalid_ids = [cid for cid in cluster_ids if str(cid) not in valid_ids]

        if invalid_ids:
            for cid in invalid_ids:
                errors.append(
                    InvalidClusterIdError(
                        message=f"Cluster ID '{cid}' not found{' in ' + context if context else ''}",
                        cluster_id=str(cid),
                        expected="Valid cluster ID from adata.obs",
                        found=f"'{cid}' not in data",
                        valid_ids=sorted(valid_ids)[:20],
                    )
                )

        return ValidationResult(
            is_valid=len(errors) == 0,
            errors=errors,
        )

    # -------------------------------------------------------------------------
    # Adaptation methods
    # -------------------------------------------------------------------------

    def adapt_cluster_annotations(
        self,
        df: pd.DataFrame,
        validate: bool = True,
    ) -> pd.DataFrame:
        """Normalize cluster_annotations.csv to canonical column names.

        Renames columns to match canonical schema and ensures cluster_id is string.

        Parameters
        ----------
        df : pd.DataFrame
            Raw cluster_annotations DataFrame
        validate : bool
            Whether to validate first (raises on error)

        Returns
        -------
        pd.DataFrame
            DataFrame with canonical column names

        Raises
        ------
        ValueError
            If validate=True and validation fails
        """
        if validate:
            result = self.validate_cluster_annotations(df)
            result.raise_if_invalid("cluster_annotations.csv")
        else:
            result = self.validate_cluster_annotations(df)

        df_adapted = df.copy()

        # Rename columns to canonical
        rename_map = {}
        for canonical, actual in result.adapted_columns.items():
            if actual != canonical and actual in df_adapted.columns:
                rename_map[actual] = canonical

        if rename_map:
            df_adapted = df_adapted.rename(columns=rename_map)

        # Ensure cluster_id is string
        if self.schema.csv_cluster_id in df_adapted.columns:
            df_adapted[self.schema.csv_cluster_id] = df_adapted[
                self.schema.csv_cluster_id
            ].astype(str)

        return df_adapted

    def adapt_marker_scores(
        self,
        df: pd.DataFrame,
        validate: bool = True,
        normalize_paths: bool = True,
    ) -> pd.DataFrame:
        """Normalize marker_scores.csv to canonical column names and path separator.

        Renames columns and normalizes path separators to " / ".

        Parameters
        ----------
        df : pd.DataFrame
            Raw marker_scores DataFrame
        validate : bool
            Whether to validate first (raises on error)
        normalize_paths : bool
            Whether to normalize path separators

        Returns
        -------
        pd.DataFrame
            DataFrame with canonical column names and normalized paths

        Raises
        ------
        ValueError
            If validate=True and validation fails
        """
        if validate:
            result = self.validate_marker_scores(df)
            result.raise_if_invalid("marker_scores.csv")
        else:
            result = self.validate_marker_scores(df)

        df_adapted = df.copy()

        # Rename columns
        rename_map = {}
        for canonical, actual in result.adapted_columns.items():
            if actual != canonical and actual in df_adapted.columns:
                rename_map[actual] = canonical

        if rename_map:
            df_adapted = df_adapted.rename(columns=rename_map)

        # Ensure cluster_id is string
        if self.schema.csv_ms_cluster_id in df_adapted.columns:
            df_adapted[self.schema.csv_ms_cluster_id] = df_adapted[
                self.schema.csv_ms_cluster_id
            ].astype(str)

        # Normalize path separator
        if normalize_paths and self.schema.csv_ms_path in df_adapted.columns:
            df_adapted[self.schema.csv_ms_path] = df_adapted[
                self.schema.csv_ms_path
            ].apply(self.normalize_path_separator)

        return df_adapted

    def get_cluster_column(
        self,
        adata: "sc.AnnData",
        prefer_refined: bool = True,
    ) -> str:
        """Get the actual cluster column name from AnnData.

        Parameters
        ----------
        adata : sc.AnnData
            AnnData object
        prefer_refined : bool
            If True, return cluster_lvl1 if available, else cluster_lvl0

        Returns
        -------
        str
            Actual column name to use

        Raises
        ------
        ValueError
            If no cluster column can be found
        """
        result = self.validate_adata(adata)
        result.raise_if_invalid("adata")

        if prefer_refined:
            refined_col = result.adapted_columns.get(self.schema.cluster_key_refined)
            if refined_col and refined_col in adata.obs:
                return refined_col

        return result.adapted_columns.get(
            self.schema.cluster_key, self.schema.cluster_key
        )

    def get_label_column(
        self,
        adata: "sc.AnnData",
    ) -> Optional[str]:
        """Get the actual label column name from AnnData.

        Parameters
        ----------
        adata : sc.AnnData
            AnnData object

        Returns
        -------
        str or None
            Actual column name, or None if not found
        """
        result = self.validate_adata(adata)
        return result.adapted_columns.get(self.schema.label_key_in)

    # -------------------------------------------------------------------------
    # Utility methods
    # -------------------------------------------------------------------------

    def detect_provenance_type(self, adata: "sc.AnnData") -> str:
        """Detect provenance type from adata.uns.

        Used to determine whether input is from clustering or a previous
        refinement iteration.

        Parameters
        ----------
        adata : sc.AnnData
            AnnData object to check

        Returns
        -------
        str
            "clustering" - Direct clustering output
            "refinement" - Chained refinement output
            "unknown" - No provenance found
        """
        refine_keys = [
            k
            for k in adata.uns.keys()
            if k.startswith(self.schema.provenance_key_refine_prefix)
        ]

        if refine_keys:
            return "refinement"
        elif self.schema.provenance_key_cluster in adata.uns:
            return "clustering"
        # Check for common alternative provenance keys
        elif "stage_h_provenance" in adata.uns:
            return "clustering"
        else:
            return "unknown"

    def get_refinement_iteration(self, adata: "sc.AnnData") -> int:
        """Get the current refinement iteration number.

        Parameters
        ----------
        adata : sc.AnnData
            AnnData object to check

        Returns
        -------
        int
            Iteration number (0 if clustering, 1+ for refinements)
        """
        refine_keys = [
            k
            for k in adata.uns.keys()
            if k.startswith(self.schema.provenance_key_refine_prefix)
        ]

        if not refine_keys:
            return 0

        # Extract iteration numbers and find max
        iterations = []
        for key in refine_keys:
            suffix = key[len(self.schema.provenance_key_refine_prefix):]
            try:
                iterations.append(int(suffix))
            except ValueError:
                continue

        return max(iterations) if iterations else 0

    def normalize_path_separator(self, path: Union[str, Any]) -> Union[str, Any]:
        """Normalize path separator to canonical format.

        Converts alternative separators ("|", "::", ">") to canonical " / ".

        Parameters
        ----------
        path : str or Any
            Path string to normalize (non-strings passed through)

        Returns
        -------
        str or Any
            Normalized path string
        """
        if pd.isna(path):
            return path
        if not isinstance(path, str):
            return path

        path_str = str(path)
        for sep in self.PATH_SEPARATORS:
            if sep != self.CANONICAL_PATH_SEPARATOR and sep in path_str:
                path_str = path_str.replace(sep, self.CANONICAL_PATH_SEPARATOR)

        return path_str

    def _find_column_by_alias(
        self,
        df_or_obs: pd.DataFrame,
        canonical: str,
        aliases: List[str],
        explicit_override: Optional[str] = None,
    ) -> Optional[str]:
        """Find column by checking canonical name, then aliases.

        Parameters
        ----------
        df_or_obs : pd.DataFrame
            DataFrame or adata.obs to search
        canonical : str
            Canonical column name to look for first
        aliases : List[str]
            List of alias column names to try
        explicit_override : str, optional
            Explicit column name that takes priority

        Returns
        -------
        str or None
            Actual column name found, or None if not found
        """
        columns = set(df_or_obs.columns)

        # Explicit override takes priority
        if explicit_override:
            if explicit_override in columns:
                return explicit_override
            self.logger.warning(
                "Explicit column override '%s' not found. Falling back to auto-detection.",
                explicit_override,
            )

        # Check canonical name first
        if canonical in columns:
            return canonical

        # Check aliases
        for alias in aliases:
            if alias in columns:
                return alias

        return None

    def _detect_path_separator(self, paths: List[str]) -> Optional[str]:
        """Detect the path separator used in a list of paths.

        Parameters
        ----------
        paths : List[str]
            Sample paths to analyze

        Returns
        -------
        str or None
            Detected separator, or None if no separator found
        """
        for sep in self.PATH_SEPARATORS:
            if any(sep in str(p) for p in paths if pd.notna(p)):
                return sep
        return None


# Backward compatibility alias
StageHAdapter = InputAdapter


def normalize_inputs(
    adata: "sc.AnnData",
    cluster_annotations_path: Optional[Path] = None,
    marker_scores_path: Optional[Path] = None,
    cluster_key: Optional[str] = None,
    label_key_in: Optional[str] = None,
    logger: Optional[logging.Logger] = None,
) -> Dict[str, Any]:
    """Convenience function to normalize all inputs at once.

    This is the main entry point for normalization.

    Parameters
    ----------
    adata : sc.AnnData
        AnnData object
    cluster_annotations_path : Path, optional
        Path to cluster_annotations.csv
    marker_scores_path : Path, optional
        Path to marker_scores.csv
    cluster_key : str, optional
        Explicit cluster column name override
    label_key_in : str, optional
        Explicit label column name override
    logger : logging.Logger, optional
        Logger for warnings

    Returns
    -------
    Dict with keys:
        - "adata_result": ValidationResult for AnnData
        - "cluster_annotations": Adapted DataFrame or None
        - "marker_scores": Adapted DataFrame or None
        - "adapter": The InputAdapter instance used
        - "provenance_type": "clustering", "refinement", or "unknown"

    Raises
    ------
    ValueError
        If adata validation fails
    """
    adapter = InputAdapter(
        cluster_key=cluster_key,
        label_key_in=label_key_in,
        logger=logger,
    )

    # Validate AnnData (required)
    adata_result = adapter.validate_adata(adata)
    if logger:
        adata_result.log_warnings(logger)
    adata_result.raise_if_invalid("input AnnData")

    # Adapt cluster_annotations if provided
    cluster_annotations = None
    if cluster_annotations_path and cluster_annotations_path.exists():
        raw_df = pd.read_csv(cluster_annotations_path)
        result = adapter.validate_cluster_annotations(raw_df)
        if logger:
            result.log_warnings(logger)
        result.raise_if_invalid("cluster_annotations.csv")
        cluster_annotations = adapter.adapt_cluster_annotations(raw_df, validate=False)

    # Adapt marker_scores if provided
    marker_scores = None
    if marker_scores_path and marker_scores_path.exists():
        raw_df = pd.read_csv(marker_scores_path)
        result = adapter.validate_marker_scores(raw_df)
        if logger:
            result.log_warnings(logger)
        result.raise_if_invalid("marker_scores.csv")
        marker_scores = adapter.adapt_marker_scores(raw_df, validate=False)

    return {
        "adata_result": adata_result,
        "cluster_annotations": cluster_annotations,
        "marker_scores": marker_scores,
        "adapter": adapter,
        "provenance_type": adapter.detect_provenance_type(adata),
    }
