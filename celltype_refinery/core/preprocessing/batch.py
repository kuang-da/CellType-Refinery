"""Batch effect correction (Stage E).

Provides linear model-based batch correction to remove
technical variation while preserving biological signal.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from .config import BatchCorrectionConfig


@dataclass
class BatchCorrectionResult:
    """Result from batch correction.

    Attributes
    ----------
    sample_id : str
        Sample identifier
    corrected_matrix : pd.DataFrame
        Batch-corrected intensity matrix
    shifts : Dict[str, float]
        Per-marker correction shifts applied
    """

    sample_id: str
    corrected_matrix: Optional[pd.DataFrame] = None
    shifts: Dict[str, float] = field(default_factory=dict)


@dataclass
class BatchDiagnostics:
    """Batch effect diagnostics.

    Attributes
    ----------
    marker_stats : pd.DataFrame
        Per-marker statistics (CV, eta2_batch, p_batch, stability)
    model_summary : str
        Linear model summary text
    """

    marker_stats: Optional[pd.DataFrame] = None
    model_summary: str = ""


class BatchCorrector:
    """Batch effect corrector using linear models.

    Uses weighted least squares to model technical batch effects
    (donor, imaging color, imaging cycle) while preserving biological
    variation (marker, region, marker x region interaction).

    Parameters
    ----------
    config : BatchCorrectionConfig
        Batch correction configuration

    Example
    -------
    >>> from celltype_refinery.core.preprocessing import BatchCorrector, BatchCorrectionConfig
    >>> config = BatchCorrectionConfig(batch_col="donor", biovar_col="region")
    >>> corrector = BatchCorrector(config)
    >>> shifts = corrector.compute_batch_shifts(summary_df, metadata)
    >>> result = corrector.correct_sample(matrix, "sample_01", shifts)
    """

    def __init__(self, config: Optional[BatchCorrectionConfig] = None):
        self.config = config or BatchCorrectionConfig()
        self._check_dependencies()

    def _check_dependencies(self) -> None:
        """Check for required dependencies."""
        try:
            import statsmodels.api
            import statsmodels.formula.api
        except ImportError:
            raise RuntimeError(
                "Batch correction requires statsmodels. "
                "Install with: pip install statsmodels"
            )

    def compute_sample_summary(
        self,
        sample_matrices: Dict[str, pd.DataFrame],
        sample_metadata: Dict[str, Dict[str, Any]],
        cell_id_col: str = "cell_mask_id",
    ) -> pd.DataFrame:
        """Compute per-sample, per-marker intensity summary.

        Parameters
        ----------
        sample_matrices : Dict[str, pd.DataFrame]
            Map of sample_id to aligned intensity matrix
        sample_metadata : Dict[str, Dict[str, Any]]
            Map of sample_id to metadata dict (donor, region, etc.)
        cell_id_col : str
            Cell ID column name

        Returns
        -------
        pd.DataFrame
            Summary with columns: sample_id, marker, intensity, n_cells,
            plus metadata columns
        """
        records = []

        for sample_id, matrix in sample_matrices.items():
            meta = sample_metadata.get(sample_id, {})
            marker_cols = [c for c in matrix.columns if c != cell_id_col]
            numeric = matrix[marker_cols].apply(pd.to_numeric, errors="coerce")

            means = numeric.mean(axis=0)
            counts = numeric.notna().sum(axis=0)

            for marker in marker_cols:
                record = {
                    "sample_id": sample_id,
                    "marker": marker,
                    "intensity": float(means.get(marker, np.nan)),
                    "n_cells": int(counts.get(marker, 0)),
                    self.config.batch_col: meta.get(self.config.batch_col),
                    self.config.biovar_col: meta.get(self.config.biovar_col),
                }
                # Add additional batch factors
                for factor in self.config.additional_batch_factors:
                    record[factor] = meta.get(factor)

                records.append(record)

        return pd.DataFrame(records)

    def compute_batch_shifts(
        self,
        summary_df: pd.DataFrame,
    ) -> Tuple[Dict[Tuple[str, str], float], BatchDiagnostics]:
        """Compute batch correction shifts using linear model.

        Parameters
        ----------
        summary_df : pd.DataFrame
            Sample summary from compute_sample_summary

        Returns
        -------
        Tuple[Dict, BatchDiagnostics]
            - Shift lookup: (sample_id, marker) -> shift value
            - Diagnostics with model summary
        """
        import statsmodels.api as sm
        import statsmodels.formula.api as smf

        diagnostics = BatchDiagnostics()

        if not self.config.enabled or self.config.correction_mode == "none":
            return {}, diagnostics

        # Prepare data
        model_df = summary_df.copy()
        model_df = model_df.dropna(
            subset=[self.config.batch_col, self.config.biovar_col, "marker"]
        )

        if model_df.empty:
            return {}, diagnostics

        weights = model_df["n_cells"].clip(lower=1)

        # Build formula with biological and technical terms
        batch_col = self.config.batch_col
        biovar_col = self.config.biovar_col

        formula_parts = [
            f"C(marker)",
            f"C({biovar_col})",
            f"C(marker):C({biovar_col})",
            f"C({batch_col})",
        ]

        # Add additional batch factors if present
        for factor in self.config.additional_batch_factors:
            if factor in model_df.columns:
                model_df[factor] = model_df[factor].fillna("missing").astype(str)
                formula_parts.append(f"C({factor})")

        formula = "intensity ~ " + " + ".join(formula_parts)

        # Fit model
        try:
            lm = smf.wls(formula, data=model_df, weights=weights).fit()
            diagnostics.model_summary = lm.summary().as_text()
        except Exception as e:
            diagnostics.model_summary = f"Model fitting failed: {e}"
            return {}, diagnostics

        # Extract batch effects
        def get_effect(prefix: str, level: Optional[str]) -> float:
            if level is None or str(level) == "nan":
                return 0.0
            key = f"C({prefix})[T.{level}]"
            return float(lm.params.get(key, 0.0))

        # Build shift lookup
        shift_lookup: Dict[Tuple[str, str], float] = {}

        for _, row in model_df.iterrows():
            sample_id = row["sample_id"]
            marker = row["marker"]

            # Sum all batch effects
            total_shift = get_effect(batch_col, row.get(batch_col))

            for factor in self.config.additional_batch_factors:
                if factor in row:
                    total_shift += get_effect(factor, row.get(factor))

            shift_lookup[(sample_id, marker)] = total_shift

        return shift_lookup, diagnostics

    def correct_sample(
        self,
        cell_matrix: pd.DataFrame,
        sample_id: str,
        shift_lookup: Dict[Tuple[str, str], float],
        cell_id_col: str = "cell_mask_id",
    ) -> BatchCorrectionResult:
        """Apply batch correction to a single sample.

        Parameters
        ----------
        cell_matrix : pd.DataFrame
            Aligned intensity matrix
        sample_id : str
            Sample identifier
        shift_lookup : Dict[Tuple[str, str], float]
            Shift values from compute_batch_shifts
        cell_id_col : str
            Cell ID column name

        Returns
        -------
        BatchCorrectionResult
            Corrected matrix and applied shifts
        """
        result = BatchCorrectionResult(sample_id=sample_id)

        if not self.config.enabled or self.config.correction_mode == "none":
            result.corrected_matrix = cell_matrix.copy()
            return result

        marker_cols = [c for c in cell_matrix.columns if c != cell_id_col]
        numeric = cell_matrix[marker_cols].apply(pd.to_numeric, errors="coerce")
        corrected = numeric.copy()

        for marker in marker_cols:
            shift = shift_lookup.get((sample_id, marker), 0.0)
            corrected[marker] = corrected[marker] - shift
            result.shifts[marker] = shift

        # Rebuild matrix with cell IDs
        result.corrected_matrix = pd.DataFrame({cell_id_col: cell_matrix[cell_id_col]})
        for marker in marker_cols:
            result.corrected_matrix[marker] = corrected[marker]

        return result

    def compute_marker_statistics(
        self,
        summary_df: pd.DataFrame,
        sample_counts: Dict[str, int],
    ) -> pd.DataFrame:
        """Compute per-marker batch effect statistics.

        Parameters
        ----------
        summary_df : pd.DataFrame
            Sample summary
        sample_counts : Dict[str, int]
            Number of cells per sample

        Returns
        -------
        pd.DataFrame
            Statistics with CV, eta2_batch, p_batch, stability
        """
        import statsmodels.api as sm
        import statsmodels.formula.api as smf

        batch_col = self.config.batch_col
        biovar_col = self.config.biovar_col

        # Pivot to get sample x marker matrix
        intensity_matrix = summary_df.pivot(
            index="sample_id", columns="marker", values="intensity"
        )

        # Get metadata indexed by sample (only sample-level columns, not per-marker factors)
        # Note: additional_batch_factors like imaging_color/cycle vary per marker, so we
        # only use batch_col and biovar_col for the per-marker ANOVA here.
        sample_meta_cols = [batch_col, biovar_col]
        metadata = summary_df[["sample_id"] + sample_meta_cols].drop_duplicates(subset=["sample_id"])
        metadata = metadata.set_index("sample_id")

        records = []

        for marker in intensity_matrix.columns:
            series = pd.to_numeric(intensity_matrix[marker], errors="coerce").dropna()

            # Compute CV
            variance = float(np.var(series, ddof=1)) if len(series) > 1 else float("nan")
            mean = float(np.mean(series)) if len(series) else float("nan")
            cv = (
                float(np.sqrt(max(variance, 0.0)) / abs(mean))
                if not (np.isnan(mean) or np.isclose(mean, 0.0))
                else float("nan")
            )

            n_cells = sum(sample_counts.get(s, 0) for s in series.index)

            record = {
                "marker": marker,
                "cv": cv,
                "variance": variance,
                "n_cells": n_cells,
                "eta2_batch": float("nan"),
                "p_batch": float("nan"),
                "stability": "stable",
            }

            # Run ANOVA if enough data
            if series.empty:
                records.append(record)
                continue

            joined = metadata.reindex(series.index).copy()
            joined = joined.dropna(subset=[batch_col])

            if joined.empty or joined[batch_col].nunique() < 2:
                records.append(record)
                continue

            joined[batch_col] = joined[batch_col].astype(str)
            anova_df = joined.assign(value=series.loc[joined.index].to_numpy())

            try:
                formula = f"value ~ C({batch_col})"
                if biovar_col in joined.columns and joined[biovar_col].nunique() >= 2:
                    joined[biovar_col] = joined[biovar_col].astype(str)
                    anova_df = joined.assign(value=series.loc[joined.index].to_numpy())
                    formula = f"value ~ C({batch_col}) + C({biovar_col})"

                model = smf.ols(formula, data=anova_df).fit()
                anova_table = sm.stats.anova_lm(model, typ=2)

                total_ss = float(anova_table["sum_sq"].sum())
                if total_ss > 0:
                    batch_term = f"C({batch_col})"
                    if batch_term in anova_table.index:
                        batch_row = anova_table.loc[batch_term]
                        record["eta2_batch"] = float(batch_row["sum_sq"]) / total_ss
                        record["p_batch"] = float(batch_row.get("PR(>F)", float("nan")))

            except Exception:
                pass

            # Classify stability
            p = record["p_batch"]
            eta = record["eta2_batch"]
            if pd.notna(p) and pd.notna(eta):
                if p < 0.01 and eta >= 0.5:
                    record["stability"] = "drifted"
                elif (0.01 <= p < 0.05) or (0.3 <= eta < 0.5):
                    record["stability"] = "observe"

            records.append(record)

        return pd.DataFrame(records).sort_values("cv", ascending=False)
