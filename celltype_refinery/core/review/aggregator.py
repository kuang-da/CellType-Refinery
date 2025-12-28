"""Metrics aggregator for the Review module.

Loads and provides access to composition, spatial, and quality metrics
from upstream analysis stages.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import pandas as pd


class MetricsAggregator:
    """Aggregator for loading and accessing review metrics.

    Loads data from composition analysis, spatial analysis, and
    consolidation outputs to provide a unified interface for rules.

    Parameters
    ----------
    composition_dir : Path, optional
        Directory with composition analysis outputs
    spatial_dir : Path, optional
        Directory with spatial analysis outputs
    consolidation_dir : Path, optional
        Directory with consolidation outputs

    Attributes
    ----------
    composition_global : pd.DataFrame
        Global composition statistics
    composition_regional : pd.DataFrame
        Regional composition statistics
    composition_by_sample : pd.DataFrame
        Per-sample composition
    biology_by_region : pd.DataFrame
        Biology metrics by region
    morans_global : pd.DataFrame
        Moran's I statistics
    enrichment_pairs : pd.DataFrame
        Neighborhood enrichment pairs
    orphan_candidates : pd.DataFrame
        Orphan rescue candidates
    """

    def __init__(
        self,
        composition_dir: Optional[Path] = None,
        spatial_dir: Optional[Path] = None,
        consolidation_dir: Optional[Path] = None,
    ):
        self.composition_dir = Path(composition_dir) if composition_dir else None
        self.spatial_dir = Path(spatial_dir) if spatial_dir else None
        self.consolidation_dir = Path(consolidation_dir) if consolidation_dir else None

        # Loaded data
        self.composition_global: Optional[pd.DataFrame] = None
        self.composition_regional: Optional[pd.DataFrame] = None
        self.composition_by_sample: Optional[pd.DataFrame] = None
        self.biology_by_region: Optional[pd.DataFrame] = None
        self.morans_global: Optional[pd.DataFrame] = None
        self.enrichment_pairs: Optional[pd.DataFrame] = None
        self.orphan_candidates: Optional[pd.DataFrame] = None
        self.confidence_data: Optional[Dict[str, int]] = None

        # Cached lookups
        self._cell_types: Optional[List[str]] = None
        self._regions: Optional[List[str]] = None
        self._global_pct_cache: Dict[str, float] = {}
        self._regional_pct_cache: Dict[str, float] = {}
        self._morans_cache: Dict[str, float] = {}
        self._enrichment_cache: Dict[str, Dict[str, float]] = {}
        self._sample_dist_cache: Dict[str, Dict[str, Any]] = {}

        self._load_data()

    def _load_data(self) -> None:
        """Load all available data files."""
        # Load composition data
        if self.composition_dir and self.composition_dir.exists():
            self._load_composition()

        # Load spatial data
        if self.spatial_dir and self.spatial_dir.exists():
            self._load_spatial()

        # Load consolidation data
        if self.consolidation_dir and self.consolidation_dir.exists():
            self._load_consolidation()

    def _load_composition(self) -> None:
        """Load composition analysis outputs."""
        comp_dir = self.composition_dir

        # Global composition
        global_path = comp_dir / "composition_global.csv"
        if global_path.exists():
            self.composition_global = pd.read_csv(global_path)

        # Regional composition
        regional_path = comp_dir / "composition_by_region.csv"
        if regional_path.exists():
            self.composition_regional = pd.read_csv(regional_path)

        # Per-sample composition
        sample_path = comp_dir / "composition_by_sample.csv"
        if sample_path.exists():
            self.composition_by_sample = pd.read_csv(sample_path)

        # Biology metrics by region
        for fname in ["ft_biology_by_region.csv", "biology_by_region.csv"]:
            bio_path = comp_dir / fname
            if bio_path.exists():
                self.biology_by_region = pd.read_csv(bio_path)
                break

    def _load_spatial(self) -> None:
        """Load spatial analysis outputs."""
        spatial_dir = self.spatial_dir

        # Moran's I
        for subdir in ["heterogeneity", ""]:
            morans_path = spatial_dir / subdir / "morans_i_global.csv"
            if morans_path.exists():
                self.morans_global = pd.read_csv(morans_path)
                break

        # Enrichment pairs
        for subdir in ["neighborhood_enrichment", ""]:
            enrich_path = spatial_dir / subdir / "enrichment_pairs.csv"
            if enrich_path.exists():
                self.enrichment_pairs = pd.read_csv(enrich_path)
                break

    def _load_consolidation(self) -> None:
        """Load consolidation outputs."""
        cons_dir = self.consolidation_dir

        # Orphan candidates
        orphan_path = cons_dir / "orphan_candidates.csv"
        if orphan_path.exists():
            self.orphan_candidates = pd.read_csv(orphan_path)

    def load_from_dataframes(
        self,
        composition_global: Optional[pd.DataFrame] = None,
        composition_regional: Optional[pd.DataFrame] = None,
        composition_by_sample: Optional[pd.DataFrame] = None,
        biology_by_region: Optional[pd.DataFrame] = None,
        morans_global: Optional[pd.DataFrame] = None,
        enrichment_pairs: Optional[pd.DataFrame] = None,
        orphan_candidates: Optional[pd.DataFrame] = None,
    ) -> None:
        """Load data directly from DataFrames.

        Useful for testing or when data is already in memory.
        """
        if composition_global is not None:
            self.composition_global = composition_global
        if composition_regional is not None:
            self.composition_regional = composition_regional
        if composition_by_sample is not None:
            self.composition_by_sample = composition_by_sample
        if biology_by_region is not None:
            self.biology_by_region = biology_by_region
        if morans_global is not None:
            self.morans_global = morans_global
        if enrichment_pairs is not None:
            self.enrichment_pairs = enrichment_pairs
        if orphan_candidates is not None:
            self.orphan_candidates = orphan_candidates

        # Clear caches
        self._cell_types = None
        self._regions = None
        self._global_pct_cache.clear()
        self._regional_pct_cache.clear()
        self._morans_cache.clear()
        self._enrichment_cache.clear()
        self._sample_dist_cache.clear()

    def get_cell_types(self) -> List[str]:
        """Get list of all cell types in the data."""
        if self._cell_types is not None:
            return self._cell_types

        cell_types = set()

        if self.composition_global is not None:
            col = "cell_type" if "cell_type" in self.composition_global.columns else None
            if col:
                cell_types.update(self.composition_global[col].unique())

        self._cell_types = sorted(cell_types)
        return self._cell_types

    def get_regions(self) -> List[str]:
        """Get list of all regions in the data."""
        if self._regions is not None:
            return self._regions

        regions = set()

        if self.composition_regional is not None:
            if "region" in self.composition_regional.columns:
                regions.update(self.composition_regional["region"].unique())

        self._regions = sorted(regions)
        return self._regions

    def get_global_pct(self, cell_type: str) -> Optional[float]:
        """Get global percentage for a cell type."""
        if cell_type in self._global_pct_cache:
            return self._global_pct_cache[cell_type]

        if self.composition_global is None:
            return None

        ct_col = "cell_type" if "cell_type" in self.composition_global.columns else None
        pct_col = None
        for col in ["pct", "percentage", "proportion", "pct_of_total"]:
            if col in self.composition_global.columns:
                pct_col = col
                break

        if ct_col is None or pct_col is None:
            return None

        mask = self.composition_global[ct_col] == cell_type
        if mask.sum() == 0:
            return None

        pct = float(self.composition_global.loc[mask, pct_col].iloc[0])
        self._global_pct_cache[cell_type] = pct
        return pct

    def get_regional_pct(self, cell_type: str, region: str) -> Optional[float]:
        """Get regional percentage for a cell type."""
        cache_key = f"{cell_type}::{region}"
        if cache_key in self._regional_pct_cache:
            return self._regional_pct_cache[cache_key]

        if self.composition_regional is None:
            return None

        ct_col = "cell_type" if "cell_type" in self.composition_regional.columns else None
        pct_col = None
        for col in ["mean_pct", "pct", "percentage", "proportion"]:
            if col in self.composition_regional.columns:
                pct_col = col
                break

        if ct_col is None or pct_col is None:
            return None

        mask = (
            (self.composition_regional[ct_col] == cell_type)
            & (self.composition_regional["region"].str.lower() == region.lower())
        )
        if mask.sum() == 0:
            return None

        pct = float(self.composition_regional.loc[mask, pct_col].iloc[0])
        self._regional_pct_cache[cache_key] = pct
        return pct

    def get_sample_distribution(self, cell_type: str) -> Optional[Dict[str, Any]]:
        """Get sample distribution statistics for a cell type."""
        if cell_type in self._sample_dist_cache:
            return self._sample_dist_cache[cell_type]

        if self.composition_by_sample is None:
            return None

        ct_col = "cell_type" if "cell_type" in self.composition_by_sample.columns else None
        if ct_col is None:
            return None

        ct_data = self.composition_by_sample[
            self.composition_by_sample[ct_col] == cell_type
        ]

        if len(ct_data) == 0:
            return None

        # Calculate distribution stats
        count_col = "count" if "count" in ct_data.columns else "n_cells"
        if count_col not in ct_data.columns:
            return None

        total_count = ct_data[count_col].sum()
        n_samples = len(ct_data)
        max_count = ct_data[count_col].max()

        result = {
            "n_samples": n_samples,
            "total_count": total_count,
            "max_sample_count": max_count,
            "max_sample_pct": max_count / total_count if total_count > 0 else 0,
        }

        self._sample_dist_cache[cell_type] = result
        return result

    def get_moran_i(self, cell_type: str) -> Optional[float]:
        """Get Moran's I value for a cell type."""
        if cell_type in self._morans_cache:
            return self._morans_cache[cell_type]

        if self.morans_global is None:
            return None

        ct_col = "cell_type" if "cell_type" in self.morans_global.columns else None
        moran_col = None
        for col in ["moran_i", "morans_i", "I"]:
            if col in self.morans_global.columns:
                moran_col = col
                break

        if ct_col is None or moran_col is None:
            return None

        mask = self.morans_global[ct_col] == cell_type
        if mask.sum() == 0:
            return None

        moran_i = float(self.morans_global.loc[mask, moran_col].iloc[0])
        self._morans_cache[cell_type] = moran_i
        return moran_i

    def get_morans_i_data(self) -> Optional[pd.DataFrame]:
        """Get full Moran's I DataFrame."""
        return self.morans_global

    def get_enrichment_pairs(self) -> Optional[pd.DataFrame]:
        """Get enrichment pairs DataFrame."""
        return self.enrichment_pairs

    def get_enrichments_for_type(self, cell_type: str) -> Optional[Dict[str, float]]:
        """Get enrichment z-scores for a cell type with all other types."""
        if cell_type in self._enrichment_cache:
            return self._enrichment_cache[cell_type]

        if self.enrichment_pairs is None:
            return None

        # Find column names
        ct1_col = None
        ct2_col = None
        z_col = None

        for col in ["cell_type_1", "source", "from_type"]:
            if col in self.enrichment_pairs.columns:
                ct1_col = col
                break

        for col in ["cell_type_2", "target", "to_type"]:
            if col in self.enrichment_pairs.columns:
                ct2_col = col
                break

        for col in ["z_score", "zscore", "z"]:
            if col in self.enrichment_pairs.columns:
                z_col = col
                break

        if ct1_col is None or ct2_col is None or z_col is None:
            return None

        # Get enrichments where cell_type is either source or target
        mask = (self.enrichment_pairs[ct1_col] == cell_type) | (
            self.enrichment_pairs[ct2_col] == cell_type
        )
        subset = self.enrichment_pairs[mask]

        enrichments = {}
        for _, row in subset.iterrows():
            other_type = (
                row[ct2_col] if row[ct1_col] == cell_type else row[ct1_col]
            )
            enrichments[other_type] = row[z_col]

        self._enrichment_cache[cell_type] = enrichments
        return enrichments

    def get_regional_composition(self) -> Optional[pd.DataFrame]:
        """Get regional composition DataFrame."""
        return self.composition_regional

    def get_biology_by_region(self) -> Optional[pd.DataFrame]:
        """Get biology metrics by region DataFrame."""
        return self.biology_by_region

    def get_orphan_candidates(self) -> Optional[pd.DataFrame]:
        """Get orphan candidates DataFrame."""
        return self.orphan_candidates

    def get_iel_candidates(self) -> Optional[pd.DataFrame]:
        """Get IEL (intraepithelial lymphocyte) candidates DataFrame."""
        # IEL data is typically part of orphan candidates
        if self.orphan_candidates is None:
            return None

        if "subtype" in self.orphan_candidates.columns:
            iel_mask = self.orphan_candidates["subtype"].str.contains(
                "intraepithelial|IEL", case=False, na=False
            )
            return self.orphan_candidates[iel_mask]

        return None

    def get_confidence_breakdown(self) -> Optional[Dict[str, int]]:
        """Get confidence level breakdown."""
        return self.confidence_data

    def set_confidence_breakdown(self, data: Dict[str, int]) -> None:
        """Set confidence breakdown data."""
        self.confidence_data = data

    def summary(self) -> Dict[str, Any]:
        """Get summary of loaded data."""
        return {
            "composition_global": self.composition_global is not None,
            "composition_regional": self.composition_regional is not None,
            "composition_by_sample": self.composition_by_sample is not None,
            "biology_by_region": self.biology_by_region is not None,
            "morans_global": self.morans_global is not None,
            "enrichment_pairs": self.enrichment_pairs is not None,
            "orphan_candidates": self.orphan_candidates is not None,
            "n_cell_types": len(self.get_cell_types()),
            "n_regions": len(self.get_regions()),
        }
