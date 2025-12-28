"""Statistical utilities for CellType-Refinery.

Provides functions for percentile computation, robust statistics,
and distribution comparison metrics.
"""

from __future__ import annotations

from typing import Iterable, Sequence, Union

import numpy as np
from scipy.stats import ks_2samp, wasserstein_distance

ArrayLike = Union[Iterable[float], np.ndarray]


def _to_clean_array(values: ArrayLike) -> np.ndarray:
    """Convert input to clean numpy array, removing non-finite values.

    Parameters
    ----------
    values : ArrayLike
        Input values (list, iterable, or array).

    Returns
    -------
    np.ndarray
        Clean array with only finite values.
    """
    arr = np.asarray(list(values), dtype=float)
    if arr.size == 0:
        return arr
    return arr[np.isfinite(arr)]


def compute_percentiles(values: ArrayLike, percentiles: Sequence[float]) -> np.ndarray:
    """Compute percentile values ignoring NaNs.

    Parameters
    ----------
    values : ArrayLike
        Input values.
    percentiles : Sequence[float]
        Percentiles to compute (0-100).

    Returns
    -------
    np.ndarray
        Computed percentile values. Returns NaN array if input is empty.
    """
    arr = _to_clean_array(values)
    if arr.size == 0:
        return np.full(len(percentiles), np.nan)
    return np.percentile(arr, percentiles)


def robust_zscore(
    values: ArrayLike,
    *,
    median: float | None = None,
    mad: float | None = None,
) -> np.ndarray:
    """Compute a robust z-score using the median absolute deviation (MAD).

    The MAD-based z-score is more robust to outliers than standard z-scores.
    Uses the standard conversion factor of 1.4826 to make MAD comparable
    to standard deviation for normally distributed data.

    Parameters
    ----------
    values : ArrayLike
        Input values.
    median : float, optional
        Pre-computed median. If None, computed from data.
    mad : float, optional
        Pre-computed MAD. If None, computed from data.

    Returns
    -------
    np.ndarray
        Robust z-scores. Non-finite inputs become NaN in output.
    """
    arr = np.asarray(list(values), dtype=float)
    if arr.size == 0:
        return arr

    mask = np.isfinite(arr)
    clean = arr[mask]
    if clean.size == 0:
        return np.full_like(arr, np.nan, dtype=float)

    if median is None:
        median = np.median(clean)
    if mad is None:
        mad = np.median(np.abs(clean - median))

    # Scale MAD to be comparable to standard deviation
    scale = mad * 1.4826 if mad else np.nan

    if not np.isfinite(scale) or scale == 0:
        z = np.zeros_like(clean)
    else:
        z = (clean - median) / scale

    result = np.full_like(arr, np.nan, dtype=float)
    result[mask] = z
    return result


def ks_distance(a: ArrayLike, b: ArrayLike) -> float:
    """Compute the two-sample Kolmogorov-Smirnov statistic.

    The KS statistic measures the maximum difference between two
    empirical distribution functions.

    Parameters
    ----------
    a : ArrayLike
        First sample.
    b : ArrayLike
        Second sample.

    Returns
    -------
    float
        KS statistic (0 to 1). Returns NaN if either sample has <2 values.
    """
    arr_a = _to_clean_array(a)
    arr_b = _to_clean_array(b)
    if arr_a.size < 2 or arr_b.size < 2:
        return float("nan")
    return float(ks_2samp(arr_a, arr_b, alternative="two-sided").statistic)


def earth_movers_distance(a: ArrayLike, b: ArrayLike) -> float:
    """Compute the Earth Mover's Distance (first Wasserstein distance).

    The EMD measures the minimum "work" required to transform one
    distribution into another, where work is defined as the amount
    of distribution weight times the distance it must be moved.

    Parameters
    ----------
    a : ArrayLike
        First sample.
    b : ArrayLike
        Second sample.

    Returns
    -------
    float
        Earth Mover's Distance. Returns NaN if either sample is empty.
    """
    arr_a = _to_clean_array(a)
    arr_b = _to_clean_array(b)
    if arr_a.size == 0 or arr_b.size == 0:
        return float("nan")
    return float(wasserstein_distance(arr_a, arr_b))
