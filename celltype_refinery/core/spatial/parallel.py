"""Parallel execution strategies for spatial analysis.

Provides efficient parallelization for:
- Permutation tests (parallel permutations within samples)
- Sample processing (parallel samples with adaptive thread allocation)

Optimized for high-core systems (128+ threads, 256+ GB RAM).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Dict, List, Optional, Tuple
import logging

import numpy as np
from scipy import sparse
from joblib import Parallel, delayed

from .interaction import compute_neighbor_counts_vectorized

logger = logging.getLogger(__name__)


@dataclass
class SampleInput:
    """Input data for processing a single sample."""

    sample_id: str
    graph: sparse.csr_matrix
    cell_types: np.ndarray
    region: Optional[str] = None

    @property
    def n_cells(self) -> int:
        return len(self.cell_types)


def compute_single_permutation(
    graph: sparse.csr_matrix,
    cell_types: np.ndarray,
    unique_types: List[str],
    seed: int,
) -> np.ndarray:
    """Compute neighbor counts for a single permutation.

    This function is designed to be called in parallel.

    Parameters
    ----------
    graph : sparse.csr_matrix
        Sparse adjacency matrix
    cell_types : np.ndarray
        Original cell type labels (will be shuffled)
    unique_types : List[str]
        List of unique cell types
    seed : int
        Random seed for this permutation

    Returns
    -------
    np.ndarray
        Neighbor count matrix for this permutation
    """
    rng = np.random.default_rng(seed)
    shuffled = rng.permutation(cell_types)
    return compute_neighbor_counts_vectorized(graph, shuffled, unique_types)


def batch_permutations_parallel(
    graph: sparse.csr_matrix,
    cell_types: np.ndarray,
    unique_types: List[str],
    n_permutations: int,
    n_jobs: int = 64,
    batch_size: int = 16,
    base_seed: int = 42,
) -> np.ndarray:
    """Run permutations in parallel batches.

    Parallelizes across permutations within a single sample.
    Optimal for large samples (>50K cells) where each permutation
    takes significant time.

    Parameters
    ----------
    graph : sparse.csr_matrix
        Sparse adjacency matrix
    cell_types : np.ndarray
        Cell type labels
    unique_types : List[str]
        List of unique cell types
    n_permutations : int
        Number of permutations to run
    n_jobs : int
        Number of parallel jobs (default: 64)
    batch_size : int
        Batch size for joblib (default: 16)
    base_seed : int
        Base random seed (each permutation uses base_seed + perm_index)

    Returns
    -------
    np.ndarray
        Null distribution of shape (n_permutations, n_types, n_types)
    """
    # Generate unique seeds for each permutation
    seeds = list(range(base_seed, base_seed + n_permutations))

    # Run permutations in parallel
    results = Parallel(n_jobs=n_jobs, batch_size=batch_size, verbose=0)(
        delayed(compute_single_permutation)(graph, cell_types, unique_types, seed)
        for seed in seeds
    )

    return np.array(results)


def run_permutations_sequential(
    graph: sparse.csr_matrix,
    cell_types: np.ndarray,
    unique_types: List[str],
    n_permutations: int,
    base_seed: int = 42,
    progress_callback: Optional[Callable[[int, int], None]] = None,
) -> np.ndarray:
    """Run permutations sequentially (for small samples or debugging).

    Parameters
    ----------
    graph : sparse.csr_matrix
        Sparse adjacency matrix
    cell_types : np.ndarray
        Cell type labels
    unique_types : List[str]
        List of unique cell types
    n_permutations : int
        Number of permutations
    base_seed : int
        Base random seed
    progress_callback : Callable, optional
        Callback(current_perm, total_perms) for progress tracking

    Returns
    -------
    np.ndarray
        Null distribution
    """
    n_types = len(unique_types)
    null_counts = np.zeros((n_permutations, n_types, n_types))

    for perm in range(n_permutations):
        seed = base_seed + perm
        null_counts[perm] = compute_single_permutation(
            graph, cell_types, unique_types, seed
        )

        if progress_callback:
            progress_callback(perm + 1, n_permutations)

    return null_counts


def process_sample_parallel_perms(
    sample: SampleInput,
    valid_types: List[str],
    n_permutations: int,
    n_perm_threads: int,
    base_seed: int,
    min_cells_per_type: int,
    process_func: Callable,
) -> Optional[Dict[str, Any]]:
    """Process a sample with parallel permutations.

    Parameters
    ----------
    sample : SampleInput
        Sample data
    valid_types : List[str]
        Cell types to analyze
    n_permutations : int
        Number of permutations
    n_perm_threads : int
        Threads for parallel permutations
    base_seed : int
        Base random seed
    min_cells_per_type : int
        Minimum cells per type
    process_func : Callable
        Function to compute all metrics for this sample

    Returns
    -------
    Dict or None
        Sample results, or None if sample should be skipped
    """
    if sample.n_cells < 100:  # Skip tiny samples
        logger.debug(f"Skipping {sample.sample_id}: too few cells ({sample.n_cells})")
        return None

    return process_func(
        sample=sample,
        valid_types=valid_types,
        n_permutations=n_permutations,
        n_perm_threads=n_perm_threads,
        base_seed=base_seed,
        min_cells_per_type=min_cells_per_type,
    )


def process_samples_adaptive(
    samples: List[SampleInput],
    valid_types: List[str],
    n_permutations: int,
    n_threads_total: int,
    base_seed: int,
    min_cells_per_type: int,
    process_func: Callable,
    large_sample_threshold: int = 50000,
) -> List[Dict[str, Any]]:
    """Adaptive parallel processing based on sample size.

    Strategy:
    - Large samples (>threshold cells): Process one at a time with
      many threads for parallel permutations
    - Small samples (<threshold cells): Process multiple samples in
      parallel with fewer threads each

    Parameters
    ----------
    samples : List[SampleInput]
        Sample data
    valid_types : List[str]
        Cell types to analyze
    n_permutations : int
        Number of permutations
    n_threads_total : int
        Total available threads
    base_seed : int
        Base random seed
    min_cells_per_type : int
        Minimum cells per type
    process_func : Callable
        Function to process each sample
    large_sample_threshold : int
        Threshold for large vs small samples

    Returns
    -------
    List[Dict]
        Results from all samples (filtered for None)
    """
    large_samples = [s for s in samples if s.n_cells > large_sample_threshold]
    small_samples = [s for s in samples if s.n_cells <= large_sample_threshold]

    results = []

    # Large samples: process sequentially, each with parallel permutations
    logger.info(
        f"Processing {len(large_samples)} large samples "
        f"(>{large_sample_threshold} cells) with {n_threads_total} threads each"
    )
    for sample in large_samples:
        result = process_func(
            sample=sample,
            valid_types=valid_types,
            n_permutations=n_permutations,
            n_perm_threads=n_threads_total,
            base_seed=base_seed,
            min_cells_per_type=min_cells_per_type,
        )
        if result is not None:
            results.append(result)

    # Small samples: process in parallel, each with fewer threads
    if small_samples:
        n_parallel_samples = min(8, len(small_samples))
        threads_per_sample = max(1, n_threads_total // n_parallel_samples)

        logger.info(
            f"Processing {len(small_samples)} small samples "
            f"({n_parallel_samples} parallel, {threads_per_sample} threads each)"
        )

        small_results = Parallel(n_jobs=n_parallel_samples, verbose=10)(
            delayed(process_func)(
                sample=sample,
                valid_types=valid_types,
                n_permutations=n_permutations,
                n_perm_threads=threads_per_sample,
                base_seed=base_seed + i * n_permutations,  # Different seeds per sample
                min_cells_per_type=min_cells_per_type,
            )
            for i, sample in enumerate(small_samples)
        )

        results.extend([r for r in small_results if r is not None])

    return results


def get_optimal_thread_allocation(
    n_samples: int,
    sample_sizes: List[int],
    n_threads_total: int = 128,
    large_threshold: int = 50000,
) -> Tuple[int, int]:
    """Get optimal thread allocation for the workload.

    Parameters
    ----------
    n_samples : int
        Number of samples
    sample_sizes : List[int]
        Cell counts per sample
    n_threads_total : int
        Total available threads
    large_threshold : int
        Threshold for large samples

    Returns
    -------
    Tuple[int, int]
        (n_parallel_samples, threads_per_perm_batch)
    """
    n_large = sum(1 for s in sample_sizes if s > large_threshold)
    n_small = n_samples - n_large

    if n_large > n_small:
        # Mostly large samples: maximize threads per sample
        return 1, n_threads_total
    else:
        # Mostly small samples: parallelize across samples
        n_parallel = min(8, n_small)
        threads_per = max(1, n_threads_total // n_parallel)
        return n_parallel, threads_per
