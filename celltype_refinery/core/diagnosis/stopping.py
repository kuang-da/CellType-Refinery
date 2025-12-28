"""Stopping criteria for iterative refinement.

This module provides classes for determining when iterative refinement should stop.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Optional

import pandas as pd


@dataclass
class StoppingStatus:
    """Result of stopping criteria evaluation.

    Attributes
    ----------
    should_stop : bool
        Whether iteration should stop.
    reason : str
        Reason for stopping (or empty if continuing).
        Values: "no_more_subcluster", "max_iterations", "below_threshold", ""
    iteration : int
        Current iteration number.
    n_subcluster : int
        Number of clusters recommended for SUBCLUSTER.
    n_relabel : int
        Number of clusters recommended for RELABEL.
    n_skip : int
        Number of clusters recommended for SKIP.
    improvement_pct : Optional[float]
        Percentage improvement from previous iteration (if applicable).
    """

    should_stop: bool
    reason: str
    iteration: int
    n_subcluster: int
    n_relabel: int
    n_skip: int
    improvement_pct: Optional[float] = None

    def __str__(self) -> str:
        """Human-readable status summary."""
        status = "STOP" if self.should_stop else "CONTINUE"
        parts = [
            f"Iteration {self.iteration}: {status}",
            f"SUBCLUSTER={self.n_subcluster}, RELABEL={self.n_relabel}, SKIP={self.n_skip}",
        ]
        if self.reason:
            parts.append(f"Reason: {self.reason}")
        if self.improvement_pct is not None:
            parts.append(f"Improvement: {self.improvement_pct:.1f}%")
        return " | ".join(parts)


class StoppingCriteria:
    """Evaluates whether iterative refinement should stop.

    Stopping conditions (checked in order):
    1. No more SUBCLUSTER recommendations (primary)
    2. Max iterations reached (safety)
    3. Improvement below threshold (optional)

    Example:
        >>> criteria = StoppingCriteria(max_iterations=10)
        >>> status = criteria.evaluate(report, iteration=1)
        >>> if status.should_stop:
        ...     print(f"Stopping: {status.reason}")
    """

    def __init__(
        self,
        max_iterations: int = 10,
        min_improvement_pct: Optional[float] = None,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize stopping criteria.

        Parameters
        ----------
        max_iterations : int
            Maximum number of iterations before forced stop.
        min_improvement_pct : float, optional
            Minimum percentage improvement to continue (0-100).
            If provided, stops when improvement drops below this.
        logger : logging.Logger, optional
            Logger for debug messages.
        """
        self.max_iterations = max_iterations
        self.min_improvement_pct = min_improvement_pct
        self.logger = logger or logging.getLogger(__name__)
        self._prev_n_subcluster: Optional[int] = None

    def evaluate(
        self,
        report: pd.DataFrame,
        iteration: int,
    ) -> StoppingStatus:
        """Evaluate stopping criteria for current iteration.

        Parameters
        ----------
        report : pd.DataFrame
            Diagnostic report with 'recommendation' column.
        iteration : int
            Current iteration number (1-indexed).

        Returns
        -------
        StoppingStatus
            Status indicating whether to stop and why.
        """
        # Count recommendations
        n_subcluster = (report["recommendation"] == "SUBCLUSTER").sum()
        n_relabel = (report["recommendation"] == "RELABEL").sum()
        n_skip = (report["recommendation"] == "SKIP").sum()

        # Calculate improvement if we have previous data
        improvement_pct = None
        if self._prev_n_subcluster is not None and self._prev_n_subcluster > 0:
            reduction = self._prev_n_subcluster - n_subcluster
            improvement_pct = (reduction / self._prev_n_subcluster) * 100

        # Check stopping conditions
        should_stop = False
        reason = ""

        # Condition 1: No more SUBCLUSTER recommendations
        if n_subcluster == 0:
            should_stop = True
            reason = "no_more_subcluster"
            self.logger.info(
                "Stopping: No more clusters need subclustering (iteration %d)",
                iteration,
            )

        # Condition 2: Max iterations reached
        elif iteration >= self.max_iterations:
            should_stop = True
            reason = "max_iterations"
            self.logger.info(
                "Stopping: Max iterations reached (%d)", self.max_iterations
            )

        # Condition 3: Improvement below threshold (optional)
        elif (
            self.min_improvement_pct is not None
            and improvement_pct is not None
            and improvement_pct < self.min_improvement_pct
        ):
            should_stop = True
            reason = "below_threshold"
            self.logger.info(
                "Stopping: Improvement %.1f%% below threshold %.1f%% (iteration %d)",
                improvement_pct,
                self.min_improvement_pct,
                iteration,
            )

        # Update state for next iteration
        self._prev_n_subcluster = n_subcluster

        status = StoppingStatus(
            should_stop=should_stop,
            reason=reason,
            iteration=iteration,
            n_subcluster=n_subcluster,
            n_relabel=n_relabel,
            n_skip=n_skip,
            improvement_pct=improvement_pct,
        )

        self.logger.debug("%s", status)
        return status

    def reset(self) -> None:
        """Reset internal state for a new iteration sequence."""
        self._prev_n_subcluster = None


def check_convergence(
    current_report: pd.DataFrame,
    previous_report: Optional[pd.DataFrame] = None,
) -> bool:
    """Check if refinement has converged (no changes between iterations).

    Parameters
    ----------
    current_report : pd.DataFrame
        Current iteration's diagnostic report.
    previous_report : pd.DataFrame, optional
        Previous iteration's diagnostic report.

    Returns
    -------
    bool
        True if converged (no meaningful changes), False otherwise.
    """
    if previous_report is None:
        return False

    # Compare SUBCLUSTER recommendations
    current_sub = set(
        current_report[current_report["recommendation"] == "SUBCLUSTER"][
            "cluster_id"
        ].tolist()
    )
    prev_sub = set(
        previous_report[previous_report["recommendation"] == "SUBCLUSTER"][
            "cluster_id"
        ].tolist()
    )

    # Converged if same clusters need subclustering
    return current_sub == prev_sub
