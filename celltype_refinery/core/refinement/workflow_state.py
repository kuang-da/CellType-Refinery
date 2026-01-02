"""Workflow state management for Stage I refinement iterations.

This module implements the WorkflowState pattern to enable diagnostic report
reuse across Stage I iterations. When Stage I runs after Stage J (or a previous
Stage I iteration), it detects and reuses the previous stage's diagnostic_report.csv
instead of regenerating it.

Key features:
- Detects previous stage diagnostic artifacts (diagnostic_report.csv, groups_derived.yaml)
- Reuses diagnostic reports to preserve Pool cluster recommendations (SUBCLUSTER)
- Resets diagnose_complete after execute so next iteration regenerates

Reference: ft/src/workflow/state.py (WorkflowState.load_or_create, mark_execute_complete)
"""

from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Optional
import logging
import shutil

import yaml

STATE_FILENAME = "workflow_state.yaml"

logger = logging.getLogger(__name__)


@dataclass
class RefinementWorkflowState:
    """Workflow state for tracking diagnostic and execution status.

    This class tracks the state of Stage I refinement workflows to enable:
    1. Diagnostic report reuse (skip regeneration if previous stage has valid diagnostic)
    2. State persistence across workflow steps (execute, review, diagnose)
    3. Proper reset after execution to enable iterative refinement

    Attributes
    ----------
    version : str
        State schema version
    out_dir : str
        Output directory path
    created_at : str
        ISO timestamp when state was created
    updated_at : str
        ISO timestamp when state was last updated
    diagnose_complete : bool
        True if diagnostic report is available and can be reused
    diagnose_timestamp : str, optional
        ISO timestamp when diagnostic was completed
    diagnostic_report_path : str, optional
        Path to diagnostic_report.csv (may be from parent stage)
    diagnostic_summary : dict
        Summary counts {"SUBCLUSTER": n, "RELABEL": n, "SKIP": n}
    groups_derived_path : str, optional
        Path to groups_derived.yaml in output directory
    execute_complete : bool
        True if execution step completed
    execute_timestamp : str, optional
        ISO timestamp when execution completed
    parent_dir : str, optional
        Path to parent stage directory (Stage J or previous Stage I)
    """

    # Core fields
    version: str = "1.0"
    out_dir: str = ""
    created_at: str = ""
    updated_at: str = ""

    # Diagnostic state
    diagnose_complete: bool = False
    diagnose_timestamp: Optional[str] = None
    diagnostic_report_path: Optional[str] = None
    diagnostic_summary: Dict[str, int] = field(default_factory=dict)
    groups_derived_path: Optional[str] = None

    # Execution state
    execute_complete: bool = False
    execute_timestamp: Optional[str] = None

    # Parent tracking
    parent_dir: Optional[str] = None

    @classmethod
    def load_or_create(
        cls,
        out_dir: Path,
        stage_h_dir: Optional[Path] = None,
    ) -> "RefinementWorkflowState":
        """Load existing state or create new with parent detection.

        This method implements the key WorkflowState pattern:
        1. If workflow_state.yaml exists in out_dir -> load it
        2. If creating new AND stage_h_dir has diagnostic artifacts:
           - Set diagnose_complete=True
           - Set diagnostic_report_path to parent's diagnostic
           - Copy groups_derived.yaml to out_dir

        Parameters
        ----------
        out_dir : Path
            Output directory for this Stage I run
        stage_h_dir : Path, optional
            Parent stage directory (Stage J or previous Stage I) to check
            for diagnostic artifacts

        Returns
        -------
        RefinementWorkflowState
            Loaded or newly created state
        """
        out_dir = Path(out_dir)
        state_path = out_dir / STATE_FILENAME

        if state_path.exists():
            logger.debug("Loading existing workflow state from: %s", state_path)
            return cls.load(state_path)

        # Create new state
        state = cls(
            out_dir=str(out_dir),
            created_at=datetime.now().isoformat(),
        )

        # Detect parent artifacts (Stage J or previous Stage I)
        if stage_h_dir and not state.execute_complete:
            h_dir = Path(stage_h_dir)
            h_diag = h_dir / "diagnostic_report.csv"
            h_groups = h_dir / "groups_derived.yaml"

            if h_diag.exists() and h_groups.exists():
                logger.info("Detected parent diagnostic artifacts in: %s", h_dir)
                logger.info("  diagnostic_report.csv: %s", h_diag)
                logger.info("  groups_derived.yaml: %s", h_groups)

                # Parent has diagnostic artifacts - reuse them
                state.diagnose_complete = True
                state.diagnose_timestamp = datetime.now().isoformat()
                state.diagnostic_report_path = str(h_diag)
                state.parent_dir = str(h_dir)

                # Load diagnostic summary from parent's workflow_state if available
                parent_state_path = h_dir / STATE_FILENAME
                if parent_state_path.exists():
                    try:
                        with open(parent_state_path) as f:
                            parent_state = yaml.safe_load(f) or {}
                        if "diagnostic_summary" in parent_state:
                            state.diagnostic_summary = parent_state["diagnostic_summary"]
                            logger.info(
                                "  Loaded diagnostic_summary: %s",
                                state.diagnostic_summary,
                            )
                    except Exception as e:
                        logger.warning("Failed to load parent state: %s", e)

                # Copy groups_derived.yaml to output directory
                out_dir.mkdir(parents=True, exist_ok=True)
                out_groups = out_dir / "groups_derived.yaml"
                if not out_groups.exists():
                    shutil.copy(h_groups, out_groups)
                    logger.info("  Copied groups_derived.yaml to: %s", out_groups)
                state.groups_derived_path = str(out_groups)

        return state

    @classmethod
    def load(cls, path: Path) -> "RefinementWorkflowState":
        """Load state from YAML file.

        Parameters
        ----------
        path : Path
            Path to workflow_state.yaml or directory containing it

        Returns
        -------
        RefinementWorkflowState
            Loaded state
        """
        path = Path(path)
        if path.is_dir():
            path = path / STATE_FILENAME
        with open(path) as f:
            data = yaml.safe_load(f) or {}
        # Filter to only known fields to avoid errors from extra fields
        known_fields = cls.__dataclass_fields__.keys()
        filtered_data = {k: v for k, v in data.items() if k in known_fields}
        return cls(**filtered_data)

    def save(self, path: Optional[Path] = None) -> Path:
        """Save state to YAML file.

        Parameters
        ----------
        path : Path, optional
            Path to save to. If None, saves to out_dir/workflow_state.yaml

        Returns
        -------
        Path
            Path where state was saved
        """
        if path is None:
            path = Path(self.out_dir) / STATE_FILENAME
        else:
            path = Path(path)

        path.parent.mkdir(parents=True, exist_ok=True)
        self.updated_at = datetime.now().isoformat()

        data = {
            "version": self.version,
            "out_dir": self.out_dir,
            "created_at": self.created_at,
            "updated_at": self.updated_at,
            "diagnose_complete": self.diagnose_complete,
            "diagnose_timestamp": self.diagnose_timestamp,
            "diagnostic_report_path": self.diagnostic_report_path,
            "diagnostic_summary": self.diagnostic_summary,
            "groups_derived_path": self.groups_derived_path,
            "execute_complete": self.execute_complete,
            "execute_timestamp": self.execute_timestamp,
            "parent_dir": self.parent_dir,
        }
        with open(path, "w") as f:
            yaml.safe_dump(data, f, default_flow_style=False, sort_keys=False)

        logger.debug("Saved workflow state to: %s", path)
        return path

    def should_skip_diagnostic(self) -> bool:
        """Check if diagnostic regeneration can be skipped.

        Returns True if:
        1. diagnose_complete is True
        2. diagnostic_report_path is set
        3. The file at diagnostic_report_path exists

        Returns
        -------
        bool
            True if diagnostic can be reused, False if regeneration needed
        """
        if not self.diagnose_complete:
            return False
        if not self.diagnostic_report_path:
            return False
        return Path(self.diagnostic_report_path).exists()

    def get_diagnostic_report_path(self) -> Optional[Path]:
        """Get path to diagnostic report if available.

        Returns
        -------
        Path or None
            Path to diagnostic_report.csv if exists, None otherwise
        """
        if self.diagnostic_report_path:
            p = Path(self.diagnostic_report_path)
            return p if p.exists() else None
        return None

    def mark_execute_complete(self) -> None:
        """Mark execution complete and RESET diagnostic for next iteration.

        This is the critical reset logic that enables iterative refinement:
        - After execute creates new subclusters, the diagnostic must be regenerated
        - Setting diagnose_complete=False ensures the next Stage I iteration
          will regenerate the diagnostic (via the diagnose step)

        This matches the reference behavior in:
        ft/src/workflow/state.py:WorkflowState.mark_execute_complete()
        """
        self.execute_complete = True
        self.execute_timestamp = datetime.now().isoformat()

        # CRITICAL: Reset diagnostic for next iteration
        # New subclusters need fresh diagnostic recommendations
        self.diagnose_complete = False
        self.diagnostic_report_path = None
        self.diagnostic_summary = {}
        self.groups_derived_path = None

        logger.info(
            "Marked execute_complete=True, reset diagnose_complete=False for next iteration"
        )

    def get_status_summary(self) -> str:
        """Get human-readable status summary.

        Returns
        -------
        str
            Formatted status summary
        """
        lines = [
            "Workflow State:",
            f"  Version: {self.version}",
            f"  Output: {self.out_dir}",
            f"  Created: {self.created_at}",
            f"  Updated: {self.updated_at}",
            "",
            "Steps:",
        ]

        # Diagnose status
        if self.diagnose_complete:
            summary_str = ", ".join(
                f"{k}={v}" for k, v in self.diagnostic_summary.items()
            )
            lines.append(f"  [x] diagnose - complete ({summary_str})")
            if self.parent_dir:
                lines.append(f"      (reusing from: {self.parent_dir})")
        else:
            lines.append("  [ ] diagnose - pending")

        # Execute status
        if self.execute_complete:
            lines.append(f"  [x] execute - complete at {self.execute_timestamp}")
        else:
            lines.append("  [ ] execute - pending")

        return "\n".join(lines)
