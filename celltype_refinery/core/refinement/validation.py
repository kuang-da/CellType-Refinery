"""
Validation errors with actionable diagnostics for refinement.

These error classes provide structured information about what went wrong,
what was expected, and how to fix it. Error codes enable programmatic handling.

Error Codes:
    E001_MISSING_COLUMN: Required column not found in DataFrame/adata.obs
    E002_INVALID_CLUSTER_ID: Cluster ID not found in AnnData
    E003_MISSING_PROVENANCE: Expected provenance keys not in adata.uns
    E004_PATH_SEPARATOR: Unrecognized hierarchical path separator
    E005_TYPE_MISMATCH: Unexpected data type for column
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from difflib import get_close_matches
from typing import Any, Dict, List, Optional


@dataclass
class ValidationError:
    """Base class for validation errors with actionable diagnostics.

    Attributes
    ----------
    message : str
        Human-readable error description
    error_code : str
        Machine-readable error code for programmatic handling
    expected : Any
        What the validator expected to find
    found : Any
        What was actually found
    suggestion : str
        Actionable suggestion for fixing the error
    context : Dict[str, Any]
        Additional context for debugging
    """

    message: str
    error_code: str = "E000_UNKNOWN"
    expected: Any = None
    found: Any = None
    suggestion: str = ""
    context: Dict[str, Any] = field(default_factory=dict)

    def __str__(self) -> str:
        """Format error as human-readable multi-line string."""
        parts = [f"[{self.error_code}] {self.message}"]
        if self.expected is not None:
            parts.append(f"  Expected: {self.expected}")
        if self.found is not None:
            parts.append(f"  Found: {self.found}")
        if self.suggestion:
            parts.append(f"  Suggestion: {self.suggestion}")
        return "\n".join(parts)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            "error_code": self.error_code,
            "message": self.message,
            "expected": str(self.expected) if self.expected is not None else None,
            "found": str(self.found) if self.found is not None else None,
            "suggestion": self.suggestion,
            "context": self.context,
        }


@dataclass
class MissingColumnError(ValidationError):
    """Error for missing required columns.

    Automatically suggests close matches from available columns using difflib.
    """

    error_code: str = "E001_MISSING_COLUMN"
    available_columns: List[str] = field(default_factory=list)

    def __post_init__(self):
        """Generate suggestion from available columns if not provided."""
        if not self.suggestion and self.available_columns and self.expected:
            # Find closest matching column names
            matches = get_close_matches(
                str(self.expected), self.available_columns, n=3, cutoff=0.4
            )
            if matches:
                self.suggestion = f"Did you mean: {', '.join(matches)}?"


@dataclass
class InvalidClusterIdError(ValidationError):
    """Error for cluster IDs not found in AnnData.

    Provides suggestions for valid cluster IDs when available.
    """

    error_code: str = "E002_INVALID_CLUSTER_ID"
    cluster_id: str = ""
    valid_ids: List[str] = field(default_factory=list)

    def __post_init__(self):
        """Generate suggestion from valid cluster IDs if not provided."""
        if not self.suggestion and self.valid_ids:
            # Try to find close matches
            matches = get_close_matches(
                self.cluster_id, self.valid_ids, n=5, cutoff=0.4
            )
            if matches:
                self.suggestion = f"Valid cluster IDs include: {', '.join(matches[:5])}"
            elif len(self.valid_ids) <= 10:
                self.suggestion = f"Valid cluster IDs are: {', '.join(sorted(self.valid_ids))}"
            else:
                # Show sample of valid IDs
                sample = sorted(self.valid_ids)[:5]
                self.suggestion = (
                    f"Valid cluster IDs include: {', '.join(sample)}... "
                    f"({len(self.valid_ids)} total)"
                )


@dataclass
class MissingProvenanceError(ValidationError):
    """Error when provenance keys are missing from adata.uns."""

    error_code: str = "E003_MISSING_PROVENANCE"
    expected_keys: List[str] = field(default_factory=list)
    found_keys: List[str] = field(default_factory=list)

    def __post_init__(self):
        """Generate suggestion for provenance issues."""
        if not self.suggestion:
            if self.expected_keys:
                self.suggestion = (
                    f"Expected one of: {', '.join(self.expected_keys)}. "
                    "Ensure input is from clustering stage or a previous refinement iteration."
                )


@dataclass
class PathSeparatorError(ValidationError):
    """Error for unrecognized path separators in hierarchical paths."""

    error_code: str = "E004_PATH_SEPARATOR"
    detected_separator: str = ""
    supported_separators: List[str] = field(default_factory=list)

    def __post_init__(self):
        """Generate suggestion for path separator issues."""
        if not self.suggestion and self.supported_separators:
            self.suggestion = (
                f"Supported separators: {self.supported_separators}. "
                "Use ' / ' (space-slash-space) as the canonical separator."
            )


@dataclass
class TypeMismatchError(ValidationError):
    """Error for unexpected data types."""

    error_code: str = "E005_TYPE_MISMATCH"
    column_name: str = ""
    expected_type: str = ""
    found_type: str = ""

    def __post_init__(self):
        """Generate type conversion suggestion."""
        if not self.suggestion and self.expected_type and self.found_type:
            self.suggestion = (
                f"Convert column '{self.column_name}' from {self.found_type} "
                f"to {self.expected_type}."
            )


@dataclass
class ValidationResult:
    """Result of validation with errors, warnings, and adaptation info.

    Attributes
    ----------
    is_valid : bool
        True if validation passed (no errors)
    errors : List[ValidationError]
        List of validation errors (empty if valid)
    warnings : List[str]
        Non-fatal warnings (e.g., column remapping)
    adapted_columns : Dict[str, str]
        Mapping from canonical column name to actual column name
    missing_columns : List[str]
        List of canonical column names that could not be found
    suggested_fixes : List[str]
        Aggregated suggestions from all errors

    Example
    -------
    >>> result = adapter.validate_adata(adata)
    >>> if result.is_valid:
    ...     cluster_col = result.adapted_columns["cluster_lvl0"]
    ... else:
    ...     result.raise_if_invalid("input AnnData")
    """

    is_valid: bool
    errors: List[ValidationError] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    adapted_columns: Dict[str, str] = field(default_factory=dict)
    missing_columns: List[str] = field(default_factory=list)
    suggested_fixes: List[str] = field(default_factory=list)

    def raise_if_invalid(self, context: str = "") -> None:
        """Raise ValueError if validation failed.

        Parameters
        ----------
        context : str
            Additional context for the error message (e.g., "cluster_annotations.csv")

        Raises
        ------
        ValueError
            If validation failed, with formatted error details
        """
        if not self.is_valid:
            error_msgs = [str(e) for e in self.errors]
            context_str = f" for {context}" if context else ""
            msg = f"Validation failed{context_str}:\n\n" + "\n\n".join(error_msgs)
            raise ValueError(msg)

    def log_warnings(self, logger: logging.Logger) -> None:
        """Log all warnings using the provided logger.

        Parameters
        ----------
        logger : logging.Logger
            Logger to use for warning output
        """
        for warning in self.warnings:
            logger.warning(warning)

    def merge(self, other: "ValidationResult") -> "ValidationResult":
        """Merge another ValidationResult into this one.

        Useful for combining results from multiple validation steps.

        Parameters
        ----------
        other : ValidationResult
            Another validation result to merge

        Returns
        -------
        ValidationResult
            New result combining both
        """
        return ValidationResult(
            is_valid=self.is_valid and other.is_valid,
            errors=self.errors + other.errors,
            warnings=self.warnings + other.warnings,
            adapted_columns={**self.adapted_columns, **other.adapted_columns},
            missing_columns=self.missing_columns + other.missing_columns,
            suggested_fixes=self.suggested_fixes + other.suggested_fixes,
        )

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            "is_valid": self.is_valid,
            "errors": [e.to_dict() for e in self.errors],
            "warnings": self.warnings,
            "adapted_columns": self.adapted_columns,
            "missing_columns": self.missing_columns,
            "suggested_fixes": self.suggested_fixes,
        }
