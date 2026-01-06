"""Configuration for audit card generation."""

from dataclasses import dataclass


@dataclass
class AuditConfig:
    """Configuration for audit card generation.

    Attributes:
        stage_name: Name to display in HTML title/header.
        include_marker_evidence: Whether to include marker evidence section.
        include_decision_steps: Whether to include gate summary section.
        high_confidence_threshold: Score threshold for high confidence.
        medium_confidence_threshold: Score threshold for medium confidence.
        low_confidence_threshold: Score threshold for low confidence.
    """

    stage_name: str = "Stage H"
    include_marker_evidence: bool = True
    include_decision_steps: bool = True
    high_confidence_threshold: float = 2.0
    medium_confidence_threshold: float = 1.0
    low_confidence_threshold: float = 0.5
