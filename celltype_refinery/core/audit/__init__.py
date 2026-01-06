"""Audit cards module for cell-type annotation review.

This module generates HTML audit cards for each cluster showing:
- Decision path trace through the hierarchy
- Top competing cell types and their scores
- Marker evidence
- QC flags and warnings
"""

from .config import AuditConfig
from .engine import AuditEngine, AuditResult
from .card import generate_audit_card_html, generate_all_audit_cards

__all__ = [
    "AuditConfig",
    "AuditEngine",
    "AuditResult",
    "generate_audit_card_html",
    "generate_all_audit_cards",
]
