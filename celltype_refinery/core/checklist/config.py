"""Configuration for review checklist generation."""

from dataclasses import dataclass


@dataclass
class ChecklistConfig:
    """Configuration for review checklist generation.

    Attributes:
        stage_name: Name to display in checklist header.
        high_priority_threshold: Priority score threshold for high priority.
        medium_priority_threshold: Priority score threshold for medium priority.
        low_priority_max_items: Max low priority items to show in checklist.

        Priority weights for computing review priority:
        - low_conf_weight: Weight for low confidence (score < 1.0)
        - medium_low_conf_weight: Weight for medium-low confidence (score < 1.5)
        - large_cluster_weight: Weight for large clusters (> 50k cells)
        - ambiguous_siblings_weight: Weight for ambiguous siblings
        - low_margin_weight: Weight for low margin (< 0.3)
        - low_coverage_weight: Weight for low coverage (< 0.75)
        - stopped_before_leaf_weight: Weight for needs refinement
    """

    stage_name: str = "Stage H"
    high_priority_threshold: int = 100
    medium_priority_threshold: int = 50
    low_priority_max_items: int = 10

    # Priority weights
    low_conf_weight: int = 100
    medium_low_conf_weight: int = 50
    large_cluster_weight: int = 30
    medium_cluster_weight: int = 20
    small_cluster_weight: int = 10
    ambiguous_siblings_weight: int = 40
    low_margin_weight: int = 35
    low_coverage_weight: int = 25
    medium_coverage_weight: int = 10
    stopped_before_leaf_weight: int = 30
