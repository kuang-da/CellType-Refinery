"""Centralized configuration for CellType-Refinery.

This module provides organ-specific configuration (region ordering, colors)
that can be used by all core modules (composition, annotation, preprocessing,
spatial).

Example
-------
>>> from celltype_refinery.config import get_organ_config, list_available_organs
>>>
>>> # List available organs
>>> print(list_available_organs())
['fallopian_tube', 'generic', 'uterus']
>>>
>>> # Get config for fallopian tube (with alias)
>>> config = get_organ_config("ft")
>>> print(config.region_order)
['fimbriae', 'infundibulum', 'ampulla', 'isthmus']
>>>
>>> # Sort regions by canonical order
>>> sorted_regions = config.sort_regions(["isthmus", "ampulla", "fimbriae"])
>>> print(sorted_regions)
['fimbriae', 'ampulla', 'isthmus']
>>>
>>> # Get region color
>>> print(config.get_region_color("fimbriae"))
'#9b59b6'
"""

from .organ import (
    OrganConfig,
    get_organ_config,
    list_available_organs,
    list_organ_aliases,
    register_organ_config,
)

__all__ = [
    "OrganConfig",
    "get_organ_config",
    "list_available_organs",
    "list_organ_aliases",
    "register_organ_config",
]
