"""Rule registry for dynamic rule discovery and instantiation.

Provides decorator-based registration of rule classes.
"""

from typing import Dict, List, Optional, Type, TYPE_CHECKING

if TYPE_CHECKING:
    from .base import BaseRule
    from ..config import TissueTemplate, RuleConfig


class RuleRegistry:
    """Registry for flagging rules.

    Provides:
    - Decorator-based registration: @RuleRegistry.register
    - Rule discovery by ID or category
    - Instantiation with template and config
    """

    _rules: Dict[str, Type["BaseRule"]] = {}

    @classmethod
    def register(cls, rule_class: Type["BaseRule"]) -> Type["BaseRule"]:
        """Register a rule class.

        Use as decorator:
            @RuleRegistry.register
            class MyRule(BaseRule):
                rule_id = "MY_RULE"
                ...

        Parameters
        ----------
        rule_class : Type[BaseRule]
            Rule class to register

        Returns
        -------
        Type[BaseRule]
            The registered class (unchanged)
        """
        cls._rules[rule_class.rule_id] = rule_class
        return rule_class

    @classmethod
    def get_rule(cls, rule_id: str) -> Optional[Type["BaseRule"]]:
        """Get rule class by ID.

        Parameters
        ----------
        rule_id : str
            Rule identifier

        Returns
        -------
        Type[BaseRule] or None
            Rule class if found
        """
        return cls._rules.get(rule_id)

    @classmethod
    def get_all_rules(cls) -> Dict[str, Type["BaseRule"]]:
        """Get all registered rules.

        Returns
        -------
        Dict[str, Type[BaseRule]]
            Map of rule_id to rule class
        """
        return cls._rules.copy()

    @classmethod
    def get_rules_by_category(cls, category: str) -> List[Type["BaseRule"]]:
        """Get all rules in a category.

        Parameters
        ----------
        category : str
            Category name (proportion, spatial, biology, quality)

        Returns
        -------
        List[Type[BaseRule]]
            List of rule classes in category
        """
        return [r for r in cls._rules.values() if r.category == category]

    @classmethod
    def list_rule_ids(cls) -> List[str]:
        """Get list of all registered rule IDs.

        Returns
        -------
        List[str]
            Sorted list of rule IDs
        """
        return sorted(cls._rules.keys())

    @classmethod
    def instantiate(
        cls,
        rule_id: str,
        template: "TissueTemplate",
        config: Optional["RuleConfig"] = None,
    ) -> Optional["BaseRule"]:
        """Instantiate a rule by ID.

        Parameters
        ----------
        rule_id : str
            Rule identifier
        template : TissueTemplate
            Tissue template
        config : RuleConfig, optional
            Rule configuration

        Returns
        -------
        BaseRule or None
            Instantiated rule if found
        """
        rule_class = cls._rules.get(rule_id)
        if rule_class is None:
            return None
        return rule_class(template=template, config=config)

    @classmethod
    def instantiate_all(
        cls,
        template: "TissueTemplate",
        configs: Optional[Dict[str, "RuleConfig"]] = None,
        skip_rules: Optional[List[str]] = None,
    ) -> List["BaseRule"]:
        """Instantiate all registered rules.

        Parameters
        ----------
        template : TissueTemplate
            Tissue template
        configs : Dict[str, RuleConfig], optional
            Per-rule configurations
        skip_rules : List[str], optional
            Rule IDs to skip

        Returns
        -------
        List[BaseRule]
            List of instantiated rules
        """
        configs = configs or {}
        skip_rules = skip_rules or []

        rules = []
        for rule_id, rule_class in cls._rules.items():
            if rule_id in skip_rules:
                continue

            config = configs.get(rule_id)

            # Skip if config explicitly disables
            if config and not config.enabled:
                continue

            rules.append(rule_class(template=template, config=config))

        return rules

    @classmethod
    def clear(cls) -> None:
        """Clear all registered rules (for testing)."""
        cls._rules.clear()

    @classmethod
    def summary(cls) -> Dict[str, List[str]]:
        """Get summary of rules by category.

        Returns
        -------
        Dict[str, List[str]]
            Map of category to list of rule IDs
        """
        summary: Dict[str, List[str]] = {}
        for rule_id, rule_class in cls._rules.items():
            category = rule_class.category
            if category not in summary:
                summary[category] = []
            summary[category].append(rule_id)

        # Sort within categories
        for cat in summary:
            summary[cat].sort()

        return summary
