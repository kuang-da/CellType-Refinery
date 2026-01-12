"""Label normalization utilities."""


def to_title_case(label: str) -> str:
    """
    Normalize label to Title Case.

    Examples:
        "helper t cells" -> "Helper T Cells"
        "T cells" -> "T Cells"
        "FIBROBLASTS" -> "Fibroblasts"
    """
    return " ".join(word.capitalize() for word in label.strip().split())
