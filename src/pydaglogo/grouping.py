from __future__ import annotations

import pandas as pd

from .models import GroupingScheme

NO_GROUPING = GroupingScheme(
    name="no",
    color={
        "A": "#CCFF00",
        "C": "#FFFF00",
        "D": "#FF0000",
        "E": "#FF0066",
        "F": "#00FF66",
        "G": "#FF9900",
        "H": "#0066FF",
        "I": "#66FF00",
        "K": "#6600FF",
        "L": "#33FF00",
        "M": "#00FF00",
        "N": "#CC00FF",
        "P": "#FFCC00",
        "Q": "#FF00CC",
        "R": "#0000FF",
        "S": "#FF3300",
        "T": "#FF6600",
        "V": "#99FF00",
        "W": "#00CCFF",
        "Y": "#00FFCC",
    },
    symbols={
        "Ala": "A",
        "Arg": "R",
        "Asn": "N",
        "Asp": "D",
        "Cys": "C",
        "Glu": "E",
        "Gln": "Q",
        "Gly": "G",
        "His": "H",
        "Ile": "I",
        "Leu": "L",
        "Lys": "K",
        "Met": "M",
        "Phe": "F",
        "Pro": "P",
        "Ser": "S",
        "Thr": "T",
        "Trp": "W",
        "Tyr": "Y",
        "Val": "V",
    },
)

CHARGE_GROUP = GroupingScheme(
    name="charge_group",
    color={"positive": "#FFB32C", "neutral": "#2000C7", "negative": "#CCCCCC"},
    symbols={"positive": "P", "neutral": "U", "negative": "N"},
    groups={
        "positive": ("H", "K", "R"),
        "neutral": ("A", "C", "F", "G", "I", "L", "M", "N", "P", "Q", "S", "T", "V", "W", "Y"),
        "negative": ("D", "E"),
    },
)

AVAILABLE_GROUPING_SCHEMES: dict[str, GroupingScheme] = {
    "no": NO_GROUPING,
    "charge_group": CHARGE_GROUP,
}


def resolve_grouping_scheme(grouping_scheme: str | GroupingScheme | None) -> GroupingScheme:
    """Resolve a grouping-scheme name or object into a concrete scheme."""
    if grouping_scheme is None:
        return NO_GROUPING
    if isinstance(grouping_scheme, GroupingScheme):
        return grouping_scheme
    if grouping_scheme not in AVAILABLE_GROUPING_SCHEMES:
        raise KeyError(f"Unknown grouping scheme: {grouping_scheme}")
    return AVAILABLE_GROUPING_SCHEMES[grouping_scheme]


def apply_grouping_matrix(matrix: pd.DataFrame, grouping_scheme: str | GroupingScheme | None) -> pd.DataFrame:
    """Collapse an aligned residue matrix according to a grouping scheme."""
    scheme = resolve_grouping_scheme(grouping_scheme)
    if scheme.groups is None:
        return matrix.copy()
    mapping: dict[str, str] = {}
    for group_name, members in scheme.groups.items():
        for member in members:
            mapping[member] = group_name
    return matrix.replace(mapping)


def ordered_index_for_scheme(grouping_scheme: str | GroupingScheme | None) -> list[str]:
    """Return the display/counting order for one grouping scheme."""
    return resolve_grouping_scheme(grouping_scheme).order


def plot_symbols_for_scheme(grouping_scheme: str | GroupingScheme | None) -> dict[str, str]:
    """Return one-character plotting symbols for one grouping scheme."""
    return resolve_grouping_scheme(grouping_scheme).plot_symbols
