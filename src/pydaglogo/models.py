from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd


@dataclass(slots=True)
class Proteome:
    """Proteome sequences plus caches used by the analysis workflow.

    The ``proteome`` dataframe stores parsed protein IDs and sequences. The
    cache fields are internal performance helpers used by mapping and
    background-generation code.
    """

    proteome: pd.DataFrame
    type: str = "fasta"
    species: str = "unknown"
    window_match_cache: dict[str, tuple[str, int, str] | None] = field(default_factory=dict)
    window_sampling_cache: dict[tuple[int, str], dict[str, object]] = field(default_factory=dict)
    background_match_cache: dict[tuple[object, ...], list[str]] = field(default_factory=dict)
    background_model_cache: dict[tuple[object, ...], "BackgroundModel"] = field(default_factory=dict)


@dataclass(slots=True)
class DagPeptides:
    """Formatted sequence windows ready for background generation and testing.

    ``data`` stores metadata such as anchors and mapped IDs, while ``peptides``
    stores the aligned residue matrix in ``V1..Vn`` columns.
    """

    data: pd.DataFrame
    peptides: pd.DataFrame
    upstream_offset: int
    downstream_offset: int
    type: str = "fasta"

    @property
    def width(self) -> int:
        """Total window width derived from upstream and downstream offsets."""
        return self.upstream_offset + self.downstream_offset + 1


@dataclass(slots=True)
class BackgroundModel:
    """Background samples and metadata for ztest or Fisher comparisons.

    For ``ztest``, ``background`` contains multiple sampled matrices. For
    ``fisher``, it contains a single deterministic background matrix.
    """

    background: list[pd.DataFrame]
    num_subsamples: int
    test_type: str
    source_background: str
    model: str
    target_position: str
    unique_seq: bool
    rand_seed: int
    replacement: bool


@dataclass(slots=True)
class GroupingScheme:
    """Residue grouping definition used to collapse amino-acid alphabets.

    Grouping schemes can preserve the full amino-acid alphabet or collapse it
    into biochemical groups for testing and plotting.
    """

    name: str
    color: dict[str, str]
    symbols: dict[str, str]
    groups: dict[str, tuple[str, ...]] | None = None

    @property
    def order(self) -> list[str]:
        """Display order for groups or individual residues."""
        if self.groups is None:
            return list(self.symbols.values())
        return list(self.groups.keys())

    @property
    def plot_symbols(self) -> dict[str, str]:
        """Single-character plot symbols derived from the grouping definition."""
        if self.groups is None:
            return {letter: letter for letter in self.order}
        assigned: dict[str, str] = {}
        used: set[str] = set()
        fallback = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"
        for key in self.order:
            desired = self.symbols.get(key, "")
            candidate = desired[:1] if desired else ""
            if not candidate or candidate in used:
                for char in fallback:
                    if char not in used:
                        candidate = char
                        break
            assigned[key] = candidate
            used.add(candidate)
        return assigned


@dataclass(slots=True)
class DAUResult:
    """Differential amino-acid usage result matrices for one analysis run.

    The result keeps the raw inferential outputs (`difference`, `statistics`,
    `pvalue`) alongside descriptive matrices (`background`, `motif`) and the
    metadata needed to plot them consistently.
    """

    difference: pd.DataFrame
    statistics: pd.DataFrame
    pvalue: pd.DataFrame
    background: pd.DataFrame
    motif: pd.DataFrame
    test_type: str
    upstream_offset: int
    downstream_offset: int
    grouping_scheme: str
    plot_symbols: dict[str, str] = field(default_factory=dict)
