from __future__ import annotations

from pathlib import Path
from typing import Iterable

import pandas as pd

from .models import DagPeptides, Proteome

AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
DEFAULT_POSITION_COLUMNS = ("P4", "P3", "P2", "P1", "P1'", "P2'", "P3'", "P4'")


def prepare_proteome(fasta_path: str | Path, species: str = "unknown") -> Proteome:
    """Parse a FASTA file into a :class:`Proteome` object.

    Parameters
    ----------
    fasta_path
        Path to a FASTA file containing protein sequences.
    species
        Short species label stored on the returned proteome object.

    Returns
    -------
    Proteome
        Parsed protein identifiers and uppercase sequences, plus lightweight
        caches used by mapping and background generation.
    """

    path = Path(fasta_path)
    identifiers: list[str] = []
    sequences: list[str] = []
    current_id: str | None = None
    current_sequence: list[str] = []

    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    identifiers.append(current_id)
                    sequences.append("".join(current_sequence).upper())
                current_id = line[1:]
                current_sequence = []
            else:
                current_sequence.append(line)

    if current_id is not None:
        identifiers.append(current_id)
        sequences.append("".join(current_sequence).upper())

    proteome = pd.DataFrame({"SEQUENCE": sequences, "ID": identifiers})
    proteome["LEN"] = proteome["SEQUENCE"].str.len()
    return Proteome(proteome=proteome, type="fasta", species=species)


def _clean_sequence(sequence: str) -> str:
    """Normalize and validate a fixed-width amino-acid window."""
    normalized = sequence.upper().strip()
    invalid = sorted({char for char in normalized if char not in AMINO_ACIDS})
    if invalid:
        invalid_text = "".join(invalid)
        raise ValueError(
            f"Sequence window contains non-canonical amino-acid characters: {invalid_text!r}."
        )
    return normalized


def _normalize_sequence_value(value: object, *, canonical_only: bool = True) -> str:
    """Normalize one table cell into an uppercase amino-acid substring."""
    if pd.isna(value):
        return ""
    sequence = str(value).upper().strip()
    if canonical_only:
        return "".join(char for char in sequence if char in AMINO_ACIDS)
    return sequence


def build_windows_from_flanks(
    frame: pd.DataFrame,
    before_col: str = "five_res_before",
    after_col: str = "five_res_after",
    left: int = 4,
    right: int = 4,
    output_col: str = "window",
    canonical_only: bool = True,
    drop_invalid: bool = True,
) -> pd.DataFrame:
    """Build fixed-width windows from upstream and downstream flank columns.

    Parameters
    ----------
    frame
        Input table containing flank columns.
    before_col, after_col
        Column names holding residues before and after the cleavage site.
    left, right
        Number of residues to take from the end of ``before_col`` and the start
        of ``after_col``.
    output_col
        Name of the generated window column.
    canonical_only
        When ``True``, retain only canonical amino-acid letters.
    drop_invalid
        When ``True``, drop rows whose constructed windows do not match the
        expected fixed width.

    Returns
    -------
    pandas.DataFrame
        Copy of the input frame with a new fixed-width sequence window column.

    Raises
    ------
    KeyError
        If either flank column is missing.
    ValueError
        If ``left`` or ``right`` is negative.
    """

    if before_col not in frame.columns:
        raise KeyError(f"Missing before_col in frame: {before_col}")
    if after_col not in frame.columns:
        raise KeyError(f"Missing after_col in frame: {after_col}")
    if left < 0 or right < 0:
        raise ValueError("left and right must be non-negative.")

    prepared = frame.copy()
    before_values = [
        _normalize_sequence_value(value, canonical_only=canonical_only)
        for value in prepared[before_col]
    ]
    after_values = [
        _normalize_sequence_value(value, canonical_only=canonical_only)
        for value in prepared[after_col]
    ]

    def _build_window(before: str, after: str) -> str:
        before_part = before[-left:] if left else ""
        after_part = after[:right] if right else ""
        return before_part + after_part

    prepared[output_col] = [
        _build_window(before, after)
        for before, after in zip(before_values, after_values)
    ]

    if drop_invalid:
        expected_width = left + right
        prepared = prepared[prepared[output_col].str.len() == expected_width].copy()

    return prepared


def build_windows_from_positions(
    frame: pd.DataFrame,
    position_cols: tuple[str, ...] = DEFAULT_POSITION_COLUMNS,
    output_col: str = "window",
    canonical_only: bool = True,
    drop_invalid: bool = True,
) -> pd.DataFrame:
    """Build fixed-width sequence windows from per-position residue columns.

    Parameters
    ----------
    frame
        Input table containing one residue column per position.
    position_cols
        Ordered residue columns to concatenate into the output window.
    output_col
        Name of the generated window column.
    canonical_only
        When ``True``, retain only canonical amino-acid letters in each field.
    drop_invalid
        When ``True``, drop rows whose constructed windows do not match the
        expected fixed width.

    Returns
    -------
    pandas.DataFrame
        Copy of the input frame with a new fixed-width sequence window column.
    """

    if not position_cols:
        raise ValueError("position_cols must not be empty.")

    missing = [column for column in position_cols if column not in frame.columns]
    if missing:
        raise KeyError(f"Missing position columns in frame: {missing}")

    prepared = frame.copy()
    expected_width = len(position_cols)
    prepared[output_col] = prepared.apply(
        lambda row: "".join(
            _normalize_sequence_value(row[column], canonical_only=canonical_only)
            for column in position_cols
        ),
        axis=1,
    )

    if drop_invalid:
        prepared = prepared[prepared[output_col].str.len() == expected_width].copy()

    return prepared


def _find_window_match(proteome: Proteome, window: str) -> tuple[str, int, str] | None:
    """Return the first cached or discovered proteome match for a window."""
    if window in proteome.window_match_cache:
        return proteome.window_match_cache[window]

    proteome_sequences = proteome.proteome["SEQUENCE"].tolist()
    proteome_ids = proteome.proteome["ID"].tolist()
    match: tuple[str, int, str] | None = None

    for idx, peptide in enumerate(proteome_sequences):
        found = peptide.find(window)
        if found != -1:
            match = (proteome_ids[idx], found + 1, peptide)
            break

    proteome.window_match_cache[window] = match
    return match


def format_sequence(
    seq: Iterable[str],
    proteome: Proteome | None,
    upstream_offset: int | None = None,
    downstream_offset: int | None = None,
    map_to_proteome: bool = True,
) -> DagPeptides:
    """Convert fixed-width windows into a :class:`DagPeptides` object.

    Parameters
    ----------
    seq
        Iterable of sequence windows with equal width.
    proteome
        Proteome used for sequence-to-proteome mapping. Required only when
        ``map_to_proteome=True``.
    upstream_offset, downstream_offset
        Number of residues left and right of the anchor position. If both are
        omitted, the anchor defaults to the center residue for odd-width
        windows.
    map_to_proteome
        When ``True``, each input window is located in the supplied proteome
        and annotated with protein identifiers and positions. When ``False``,
        windows are kept as unmapped inputs.

    Returns
    -------
    DagPeptides
        Structured peptide representation used by background generation and
        differential amino-acid usage testing.

    Raises
    ------
    ValueError
        If the input is empty, contains mixed widths, has inconsistent offsets,
        contains non-canonical amino-acid characters, or requires proteome
        mapping without a proteome.

    Notes
    -----
    For even-width cleavage windows, higher-level helpers such as
    :func:`pydaglogo.run_daglogo` choose the residue just left of center as the
    default anchor. This function only uses that behavior when offsets are
    explicitly supplied by the caller.
    """

    sequences = [_clean_sequence(str(item)) for item in seq]
    if not sequences:
        raise ValueError("seq is required and must not be empty.")

    widths = {len(sequence) for sequence in sequences}
    if len(widths) != 1:
        raise ValueError("seq must contain characters of the same length.")
    width = widths.pop()

    if upstream_offset is None and downstream_offset is None:
        upstream_offset = width // 2
        downstream_offset = width - upstream_offset - 1
    elif upstream_offset is None:
        upstream_offset = width - downstream_offset - 1
    elif downstream_offset is None:
        downstream_offset = width - upstream_offset - 1

    if upstream_offset < 0 or downstream_offset < 0:
        raise ValueError("upstream_offset and downstream_offset must be non-negative.")
    if upstream_offset + downstream_offset + 1 != width:
        raise ValueError("upstream_offset and downstream_offset do not match the sequence width.")

    center_index = upstream_offset
    ids: list[str] = []
    anchor_positions: list[int] = []
    matched_peptides: list[str] = []
    anchors: list[str] = []
    upstream_sequences: list[str] = []
    downstream_sequences: list[str] = []
    matrix_rows: list[list[str]] = []

    for window in sequences:
        if map_to_proteome:
            if proteome is None:
                raise ValueError("proteome is required when map_to_proteome=True.")
            match = _find_window_match(proteome, window)
            if match is None:
                raise ValueError(f"Could not map sequence window to proteome: {window}")
            matched_id, match_pos, matched_peptide = match
        else:
            matched_id = "UNMAPPED"
            match_pos = center_index + 1
            matched_peptide = window

        anchor = window[center_index]
        anchor_pos = match_pos + center_index if map_to_proteome else match_pos
        ids.append(matched_id)
        anchor_positions.append(anchor_pos)
        matched_peptides.append(matched_peptide)
        anchors.append(anchor)
        upstream_sequences.append(window[:upstream_offset])
        downstream_sequences.append(window[upstream_offset + 1 :])
        matrix_rows.append(list(window))

    data = pd.DataFrame(
        {
            "IDs": ids,
            "anchorAA": anchors,
            "anchorPos": anchor_positions,
            "peptide": matched_peptides,
            "anchor": anchors,
            "upstream": upstream_sequences,
            "downstream": downstream_sequences,
        }
    )
    columns = [f"V{i}" for i in range(1, width + 1)]
    peptides = pd.DataFrame(matrix_rows, columns=columns, index=range(1, len(matrix_rows) + 1))
    result_type = proteome.type if proteome is not None else "unmapped"
    return DagPeptides(
        data=data,
        peptides=peptides,
        upstream_offset=upstream_offset,
        downstream_offset=downstream_offset,
        type=result_type,
    )
