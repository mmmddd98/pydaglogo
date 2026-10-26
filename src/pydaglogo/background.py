from __future__ import annotations

import re
import warnings
from collections.abc import Iterable
from pathlib import Path

import numpy as np
import pandas as pd

from .models import BackgroundModel, DagPeptides, Proteome

BACKGROUND_ALIASES = {
    "wholeProteome": "wholeProteome",
    "whole_proteome": "wholeProteome",
    "inputSet": "inputSet",
    "input_set": "inputSet",
    "nonInputSet": "nonInputSet",
    "non_input_set": "nonInputSet",
}

MODEL_ALIASES = {"any": "any", "anchored": "anchored"}
TARGET_ALIASES = {"any": "any", "Nterminus": "Nterminus", "Cterminus": "Cterminus"}
TEST_TYPE_ALIASES = {"ztest": "ztest", "fisher": "fisher"}


def _normalize(value: str, aliases: dict[str, str], label: str) -> str:
    """Normalize a user-facing option string through an alias table."""
    if value not in aliases:
        raise ValueError(f"Unknown {label}: {value}")
    return aliases[value]


def _regex_matches(sequence: str, pattern: str) -> list[str]:
    """Collect all captured regex window matches from a sequence."""
    return [match.group(1) for match in re.finditer(pattern, sequence)]


def _sequence_window_count(sequence_length: int, width: int, target_position: str) -> int:
    """Return how many windows of one width can be drawn from a sequence."""
    if sequence_length < width:
        return 0
    if target_position in {"Nterminus", "Cterminus"}:
        return 1
    return sequence_length - width + 1


def _sequence_start_from_offset(sequence_length: int, width: int, target_position: str, offset: int) -> int:
    """Convert a per-sequence offset into an absolute window start."""
    if target_position == "Nterminus":
        return 0
    if target_position == "Cterminus":
        return sequence_length - width
    return offset


def _build_sampling_pool(proteome: Proteome, width: int, target_position: str) -> dict[str, object]:
    """Build or reuse the cached global sampling pool for whole-proteome ztest backgrounds."""
    cache_key = (width, target_position)
    if cache_key in proteome.window_sampling_cache:
        return proteome.window_sampling_cache[cache_key]

    sequences = proteome.proteome["SEQUENCE"].tolist()
    lengths = np.asarray([len(sequence) for sequence in sequences], dtype=np.int64)
    counts = np.asarray(
        [_sequence_window_count(length, width, target_position) for length in lengths],
        dtype=np.int64,
    )
    valid_mask = counts > 0
    valid_sequences = [sequence for sequence, keep in zip(sequences, valid_mask.tolist()) if keep]
    valid_lengths = lengths[valid_mask]
    valid_counts = counts[valid_mask]
    if valid_counts.size == 0:
        raise ValueError("No background matches were generated.")

    pool = {
        "sequences": valid_sequences,
        "lengths": valid_lengths,
        "counts": valid_counts,
        "cum_counts": np.cumsum(valid_counts, dtype=np.int64),
        "starts_before": np.concatenate(([0], np.cumsum(valid_counts[:-1], dtype=np.int64))),
        "total": int(valid_counts.sum()),
        "target_position": target_position,
        "width": width,
    }
    proteome.window_sampling_cache[cache_key] = pool
    return pool


def _sample_global_indices(
    total: int,
    sample_size: int,
    random_state: np.random.RandomState,
    replacement: bool,
) -> np.ndarray:
    """Sample global window indices from the cached whole-proteome pool."""
    if replacement:
        return random_state.randint(0, total, size=sample_size)

    if sample_size > total:
        raise ValueError("Too few matches in the background. Please try different parameters.")

    # For tiny samples drawn from a huge global window pool, np.random.choice(..., replace=False)
    # can be disproportionately slow. A rejection-based sampler is much faster here.
    if sample_size * 32 < total:
        sampled = np.empty(sample_size, dtype=np.int64)
        seen: set[int] = set()
        filled = 0
        while filled < sample_size:
            remaining = sample_size - filled
            draw_count = max(remaining * 2, 8)
            for draw in random_state.randint(0, total, size=draw_count).tolist():
                if draw in seen:
                    continue
                seen.add(draw)
                sampled[filled] = draw
                filled += 1
                if filled == sample_size:
                    break
        return sampled

    return random_state.choice(total, size=sample_size, replace=False)


def _sample_windows_from_pool(
    pool: dict[str, object],
    sample_size: int,
    random_state: np.random.RandomState,
    replacement: bool,
) -> pd.DataFrame:
    """Materialize sampled sequence windows from a cached global pool."""
    total = int(pool["total"])
    sampled_indices = _sample_global_indices(total, sample_size, random_state, replacement)
    cum_counts = pool["cum_counts"]
    seq_indices = np.searchsorted(cum_counts, sampled_indices, side="right")
    starts_before = pool["starts_before"]
    offsets = sampled_indices - starts_before[seq_indices]
    sequences = pool["sequences"]
    lengths = pool["lengths"]
    width = int(pool["width"])
    target_position = str(pool["target_position"])

    sampled_windows = [
        sequences[seq_idx][
            start : start + width
        ]
        for seq_idx, start in zip(
            seq_indices.tolist(),
            [
                _sequence_start_from_offset(int(lengths[seq_idx]), width, target_position, int(offset))
                for seq_idx, offset in zip(seq_indices.tolist(), offsets.tolist())
            ],
        )
    ]
    return _split_matches(sampled_windows)


def _clone_background_model(background_model: BackgroundModel) -> BackgroundModel:
    """Return a deep-copy-style clone suitable for cache reuse."""
    return BackgroundModel(
        background=[frame.copy(deep=True) for frame in background_model.background],
        num_subsamples=background_model.num_subsamples,
        test_type=background_model.test_type,
        source_background=background_model.source_background,
        model=background_model.model,
        target_position=background_model.target_position,
        unique_seq=background_model.unique_seq,
        rand_seed=background_model.rand_seed,
        replacement=background_model.replacement,
    )


def _collect_matches(
    sequences: Iterable[str],
    width: int,
    model: str,
    target_position: str,
    anchor_letters: list[str],
    upstream_offset: int,
    downstream_offset: int,
) -> list[str]:
    """Collect candidate background windows under one background-extraction rule."""
    matches: list[str] = []
    anchored_mode = model == "anchored" and len(set(anchor_letters)) == 1
    if model == "anchored" and not anchored_mode:
        warnings.warn("Anchor amino acid is not unique. Model is set to 'any'.", stacklevel=2)

    if anchored_mode:
        anchor = anchor_letters[0]

    for sequence in sequences:
        sequence_length = len(sequence)
        count = _sequence_window_count(sequence_length, width, target_position)
        if count == 0:
            continue
        for offset in range(count):
            start = _sequence_start_from_offset(sequence_length, width, target_position, offset)
            window = sequence[start : start + width]
            if anchored_mode and window[upstream_offset] != anchor:
                continue
            matches.append(window)
    return matches


def _split_matches(matches: Iterable[str]) -> pd.DataFrame:
    """Convert raw window strings into a V1..Vn residue matrix."""
    rows = [list(match) for match in matches]
    if not rows:
        raise ValueError("No background matches were generated.")
    width = len(rows[0])
    return pd.DataFrame(rows, columns=[f"V{i}" for i in range(1, width + 1)])


def _cached_whole_proteome_matches(
    proteome: Proteome,
    dag_peptides: DagPeptides,
    target_position: str,
) -> list[str]:
    """Return cached deterministic whole-proteome matches for Fisher backgrounds."""
    cache_key = ("wholeProteome", dag_peptides.width, target_position, "any")
    if cache_key in proteome.background_match_cache:
        return proteome.background_match_cache[cache_key]

    matches = _collect_matches(
        sequences=proteome.proteome["SEQUENCE"].tolist(),
        width=dag_peptides.width,
        model="any",
        target_position=target_position,
        anchor_letters=dag_peptides.data["anchorAA"].tolist(),
        upstream_offset=dag_peptides.upstream_offset,
        downstream_offset=dag_peptides.downstream_offset,
    )
    proteome.background_match_cache[cache_key] = matches
    return matches


def load_background_samples(path: str | Path) -> list[pd.DataFrame]:
    """Load committed background samples from a CSV or gzipped CSV bundle.

    Parameters
    ----------
    path
        CSV file, including ``.csv.gz`` bundles, containing a ``sample_id``
        column plus sequence columns named ``V1`` through ``Vn``.

    Returns
    -------
    list[pandas.DataFrame]
        One dataframe per background sample, each containing only the ``V``
        sequence columns.

    Raises
    ------
    FileNotFoundError
        If the provided file does not exist.
    ValueError
        If the CSV does not contain the expected columns or no samples can be
        loaded.
    """

    sample_path = Path(path)
    if not sample_path.exists():
        compressed_fallback = sample_path.with_name(sample_path.name + ".gz")
        if compressed_fallback.exists():
            sample_path = compressed_fallback
        else:
            raise FileNotFoundError(f"Background samples file not found: {sample_path}")

    frame = pd.read_csv(sample_path)
    required = {"sample_id"}
    if not required.issubset(frame.columns):
        raise ValueError("background samples file must contain a 'sample_id' column.")

    sequence_columns = [column for column in frame.columns if re.fullmatch(r"V\d+", str(column))]
    if not sequence_columns:
        raise ValueError("background samples file must contain V1..Vn sequence columns.")

    def _sort_key(value: str) -> int:
        return int(str(value)[1:])

    sequence_columns = sorted(sequence_columns, key=_sort_key)
    sample_frames: list[pd.DataFrame] = []
    for _, sample_frame in frame.groupby("sample_id", sort=True):
        sample_frames.append(sample_frame.loc[:, sequence_columns].reset_index(drop=True))

    if not sample_frames:
        raise ValueError("No background samples were loaded from the provided file.")
    return sample_frames


def build_background_model_from_samples(
    background_samples: Iterable[pd.DataFrame],
    *,
    test_type: str = "ztest",
    source_background: str = "externalSamples",
    model: str = "any",
    target_position: str = "any",
    unique_seq: bool = False,
    rand_seed: int = 1,
    replacement: bool = False,
) -> BackgroundModel:
    """Wrap precomputed background samples in a :class:`BackgroundModel`.

    Parameters
    ----------
    background_samples
        Iterable of sample matrices with sequence columns ``V1`` through ``Vn``.
    test_type
        Statistical mode associated with the samples: ``"ztest"`` or
        ``"fisher"``.
    source_background
        Free-text label describing the source of the committed samples.
    model, target_position, unique_seq, rand_seed, replacement
        Metadata copied into the returned background model.

    Returns
    -------
    BackgroundModel
        Background wrapper suitable for :func:`pydaglogo.test_dau`.

    Raises
    ------
    ValueError
        If no samples are provided, sample columns do not match, or Fisher is
        supplied with more than one background matrix.
    """

    normalized_model = _normalize(model, MODEL_ALIASES, "model")
    normalized_target = _normalize(target_position, TARGET_ALIASES, "target_position")
    normalized_test_type = _normalize(test_type, TEST_TYPE_ALIASES, "test_type")

    frames = [sample.copy(deep=True).reset_index(drop=True) for sample in background_samples]
    if not frames:
        raise ValueError("background_samples must not be empty.")
    if normalized_test_type == "fisher" and len(frames) != 1:
        raise ValueError("Fisher background samples must contain exactly one background matrix.")

    first_columns = list(frames[0].columns)
    for sample in frames:
        if list(sample.columns) != first_columns:
            raise ValueError("All background samples must have the same sequence columns.")

    return BackgroundModel(
        background=frames,
        num_subsamples=len(frames),
        test_type=normalized_test_type,
        source_background=source_background,
        model=normalized_model,
        target_position=normalized_target,
        unique_seq=unique_seq,
        rand_seed=rand_seed,
        replacement=replacement,
    )


def build_background_model(
    dag_peptides: DagPeptides,
    proteome: Proteome | None = None,
    background: str = "whole_proteome",
    model: str = "any",
    target_position: str = "any",
    unique_seq: bool = False,
    num_subsamples: int = 300,
    rand_seed: int = 1,
    replacement: bool = False,
    test_type: str = "ztest",
) -> BackgroundModel:
    """Generate a background model from a proteome or the input set itself.

    Parameters
    ----------
    dag_peptides
        Formatted input windows to compare against the background.
    proteome
        Proteome source used for ``whole_proteome`` or ``non_input_set``
        backgrounds.
    background
        Background source: ``whole_proteome``, ``input_set``, or
        ``non_input_set``.
    model
        Background extraction mode. ``"any"`` scans all matching windows, while
        ``"anchored"`` constrains windows to a shared center residue when that
        residue is unique.
    target_position
        Whether to sample any matching window, only N-termini, or only
        C-termini.
    unique_seq
        When ``True``, de-duplicate matched background windows before sampling.
    num_subsamples, rand_seed, replacement
        Sampling controls for Monte Carlo ``ztest`` backgrounds.
    test_type
        Statistical mode: ``"ztest"`` or ``"fisher"``.

    Returns
    -------
    BackgroundModel
        Background samples or a single Fisher background matrix.

    Raises
    ------
    ValueError
        If the requested background cannot be constructed, if a proteome is
        required but missing, or if too few background matches exist.

    Notes
    -----
    ``ztest`` backgrounds are sampled repeatedly and cached for the optimized
    whole-proteome path. ``fisher`` backgrounds are deterministic and contain a
    single background matrix.
    """

    normalized_background = _normalize(background, BACKGROUND_ALIASES, "background")
    normalized_model = _normalize(model, MODEL_ALIASES, "model")
    normalized_target = _normalize(target_position, TARGET_ALIASES, "target_position")
    normalized_test_type = _normalize(test_type, TEST_TYPE_ALIASES, "test_type")

    if normalized_test_type != "fisher" and num_subsamples < 2:
        raise ValueError("num_subsamples must be greater than 1.")

    n = len(dag_peptides.data)

    if (
        proteome is not None
        and normalized_background == "wholeProteome"
        and normalized_model == "any"
        and not unique_seq
        and normalized_test_type == "ztest"
    ):
        cache_key = (
            "wholeProteome",
            dag_peptides.width,
            normalized_target,
            n,
            int(num_subsamples),
            int(rand_seed),
            bool(replacement),
            normalized_test_type,
        )
        if cache_key in proteome.background_model_cache:
            return _clone_background_model(proteome.background_model_cache[cache_key])

        pool = _build_sampling_pool(proteome, dag_peptides.width, normalized_target)
        random_state = np.random.RandomState(rand_seed)
        background_samples = [
            _sample_windows_from_pool(pool, n, random_state, replacement)
            for _ in range(int(num_subsamples))
        ]
        background_model = BackgroundModel(
            background=background_samples,
            num_subsamples=int(num_subsamples),
            test_type=normalized_test_type,
            source_background=normalized_background,
            model=normalized_model,
            target_position=normalized_target,
            unique_seq=unique_seq,
            rand_seed=rand_seed,
            replacement=replacement,
        )
        proteome.background_model_cache[cache_key] = _clone_background_model(background_model)
        return background_model

    if (
        proteome is not None
        and normalized_background == "wholeProteome"
        and normalized_model == "any"
        and not unique_seq
        and normalized_test_type == "fisher"
    ):
        matches = _cached_whole_proteome_matches(proteome, dag_peptides, normalized_target)
        return BackgroundModel(
            background=[_split_matches(matches)],
            num_subsamples=1,
            test_type=normalized_test_type,
            source_background=normalized_background,
            model=normalized_model,
            target_position=normalized_target,
            unique_seq=unique_seq,
            rand_seed=rand_seed,
            replacement=replacement,
        )

    if normalized_background != "inputSet":
        if proteome is None:
            raise ValueError("proteome is required unless background='inputSet'.")
        if normalized_background == "wholeProteome":
            sequence_pool = proteome.proteome["SEQUENCE"].tolist()
        else:
            excluded = set(dag_peptides.data["peptide"].tolist())
            sequence_pool = [sequence for sequence in proteome.proteome["SEQUENCE"].tolist() if sequence not in excluded]
    else:
        sequence_pool = dag_peptides.data["peptide"].tolist()

    matches = _collect_matches(
        sequences=sequence_pool,
        width=dag_peptides.width,
        model=normalized_model,
        target_position=normalized_target,
        anchor_letters=dag_peptides.data["anchorAA"].tolist(),
        upstream_offset=dag_peptides.upstream_offset,
        downstream_offset=dag_peptides.downstream_offset,
    )
    if unique_seq:
        matches = list(dict.fromkeys(matches))

    if normalized_test_type == "fisher":
        return BackgroundModel(
            background=[_split_matches(matches)],
            num_subsamples=1,
            test_type=normalized_test_type,
            source_background=normalized_background,
            model=normalized_model,
            target_position=normalized_target,
            unique_seq=unique_seq,
            rand_seed=rand_seed,
            replacement=replacement,
        )

    if len(matches) < n:
        raise ValueError("Too few matches in the background. Please try different parameters.")

    random_state = np.random.RandomState(rand_seed)
    background_samples: list[pd.DataFrame] = []
    matches_array = np.asarray(matches, dtype=object)
    for _ in range(int(num_subsamples)):
        sampled = random_state.choice(matches_array, size=n, replace=replacement)
        background_samples.append(_split_matches(sampled.tolist()))

    return BackgroundModel(
        background=background_samples,
        num_subsamples=int(num_subsamples),
        test_type=normalized_test_type,
        source_background=normalized_background,
        model=normalized_model,
        target_position=normalized_target,
        unique_seq=unique_seq,
        rand_seed=rand_seed,
        replacement=replacement,
    )
