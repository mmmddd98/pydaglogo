from __future__ import annotations

from collections.abc import Iterable
from typing import Literal

import pandas as pd

from .background import build_background_model
from .io import build_windows_from_flanks, build_windows_from_positions, format_sequence
from .models import BackgroundModel, DAUResult, Proteome
from .plotting import DEFAULT_Y_LABEL, plot_logo
from .statistics import test_dau


def _prepare_windows(
    input_data: pd.DataFrame | Iterable[str],
    *,
    input_mode: Literal["windows", "flanks", "positions"],
    window_col: str,
    before_col: str,
    after_col: str,
    left: int,
    right: int,
    position_cols: tuple[str, ...],
    canonical_only: bool,
    drop_invalid: bool,
) -> list[str]:
    """Normalize supported input styles into a plain list of fixed-width windows."""
    if isinstance(input_data, pd.DataFrame):
        if input_mode == "windows":
            if window_col not in input_data.columns:
                raise KeyError(f"Missing window column in frame: {window_col}")
            return input_data[window_col].astype(str).tolist()
        if input_mode == "flanks":
            prepared = build_windows_from_flanks(
                input_data,
                before_col=before_col,
                after_col=after_col,
                left=left,
                right=right,
                output_col=window_col,
                canonical_only=canonical_only,
                drop_invalid=drop_invalid,
            )
            return prepared[window_col].astype(str).tolist()
        prepared = build_windows_from_positions(
            input_data,
            position_cols=position_cols,
            output_col=window_col,
            canonical_only=canonical_only,
            drop_invalid=drop_invalid,
        )
        return prepared[window_col].astype(str).tolist()

    if input_mode != "windows":
        raise TypeError("DataFrame input is required for input_mode='flanks' or input_mode='positions'.")

    return [str(value) for value in input_data]


def _infer_offsets(
    windows: list[str],
    upstream_offset: int | None,
    downstream_offset: int | None,
) -> tuple[int | None, int | None]:
    """Infer anchor offsets from window width when the caller leaves them unset."""
    if upstream_offset is not None or downstream_offset is not None:
        return upstream_offset, downstream_offset
    if not windows:
        return upstream_offset, downstream_offset

    width = len(windows[0])
    if width % 2 == 0:
        # Cleavage-site windows are usually even-width and anchored on the
        # residue just left of center (for example P1 in an 8-mer P4..P4').
        inferred_upstream = max((width // 2) - 1, 0)
        inferred_downstream = width - inferred_upstream - 1
        return inferred_upstream, inferred_downstream

    inferred_upstream = width // 2
    inferred_downstream = width - inferred_upstream - 1
    return inferred_upstream, inferred_downstream


def run_daglogo(
    input_data: pd.DataFrame | Iterable[str],
    *,
    proteome: Proteome | None = None,
    background_model: BackgroundModel | None = None,
    input_mode: Literal["windows", "flanks", "positions"] = "windows",
    window_col: str = "window",
    before_col: str = "five_res_before",
    after_col: str = "five_res_after",
    left: int = 4,
    right: int = 4,
    position_cols: tuple[str, ...] = ("P4", "P3", "P2", "P1", "P1'", "P2'", "P3'", "P4'"),
    canonical_only: bool = True,
    drop_invalid: bool = True,
    upstream_offset: int | None = None,
    downstream_offset: int | None = None,
    map_to_proteome: bool = False,
    background: str = "whole_proteome",
    model: str = "any",
    target_position: str = "any",
    unique_seq: bool = False,
    num_subsamples: int = 300,
    rand_seed: int = 1,
    replacement: bool = False,
    test_type: str = "ztest",
    grouping_scheme: str | None = "no",
    bg_noise: float | None = None,
    p_adjust_method: str = "none",
) -> DAUResult:
    """Run the core dagLogo workflow in one step.

    Parameters
    ----------
    input_data
        Either an iterable of prebuilt windows or a dataframe containing window,
        flank, or position columns.
    proteome
        Proteome used for mapping or on-the-fly background generation.
    background_model
        Optional committed or precomputed background model. When supplied, it is
        used directly instead of generating a new background.
    input_mode
        How to interpret ``input_data``: prebuilt windows, flank columns, or
        per-position columns such as ``P4..P4'``.
    window_col, before_col, after_col, left, right, position_cols
        Input-preparation controls for the supported input modes.
    canonical_only, drop_invalid
        Input-cleaning controls applied before formatting windows.
    upstream_offset, downstream_offset
        Explicit anchor offsets. If omitted, even-width windows default to the
        residue just left of center.
    map_to_proteome
        Whether to map input windows back onto the provided proteome.
    background, model, target_position, unique_seq, num_subsamples, rand_seed, replacement, test_type
        Background-generation controls used only when ``background_model`` is
        not supplied.
    grouping_scheme, bg_noise, p_adjust_method
        Statistical options forwarded into :func:`pydaglogo.test_dau`.

    Returns
    -------
    DAUResult
        Full analysis result for one dagLogo run.

    Notes
    -----
    The most important behavior difference is whether you pass
    ``background_model=``. If you do, the wrapper uses that committed or custom
    background directly. If you do not, it builds a new background from the
    requested proteome/background settings.
    """

    windows = _prepare_windows(
        input_data,
        input_mode=input_mode,
        window_col=window_col,
        before_col=before_col,
        after_col=after_col,
        left=left,
        right=right,
        position_cols=position_cols,
        canonical_only=canonical_only,
        drop_invalid=drop_invalid,
    )
    if not windows:
        raise ValueError("No valid sequence windows were available for analysis.")

    resolved_upstream, resolved_downstream = _infer_offsets(
        windows,
        upstream_offset,
        downstream_offset,
    )
    dag_peptides = format_sequence(
        windows,
        proteome=proteome,
        upstream_offset=resolved_upstream,
        downstream_offset=resolved_downstream,
        map_to_proteome=map_to_proteome,
    )

    resolved_background = background_model
    if resolved_background is None:
        resolved_background = build_background_model(
            dag_peptides,
            proteome=proteome,
            background=background,
            model=model,
            target_position=target_position,
            unique_seq=unique_seq,
            num_subsamples=num_subsamples,
            rand_seed=rand_seed,
            replacement=replacement,
            test_type=test_type,
        )

    return test_dau(
        dag_peptides,
        resolved_background,
        grouping_scheme=grouping_scheme,
        bg_noise=bg_noise,
        p_adjust_method=p_adjust_method,
    )


def make_daglogo(
    input_data: pd.DataFrame | Iterable[str],
    *,
    metric: str = "difference",
    plot_mode: str = "positive_only",
    labels: list[str] | None = None,
    color_scheme: str = "chemistry",
    y_label: str | None = DEFAULT_Y_LABEL,
    pvalue_cutoff: float | None = 0.05,
    plot_p_adjust_method: str = "none",
    use_global_bh: bool = False,
    ax=None,
    analysis_p_adjust_method: str = "none",
    **run_kwargs,
) -> tuple[DAUResult, object, pd.DataFrame]:
    """Run the workflow and immediately create a logo plot.

    Parameters
    ----------
    input_data
        Input windows or dataframe accepted by :func:`pydaglogo.run_daglogo`.
    metric, plot_mode, labels, color_scheme, y_label, pvalue_cutoff, plot_p_adjust_method, use_global_bh, ax
        Plotting controls forwarded to :func:`pydaglogo.plot_logo`.
    analysis_p_adjust_method
        P-value adjustment applied during the statistical analysis step.
    **run_kwargs
        Additional arguments forwarded to :func:`pydaglogo.run_daglogo`.

    Returns
    -------
    tuple
        ``(result, ax, matrix)`` where ``result`` is the full analysis result,
        ``ax`` is the matplotlib axis, and ``matrix`` is the plotted logo
        matrix.
    """

    result = run_daglogo(
        input_data,
        p_adjust_method=analysis_p_adjust_method,
        **run_kwargs,
    )
    axis, matrix = plot_logo(
        result,
        metric=metric,
        plot_mode=plot_mode,
        labels=labels,
        color_scheme=color_scheme,
        y_label=y_label,
        pvalue_cutoff=pvalue_cutoff,
        p_adjust_method=plot_p_adjust_method,
        use_global_bh=use_global_bh,
        ax=ax,
    )
    return result, axis, matrix
