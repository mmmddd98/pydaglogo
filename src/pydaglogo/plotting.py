from __future__ import annotations

import matplotlib.pyplot as plt
import logomaker
import pandas as pd

from .models import DAUResult
from .statistics import to_logo_matrix

DEFAULT_Y_LABEL = "Δ residue frequency"


def plot_logo(
    result_or_matrix: DAUResult | pd.DataFrame,
    metric: str = "difference",
    plot_mode: str = "signed",
    labels: list[str] | None = None,
    color_scheme: str = "chemistry",
    y_label: str | None = DEFAULT_Y_LABEL,
    pvalue_cutoff: float | None = None,
    p_adjust_method: str = "none",
    use_global_bh: bool = False,
    ax=None,
):
    """Plot a logo from a :class:`DAUResult` or a precomputed matrix.

    Parameters
    ----------
    result_or_matrix
        Either a full analysis result or a matrix already arranged for
        plotting.
    metric
        Matrix to derive from a ``DAUResult``. ``difference`` and
        ``statistics`` are inferential; ``motif`` and ``background`` are
        descriptive.
    plot_mode
        Whether to show signed values, only enriched values, or only depleted
        values.
    labels
        Optional x-axis labels matching the number of positions.
    color_scheme
        Logomaker color scheme name.
    y_label
        Optional y-axis label. Pass ``None`` to suppress the label.
    pvalue_cutoff, p_adjust_method, use_global_bh
        Significance masking controls used only when plotting from a
        ``DAUResult``.
    ax
        Existing matplotlib axis to draw into.

    Returns
    -------
    tuple
        ``(ax, matrix)`` where ``ax`` is the matplotlib axis and ``matrix`` is
        the actual plotted logo matrix.
    """

    if isinstance(result_or_matrix, DAUResult):
        matrix = to_logo_matrix(
            result_or_matrix,
            metric=metric,
            plot_mode=plot_mode,
            pvalue_cutoff=pvalue_cutoff,
            p_adjust_method=p_adjust_method,
            use_global_bh=use_global_bh,
        )
    else:
        matrix = result_or_matrix.copy()

    plot_matrix = matrix.copy()
    plot_matrix.index = range(len(plot_matrix.index))

    if ax is None:
        _, ax = plt.subplots(figsize=(6.0, 2.5))

    logomaker.Logo(
        plot_matrix,
        ax=ax,
        color_scheme=color_scheme,
        baseline_width=0.8,
        stack_order="big_on_top",
    )
    ax.axhline(0.0, color="black", linewidth=0.8)
    if labels is not None:
        if len(labels) != len(matrix.index):
            raise ValueError("labels length must match the number of positions.")
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels)
    ax.set_ylabel("" if y_label is None else y_label)
    if plot_mode == "positive_only":
        ax.set_ylim(bottom=0.0)
    return ax, matrix
