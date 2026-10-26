from __future__ import annotations

import io
import os
import subprocess
import sys
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, norm
from scipy.stats.contingency import odds_ratio

from .grouping import apply_grouping_matrix, ordered_index_for_scheme, plot_symbols_for_scheme, resolve_grouping_scheme
from .models import BackgroundModel, DAUResult, DagPeptides


def _resolve_fisher_backend() -> str:
    """Resolve the requested Fisher backend from the environment."""
    backend = os.environ.get("PYDAGLOGO_FISHER_BACKEND", "auto").strip().lower()
    if backend not in {"auto", "scipy", "r"}:
        raise ValueError(
            "PYDAGLOGO_FISHER_BACKEND must be one of: auto, scipy, r."
        )
    return backend


def _candidate_rscript_paths() -> list[Path]:
    """Return plausible local ``Rscript`` paths for the exact Fisher backend."""
    candidates: list[Path] = []
    override = os.environ.get("PYDAGLOGO_RSCRIPT")
    if override:
        candidates.append(Path(override))

    python_executable = Path(sys.executable).resolve()
    base_prefixes = [python_executable.parent]
    if python_executable.parent.parent != python_executable.parent:
        base_prefixes.append(python_executable.parent.parent)
    for base_prefix in base_prefixes:
        candidates.extend(
            [
                base_prefix / "Scripts" / "Rscript.exe",
                base_prefix / "bin" / "Rscript",
            ]
        )

    seen: set[str] = set()
    existing: list[Path] = []
    for candidate in candidates:
        key = str(candidate).lower()
        if key in seen:
            continue
        seen.add(key)
        if candidate.exists():
            existing.append(candidate)
    return existing


def _build_r_subprocess_env(rscript_path: Path) -> dict[str, str]:
    """Build a subprocess environment that can launch the detected R runtime."""
    env = os.environ.copy()
    base_prefix = rscript_path.parent.parent if rscript_path.parent.name.lower() == "scripts" else None

    if "R_HOME" not in env or not env["R_HOME"]:
        if base_prefix is not None:
            r_home_candidate = base_prefix / "Lib" / "R"
            if r_home_candidate.exists():
                env["R_HOME"] = str(r_home_candidate)

    extra_paths: list[str] = []
    if base_prefix is not None:
        extra_paths.extend(
            [
                str(base_prefix / "Library" / "mingw-w64" / "bin"),
                str(base_prefix / "Library" / "usr" / "bin"),
                str(base_prefix / "Library" / "bin"),
                str(base_prefix / "Scripts"),
            ]
        )
    if env.get("R_HOME"):
        r_home = Path(env["R_HOME"])
        extra_paths.extend(
            [
                str(r_home / "bin" / "x64"),
                str(r_home / "bin"),
            ]
        )
    current_path = env.get("PATH", "")
    env["PATH"] = os.pathsep.join([path for path in extra_paths if path] + [current_path])
    return env


def _run_r_fisher_backend(
    exp_freq: pd.DataFrame,
    bg_freq: pd.DataFrame,
    *,
    required: bool,
) -> tuple[np.ndarray, np.ndarray] | None:
    """Execute R ``fisher.test`` over every residue-position contingency table."""
    exp_values = exp_freq.to_numpy(dtype=np.int64)
    bg_values = bg_freq.to_numpy(dtype=np.int64)
    exp_totals = exp_values.sum(axis=0)
    bg_totals = bg_values.sum(axis=0)

    x_values = exp_values.ravel(order="F")
    X_values = (exp_totals[np.newaxis, :] - exp_values).ravel(order="F")
    y_values = bg_values.ravel(order="F")
    Y_values = (bg_totals[np.newaxis, :] - bg_values).ravel(order="F")

    script = f"""
x <- c({",".join(map(str, x_values.tolist()))})
X <- c({",".join(map(str, X_values.tolist()))})
y <- c({",".join(map(str, y_values.tolist()))})
Y <- c({",".join(map(str, Y_values.tolist()))})
out <- mapply(function(x, X, y, Y) {{
  ft <- fisher.test(matrix(c(x, X, y, Y), nrow = 2, byrow = TRUE), alternative = "two.sided")
  c(statistics = unname(ft$estimate), pvalues = ft$p.value)
}}, x, X, y, Y)
write.table(t(out), file = "", sep = "\\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
""".strip()

    last_error: Exception | None = None
    for rscript_path in _candidate_rscript_paths():
        try:
            result = subprocess.run(
                [str(rscript_path), "-e", script],
                capture_output=True,
                text=True,
                env=_build_r_subprocess_env(rscript_path),
                check=True,
            )
            parsed = pd.read_csv(
                io.StringIO(result.stdout),
                sep="\t",
                header=None,
                names=["statistics", "pvalues"],
            )
            statistics = parsed["statistics"].to_numpy(dtype=float).reshape(exp_freq.shape, order="F")
            pvalues = parsed["pvalues"].to_numpy(dtype=float).reshape(exp_freq.shape, order="F")
            return statistics, pvalues
        except Exception as exc:  # pragma: no cover - backend discovery is environment-specific
            last_error = exc

    if required:
        raise RuntimeError("Unable to execute an R Fisher backend for exact odds-ratio estimates.") from last_error
    return None


def _position_labels(upstream_offset: int, downstream_offset: int) -> list[str]:
    """Build human-readable positional labels such as ``AA-4`` and ``AA1``."""
    labels: list[str] = []
    if upstream_offset > 0:
        labels.extend([f"AA{-value}" for value in range(upstream_offset, 0, -1)])
    labels.append("AA0")
    if downstream_offset > 0:
        labels.extend([f"AA{value}" for value in range(1, downstream_offset + 1)])
    return labels


def _count_matrix(matrix: pd.DataFrame, order: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Count and normalize residue frequencies for each aligned position."""
    freq_columns: list[pd.Series] = []
    percent_columns: list[pd.Series] = []
    for column in matrix.columns:
        counts = matrix[column].value_counts().reindex(order, fill_value=0).astype(float)
        total = counts.sum()
        freq_columns.append(counts)
        percent_columns.append(counts / total if total else counts * np.nan)
    freq = pd.concat(freq_columns, axis=1)
    percent = pd.concat(percent_columns, axis=1)
    freq.columns = list(matrix.columns)
    percent.columns = list(matrix.columns)
    freq.index = order
    percent.index = order
    return freq, percent


def _dirichlet_noise(num_rows: int, num_cols: int, random_state: np.random.RandomState) -> np.ndarray:
    """Generate Dirichlet-style noise used for optional ztest background smoothing."""
    samples = random_state.gamma(shape=1.0, scale=1.0, size=(num_rows, num_cols))
    samples /= samples.sum(axis=0, keepdims=True)
    return samples


def _adjust_1d(values: np.ndarray, method: str) -> np.ndarray:
    """Adjust one vector of p-values with Bonferroni or Benjamini-Hochberg."""
    array = values.astype(float, copy=True)
    valid_mask = ~np.isnan(array)
    valid = array[valid_mask]
    if valid.size == 0:
        return array
    if method == "bonferroni":
        adjusted = np.clip(valid * valid.size, 0.0, 1.0)
    elif method in {"BH", "fdr_bh"}:
        order = np.argsort(valid)
        ranked = valid[order]
        adjusted_ranked = ranked * valid.size / np.arange(1, valid.size + 1)
        adjusted_ranked = np.minimum.accumulate(adjusted_ranked[::-1])[::-1]
        adjusted_ranked = np.clip(adjusted_ranked, 0.0, 1.0)
        adjusted = np.empty_like(valid)
        adjusted[order] = adjusted_ranked
    else:
        raise ValueError(f"Unsupported p-value adjustment method: {method}")
    array[valid_mask] = adjusted
    return array


def adjust_pvalues(frame: pd.DataFrame, method: str = "BH", global_adjust: bool = False) -> pd.DataFrame:
    """Adjust p-values within one result matrix.

    Parameters
    ----------
    frame
        P-value matrix to adjust.
    method
        Adjustment method: ``"BH"``, ``"fdr_bh"``, ``"bonferroni"``, or
        ``"none"``.
    global_adjust
        When ``True``, adjust across the full flattened matrix instead of
        column by column.

    Returns
    -------
    pandas.DataFrame
        Adjusted p-value matrix with the same shape and labels.
    """
    if method == "none":
        return frame.copy()
    if global_adjust:
        flattened = _adjust_1d(frame.to_numpy().reshape(-1), method)
        adjusted = flattened.reshape(frame.shape)
        return pd.DataFrame(adjusted, index=frame.index, columns=frame.columns)
    adjusted_columns = [_adjust_1d(frame[column].to_numpy(), method) for column in frame.columns]
    adjusted = np.column_stack(adjusted_columns)
    return pd.DataFrame(adjusted, index=frame.index, columns=frame.columns)


def test_dau(
    dag_peptides: DagPeptides,
    background_model: BackgroundModel,
    grouping_scheme: str | None = "no",
    bg_noise: float | None = None,
    p_adjust_method: str = "none",
) -> DAUResult:
    """Run differential amino-acid usage testing against a background model.

    Parameters
    ----------
    dag_peptides
        Formatted input windows.
    background_model
        Background samples or a single Fisher background matrix.
    grouping_scheme
        Residue grouping to apply before counting frequencies.
    bg_noise
        Optional Dirichlet-style smoothing applied only to ``ztest``
        background frequencies.
    p_adjust_method
        Column-wise p-value adjustment method applied to the returned result.

    Returns
    -------
    DAUResult
        Difference, statistic, p-value, motif, and background matrices for one
        analysis run.

    Raises
    ------
    ValueError
        If input and background widths do not match or if Fisher receives more
        than one background matrix.

    Notes
    -----
    Fisher exact testing uses an exact R-backed odds-ratio path when available.
    If no working ``Rscript`` is found and the backend is not forced to ``r``,
    the implementation falls back to SciPy's conditional odds ratio.
    """

    scheme = resolve_grouping_scheme(grouping_scheme)
    if bg_noise is not None and not 0.0 <= bg_noise <= 1.0:
        raise ValueError("bg_noise must be between 0 and 1.")

    exp = apply_grouping_matrix(dag_peptides.peptides, scheme)
    bg = [apply_grouping_matrix(sample, scheme) for sample in background_model.background]

    if exp.shape[1] != bg[0].shape[1]:
        raise ValueError("The length of background is different from inputs.")

    order = ordered_index_for_scheme(scheme)
    exp_freq, exp_percent = _count_matrix(exp, order)
    raw_columns = list(exp.columns)
    labels = _position_labels(dag_peptides.upstream_offset, dag_peptides.downstream_offset)
    exp_freq.columns = labels

    if background_model.test_type == "ztest":
        bg_counts = [_count_matrix(sample, order)[1] for sample in bg]
        bg_stack = np.stack([frame.to_numpy() for frame in bg_counts], axis=2)
        if bg_noise is not None:
            random_state = np.random.RandomState(background_model.rand_seed)
            noisy_slices: list[np.ndarray] = []
            for pos in range(bg_stack.shape[1]):
                current = bg_stack[:, pos, :]
                noise = _dirichlet_noise(current.shape[0], current.shape[1], random_state)
                noisy_slices.append(((1 - bg_noise) * current) + (bg_noise * noise))
            bg_stack = np.stack(noisy_slices, axis=1)
        mu = np.nanmean(bg_stack, axis=2)
        difference = exp_percent.to_numpy() - mu
        std = np.sqrt(mu * (1.0 - mu) / background_model.num_subsamples)
        with np.errstate(divide="ignore", invalid="ignore"):
            statistics = difference / std
        pvalue = 2.0 * norm.sf(np.abs(statistics))
        background = pd.DataFrame(mu, index=order, columns=raw_columns)
        difference_frame = pd.DataFrame(difference, index=order, columns=labels)
        statistics_frame = pd.DataFrame(statistics, index=order, columns=labels)
        pvalue_frame = pd.DataFrame(pvalue, index=order, columns=labels)
    else:
        if len(bg) != 1:
            raise ValueError("Fisher requires exactly one background matrix.")
        bg_freq, bg_percent = _count_matrix(bg[0], order)
        bg_freq.columns = labels
        difference_frame = exp_percent - bg_percent
        odds = np.zeros_like(difference_frame.to_numpy(), dtype=float)
        pvalues = np.zeros_like(difference_frame.to_numpy(), dtype=float)
        exp_values = exp_freq.to_numpy()
        bg_values = bg_freq.to_numpy()
        exp_totals = exp_values.sum(axis=0)
        bg_totals = bg_values.sum(axis=0)
        fisher_backend = _resolve_fisher_backend()
        if fisher_backend != "scipy":
            r_backend_result = _run_r_fisher_backend(
                exp_freq,
                bg_freq,
                required=fisher_backend == "r",
            )
            if r_backend_result is not None:
                odds, pvalues = r_backend_result
            else:
                fisher_backend = "scipy"

        if fisher_backend == "scipy":
            for row_idx in range(exp_values.shape[0]):
                for col_idx in range(exp_values.shape[1]):
                    x = exp_values[row_idx, col_idx]
                    X = exp_totals[col_idx] - x
                    y = bg_values[row_idx, col_idx]
                    Y = bg_totals[col_idx] - y
                    table = np.array(
                        [[int(x), int(X)], [int(y), int(Y)]],
                        dtype=np.int64,
                    )
                    stat = odds_ratio(table, kind="conditional").statistic
                    _, pval = fisher_exact(table, alternative="two-sided")
                    odds[row_idx, col_idx] = stat
                    pvalues[row_idx, col_idx] = pval
        difference_frame.columns = labels
        background = bg_percent
        statistics_frame = pd.DataFrame(odds, index=order, columns=labels)
        pvalue_frame = pd.DataFrame(pvalues, index=order, columns=labels)

    if p_adjust_method != "none":
        pvalue_frame = adjust_pvalues(pvalue_frame, method=p_adjust_method, global_adjust=False)

    return DAUResult(
        difference=difference_frame,
        statistics=statistics_frame,
        pvalue=pvalue_frame,
        background=background,
        motif=exp_percent,
        test_type=background_model.test_type,
        upstream_offset=dag_peptides.upstream_offset,
        downstream_offset=dag_peptides.downstream_offset,
        grouping_scheme=scheme.name,
        plot_symbols=plot_symbols_for_scheme(scheme),
    )


test_dau.__test__ = False


def to_logo_matrix(
    result: DAUResult,
    metric: Literal["difference", "statistics", "motif", "background"] = "difference",
    plot_mode: Literal["signed", "positive_only", "negative_only"] = "signed",
    pvalue_cutoff: float | None = None,
    p_adjust_method: str = "none",
    use_global_bh: bool = False,
) -> pd.DataFrame:
    """Convert a result matrix into the format expected by logo plotting.

    Parameters
    ----------
    result
        Differential amino-acid usage result.
    metric
        Matrix to transform. ``difference`` and ``statistics`` are inferential
        metrics that can be significance-masked. ``motif`` and ``background``
        are descriptive frequency matrices.
    plot_mode
        Whether to keep signed values, only positive values, or only negative
        values.
    pvalue_cutoff, p_adjust_method, use_global_bh
        Significance masking controls applied only for inferential metrics.

    Returns
    -------
    pandas.DataFrame
        Position-by-symbol matrix ready for :func:`pydaglogo.plot_logo`.
    """

    values = getattr(result, metric).copy()
    if pvalue_cutoff is not None and metric in {"difference", "statistics"}:
        pvalues = result.pvalue.copy()
        if use_global_bh:
            pvalues = adjust_pvalues(pvalues, method="BH", global_adjust=True)
        elif p_adjust_method != "none":
            pvalues = adjust_pvalues(pvalues, method=p_adjust_method, global_adjust=False)
        values = values.where(pvalues < pvalue_cutoff, 0.0)

    if plot_mode == "positive_only":
        values = values.clip(lower=0.0)
    elif plot_mode == "negative_only":
        values = values.clip(upper=0.0)

    matrix = values.T
    matrix.columns = [result.plot_symbols.get(column, column) for column in matrix.columns]
    return matrix
