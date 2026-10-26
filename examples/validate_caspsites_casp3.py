from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = PACKAGE_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pydaglogo import (  # noqa: E402
    build_background_model_from_samples,
    build_windows_from_positions,
    load_background_samples,
    plot_logo,
    run_daglogo,
    to_logo_matrix,
)

VALIDATION_ROOT = PACKAGE_ROOT / "validation_data" / "caspsites_casp3"
DATASET_PATH = VALIDATION_ROOT / "caspsites_casp3.csv"
WINDOW_COLUMNS = ["P4", "P3", "P2", "P1", "P1'", "P2'", "P3'", "P4'"]
POSITION_LABELS = WINDOW_COLUMNS
SUPPORTED_TEST_TYPES = ("ztest", "fisher")


def reference_root_for_test_type(test_type: str) -> Path:
    if test_type == "ztest":
        return VALIDATION_ROOT / "reference_bundle"
    if test_type == "fisher":
        return VALIDATION_ROOT / "reference_bundle_fisher"
    raise ValueError(f"Unsupported test_type: {test_type}")


def expected_summary_path_for_test_type(test_type: str) -> Path:
    if test_type == "ztest":
        return VALIDATION_ROOT / "expected_parity_summary.tsv"
    if test_type == "fisher":
        return VALIDATION_ROOT / "expected_parity_summary_fisher.tsv"
    raise ValueError(f"Unsupported test_type: {test_type}")


def generated_output_root_for_test_type(test_type: str) -> Path:
    if test_type == "ztest":
        return VALIDATION_ROOT / "generated_outputs_ztest"
    if test_type == "fisher":
        return VALIDATION_ROOT / "generated_outputs_fisher"
    raise ValueError(f"Unsupported test_type: {test_type}")


def background_samples_path_for_test_type(test_type: str) -> Path:
    if test_type == "ztest":
        return reference_root_for_test_type(test_type) / "background_samples.csv"
    if test_type == "fisher":
        return reference_root_for_test_type(test_type) / "background_samples.csv.gz"
    raise ValueError(f"Unsupported test_type: {test_type}")


def build_public_windows(frame: pd.DataFrame) -> pd.DataFrame:
    return build_windows_from_positions(
        frame,
        position_cols=tuple(WINDOW_COLUMNS),
        output_col="window",
        canonical_only=True,
        drop_invalid=True,
    )


def _load_reference_frame(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, index_col=0)


def _safe_correlation(python_values, reference_values) -> float:
    finite_mask = np.isfinite(python_values) & np.isfinite(reference_values)
    python_finite = python_values[finite_mask]
    reference_finite = reference_values[finite_mask]
    if python_finite.size < 2:
        return float("nan")
    if (
        (python_finite == python_finite[0]).all()
        or (reference_finite == reference_finite[0]).all()
    ):
        return 1.0 if (python_finite == reference_finite).all() else float("nan")
    return float(pd.Series(python_finite).corr(pd.Series(reference_finite)))


def summarize_table_comparison(
    name: str,
    python_frame: pd.DataFrame,
    reference_frame: pd.DataFrame,
) -> dict[str, object]:
    python_aligned = python_frame.copy()
    reference_aligned = reference_frame.copy()
    python_aligned.index = python_aligned.index.map(str)
    reference_aligned.index = reference_aligned.index.map(str)
    python_aligned.columns = python_aligned.columns.map(str)
    reference_aligned.columns = reference_aligned.columns.map(str)

    if set(python_aligned.index) == set(reference_aligned.index):
        python_aligned = python_aligned.reindex(reference_aligned.index)
    if set(python_aligned.columns) == set(reference_aligned.columns):
        python_aligned = python_aligned.reindex(columns=reference_aligned.columns)

    shape_match = python_aligned.shape == reference_aligned.shape
    index_match = list(python_aligned.index) == list(reference_aligned.index)
    column_match = list(python_aligned.columns) == list(reference_aligned.columns)
    if not shape_match:
        return {
            "table": name,
            "shape_match": shape_match,
            "index_match": index_match,
            "column_match": column_match,
            "exact_allclose": False,
            "max_abs_error": float("nan"),
            "mean_abs_error": float("nan"),
            "correlation": float("nan"),
            "masking_agreement": float("nan"),
        }

    python_values = python_aligned.to_numpy(dtype=float)
    reference_values = reference_aligned.to_numpy(dtype=float)
    finite_mask = np.isfinite(python_values) & np.isfinite(reference_values)
    abs_error = abs(python_values[finite_mask] - reference_values[finite_mask])
    exact_allclose = bool(
        pd.DataFrame(python_values).equals(pd.DataFrame(reference_values))
        or np.allclose(
            python_values,
            reference_values,
            rtol=1e-9,
            atol=1e-12,
            equal_nan=True,
        )
    )
    masking_agreement = float("nan")
    if name == "logo_matrix":
        masking_agreement = float(((python_values != 0) == (reference_values != 0)).mean())
    return {
        "table": name,
        "shape_match": shape_match,
        "index_match": index_match,
        "column_match": column_match,
        "exact_allclose": exact_allclose,
        "max_abs_error": float(abs_error.max()) if abs_error.size else float("nan"),
        "mean_abs_error": float(abs_error.mean()) if abs_error.size else float("nan"),
        "correlation": _safe_correlation(
            python_values.reshape(-1),
            reference_values.reshape(-1),
        ),
        "masking_agreement": masking_agreement,
    }


def save_logo_figure(
    result,
    output_path: Path,
    *,
    title: str,
    metric: str = "difference",
    plot_mode: str = "positive_only",
    logo_matrix: pd.DataFrame | None = None,
    pvalue_cutoff: float | None = 0.05,
    use_global_bh: bool = True,
) -> None:
    plot_input = logo_matrix if logo_matrix is not None else result
    ax, _ = plot_logo(
        plot_input,
        metric=metric,
        plot_mode=plot_mode,
        labels=POSITION_LABELS,
        color_scheme="chemistry",
        pvalue_cutoff=None if logo_matrix is not None else pvalue_cutoff,
        use_global_bh=False if logo_matrix is not None else use_global_bh,
    )
    ax.set_title(title, fontsize=10, fontweight="bold")
    figure = ax.figure
    figure.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)


def run_validation(
    output_dir: Path | None = None,
    *,
    test_type: str = "ztest",
) -> tuple[pd.DataFrame, dict[str, Path]]:
    if test_type not in SUPPORTED_TEST_TYPES:
        raise ValueError(f"Unsupported test_type: {test_type}")
    reference_root = reference_root_for_test_type(test_type)
    if output_dir is None:
        output_dir = generated_output_root_for_test_type(test_type)
    output_dir.mkdir(parents=True, exist_ok=True)

    input_frame = pd.read_csv(DATASET_PATH)
    prepared = build_public_windows(input_frame)
    background_samples = load_background_samples(background_samples_path_for_test_type(test_type))
    background_model = build_background_model_from_samples(
        background_samples,
        test_type=test_type,
        source_background=f"committed_r_{test_type}_reference",
        model="any",
        target_position="any",
        unique_seq=False,
        rand_seed=404,
        replacement=False,
    )
    result = run_daglogo(
        prepared,
        background_model=background_model,
        input_mode="windows",
        window_col="window",
        map_to_proteome=False,
        grouping_scheme="no",
    )
    logo_matrix = to_logo_matrix(
        result,
        metric="difference",
        plot_mode="positive_only",
        pvalue_cutoff=0.05,
        use_global_bh=True,
    )
    logo_matrix.index = POSITION_LABELS

    table_map = {
        "difference_raw": result.difference,
        "statistics_raw": result.statistics,
        "pvalue_raw": result.pvalue,
        "background_mean": result.background,
        "motif_observed": result.motif,
        "logo_matrix": logo_matrix,
    }
    reference_tables = {
        name: _load_reference_frame(reference_root / f"{name}.csv")
        for name in table_map
    }
    summary = pd.DataFrame(
        [
            summarize_table_comparison(name, table_map[name], reference_tables[name])
            for name in table_map
        ]
    )
    summary_path = output_dir / f"caspsites_casp3_{test_type}_parity_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)

    for name, frame in table_map.items():
        frame.to_csv(output_dir / f"caspsites_casp3_{test_type}_{name}.tsv", sep="\t")

    save_logo_figure(
        result,
        output_dir / f"caspsites_casp3_{test_type}_diff_logo.png",
        title=f"Caspase-3 public validation ({test_type})\nDiff logo (n={len(prepared)})",
        logo_matrix=logo_matrix,
    )
    save_logo_figure(
        result,
        output_dir / f"caspsites_casp3_{test_type}_motif_logo.png",
        title=f"Caspase-3 public validation ({test_type})\nMotif occupancy (n={len(prepared)})",
        metric="motif",
        plot_mode="positive_only",
        pvalue_cutoff=None,
        use_global_bh=False,
    )

    outputs = {
        "summary": summary_path,
        "diff_logo": output_dir / f"caspsites_casp3_{test_type}_diff_logo.png",
        "motif_logo": output_dir / f"caspsites_casp3_{test_type}_motif_logo.png",
    }
    return summary, outputs


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Validate pydaglogo against committed CaspSites Caspase-3 R references."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for generated validation outputs.",
    )
    parser.add_argument(
        "--test-type",
        choices=SUPPORTED_TEST_TYPES,
        default="ztest",
        help="Statistical mode to validate against the committed R reference bundle.",
    )
    args = parser.parse_args()

    summary, outputs = run_validation(args.output_dir, test_type=args.test_type)
    expected = pd.read_csv(expected_summary_path_for_test_type(args.test_type), sep="\t")
    merged = expected.merge(summary[["table", "exact_allclose"]], on="table", how="left")
    exact_match = bool(
        (merged["expected_exact_allclose"].astype(bool) == merged["exact_allclose"].astype(bool)).all()
    )

    print("pydaglogo public Caspase-3 validation")
    print(f"Test type: {args.test_type}")
    print(f"Dataset: {DATASET_PATH}")
    print(f"Reference bundle: {reference_root_for_test_type(args.test_type)}")
    print(f"Saved parity summary: {outputs['summary']}")
    print(f"Saved diff logo: {outputs['diff_logo']}")
    print(f"Saved motif logo: {outputs['motif_logo']}")
    print(summary.to_string(index=False))

    if exact_match:
        print("Validation result: PASS")
        return 0

    print("Validation result: FAIL")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
