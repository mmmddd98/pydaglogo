"""Minimal one-shot Fisher demo using flank columns."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = PACKAGE_ROOT / "src"
# Allow direct execution from a repo checkout without installing the package first.
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pydaglogo import (  # noqa: E402
    build_background_model_from_samples,
    load_background_samples,
    make_daglogo,
    plot_logo,
)

VALIDATION_ROOT = PACKAGE_ROOT / "validation_data" / "demo_flanks"
DATASET_PATH = VALIDATION_ROOT / "demo_flanks.csv"
BACKGROUND_SAMPLES_PATH = VALIDATION_ROOT / "reference_bundle_fisher" / "background_samples.csv"
OUTPUT_DIR = VALIDATION_ROOT / "generated_outputs_demo_oneshot_fisher"
POSITION_LABELS = ["P4", "P3", "P2", "P1", "P1'", "P2'", "P3'", "P4'"]


def save_figure(ax, output_path: Path, title: str) -> None:
    ax.set_title(title, fontsize=10, fontweight="bold")
    figure = ax.figure
    figure.tight_layout()
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)


def main() -> int:
    if not DATASET_PATH.exists():
        print(
            f"Error: data file not found:\n  {DATASET_PATH}\n\n"
            "These example scripts require a local clone of the pydaglogo repository.\n"
            "Clone the repo first:\n"
            "  git clone https://github.com/mmmddd98/pydaglogo.git\n"
            "  cd pydaglogo",
            file=sys.stderr,
        )
        return 1
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    frame = pd.read_csv(DATASET_PATH)
    # Reuse the committed demo background so the example stays fast and deterministic.
    background_model = build_background_model_from_samples(
        load_background_samples(BACKGROUND_SAMPLES_PATH),
        test_type="fisher",
        source_background="demo_flanks_fisher_reference",
        model="any",
        target_position="any",
        unique_seq=False,
        rand_seed=404,
        replacement=False,
    )

    # This is the shortest real analysis path for flank-style input.
    result, ax, _ = make_daglogo(
        frame,
        background_model=background_model,
        input_mode="flanks",
        map_to_proteome=False,
        metric="difference",
        plot_mode="positive_only",
        labels=POSITION_LABELS,
        pvalue_cutoff=0.05,
        use_global_bh=True,
    )
    diff_logo_path = OUTPUT_DIR / "demo_flanks_oneshot_fisher_diff_logo.png"
    save_figure(ax, diff_logo_path, "Flank input demo\nOne-shot Fisher")

    motif_ax, _ = plot_logo(
        result,
        metric="motif",
        plot_mode="positive_only",
        labels=POSITION_LABELS,
        pvalue_cutoff=None,
        use_global_bh=False,
    )
    motif_logo_path = OUTPUT_DIR / "demo_flanks_oneshot_fisher_motif_logo.png"
    save_figure(motif_ax, motif_logo_path, "Flank input demo\nOne-shot Fisher motif occupancy")

    motif_columns = dict(zip(POSITION_LABELS, result.motif.columns))
    observed = {
        "P4:D": float(result.motif.loc["D", motif_columns["P4"]]),
        "P3:E": float(result.motif.loc["E", motif_columns["P3"]]),
        "P2:V": float(result.motif.loc["V", motif_columns["P2"]]),
        "P1:D": float(result.motif.loc["D", motif_columns["P1"]]),
    }

    print("pydaglogo one-shot flank Fisher demo")
    print(f"Dataset: {DATASET_PATH}")
    print(f"Background samples: {BACKGROUND_SAMPLES_PATH}")
    print(f"Saved diff logo: {diff_logo_path}")
    print(f"Saved motif logo: {motif_logo_path}")
    print(
        "Observed motif occupancy: "
        + ", ".join(f"{key}={value:.3f}" for key, value in observed.items())
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
