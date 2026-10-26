# pydaglogo

`pydaglogo` is a Python implementation of the core `dagLogo` workflow for
differential amino-acid usage analysis and sequence logo plotting.

It is built for cleavage-site style inputs, supports both `ztest` and
Fisher exact testing, and now has a short high-level API for the common case:
prepare input windows, run the analysis, and plot the logo.

`pydaglogo` focuses on background-aware differential amino-acid logos for
proteomics workflows, with exact validation against R `dagLogo`.

## At a Glance

| ztest | Fisher |
| --- | --- |
| ![ztest diff logo](assets/hero_ztest_diff_logo.png) | ![fisher diff logo](assets/hero_fisher_diff_logo.png) |

- exact public R/Python validation for `ztest` and `fisher`
- local-first workflow: input sequences and backgrounds are processed locally in Python (and local R when used for exact Fisher parity), without requiring submission to a remote logo server
- position-based and flank-based input helpers
- one-step convenience wrappers via `run_daglogo(...)` and `make_daglogo(...)`
- committed public benchmark backgrounds so demos and validation stay fast

## Installation

`pydaglogo` requires Python 3.10 or newer.

Install the latest version directly from GitHub with:

```bash
pip install "git+https://github.com/mmmddd98/pydaglogo.git"
```

The package is currently validated in GitHub Actions with Python 3.11, and the
local verification workflow in this repository was also run with Python 3.13.9.

## Quickest Example

If you already have cleavage-site residues in `P4..P4'` columns:

```python
import pandas as pd

from pydaglogo import make_daglogo, prepare_proteome

frame = pd.read_csv("my_cleavage_sites.csv")
proteome = prepare_proteome("/path/to/proteome.fasta")

result, ax, matrix = make_daglogo(
    frame,
    proteome=proteome,
    input_mode="positions",
    test_type="ztest",
    labels=["P4", "P3", "P2", "P1", "P1'", "P2'", "P3'", "P4'"],
)
```

This assumes `frame` contains one residue column per position from `P4` through
`P4'`. Switch to `test_type="fisher"` if you want Fisher testing instead.

If your table has flank columns instead, switch to `input_mode="flanks"`.

If you want the analysis result without plotting, use `run_daglogo(...)`
instead. For even-width cleavage windows, these workflow helpers default to the
residue just left of center as the anchor, which matches windows such as
`P4..P1 + P1'..P4'`.

Logo plots use `Δ residue frequency` as the default y-axis label. You can
customize it with `y_label="..."` or suppress it with `y_label=None`.

## Input Styles

> **Note:** The example scripts below use committed reference data from
> `validation_data/` inside the repository. They require a local clone of the
> repo — they cannot be run after a plain `pip install` without the repo
> present. Clone first: `git clone https://github.com/mmmddd98/pydaglogo.git`

### Position-based one-shot demos

- `python examples/demo_caspsites_casp3_oneshot_ztest.py`
- `python examples/demo_caspsites_casp3_oneshot_fisher.py`

These use the public CaspSites Caspase-3 dataset with committed background
bundles and produce one diff logo plus one motif occupancy logo.

### Flank-based one-shot demos

- `python examples/demo_flanks_oneshot_ztest.py`
- `python examples/demo_flanks_oneshot_fisher.py`

These use a tiny committed public flank fixture under
[`validation_data/demo_flanks`](../validation_data/demo_flanks) and are intended as
copy/paste integration examples for pipelines that already work with
`five_res_before` / `five_res_after`.

## Why Trust It?

`pydaglogo` is validated exactly against committed R `dagLogo` reference bundles
for the shipped public benchmark:

- foreground motif counts
- background mean matrix
- differential amino-acid usage matrix
- test statistics
- p-values
- final positive-only diff-logo matrix

`pydaglogo` does not require a remote server for its core analysis workflow.
Input tables, sequence windows, background models, and logo outputs are handled
locally in Python, with local `Rscript` used only when requested for exact
Fisher parity.

The official automated validation path is a Windows GitHub Actions run using
Python 3.11 and an exact R-backed Fisher backend.

For Fisher, `pydaglogo` uses an exact R-backed odds-ratio path when a working
`Rscript` is available. If no working `Rscript` is detected, it falls back to
SciPy's conditional odds ratio. The SciPy fallback remains supported for local
use, but the official exact-parity workflow in this repository uses `Rscript`
via `PYDAGLOGO_RSCRIPT` together with `PYDAGLOGO_FISHER_BACKEND=r`.
On the shipped public benchmark, that SciPy fallback is already extremely close
to the R reference, with only small residual differences in raw odds-ratio
estimates.

## Public Benchmark

The repository includes a public, reproducible benchmark based on the CaspSites
Caspase-3 dataset:

- dataset: [validation_data/caspsites_casp3/caspsites_casp3.csv](../validation_data/caspsites_casp3/caspsites_casp3.csv)
- committed R ztest reference bundle: [validation_data/caspsites_casp3/reference_bundle](../validation_data/caspsites_casp3/reference_bundle)
- committed R Fisher reference bundle: [validation_data/caspsites_casp3/reference_bundle_fisher](../validation_data/caspsites_casp3/reference_bundle_fisher)

Most users only need the Python validator:

```bash
python examples/validate_caspsites_casp3.py
python examples/validate_caspsites_casp3.py --test-type fisher
```

These commands write a parity summary plus one diff logo and one motif-frequency
QC logo.

Maintainers can regenerate the committed public R references with:

```bash
Rscript examples/generate_caspsites_casp3_references.R
Rscript examples/generate_caspsites_casp3_references.R fisher
```

## Current Scope

- amino-acid windows only
- FASTA-based proteome input
- `wholeProteome`, `inputSet`, and `nonInputSet` background modes
- `any` and `anchored` background extraction models
- Fisher exact test and `ztest` support
- matrix export for plotting
- Logomaker-based plotting with signed and positive-only modes
- reusable preprocessing via `build_windows_from_flanks(...)` and
  `build_windows_from_positions(...)`
- one-step workflow helpers via `run_daglogo(...)` and `make_daglogo(...)`

## Acknowledgements

`pydaglogo` is a Python reimplementation and adaptation of the R/Bioconductor
package [dagLogo](https://bioconductor.org/packages/dagLogo/), which is
licensed under GPL (>= 2).

If you use `pydaglogo` in your research, please cite the original dagLogo paper:

> Ou J, Liu H, Nirala KN, Stukalov A, Acharya U, Green RM, Zhu LJ (2020).
> "dagLogo: An R/Bioconductor package for identifying and visualizing
> differential amino acid group usage in proteomics data."
> *PLoS ONE*, **15**(11), e0242030.
> <https://doi.org/10.1371/journal.pone.0242030>

The original dagLogo authors are: Jianhong Ou, Haibo Liu, Alexey Stukalov,
Niraj Nirala, Usha Acharya, and Lihua Julie Zhu.
Maintainer: Jianhong Ou (<jou@morgridge.org>).
