# Changelog

All notable changes to this project will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).
This project uses [Semantic Versioning](https://semver.org/).

## [1.0.0] - 2026-05-05

### Added

- First public release of `pydaglogo`.
- Core differential amino-acid usage (DAU) statistics: `ztest` and Fisher exact
  test, with an exact R-backed Fisher path and a SciPy fallback.
- Three input modes: prebuilt windows (`windows`), flank columns (`flanks`), and
  per-position residue columns (`positions`).
- One-step workflow helpers `run_daglogo()` and `make_daglogo()`.
- Background generation from a FASTA proteome (`wholeProteome`, `inputSet`,
  `nonInputSet`) with `any` and `anchored` extraction models.
- Committed background loading via `load_background_samples()` and
  `build_background_model_from_samples()`.
- Logomaker-based logo plotting (`plot_logo()`) with signed, positive-only, and
  negative-only modes.
- P-value masking with optional Bonferroni or Benjamini-Hochberg adjustment.
- Grouping schemes for collapsing the amino-acid alphabet (`no`, `charge_group`).
- Reusable preprocessing helpers `build_windows_from_flanks()` and
  `build_windows_from_positions()`.
- Public benchmark: CaspSites Caspase-3 dataset with committed R reference
  bundles for both `ztest` and Fisher, validated via `validate_caspsites_casp3.py`.
- GitHub Actions CI workflow (`pydaglogo-validation.yml`) running on
  `windows-latest` with Python 3.11 and R for exact Fisher parity.
- Four copy-paste demo scripts covering position-based and flank-based inputs for
  both test types.
- Sphinx documentation skeleton under `docs/`.
