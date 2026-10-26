# Public Caspase-3 Validation Data

This folder contains a public, repo-safe validation benchmark for `pydaglogo`.

Files:
- `caspsites_casp3.csv`: public CaspSites Caspase-3 cleavage dataset used as the foreground set
- `reference_bundle/`: committed R `ztest` `dagLogo` reference outputs generated from the same CSV
- `reference_bundle_fisher/`: Fisher `dagLogo` reference outputs generated from the same CSV
- `expected_parity_summary.tsv`: expected exact-table parity flags for the Python validator
- `expected_parity_summary_fisher.tsv`: expected exact-table parity flags for the Fisher Python validator

`reference_bundle/background_samples.csv` contains the exact committed R
`ztest` background samples that the Python public validator reuses for exact parity
checks. That keeps the public benchmark lightweight and avoids requiring a full
UniProt proteome background build just to validate the package.

`reference_bundle_fisher/background_samples.csv.gz` contains the exact deterministic
Fisher background matrix that the Python public validator reuses for exact
parity checks.
It is stored as ``.csv.gz`` to keep the committed public bundle below GitHub's
file-size limit while remaining directly readable by the Python loader.

`reference_bundle/background_provenance.txt` and
`reference_bundle_fisher/background_provenance.txt` record the proteome source
used to generate those committed background matrices. The current public bundles
were built from:
- `human_uniprot_2026_01/UniProt_Proteome_Human_2026_01.fasta`

This benchmark is intended for public reproducibility and trust-building. It is
the official public validation surface used by the repository tests and CI.
