# Reproduction verification

Verification date: 2026-10-06.

The full analyses were run in a separate directory with no prior result files. ACTG 175 data were reconstructed from the pinned public CRAN source rather than copied from an author-only input directory.

## Numerical agreement

All 21 numerical output checks passed. Twenty reference CSV files were compared cell by cell. Every compared scientific value agrees exactly at the precision stored in the CSV files; the largest absolute difference was 0. The participant-level ACTG evaluation file matched its original SHA-256 fingerprint, row count, and column count exactly. That file is generated locally and is not redistributed in this repository. The checks included all 3,000 main-simulation replicates, 10,800 learner records from 5,400 sensitivity replicates, 1,000 ACTG bootstrap replicates, 200 validation splits, truth summaries, descriptive cross-fitted outputs, and both supplementary illustrative replicates.

The sole schema normalization maps the documented legacy Stage 1 rejection column names. No reference result was overwritten by a new result. See `verification/numerical_comparison.json` for row counts, hashes, and individual comparisons.

## Environment and runtime

R 4.5.2; Python 3.12.14; Pillow 12.3.0; arm64 macOS on an Apple M4 Max with 16 CPU cores and 128 GiB RAM. Eight workers were used where specified. The final run used a restored repository-local R library plus the R framework library, excluding the global user package library. All 73 locked analytical package versions matched. Detailed R versions are in `renv.lock` and `verification/session_info.txt`.

An initial restore with a different binary build of grf 2.4.0 did not reproduce every numerical output, despite matching package version strings. Single-package crossover checks isolated the difference to the compiled grf library. The final installer builds the unchanged CRAN grf 2.4.0 source with `-ffp-contract=off`. The source archive is checked by SHA-256, and the build uses the selected R 4.5.2 runtime. No statistical method, seed, reference value, or verification tolerance was changed to resolve this issue. Compiler and build evidence are recorded in `verification/grf_build.json` and `verification/package_builds.csv`.

Measured stage times exclude package installation and runner preflight. Stages ran in the documented order. Timings will vary with hardware and worker settings.

| Stage | Seconds |
| --- | ---: |
| environment | 1.064 |
| data | 1.256 |
| main | 599.219 |
| truth | 0.582 |
| sensitivity | 405.529 |
| actg_crossfit | 3.121 |
| actg_validation | 43.021 |
| tables | 0.036 |
| figure2 | 0.699 |
| supplementary_figures | 2.161 |
| figure1 | 0.236 |

## Figure checks

| Figure | Matches embedded manuscript pixels |
| --- | --- |
| Figure 1 | Yes |
| Figure 2 | Yes |
| Figure 3 | Yes |
| Supplementary Figure S1 | Yes |
| Supplementary Figure S2 | Yes |

## Scope

Quick tests check execution only. They are stored separately from this full run. The default comparison tolerance is 1e-10, but no nonzero numerical difference was observed here. The public data reconstruction also passed byte-level comparison against the original input. Runtime worker-count metadata is checked against the run settings separately from scientific values.

The verifier additionally requires a completed same-run execution record, unchanged source and output hashes, all required artifacts, and five decodable nonblank figure PNGs. The manuscript pixel comparison above is a separate author-side check. Local filesystem paths are redacted from the published execution record; checksums are retained.

These checks establish reproduction in the recorded environment. They are not the journal's independent reproducibility assessment. Linux and Windows were not executed for this release, and graphics can depend on system fonts and rendering libraries.
