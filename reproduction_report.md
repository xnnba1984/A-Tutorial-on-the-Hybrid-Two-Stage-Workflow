# Reproduction verification

Verification date: 2026-10-26.

The complete analyses were rerun in a fresh directory after the simulation scripts and ACTG case study scripts were placed in separate subfolders. There were no prior result files in that directory. ACTG 175 data were reconstructed from the pinned public CRAN source. The directory changes did not alter statistical methods, parameters, seeds, or reference results.

## Numerical agreement

All 21 numerical output checks passed. Twenty reference CSV files were compared cell by cell, covering 320,151 data cells. Every compared scientific value agrees exactly at the precision stored in the CSV files. The largest absolute difference was 0. The participant-level ACTG evaluation file matched its original SHA-256 fingerprint, row count, and column count exactly. That file is generated locally and is not redistributed.

The checks covered all 3,000 main-simulation replicates, 10,800 learner records from 5,400 sensitivity replicates, 1,000 ACTG bootstrap replicates, 200 validation splits, truth summaries, descriptive cross-fitted outputs, and both supplementary illustrative replicates.

The sole schema normalization maps the documented legacy Stage 1 rejection column names. No reference result was overwritten. See verification/numerical_comparison.json for individual checks and verification/run_manifest.json for the new execution record.

## Environment and runtime

R 4.5.2, Python 3.12.14, and Pillow 12.3.0 were used on arm64 macOS. The run reused the verified R library restored for the 2026-10-06 release, together with the R framework library. The global user package library was excluded. All 73 locked R package versions matched. Eight workers were used where specified. This check validates the reorganized workflow in that existing environment; it does not represent a new dependency installation.

The previous release identified numerical differences between two macOS builds of grf 2.4.0 despite identical package version strings. The documented installer builds the unchanged CRAN source with -ffp-contract=off, which disables fused floating-point multiply-add contraction. The source checksum and build evidence from that installation are retained in verification/grf_build.json and verification/package_builds.csv. No statistical method, seed, reference value, or comparison tolerance was changed.

Stage times below exclude package installation and preflight. Stages ran in the documented order. Timings vary with hardware and worker settings.

| Stage | Seconds |
| --- | ---: |
| environment | 1.235 |
| data | 1.108 |
| main | 663.461 |
| truth | 0.687 |
| sensitivity | 448.069 |
| actg_crossfit | 8.224 |
| actg_validation | 53.487 |
| tables | 2.28 |
| figure2 | 1.941 |
| supplementary_figures | 3.547 |
| figure1 | 3.145 |

## Figure checks

All five regenerated figures were compared with the accepted manuscript's original figure files. Their dimensions and RGB pixels were identical. These checks cover Figures 1, 2, and 3 and Supplementary Figures S1 and S2. See verification/figure_comparison.json.

## Scope

The reorganized package also passed a separate 11-stage quick execution check and 74 tooling regression tests. Quick outputs are not the manuscript results. Full-run verification retains the original absolute and relative tolerances of 1e-10, although no nonzero difference was observed. Public-data reconstruction matched the original analysis input byte for byte.

The verifier requires a completed same-run record, unchanged source and output hashes, and all expected artifacts, including five readable, nonblank figures. The pixel comparison above is a separate release check. Local filesystem paths are redacted from published records; their checksums are retained. The 2026-10-06 release remains available at commit 3059b3ad1eff3f4dc66bede244a34f10a2f7ac87.

Linux and Windows were not executed for this release. Graphics can depend on system fonts and rendering libraries.
