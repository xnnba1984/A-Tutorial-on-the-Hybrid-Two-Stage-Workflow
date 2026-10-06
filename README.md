# From Heterogeneous Treatment Effects to Treatment Policies

This repository accompanies the manuscript **From Heterogeneous Treatment Effects to Treatment Policies: A Decision Focused Framework for Treatment Personalization in Drug Development**. The paper separates a treatment effect claim from a policy claim. Evidence of heterogeneous treatment effects does not, by itself, establish that a learned treatment policy improves outcomes. The simulations and ACTG 175 example evaluate these questions separately.

The ACTG 175 example uses the publicly available trial dataset distributed in the R package [speff2trial](https://CRAN.R-project.org/package=speff2trial), version 1.0.5. The analysis compares zidovudine monotherapy with the two combination arms and evaluates event-free survival at 96 weeks. Data preparation retains all original records. The analysis scripts make the treatment-arm selection and construct the endpoint.

## Requirements

The verified environment uses R 4.5.2, Python 3.12.14, and Pillow 12.3.0 on arm64 macOS. All 73 R package versions and their dependencies are recorded in `renv.lock` and `environment_versions.csv`. Key packages include grf 2.4.0, ggplot2 4.0.0, dplyr 1.1.4, survival 3.8-3, sandwich 3.1-1, and car 3.1-3. The scripts do not require private data or access to the authors' project directory.

From the repository directory, install Python dependencies and restore the recorded R packages:

```sh
python3 -m pip install -r requirements.txt
Rscript install_dependencies.R
```

The reference installer supports R 4.5.2 on arm64 macOS and restores packages into `renv/library`; the runner uses that directory automatically. It does not install R itself. It requires the Apple command-line development tools (`clang++`, `make`, and `shasum`). Other source dependencies may require Fortran and development libraries. Installation time is separate from analysis time.

Package version numbers alone did not ensure exact agreement in the release checks. Two macOS builds of grf 2.4.0 produced different forest predictions. The installer therefore builds grf from a checksum-verified CRAN source archive with `-ffp-contract=off`, which disables fused floating-point multiply-add contraction. It preserves any previous repository-local grf installation, leaves global libraries unchanged, and records the source checksum, compiler, build information, and shared-library checksum under `renv/library/.grf-build/`. The analysis methods, seeds, and reference results are unchanged. See `reproduction_report.md` for the completed checks and their scope.

An existing environment can be used without running the installer, but matching version strings alone is insufficient evidence of numerical agreement. The environment check records actual versions and warns about differences. Every new run must be checked against the protected reference outputs. Linux and Windows have not been validated for this release; the exact-reference installer deliberately stops on those platforms rather than reporting an untested environment as restored.

For exploratory execution on another platform, the versions can be restored from R with renv installed:

```r
dir.create("renv/library", recursive = TRUE, showWarnings = FALSE)
renv::restore(project = ".", library = "renv/library",
              lockfile = "renv.lock", prompt = FALSE)
```

This alternative does not reproduce the validated grf build recipe. Run the full output comparison and investigate any differences; do not replace reference values or interpret a failed check as successful reproduction.

Figure 1 uses Times New Roman. On a system without that font, install it or use Liberation Serif. For another installation path, set `FIGURE_FONT_REGULAR` and `FIGURE_FONT_BOLD` to the two font files. The script reports the selected font and stops if none is found. Font substitution may change the conceptual figure's appearance, but not the analyses. Other figures use the R graphics device and may also differ slightly across systems.

## Run the complete analysis

Start in a fresh checkout or extracted directory with no `result/` folder:

```sh
python3 reproduce.py --mode full
python3 verify_results.py
```

The first command checks the environment, prepares ACTG data, runs all reported numerical analyses, creates numerical table outputs, and generates the five figures. It writes new outputs to `result/`. The second command checks 21 new numerical output files against the protected manuscript references in `reference_results/` and writes `result/verification.json`. Twenty files are compared cell by cell, with default absolute and relative tolerances of 1e-10. The participant-level ACTG evaluation file is checked against its SHA-256 fingerprint without redistributing its original covariates and outcomes. That comparison requires byte identity. Missing outputs or unexplained differences cause a nonzero exit. A different software or operating-system environment can produce differences; inspect them rather than replacing the references or treating the run as verified.

Verification also requires a completed run record, unchanged source and output hashes, and all expected artifacts, including five readable, nonblank figure PNGs. This checks execution and artifact completeness; pixel agreement with the manuscript figures is documented separately in `reproduction_report.md`. A failed or partial run cannot retain a successful verification status. Keep each `result/` directory with the intact repository that generated it.

Run a small execution check in a **separate fresh copy**:

```sh
python3 reproduce.py --mode quick
```

Quick mode reduces replication counts, tree counts, and selected sample sizes. Its outputs are not the manuscript results. The runner refuses to mix quick and full outputs in the same `result/` directory. The descriptive ACTG cross-fitted analysis retains its original settings in both modes.

Use `--stages` only to resume or rerun stages in an existing, recorded run. The runner validates its inputs and retains its mode, settings, and environment:

```sh
python3 reproduce.py --mode full --stages figure2
python3 reproduce.py --mode full --stages tables
```

Rerunning a stage invalidates its dependent stages. The runner lists any stages still needed before verification can pass. A fresh run refuses a nonempty `result/` directory. Use another checkout to change settings or preserve a prior run. A lock prevents two processes from writing the same result tree. After an abrupt process termination, confirm that no child analysis process remains before addressing a leftover lock.

Stages, in dependency order, are `environment`, `data`, `main`, `truth`, `sensitivity`, `actg_crossfit`, `actg_validation`, `tables`, `figure2`, `supplementary_figures`, and `figure1`. The runner sets the reported settings explicitly and disables personal R startup profiles. `R_BIN` may specify a nondefault Rscript executable or a wrapper that accepts the R script as its first argument. `--workers 1` runs sensitivity and ACTG validation with one worker; the reference configuration uses eight. Worker-count metadata is checked against the recorded run settings separately from scientific results. Windows uses sequential execution where forked workers are unavailable; Windows and Linux were not tested for this release.

For data preparation without a network connection, use the pinned archive described in `DATA_PREPARATION.md`:

```sh
python3 reproduce.py --mode full --data-source-tarball /path/to/speff2trial_1.0.5.tar.gz
```

## Data preparation

`prepare_data.R` downloads the pinned `speff2trial` 1.0.5 source archive over HTTPS, verifies its checksum, and extracts the documented ACTG 175 text data without installing or executing the package. It writes `data/ACTG175.csv` with 2,139 rows and 27 columns. No sorting, filtering, recoding, or imputation is applied at this step. Archive, data, and output fingerprints are checked; an existing mismatching CSV is never overwritten.

```sh
Rscript prepare_data.R
```

The generated CSV reproduces the original analysis input byte for byte:

```text
SHA-256: debab27bbc6f741724b79e49a8e13e1feb8db57850b17010b0635202bfc7141a
```

See `DATA_PREPARATION.md` for source URLs, the schema, and offline preparation from a previously downloaded archive. Raw participant data are not redistributed here. Upstream data retain their original terms; the repository's MIT code license does not change those terms.

## Settings and random seeds

| Analysis | Reported configuration | Random seed |
| --- | --- | --- |
| Main simulation | Six scenarios; 500 replicates per scenario; 2,000 participants; 50/25/25 training/tuning/test split; 2,000-tree causal forest | Sequential RNG stream initialized at 2025 |
| Population truth | 1,000,000 simulated individuals per scenario | 20260517 |
| Information sensitivity | Three scenarios by three sample sizes (500, 1,000, 2,000) by two dimensions (3, 20); 300 replicates per setting; 1,000 trees; two learners on shared splits | `202608230 + 100000 * setting_id + replicate`; split/forest offsets +1/+2 |
| ACTG descriptive cross-fitting | 1,578 eligible patients; five folds; inverse probability of censoring weighting | 175; model seed 175 + fold |
| ACTG Stage 1 inference | 1,000 patient bootstrap replicates, refitting censoring and outcome models | 202608230 |
| ACTG policy validation | 200 treatment-stratified 50/25/25 splits; 1,500 trees; one training model bundle reused for tuning and test | `175000 + split`; forest offset +10000 |
| Supplementary single-replicate figures | 2,000 participants per illustrative replicate | 20260516 and 20260517 |

Stage 2 is evaluated in every simulation replicate, without conditioning on Stage 1 rejection. The main script preserves its sequential RNG stream, including the seed drawn internally by grf. Parallelizing that loop or inserting random draws changes the simulated samples and is not an exact rerun of this configuration.

## Outputs and manuscript correspondence

All paths below are relative to `result/`.

| Manuscript item | Script | Output |
| --- | --- | --- |
| Table 3 and main simulation results | `sim_1.R` | `sim_main_SiM.csv`, `sim_main_SiM_summary.csv` |
| Figure 1 | `make_framework_figure1.py` | `figures/Figure_1_revised.png` |
| Figure 2 | `make_simulation_figure2.R` | `figures/Figure_2_revised.png` |
| Supplementary Table S1 | `sim_truth_summary.R` | `sim_truth_summary_step6.csv` |
| Supplementary Table S2 | `sim_sensitivity.R`, then `prepare_revision_outputs.py` | `information_sensitivity/supplementary_table_s2.csv` and its raw/summary inputs |
| ACTG final Stage 1 inference | `actg_validation_uncertainty.R` | `actg_validation/actg175_ipcw_stage1_global.csv`, `actg175_ipcw_stage1_effect_intervals.csv`, and bootstrap CSVs |
| ACTG repeated split validation | `actg_validation_uncertainty.R` | `actg_validation/actg175_ipcw_validation_splits.csv`, `actg175_ipcw_validation_summary.csv` |
| Figure 3 and descriptive ACTG curves | `real_data.R` | `actg_crossfit/figures/actg175_figure3_ipcw.png` and the corresponding `actg_crossfit/` CSVs |
| Supplementary Figures S1 and S2 | `make_supplementary_single_replicate_figures.R` | `figures/supplementary_figure_s1_constant_benefit_no_hte.png`, `figures/supplementary_figure_s2_strong_qualitative_hte.png`, `supplementary_single_replicate_summary.csv` |
| Software and execution record | `check_environment.R`, `reproduce.py` | `package_versions.csv`, `session_info.txt`, `logs/*.log`, `logs/*.json` |

Tables 1, 2, and 4 are conceptual summaries, not calculated outputs. Their text remains in the manuscript. Numerical table values come from the CSV sources above; the original Word formatting is retained in the submission source files.

The descriptive `actg_crossfit/actg175_ipcw_summary.csv` and `actg175_ipcw_optionC.csv` retain older HC3-based Stage 1 summaries for source continuity. **They are not the final Stage 1 inference reported in the manuscript.** Use the bootstrap files in `actg_validation/`. Cross-fitted policy curves are descriptive; the 200-split analysis evaluates independently selected thresholds on held-out test sets.

Some original references use `proceed`, `proceed_rate`, and `se_proceed_rate`. These mean Stage 1 rejection, not a gate determining whether Stage 2 ran. The verifier maps them to `stage1_reject`, `stage1_rejection_rate`, and `se_stage1_rejection_rate`. Main-simulation `AUQC` and `mean_AUQC` columns contain cAUQC; columns explicitly labeled `raw` contain uncentered AUQC.

## Verification and version control

`reference_results/` contains protected numerical outputs supporting the manuscript, with a SHA-256 manifest. It is separate from newly generated `result/` files. `reproduction_report.md` records completed release verification, actual timings, and platform limitations. Reported results should not be changed merely to make a comparison pass.

All scripts use repository-relative paths. `run_revision_analyses.sh` is a legacy partial runner for sensitivity, ACTG validation, and postprocessing only. Use `reproduce.py` for complete reproduction. The root-level ACTG session files document prior reference runs; each new run records its own environment in `result/`.

## License

Code retains the existing MIT license in `LICENSE`. Data attribution and terms are described in `DATA_PREPARATION.md`.
