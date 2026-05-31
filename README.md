# From Heterogeneous Treatment Effects to Treatment Policies

This repository contains code and reproducibility materials for the manuscript:

**From Heterogeneous Treatment Effects to Treatment Policies: A Decision Focused Framework for Clinical Trial Personalization**

The paper studies a practical problem in clinical trials: evidence that treatment effects vary across patients does not automatically mean that an individualized treatment policy will improve clinical decisions. The manuscript proposes a two-stage framework. Stage 1 evaluates prespecified population level evidence for heterogeneous treatment effects (HTE). Stage 2 evaluates whether estimated conditional average treatment effects (CATEs) are useful for ranking patients, improving policy value, and describing empirical surrogate benefit and harm operating metrics.

The repository reproduces the simulation study and the ACTG 175 case study reported in the paper. The ACTG 175 data are publicly available through the R package `speff2trial`; the analysis here uses a derived local CSV file, `data/ACTG175.csv`, with inverse probability of censoring weighting (IPCW) for the 96-week event-free endpoint.

## Repository Contents

Recommended repository layout:

```text
.
├── README.md
├── data/
│   └── ACTG175.csv
├── result/
│   └── figures/
├── sim_1.R
├── sim_sensitivity.R
├── sim_truth_summary.R
├── make_simulation_figure2.R
├── make_supplementary_single_replicate_figures.R
├── real_data.R
├── actg_validation_uncertainty.R
├── make_framework_figure1.py
├── actg175_ipcw_session_info.txt
└── actg175_ipcw_validation_session_info.txt
```

## Data

The ACTG 175 analysis uses publicly available deidentified clinical trial data. In the manuscript, the data source is described as:

> The ACTG 175 clinical trial data used in this study are publicly available from the R package `speff2trial`.

For the scripts in this repository, place the ACTG 175 CSV file at:

```text
data/ACTG175.csv
```

The real data scripts expect the repository root to contain the `data/` and `result/` folders.

## Software Requirements

The main analyses were run in R. The framework figure was generated in Python.

R packages used across the scripts include:

- `grf`
- `ggplot2`
- `pROC`
- `PRROC`
- `dplyr`
- `data.table`
- `survival`
- `sandwich`
- `lmtest`
- `car`
- `patchwork`
- `scales`
- `gridExtra`

Python packages used:

- `Pillow`

Session information from the final ACTG analyses is saved in:

- `actg175_ipcw_session_info.txt`
- `actg175_ipcw_validation_session_info.txt`

## How to Reproduce the Main Results

Run scripts from the repository root. If any script contains a machine-specific `setwd(...)` line, remove it or replace it with the path to the repository root before running.

### 1. Simulation Study

Main simulation:

```r
source("sim_1.R")
```

Outputs:

```text
result/sim_main_SiM.csv
result/sim_main_SiM_summary.csv
```

Sensitivity simulation:

```r
source("sim_sensitivity.R")
```

Outputs:

```text
result/sim_sensitivity_SiM.csv
result/sim_sensitivity_SiM_summary.csv
```

Large-population truth summary for the simulation scenarios:

```r
source("sim_truth_summary.R")
```

Output:

```text
result/sim_truth_summary_step6.csv
```

Figure 2:

```r
source("make_simulation_figure2.R")
```

Output:

```text
result/figures/simulation_figure2_final.png
```

Supplementary Figures S1 and S2:

```r
source("make_supplementary_single_replicate_figures.R")
```

Outputs:

```text
result/figures/supplementary_figure_s1_constant_benefit_no_hte.png
result/figures/supplementary_figure_s2_strong_qualitative_hte.png
result/supplementary_single_replicate_summary.csv
```

### 2. ACTG 175 Case Study

Main ACTG 175 IPCW analysis:

```r
source("real_data.R")
```

Outputs:

```text
result/actg175_ipcw_summary.csv
result/actg175_ipcw_optionC.csv
result/actg175_ipcw_stepp_cd4.csv
result/actg175_ipcw_stepp_karnofsky.csv
result/actg175_ipcw_uplift_curve.csv
result/actg175_ipcw_policy_curve.csv
result/actg175_ipcw_np_curve.csv
result/actg175_ipcw_eval.csv
result/actg175_ipcw_session_info.txt
result/figures/actg175_figure3_ipcw.png
```

Repeated train/tune/test validation for the ACTG 175 case study:

```r
source("actg_validation_uncertainty.R")
```

Outputs:

```text
result/actg175_ipcw_validation_splits.csv
result/actg175_ipcw_validation_summary.csv
result/actg175_ipcw_stage1_effect_intervals.csv
result/actg175_ipcw_validation_session_info.txt
```

By default, `actg_validation_uncertainty.R` uses 200 repeated splits and 1500 trees. These can be changed with environment variables:

```bash
ACTG_VALIDATION_B=50 ACTG_VALIDATION_TREES=1000 Rscript actg_validation_uncertainty.R
```

### 3. Framework Figure

Figure 1:

```bash
python3 make_framework_figure1.py
```

Output:

```text
result/figures/framework_figure1_final.png
```

## Simulation Scenarios

The final simulation study includes six scenarios:

1. No treatment effect.
2. Constant benefit without HTE.
3. Weak quantitative HTE.
4. Strong quantitative HTE.
5. Weak qualitative HTE.
6. Strong qualitative HTE.

These scenarios are used to show that population level HTE evidence, CATE ranking, policy value, and empirical surrogate harm metrics answer related but different questions. In particular, HTE evidence can be statistically real but still add little treatment assignment value when most patients share the same preferred treatment or when the signal is weak.

## ACTG 175 Analysis Summary

The ACTG 175 case study is a secondary analysis of a publicly available randomized trial. The scripts compare zidovudine (AZT) monotherapy with the two AZT-containing combination regimens, AZT plus didanosine (ddI) and AZT plus zalcitabine (ddC), after excluding ddI monotherapy. This contrast compares AZT alone with treatment intensification that retained AZT. The analysis constructs a 96-week event-free endpoint and uses IPCW to address early censoring.

The real data analysis is intended as an illustration of the two-stage framework. It should not be read as a new clinical recommendation for ACTG 175 treatment assignment. In the manuscript, the ACTG analysis shows that credible population level HTE evidence can coexist with limited evidence that a learned individualized policy improves held-out policy value.

## Reproducibility Notes

- Fixed random seeds are set inside the main scripts.
- The main simulation uses 500 replicates per scenario.
- The ACTG repeated validation script saves split-level results so that summary percentiles can be checked without rerunning the full analysis.
- Output directories are created by the scripts when possible.
- Run time depends on machine resources, especially for causal forest fitting and repeated validation.

## Citation

If you use this code, please cite the manuscript:

Xi NM, Huang X, Wang L. *From Heterogeneous Treatment Effects to Treatment Policies: A Decision Focused Framework for Clinical Trial Personalization.*
