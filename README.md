# Subgroup Identification and Individualized Treatment Policies: A Tutorial on the Hybrid Two-Stage Workflow

This repository contains the code that supports the results in the manuscript "Subgroup Identification and Individualized Treatment Policies: A Tutorial on the Hybrid Two-Stage Workflow". The ACTG 175 clinical trial data used in this study is extacted from the R package [speff2trial](https://cran.r-project.org/web/packages/speff2trial/index.html).

File | Purpose | Key outputs
-----|-----|-----|
sim_1.R | Implements the simulation study for the hybrid two-stage workflow. Performs STEPP plots, cross fitted DR policy value estimation, and NP threshold selection. Implements the binary outcome DGPs for No/Weak/Strong HTE, Stage 1 interaction testing, Stage 2 CATE learning via causal forests, and Monte Carlo evaluation across scenarios. | (i) Stage 1 p-values/gate decision, AUQC, DR value gain, NP feasibility/harm/benefit capture, AUROC/AUPRC vs true benefiter label. (ii) Summary statistics (proceed rate, mean AUQC, mean value gain, NP feasibility rate). (iii) CSV export of Monte Carlo results. (iv) STEPP curve, uplift curve, policy value curve, and NP safety dial plots. ​
Real_data.R | Analyzes the ACTG 175 clinical trial data using the hybrid workflow. Implements Stage 1 global heterogeneity testing with STEPP plots for key biomarkers. Performs Stage 2 nuisance modeling and causal forest CATE estimation, followed by uplift/AUQC computation, DR policy value curves with threshold selection, and an NP frontier under a specified harm tolerance. | (i) Clean analytic dataset. (ii) Stage 1 LRT statistic, p-values, STEPP plots along baseline CD4 and Karnofsky. (iii) Stage 2 uplift summary, AUQC, policy value, selected threshold with corresponding DR value, and NP frontier. ​​

## Prerequisites
R ≥ 4.2

R packages required across scripts: 
```data.table, mgcv, grf, speff2trial, ggplot2, dplyr```

Install them with
```r
install.packages(c(
  "data.table", "mgcv", "grf", "speff2trial", "ggplot2", "dplyr"
))
```

## Contact
For any questions or issue reports, pleasae open an issue or email nxi@ucla.edu.
