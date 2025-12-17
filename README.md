# Subgroup Identification and Individualized Treatment Policies: A Tutorial on the Hybrid Two-Stage Workflow

This repository contains the code that supports the results in the manuscript "Subgroup Identification and Individualized Treatment Policies: A Tutorial on the Hybrid Two-Stage Workflow". The ACTG 175 clinical trial data used in this study is extacted from the R package [speff2trial](https://cran.r-project.org/web/packages/speff2trial/index.html).

File | Purpose | Key outputs
-----|-----|-----|
sim_1.R | Implements the simulation study for the hybrid two-stage workflow. Performs STEPP plots, cross fitted DR policy value estimation, and NP threshold selection. Implements the binary outcome DGPs for No/Weak/Strong HTE, Stage 1 interaction testing, Stage 2 CATE learning via causal forests, and Monte Carlo evaluation across scenarios. | (i) Stage 1 p-values/gate decision, AUQC, DR value gain, NP feasibility/harm/benefit capture, AUROC/AUPRC vs true benefiter label. (ii) Summary statistics (proceed rate, mean AUQC, mean value gain, NP feasibility rate). (iii) CSV export of Monte Carlo results. (iv) STEPP curve, uplift curve, policy value curve, and NP safety dial plots. ​
zero_few_shot.R | Runs zero- and few-shot in-context learning for rare disease NER. Calculates precision, recall, and F1 score for each entity type by comparing model output against ground-truth annotations. | Predicted entities for each task and accuracy measurements. ​
RAG_error.R | Implements hybrid in-context learning for NER, combining semantically similar few-shot examples and RAG snippets for each entity type. Performs fine-grained error analysis using token overlap alignment. | Performance scores and six error rates by entity type. ​
fine_tune.R | Constructs training and validation jsonl files for fine-tuning GPT models on rare disease NER. | Structured prompts with labeled examples from RareDis corpus. ​
embedding.R | Computes 3,072-dimensional semantic embeddings for four datasets: test set, training set, validation set, and Orphar corpus. | Numerical embedding vectors. 
cost.R | Evaluates the trade-off between inference cost and F1 score across different few-shot settings. | A combined overlay plot comparing cost-efficiency across all four entity types. ​
external_DB.R | Create external RAG reference using Orphanet database. | A csv file for RAG snippet retrieval. ​
figure.R | Generates figures in the manuscript. | Figures 1 - 3. ​​

## Prerequisites
R ≥ 4.2

R packages required across scripts for evaluation, visualization, embedding, and modeling: httr, xml2, arrow, stringr, readxl, ggplot2, dplyr, tidyr, purrr, nls, tidyverse

Install them with
```r
install.packages(c(
  "httr", "xml2", "arrow", "stringr", "readxl",
  "ggplot2", "dplyr", "tidyr", "purrr",
  "nls", "tidyverse"
))
```

## Contact
For any questions or issue reports, pleasae open an issue or email nxi@ucla.edu.
