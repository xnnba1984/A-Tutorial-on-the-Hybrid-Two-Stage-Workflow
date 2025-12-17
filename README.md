# Leveraging Large Language Models for Rare Disease Named Entity Recognition

This repository contains the code that supports all results in the manuscript "Leveraging Large Language Models for Rare Disease Named Entity Recognition". The data used in this study can be downloaded [here](https://github.com/isegura/NLP4RARE-CM-UC3M) (RareDis corpus) and [here](https://www.orphadata.com/orphanet-scientific-knowledge/) (Orphanet database).

File | Purpose | Key outputs
-----|-----|-----|
data_clean_rare.R | Reads and processes .ann and .txt files from the RareDis dataset. Extracts annotated entities, filters by four target types, normalizes text, and removes duplicates. Combines text and labels into a structured list. | Named list with document-level text and label objects used in NER evaluation. ​
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
