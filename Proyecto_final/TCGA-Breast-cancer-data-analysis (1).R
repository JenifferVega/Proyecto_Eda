---
title: "TCGA Breast cancer data analysis"
author: "Andrea Bustos, Jeniffer Funez, Miguel Naranjo"
date: "2023-11-11"
output:
  html_document:
    df_print: paged
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(TCGAbiolinks)
library(biomaRt)
library(dplyr)
library(SummarizedExperiment)
```

```{r}
# a. TCGA breast cancer data download
query.exp.hg38 <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Transcriptome Profiling", 
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts"
)

# Check if the file already exists before attempting to download
if (!file.exists("exp.rda")) {
  GDCdownload(query.exp.hg38)
}

# Check for successful download
expdat <- GDCprepare(
  query = query.exp.hg38,
  save = TRUE, 
  save.filename = "exp.rda"
)

```

```{r}
# b. Data preprocessing – get gene expression

# Check available assays
available_assays <- names(assays(expdat))
print(available_assays)

# Choose the assay you want to use (replace "your_chosen_assay" with the actual assay name)
chosen_assay <- "tpm_unstrand"

# Check if the chosen assay is present
if (!(chosen_assay %in% available_assays)) {
  stop(paste("The chosen assay", chosen_assay, "is not present in the expdat object."))
}

# Extract gene expression
gene_expression <- assays(expdat)[[chosen_assay]]

```

```{r}
# c. Data preprocessing – filter gene expression
# Remove genes not expressed in at least 20% of samples
gene_expression_filtered <- gene_expression[, colMeans(gene_expression > 0) > 0.2]


```

```{r}
# d. Data preprocessing – Ensembl IDs
# Extract Ensembl IDs
ensembl_ids <- rownames(gene_expression_filtered)


```

```{r}
# e. Data preprocessing – Hugo Gene symbols
# Define the biomart dataset
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Get mapping from Ensembl IDs to Hugo Gene symbols
ensembl_to_symbol <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                            filters = "ensembl_gene_id",
                            values = ensembl_ids,
                            mart = mart)

# Merge the mapping with the gene expression data
gene_expression_mapped <- merge(ensembl_to_symbol, gene_expression_filtered,
                                by.x = "ensembl_gene_id", by.y = "row.names")

# Remove redundant column
gene_expression_mapped <- gene_expression_mapped[, -1]

```

```{r}
# f. Data preprocessing – Sample type
# Assuming the sample type information is in colData
sample_data <- DataFrame(colData(expdat))
sample_data <- as.data.frame(sample_data)  # Convert to a regular data frame

# Replace "SampleType" with the actual column containing the sample type information
sample_data <- sample_data %>% mutate(Sample_Type = ifelse(SampleType == "Normal", "Normal", "Tumor"))

```

```{r}
# g. Data preprocessing – Final data frame
# Assuming you want to keep only normal tissue and primary tumor samples
final_data <- gene_expression_filtered[, sample_data$Sample_Type %in% c("Normal", "Tumor")]

```

```{r}
# h. EDA: Exploratory Data Analysis
# (Your specific analysis will depend on the characteristics of your data)
# Perform a comparison between normal tissue and primary tumor
# Use the correct sampling technique to explore the relationships between the top 5 genes with the largest mean expression
# and the top 5 with the smallest mean expression. Explore a similar number of samples of each type.

```

```{r}
# i. DA: Perform Differential Analysis
# Use the provided PPI network to select the top 100 genes with the highest degree
# (This step requires additional information about the PPI network)
# Assuming network_data is a data frame with Ensembl IDs and network degrees
set.seed(125)
group1 <- c(rnorm(100, mean = 24, sd = 3))
group2 <- c(rnorm(100, mean = 43, sd = 2.4))

# Perform the t-test
ttest_results <- apply(network_data$Ensembl_ID, 1, function(gene) {
  ttest_res <- t.test(final_data[gene, sample_data$Sample_Type == "Normal"],
                      final_data[gene, sample_data$Sample_Type == "Tumor"])
  c(gene, p_value = ttest_res$p.value, statistic = abs(ttest_res$statistic))
})

# Convert results to data frame
ttest_results_df <- as.data.frame(t(ttest_results))
colnames(ttest_results_df) <- c("Ensembl_ID", "p_value", "statistic")

# Filter differentially expressed genes
DE_genes <- ttest_results_df[ttest_results_df$p_value < 0.05, ]

# List top 15 genes with the largest t test-statistic (absolute value)
top_15_genes <- DE_genes[order(abs(DE_genes$statistic), decreasing = TRUE), ][1:15, ]

```
