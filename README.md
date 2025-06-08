# Differential Expression Analysis Shiny App with DESeq2

This Shiny application enables interactive differential gene expression analysis using the **DESeq2** package. Users can upload their own count and condition data, run DESeq2 analysis, and visualize results through various plots including MA plot, boxplot, scatter plot, volcano plot, and heatmap.

---

## Features

- Upload count matrix and sample condition files (CSV or TXT).
- Perform DESeq2 differential expression analysis with a single click.
- Interactive previews of uploaded data.
- Visualization of:
  - MA Plot
  - Boxplot of normalized counts
  - Scatter plot comparing conditions
  - Volcano plot of differential expression results
  - Heatmap of top variable genes
- Downloadable plots as PNG files.
- Progress bar and log messages during analysis.

---

## Installation

Make sure you have the following R packages installed:

```r
install.packages(c("shiny", "ggplot2", "pheatmap", "RColorBrewer", "dplyr", "tidyverse", "matrixStats", "DT"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("DESeq2", "EnhancedVolcano"))
```

## How to Use
- Prepare Input Files

- Count Data: A CSV or tab-delimited file with genes as rows and samples as columns. The first column should be gene identifiers (used as row names).
- Example:
```
Gene	Sample1	Sample2	Sample3	Sample4
GeneA	10	15	20	5
GeneB	0	3	2	0
```
- Condition Data: A CSV or tab-delimited file describing sample conditions. It must contain two columns: Sample and Condition.
- Example:
```
Sample	Condition
Sample1	Control
Sample2	Control
Sample3	Treatment
Sample4	Treatment
```
---
- Upload Files

Upload the count data file under "Upload Count Data" and the condition data file under "Upload Condition Data".

- Run Analysis

Click the "Run DESeq2 Analysis" button. Wait for the progress bar and log message confirming completion.

- Explore Results

Use the tabs in the main panel to view data previews and plots.

- Download Plots

Each plot tab includes a download button to save the respective plot as a PNG file.

## Notes
The app automatically filters out genes with zero counts across all samples.

At least two distinct conditions are required for scatter plot and differential expression analysis.

The heatmap shows the top 50 most variable genes based on normalized counts.

Maximum upload file size is set to 100 MB.
