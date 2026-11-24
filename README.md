# rnaseq-differential-expression-r
Command-line R pipeline for exploratory data analysis and differential gene expression on RNA-seq count matrices.
# Gene Expression & Differential Analysis Pipeline in R

Command-line R pipeline developed by **Margot María Fernández Cuadros** for exploratory data analysis (EDA) and gene-by-gene differential expression testing on RNA-seq count matrices (GEO dataset: **GSE294583**).

##  Key Features
* **Command-Line Interface (CLI):** Dynamically parses input file paths, output directories, and an optional statistical significance threshold (`--pvalue`) via `commandArgs()`.
* **Data Preprocessing:** Applies a $\log_2(\text{counts} + 1)$ transformation to stabilize variance across samples.
* **Statistical Testing & Multiple Testing Correction:** Performs gene-wise two-sample Student's t-tests (with safe error handling via `tryCatch`), computes group means (`MEAN_TREATED`, `MEAN_CONTROL`), $\log_2$ fold changes ($\log_2\text{FC}$), and adjusts p-values using the **Benjamini-Hochberg (FDR)** method.
* **Comprehensive Visualization Suite:** Generates 9 publication-ready plots (300 DPI) covering sample distributions, dimensionality reduction, and differential expression patterns.
* **Reproducible Reporting:** Accompanied by an RMarkdown.

---

##  Requirements

* **R** (>= 4.0)
* **Required R Packages:** `ggplot2`, `reshape2`, `pheatmap`

Install all dependencies in R before running the pipeline:

    install.packages(c("ggplot2", "reshape2", "pheatmap"))

---

##  How to Run (Command-Line Usage)

Open your terminal, navigate to the root directory of this repository, and run the script using `Rscript`.

The program requires **3 mandatory arguments** in order:
1. Path to the counts matrix (`data/counts.tsv`)
2. Path to the metadata file (`data/metadata.tsv`)
3. Path to the output directory (`results/`)

Optionally, you can pass `--pvalue` followed by a custom significance threshold between `0` and `1` (defaults to `0.05` if omitted).

### 1. Standard Execution (Default threshold: p < 0.05)

    Rscript rnaseq_dea.R data/counts.tsv data/metadata.tsv results/

### 2. Custom P-value Threshold Execution (e.g., p < 0.01)

    Rscript rnaseq_dea.R data/counts.tsv data/metadata.tsv results/ --pvalue 0.01

---

##  Generated Outputs

Running the script populates the `results/` folder with the following files:

### Tabular Results (`.tsv`)
* **`differential_expression.tsv`**: Full gene-level statistical results (`GeneID`, `MEAN_TREATED`, `MEAN_CONTROL`, `Log2FC`, `P_VALUE`, `ADJUSTED_P_VALUE`).
* **`significant_genes.tsv`**: Filtered subset containing only genes with `P_VALUE` below the selected threshold.

### Visualizations (`.png`)
* **`Boxplot.png`**: Sample-level $\log_2(\text{counts} + 1)$ expression distribution colored by condition (`Control` vs. `Treated`).
* **`Scatterplot.png`**: Linear correlation between `MEAN_CONTROL` and `MEAN_TREATED` with regression line, Pearson's $r$, and p-value.
* **`Histogram.png`**: Overall distribution of transformed expression values.
* **`Violinplot.png`**: Condition-level expression density and summary distribution.
* **`densityplot.png`**: Overlapping kernel density curves comparing `Control` and `Treated` groups.
* **`barplot.png`**: Summary counts of `Upregulated`, `Downregulated`, and `Unchanged` genes.
* **`heatmap.png`**: Row-scaled hierarchical clustering heatmap of the **top 50 most variable genes**.
* **`PCA.png`**: Principal Component Analysis (`PC1` vs. `PC2`) computed on the top 50 most variable genes, displaying explained variance percentages.
* **`volcano.png`**: Volcano plot mapping $\log_2\text{FC}$ against $-\log_{10}(\text{P-value})$ colored by regulation status.

---

##  Repository Structure

    ├── README.md                 # Project documentation and execution instructions
    ├── data/
    │   ├── counts.tsv            # Processed gene expression count matrix
    │   └── metadata.tsv          # Sample metadata and experimental conditions
    ├── results/                  # Default output files generated with p < 0.05
    ├── prepare_data.R            # Auxiliary script for filtering and formatting raw GEO data
    ├── rnaseq_dea.R              # Main command-line differential expression pipeline
    ├── rnaseq_report.Rmd         # Reproducible RMarkdown source file
    ├── rnaseq_report.pdf         # Rendered RMarkdown report with embedded plots
    
