# ================================================
# Program: Gene Expression Analysis
# Author: Margot María Fernández Cuadros
# Date: 2026-11-16
# Description: Differential gene expression analysis
# ================================================

rm(list =ls())
library(ggplot2)
library(reshape2)
library(pheatmap)

# ===============================================
# STEP 1: READ THE FILES
# ===============================================

#Get command-line arguments
args<-commandArgs(trailingOnly = TRUE)


#The user can choose the p-value threshold (if the user )
if (length(args) != 3 && length(args) != 5) {
  stop("The program requires 3 arguments and an optional --pvalue value")
}

#Define input and output file paths
counts_file<-args[1]
metadata_file<-args[2]
output_dir<-args[3]

#check errors
if (!file.exists(counts_file)){
  stop("The counts file does not exist")
}

if (!file.exists(metadata_file)){
  stop("The metadata file does not exist")
}

# Set the default p-value threshold
pvalue_threshold <- 0.05

# Read the optional p-value threshold
if (length(args) == 5) {
  if (args[4] != "--pvalue") {
    stop("The optional argument must be --pvalue")
  }
  pvalue_threshold <- as.numeric(args[5])
}

# Check that the p-value threshold is valid
if (is.na(pvalue_threshold) || pvalue_threshold <= 0 || pvalue_threshold >= 1) {
  stop("The p-value threshold must be a number between 0 and 1")
}

cat(" Leyendo archivo:", counts_file, "\n")
cat(" Leyendo archivo:", metadata_file, "\n")
cat(" Leyendo archivo:", output_dir, "\n")

#Read the files
counts <- read.delim(counts_file,header = TRUE)
metadata<-read.delim(metadata_file,header=TRUE)

cat("\n Files uploaded successfully\n")
cat("Number of genes:", nrow(counts), "\n")
cat("Number of columns:", ncol(counts), "\n")
cat("Groups found in metadata:",unique(metadata$group),"\n")

cat("Starting gene expression analysis...\n")

# ===============================================
# STEP 2: PREPARE THE EXPRESSION MATRIX
# ===============================================

#assign a name to each row
rownames(counts)<-counts$GeneID

#Retaining the necessary data
count_matrix<-counts[,-1]

#logarithmic transformation
log_matrix<-log2(count_matrix+1)

# ===============================================
# STEP 3: DIFFERENTIAL EXPRESSION OF EACH GENE
# ===============================================

cat("Starting differential expression analysis of each gene...\n")

#T-TEST
# -----------------------------------------------
treated<-count_matrix[,c("T1","T2","T3")]
control<-count_matrix[,c("Con1","Con2","Con3")]

p_values<-numeric(nrow(count_matrix))

for (i in 1:nrow(count_matrix)){
  x<-treated[i,]
  y<-control[i,]
  test_result<-tryCatch(t.test(x,y),error=function(e) NULL)
    if (!is.null(test_result)){
      p_values[i]<-test_result$p.value
    } else {
      p_values[i]<-NA
    }
}

cat("T-test completed successfully!\n")

#MEAN
#------------------------------------------------
meanTreated<-numeric(nrow(count_matrix))
for(i in 1:nrow(treated)){
  meanT<-mean(unlist(treated[i,]))
  meanTreated[i]<-meanT
}

meanControl<-numeric(nrow(count_matrix))
for(i in 1:nrow(control)){
  meanC<-mean(unlist(control[i,]))
  meanControl[i]<-meanC
}

cat("Mean completed successfully!\n")

#LOG2 FOLD CHANGE
#-----------------------------------------------
log2FC<-log2((meanTreated+1)/(meanControl+1))

cat("Log2 fold change completed successfully!\n")

#Results of differential expression of each gene
#-------------------------------------------------
results_differential_expression<-data.frame(GeneID=counts$GeneID,MEAN_TREATED=meanTreated,MEAN_CONTROL=meanControl,Log2FC=log2FC,P_VALUE=p_values)

# Adjust p-values using the Benjamini-Hochberg method
results_differential_expression$ADJUSTED_P_VALUE <-
  p.adjust(results_differential_expression$P_VALUE, method = "BH")

cat("Benjamini-Hochberg method completed successfully!\n")

#Save differential expression results
write.table(
  results_differential_expression,
  file = file.path(output_dir, "differential_expression.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Differential expression results saved successfully!\n")

# ======================================================================================
# STEP 4: SIGNIFICANT GENES (P-VALUE THRESHOLD INDICATED BY THE USER, OR 0.05 BY DEFAULT)
# ======================================================================================

#Create the table of significant genes using the original p-value
significant_genes<-subset(results_differential_expression,P_VALUE<pvalue_threshold)

cat("P-value threshold:", pvalue_threshold, "\n")

# Create the table of significant genes using the adjusted p-value
significant_genes_BH <- subset(results_differential_expression, ADJUSTED_P_VALUE < pvalue_threshold)

cat("Number of significant genes (original p-value):", nrow(significant_genes), "\n")

cat("Number of significant genes (adjusted p-value - BH):", nrow(significant_genes_BH), "\n")

#Save significant genes
write.table(
  significant_genes,
  file = file.path(output_dir, "significant_genes.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Significant genes results saved successfully!\n")

# ===============================================
# STEP 5: BOXPLOT
# ===============================================

cat("Starting Boxplot...\n")

#Reorganise the matrix to use GGPLOT2
#------------------------------------------------

#Conversion to "long"
log_matrix_df <- as.data.frame(log_matrix)
log_matrix_df$GeneID <- rownames(log_matrix)

log_matrix_long <- melt(
  log_matrix_df,
  id.vars = "GeneID",
  variable.name = "Sample",
  value.name = "Expression"
)

log_matrix_long$Condition <- ifelse(
  log_matrix_long$Sample %in% c("T1", "T2", "T3"),
  "Treated",
  "Control"
)

#Create the boxplot
#----------------------------------------------
ggplot(log_matrix_long,
       aes(x = Sample, y = Expression, fill = Condition)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  scale_fill_manual(values = c(
    "Treated" = "#D00000",
    "Control" = "#00A83B"
  )) +
  labs(
    title = "Gene Expression Distribution by Sample",
    x = "Samples",
    y = "log2(Counts+1)",
    fill = "Condition:"
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_line(
      color = "gray50",
      linetype = "dotted",
      linewidth = 0.5
    ),
    panel.grid.minor = element_line(
      color = "gray50",
      linetype = "dotted",
      linewidth = 0.3
    ),
    axis.line = element_line(color = "black", linewidth = 0.7),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(face = "bold", size = 10),
    legend.position = "top",
    plot.title = element_text(face = "bold", size = 16, hjust=0.5)
  )
#Save the boxplot
ggsave(
  filename = file.path(output_dir, "Boxplot.png"),
  width = 10,
  height = 6,
  dpi = 300
)

cat("Boxplot saved successfully!\n")

# ===============================================
# STEP 6: SCATTER PLOT
# ===============================================

cat("Starting scatter plot...\n")

#We use "results_differential_expression", which contains the means
#Create the scatter plot
ggplot(results_differential_expression, aes(x = MEAN_CONTROL, y = MEAN_TREATED)) +
  geom_point(alpha = 0.4, color = "#4040fb") +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  labs(title = "Linear relation between Control and Treated",
       subtitle = sprintf("r = %.3f, p = %.2g",
                          cor(results_differential_expression$MEAN_CONTROL, results_differential_expression$MEAN_TREATED),
                          cor.test(results_differential_expression$MEAN_CONTROL, results_differential_expression$MEAN_TREATED)$p.value),
       x = "CONTROL", y = "TREATED") +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_line(
      color = "gray50",
      linetype = "dotted",
      linewidth = 0.5
    ),
    panel.grid.minor = element_line(
      color = "gray50",
      linetype = "dotted",
      linewidth = 0.3
    ),
    axis.line = element_line(color = "black", linewidth = 0.7),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold", size = 16, hjust=0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

#Save the Scatter plot
ggsave(
  filename = file.path(output_dir, "Scatterplot.png"),
  width = 10,
  height = 6,
  dpi = 300
)

cat("Scatter plot saved successfully!\n")

# ===============================================
# STEP 7: HISTOGRAM
# ===============================================

cat("Stating histogram...\n")

ggplot(log_matrix_long, aes(x = Expression)) +
  geom_histogram(bins = 30, fill = "#ff8000", color = "white") +
  labs(title = "General distribution of gene expression",
       x = "Expression (log2 counts + 1)",
       y = "Number de observations") +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_line(
      color = "gray50",
      linetype = "dotted",
      linewidth = 0.5
    ),
    panel.grid.minor = element_line(
      color = "gray50",
      linetype = "dotted",
      linewidth = 0.3
    ),
    axis.line = element_line(color = "black", linewidth = 0.7),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold", size = 16, hjust=0.5)
  )

#Save the histogram
ggsave(
  filename = file.path(output_dir, "Histogram.png"),
  width = 10,
  height = 6,
  dpi = 300
)

cat("Histogram saved successfully!\n")

# ===============================================
# STEP 8: VIOLIN PLOT
# ===============================================

cat("Starting violin plot...\n")

ggplot(log_matrix_long, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  scale_fill_manual(values = c(
    "Treated" = "#D00000",
    "Control" = "#00A83B"
  )) +
  labs(title = "Expression distribution",
       y = "Expression (UA)") +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_line(
      color = "gray50",
      linetype = "dotted",
      linewidth = 0.5
    ),
    panel.grid.minor = element_line(
      color = "gray50",
      linetype = "dotted",
      linewidth = 0.3
    ),
    axis.line = element_line(color = "black", linewidth = 0.7),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold", size = 16, hjust=0.5)
  )

#Save the violin plot
ggsave(
  filename = file.path(output_dir, "Violinplot.png"),
  width = 10,
  height = 6,
  dpi = 300
)

cat("Violin plot saved successfully!\n")

# ===============================================
# STEP 9: DENSITY PLOT BY CONDITION
# ===============================================

cat("Starting density plot by condition...\n")

ggplot(log_matrix_long, aes(x = Expression, fill = Condition)) +
  geom_density(alpha = 0.5, linewidth=0.8) +
  scale_fill_manual(values = c(
    "Treated" = "#D00000",
    "Control" = "#00A83B"
  )) +
  labs(title = "Expression distribution by condition",
       x = "Expression",
       y = "Density") +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_line(
      color = "gray50",
      linetype = "dotted",
      linewidth = 0.5
    ),
    panel.grid.minor = element_line(
      color = "gray50",
      linetype = "dotted",
      linewidth = 0.3
    ),
    axis.line = element_line(color = "black", linewidth = 0.7),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold", size = 16, hjust=0.5)
  )

#Save the density plot
ggsave(
  filename = file.path(output_dir, "densityplot.png"),
  width = 10,
  height = 6,
  dpi = 300
)

cat("Density plot saved successfully!\n")

# ===============================================
# STEP 10: BAR PLOT
# ===============================================

cat("Starting bar plot...\n")

#Classify the genes into: "Upregulated","Downregulated" and "Unchanged"
results_differential_expression$Regulation <- ifelse(
  is.na(results_differential_expression$P_VALUE) |
    is.na(results_differential_expression$Log2FC),
  "Not evaluated",
  ifelse(
    results_differential_expression$P_VALUE < pvalue_threshold &
      results_differential_expression$Log2FC > 0,
    "Upregulated",
    ifelse(
      results_differential_expression$P_VALUE < pvalue_threshold &
        results_differential_expression$Log2FC < 0,
      "Downregulated",
      "Unchanged"
    )
  )
)
#Create a summary table with the count for each classification
data_bar<-as.data.frame(table(results_differential_expression$Regulation)
)

colnames(data_bar) <- c("REGULATION", "COUNT")

# Exclude genes that could not be evaluated
data_bar_plot <- subset(data_bar, REGULATION != "Not evaluated")

#create the bar plot
ggplot(data_bar_plot, aes(x = REGULATION, y = COUNT, fill = REGULATION)) +
  geom_bar(stat = "identity", alpha = 0.7) +
  scale_fill_manual(values = c(
    "Upregulated" = "#D00000",
    "Downregulated" = "#00A83B",
    "Unchanged" = "#808080"
  )) +
  labs(
    title = "Distribution of Gene Regulation",
    x = "Regulation",
    y = "Number of genes",
    fill = "Regulation"
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_line(
      color = "gray75",
      linetype = "dotted",
      linewidth = 0.5
    ),
    panel.grid.minor = element_line(
      color = "gray88",
      linetype = "dotted",
      linewidth = 0.3
    ),
    axis.line = element_line(color = "black", linewidth = 0.7),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
  )

#Save the bar plot
ggsave(
  filename = file.path(output_dir, "barplot.png"),
  width = 10,
  height = 6,
  dpi = 300
)

cat("Bar plot saved successfully!\n")

# ===============================================
# STEP 11: HEATMAP
# ===============================================

cat("Starting heatmap...\n")

#Calculate the variance of each gene
variance<-apply(log_matrix,1,var,na.rm=TRUE)

#Select the 50 most variable genes
top_n<-50

top_50_genes<-names(sort(variance,decreasing=TRUE))[seq_len(top_n)]
mat_top<-log_matrix[top_50_genes,]

#Group information
grp<-factor(metadata$group)
grp_lvls<-levels(grp)

#Assign colors to group
custom_colors<-c("#00A83B", "#D00000")
grp_cols_map<-setNames(custom_colors, grp_lvls)

# Create column annotation
annotation_col <- data.frame(
  Group = grp
)

rownames(annotation_col) <- metadata$sample

# Define heatmap colors
pal <- colorRampPalette(c("blue", "white", "red"))(101)


# Create heatmap
pheatmap(
  mat_top,
  scale = "row",
  color = pal,
  annotation_col = annotation_col,
  annotation_colors = list(Group = grp_cols_map),
  fontsize_row = 8,
  main = paste0("50 MOST VARIABLE GENES"),
  filename = file.path(output_dir, "heatmap.png"),
  width = 10,
  height = 6
)

cat("Heatmap saved successfully!\n")

# ===============================================
# STEP 12: PCA
# ===============================================

cat("Starting principal component analysis (PCA)...\n")

set.seed(1313)

#transpose expression matrix
mat_pca<-t(mat_top)

#Prcomp
pca<-prcomp(mat_pca,center=TRUE,scale.=TRUE)

# Calculate the percentage of variance explained
variance <- pca$sdev^2
variance_percent <- variance / sum(variance) * 100

# Assign colors to samples
pca_colors <- grp_cols_map[grp]

# Save the PCA plot
png(file.path(output_dir, "PCA.png"), width = 800, height = 600)

#Create the PCA plot
plot(pca$x[,1], pca$x[,2],
    pch = 19, cex=1.5, col = pca_colors,
     xlab = paste0("PC1 (", round(variance_percent[1], 2), "%)"),
     ylab = paste0("PC2 (", round(variance_percent[2], 2), "%)"),
     main = "Principal Component Analysis(PCA)")
text(pca$x[,1], pca$x[,2], labels = rownames(pca$x),pos=3, offset=1, font=2, cex=0.9)
abline(h=0,v=0,lty=2,col="gray70")
legend("topright", inset = c(0.02), legend = names(grp_cols_map),
       col = grp_cols_map, pch = 19, bty = "n", title = "Experimental Group",text.font=2,title.font=2)

cat("PCA plot saved successfully!\n")


# ===============================================
# STEP 13: VOLCANO PLOT
# ===============================================

cat("Starting volvano plot...\n")

# Calculate the negative log10 of the p-value
results_differential_expression$NEG_LOG10_PVALUE <- -log10(results_differential_expression$P_VALUE)

# Exclude genes that were not evaluated
volcano_data <- subset(results_differential_expression, Regulation != "Not evaluated")

# Create the volcano plot
ggplot(volcano_data,
       aes(x = Log2FC, y = NEG_LOG10_PVALUE, color = Regulation)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c(
    "Upregulated" = "#D00000",
    "Downregulated" = "#00A83B",
    "Unchanged" = "#808080",
    "Not evaluated" = "#000000"
  )) +
  labs(
    title = "VOLCANO PLOT",
    x = "Log2 Fold Change",
    y = "-log10(P-value)",
    color = "Regulation"
  ) +
  theme_minimal()+
  theme(
    plot.title = element_text(hjust = 0.5)
  )

#Save the volcano plot
ggsave(
  filename = file.path(output_dir, "volcano.png"),
  width = 10,
  height = 6,
  dpi = 300
)

# Close the graphics device
dev.off()

cat("Volcano plot saved successfully!\n")


cat("Analysis completed successfully!\n")