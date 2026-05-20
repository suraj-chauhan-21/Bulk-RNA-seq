#!/usr/bin/env Rscript

################################################################################
# Differential Gene Expression Analysis with DESeq2
# Performs DGE analysis, creates visualizations (PCA, MA plot, volcano plot)
# Prerequisites: Install required packages - see setup section below
################################################################################

# Set random seed for reproducibility
set.seed(42)

# ============================
# PART 1: Install/Load Packages
# ============================

required_packages <- c(
    "DESeq2",
    "tidyverse",
    "pheatmap",
    "RColorBrewer",
    "ggplot2",
    "ggrepel",
    "gridExtra"
)

# Check and install packages if needed
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if (length(new_packages) > 0) {
    cat("Installing missing packages...\n")
    
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
    }
    
    # Bioconductor packages
    bioconductor_pkgs <- c("DESeq2", "pheatmap")
    cran_pkgs <- c("tidyverse", "RColorBrewer", "ggplot2", "ggrepel", "gridExtra")
    
    if (any(bioconductor_pkgs %in% new_packages)) {
        BiocManager::install(bioconductor_pkgs[bioconductor_pkgs %in% new_packages])
    }
    
    if (any(cran_pkgs %in% new_packages)) {
        install.packages(cran_pkgs[cran_pkgs %in% new_packages], repos="http://cran.r-project.org")
    }
}

# Load libraries
suppressPackageStartupMessages({
    library(DESeq2)
    library(tidyverse)
    library(pheatmap)
    library(RColorBrewer)
    library(ggplot2)
    library(ggrepel)
})

cat("All required packages loaded successfully!\n\n")

# ============================
# PART 2: Setup Paths
# ============================

# Get environment variables
fastq_dir <- Sys.getenv("FASTQ")
counts_dir <- Sys.getenv("COUNTS")
results_dir <- Sys.getenv("RESULTS")
data_dir <- Sys.getenv("DATA")

# If environment variables not set, use default paths
if (fastq_dir == "") fastq_dir <- "BulkTranscriptomics/FASTQ"
if (counts_dir == "") counts_dir <- "BulkTranscriptomics/COUNTS"
if (results_dir == "") results_dir <- "BulkTranscriptomics/RESULTS"
if (data_dir == "") data_dir <- "BulkTranscriptomics/DATA"

# Create results directory if it doesn't exist
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

cat("Analysis Paths:\n")
cat("Counts directory:", counts_dir, "\n")
cat("Results directory:", results_dir, "\n")
cat("Data directory:", data_dir, "\n\n")

# ============================
# PART 3: Load Count Data
# ============================

cat("Loading count matrix...\n")

# Read count matrix
count_file <- file.path(counts_dir, "count_matrix.txt")
if (!file.exists(count_file)) {
    stop("Count matrix file not found at: ", count_file, "\n",
         "Please run quantification.sh first")
}

# Read count data
counts <- read.delim(count_file, row.names = 1, comment.char = "#")

# Remove metadata columns (first 5 columns are: Geneid, Chr, Start, End, Strand, Length)
counts <- counts[, -(1:5)]

# Clean sample names
colnames(counts) <- gsub("\\.bam$", "", colnames(counts))
colnames(counts) <- gsub(".*\\/", "", colnames(counts))

cat("Count matrix dimensions:", nrow(counts), "genes x", ncol(counts), "samples\n\n")

# ============================
# PART 4: Load Sample Metadata
# ============================

cat("Loading sample metadata...\n")

# Read sample metadata
metadata_file <- file.path(data_dir, "sample_metadata.csv")
if (!file.exists(metadata_file)) {
    stop("Sample metadata file not found at: ", metadata_file)
}

colData <- read.csv(metadata_file, row.names = 1)

# Ensure rownames match between counts and colData
common_samples <- intersect(rownames(colData), colnames(counts))
if (length(common_samples) == 0) {
    stop("No matching sample names between count matrix and metadata!")
}

counts <- counts[, common_samples]
colData <- colData[common_samples, ]

cat("Sample metadata loaded for", nrow(colData), "samples\n")
cat("Conditions:", paste(unique(colData$Condition), collapse = ", "), "\n\n")

# ============================
# PART 5: Create DESeq2 Object
# ============================

cat("Creating DESeq2 object...\n")

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = colData,
    design = ~ Condition
)

# Set reference level (for comparison)
dds$Condition <- factor(dds$Condition, levels = c("Normal", "Cancer"))

cat("DESeq2 object created\n")
cat("Design formula:", paste(as.formula(design(dds))), "\n\n")

# ============================
# PART 6: Quality Control - Pre-filtering
# ============================

cat("Pre-filtering low-count genes...\n")

# Pre-filter: keep genes with at least 10 reads total across all samples
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

cat("Genes after filtering:", nrow(dds), "\n\n")

# ============================
# PART 7: Run DESeq2
# ============================

cat("Running DESeq2 analysis (this may take a few minutes)...\n")

dds <- DESeq(dds, quiet = FALSE)

cat("DESeq2 analysis completed!\n\n")

# ============================
# PART 8: Extract Results
# ============================

cat("Extracting and processing results...\n")

# Get results
res <- results(dds, contrast = c("Condition", "Cancer", "Normal"))

# Add gene names
res$gene_id <- rownames(res)

# Sort by adjusted p-value
res_sorted <- res[order(res$padj), ]

# Summary statistics
summary(res)

cat("\n")
cat("Significant genes (padj < 0.05):", sum(res$padj < 0.05, na.rm = TRUE), "\n")
cat("Upregulated genes (log2FC > 0.5, padj < 0.05):", 
    sum(res$padj < 0.05 & res$log2FoldChange > 0.5, na.rm = TRUE), "\n")
cat("Downregulated genes (log2FC < -0.5, padj < 0.05):", 
    sum(res$padj < 0.05 & res$log2FoldChange < -0.5, na.rm = TRUE), "\n\n")

# ============================
# PART 9: Save Results
# ============================

cat("Saving results to files...\n")

# Save all results
write.csv(as.data.frame(res_sorted), 
          file = file.path(results_dir, "DEG_results_all.csv"))

# Save filtered results (padj < 0.05)
deg_sig <- res_sorted[res_sorted$padj < 0.05, ]
write.csv(as.data.frame(deg_sig), 
          file = file.path(results_dir, "DEG_results_significant.csv"))

cat("Results saved to:\n")
cat("  - DEG_results_all.csv\n")
cat("  - DEG_results_significant.csv\n\n")

# ============================
# PART 10: Visualizations
# ============================

cat("Creating visualizations...\n")

# 10.1: PCA Plot
cat("  Generating PCA plot...\n")
vst <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vst, intgroup = "Condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = Condition)) +
    geom_point(size = 4) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "right"
    ) +
    ggtitle("PCA Plot - Sample Clustering")

ggsave(file.path(results_dir, "01_PCA_plot.pdf"), pca_plot, width = 8, height = 6)

# 10.2: MA Plot
cat("  Generating MA plot...\n")
ma_plot <- ggplot(as.data.frame(res), aes(x = baseMean, y = log2FoldChange)) +
    geom_point(aes(color = padj < 0.05 & abs(log2FoldChange) > 0.5), size = 2, alpha = 0.6) +
    scale_color_manual(
        values = c("TRUE" = "red", "FALSE" = "gray"),
        labels = c("TRUE" = "Significant (padj<0.05, |log2FC|>0.5)", "FALSE" = "Not significant")
    ) +
    scale_x_log10() +
    xlab("Mean of Normalized Counts (log10)") +
    ylab("Log2 Fold Change") +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "top"
    ) +
    ggtitle("MA Plot - Fold Change vs Expression Level")

ggsave(file.path(results_dir, "02_MA_plot.pdf"), ma_plot, width = 10, height = 6)

# 10.3: Volcano Plot
cat("  Generating Volcano plot...\n")
volcano_data <- as.data.frame(res) %>%
    mutate(
        sig = padj < 0.05 & abs(log2FoldChange) > 0.5,
        label = ifelse(sig, rownames(as.data.frame(res)), NA)
    )

volcano_plot <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = sig), size = 2, alpha = 0.6) +
    scale_color_manual(
        values = c("TRUE" = "red", "FALSE" = "gray"),
        labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
    ) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
    xlab("Log2 Fold Change") +
    ylab("-Log10 Adjusted P-value") +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "top"
    ) +
    ggtitle("Volcano Plot - Differential Expression")

ggsave(file.path(results_dir, "03_Volcano_plot.pdf"), volcano_plot, width = 10, height = 7)

# 10.4: Heatmap of Top Genes
cat("  Generating heatmap of top genes...\n")
top_genes <- rownames(deg_sig)[1:min(50, nrow(deg_sig))]
heatmap_data <- assay(vst)[top_genes, ]

pdf(file.path(results_dir, "04_Heatmap_top_genes.pdf"), width = 10, height = 12)
pheatmap(
    heatmap_data,
    scale = "row",
    annotation_col = as.data.frame(colData[, "Condition", drop = FALSE]),
    color = colorRampPalette(c("blue", "white", "red"))(50),
    main = "Heatmap of Top 50 Significant Genes",
    cluster_cols = TRUE,
    cluster_rows = TRUE
)
dev.off()

cat("  Generating expression boxplot...\n")
# 10.5: Expression Boxplot
top_3_genes <- rownames(deg_sig)[1:3]
for (gene in top_3_genes) {
    expr_data <- data.frame(
        expression = assay(vst)[gene, ],
        condition = colData$Condition,
        sample = colnames(vst)
    )
    
    box_plot <- ggplot(expr_data, aes(x = condition, y = expression, fill = condition)) +
        geom_boxplot(alpha = 0.7) +
        geom_point(size = 3, position = position_jitter(width = 0.2)) +
        xlab("Condition") +
        ylab("Expression (VST)") +
        theme_minimal() +
        theme(
            plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
            legend.position = "none"
        ) +
        ggtitle(paste("Expression of", gene))
    
    ggsave(file.path(results_dir, paste0("05_Boxplot_", gene, ".pdf")), box_plot, width = 6, height = 5)
}

cat("✓ All visualizations saved\n\n")

# ============================
# PART 11: Summary Report
# ============================

cat("========================================\n")
cat("ANALYSIS SUMMARY\n")
cat("========================================\n")
cat("Total genes analyzed:", nrow(dds), "\n")
cat("Total samples:", ncol(dds), "\n")
cat("Conditions:", paste(unique(colData$Condition), collapse = ", "), "\n")
cat("\nStatistical Results:\n")
cat("  Significant genes (padj < 0.05):", sum(res$padj < 0.05, na.rm = TRUE), "\n")
cat("  Upregulated:", sum(res$padj < 0.05 & res$log2FoldChange > 0.5, na.rm = TRUE), "\n")
cat("  Downregulated:", sum(res$padj < 0.05 & res$log2FoldChange < -0.5, na.rm = TRUE), "\n")
cat("\nOutput files:\n")
cat("  Results: DEG_results_all.csv, DEG_results_significant.csv\n")
cat("  Plots: PCA, MA, Volcano, Heatmap, Boxplots\n")
cat("========================================\n")

cat("\nDifferential expression analysis completed!\n")
cat("Next step: Run functional enrichment analysis\n")
cat("  Rscript ", file.path(Sys.getenv("SCRIPTS"), "enrichment_analysis.R"), "\n")
