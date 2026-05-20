#!/usr/bin/env Rscript

################################################################################
# Functional Enrichment Analysis
# Performs Gene Ontology (GO) and KEGG pathway enrichment analysis
# Uses clusterProfiler package
################################################################################

set.seed(42)

# ============================
# PART 1: Install/Load Packages
# ============================

required_packages <- c(
    "clusterProfiler",
    "org.Hs.eg.db",
    "tidyverse",
    "ggplot2"
)

# Check and install packages if needed
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if (length(new_packages) > 0) {
    cat("Installing missing packages...\n")
    
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
    }
    
    # Bioconductor packages
    bioconductor_pkgs <- c("clusterProfiler", "org.Hs.eg.db")
    cran_pkgs <- c("tidyverse", "ggplot2")
    
    if (any(bioconductor_pkgs %in% new_packages)) {
        BiocManager::install(bioconductor_pkgs[bioconductor_pkgs %in% new_packages])
    }
    
    if (any(cran_pkgs %in% new_packages)) {
        install.packages(cran_pkgs[cran_pkgs %in% new_packages], repos="http://cran.r-project.org")
    }
}

# Load libraries
suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(tidyverse)
    library(ggplot2)
})

cat("All required packages loaded successfully!\n\n")

# ============================
# PART 2: Setup Paths
# ============================

# Get environment variables
results_dir <- Sys.getenv("RESULTS")
data_dir <- Sys.getenv("DATA")

# If environment variables not set, use default paths
if (results_dir == "") results_dir <- "BulkTranscriptomics/RESULTS"
if (data_dir == "") data_dir <- "BulkTranscriptomics/DATA"

# Create results subdirectory for enrichment
enrichment_dir <- file.path(results_dir, "functional_enrichment")
dir.create(enrichment_dir, showWarnings = FALSE, recursive = TRUE)

cat("Analysis Paths:\n")
cat("Results directory:", results_dir, "\n")
cat("Enrichment directory:", enrichment_dir, "\n\n")

# ============================
# PART 3: Load DEG Results
# ============================

cat("Loading differential expression results...\n")

# Read DEG results
deg_file <- file.path(results_dir, "DEG_results_significant.csv")
if (!file.exists(deg_file)) {
    stop("DEG results file not found at: ", deg_file, "\n",
         "Please run dge_analysis.R first")
}

deg_results <- read.csv(deg_file, row.names = 1)

cat("Loaded", nrow(deg_results), "significant genes\n\n")

# ============================
# PART 4: Gene ID Conversion
# ============================

cat("Converting gene names to Entrez IDs...\n")

# Get unique gene IDs
gene_names <- rownames(deg_results)

# Convert gene symbols to Entrez IDs
gene_to_entrez <- bitr(
    gene_names,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
)

# Create gene lists (all significant genes and subsets by direction)
all_genes <- gene_to_entrez$ENTREZID
upregulated_genes <- gene_to_entrez$ENTREZID[
    gene_to_entrez$SYMBOL %in% rownames(deg_results)[deg_results$log2FoldChange > 0]
]
downregulated_genes <- gene_to_entrez$ENTREZID[
    gene_to_entrez$SYMBOL %in% rownames(deg_results)[deg_results$log2FoldChange < 0]
]

cat("Successfully converted", length(all_genes), "genes\n")
cat("  Upregulated:", length(upregulated_genes), "\n")
cat("  Downregulated:", length(downregulated_genes), "\n\n")

# ============================
# PART 5: GO Enrichment Analysis
# ============================

cat("========================================\n")
cat("GENE ONTOLOGY ENRICHMENT ANALYSIS\n")
cat("========================================\n\n")

# GO analysis for all significant genes
cat("Performing GO enrichment for all significant genes...\n")

go_results <- enrichGO(
    gene = all_genes,
    OrgDb = org.Hs.eg.db,
    ont = "BP",  # Biological Process
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    readable = TRUE
)

cat("Found", nrow(go_results), "significant GO terms (BP)\n")

# Save GO results
write.csv(
    as.data.frame(go_results),
    file = file.path(enrichment_dir, "GO_enrichment_BP.csv")
)

# Visualize GO results
if (nrow(go_results) > 0) {
    cat("Generating GO enrichment plots...\n")
    
    # Dotplot
    pdf(file.path(enrichment_dir, "GO_dotplot.pdf"), width = 10, height = 8)
    print(dotplot(go_results, showCategory = 20, title = "GO Enrichment - Top 20 Terms"))
    dev.off()
    
    # Barplot
    pdf(file.path(enrichment_dir, "GO_barplot.pdf"), width = 10, height = 8)
    print(barplot(go_results, showCategory = 20, title = "GO Enrichment - Bar Plot"))
    dev.off()
    
    cat("✓ GO enrichment plots saved\n\n")
} else {
    cat("No significant GO terms found\n\n")
}

# ============================
# PART 6: KEGG Pathway Analysis
# ============================

cat("========================================\n")
cat("KEGG PATHWAY ENRICHMENT ANALYSIS\n")
cat("========================================\n\n")

cat("Performing KEGG pathway enrichment...\n")

kegg_results <- enrichKEGG(
    gene = all_genes,
    organism = "hsa",  # Homo sapiens
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
)

# Convert KEGG IDs to readable pathway names
kegg_results <- setReadable(kegg_results, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

cat("Found", nrow(kegg_results), "significant KEGG pathways\n")

# Save KEGG results
write.csv(
    as.data.frame(kegg_results),
    file = file.path(enrichment_dir, "KEGG_enrichment.csv")
)

# Visualize KEGG results
if (nrow(kegg_results) > 0) {
    cat("Generating KEGG enrichment plots...\n")
    
    # Dotplot
    pdf(file.path(enrichment_dir, "KEGG_dotplot.pdf"), width = 10, height = 8)
    print(dotplot(kegg_results, showCategory = 20, title = "KEGG Pathway Enrichment - Top 20 Pathways"))
    dev.off()
    
    # Barplot
    pdf(file.path(enrichment_dir, "KEGG_barplot.pdf"), width = 10, height = 8)
    print(barplot(kegg_results, showCategory = 20, title = "KEGG Pathway Enrichment - Bar Plot"))
    dev.off()
    
    cat("✓ KEGG enrichment plots saved\n\n")
} else {
    cat("No significant KEGG pathways found\n\n")
}

# ============================
# PART 7: Comparison Analysis
# ============================

cat("========================================\n")
cat("COMPARATIVE ANALYSIS: UP vs DOWN GENES\n")
cat("========================================\n\n")

if (length(upregulated_genes) > 0) {
    cat("GO enrichment for upregulated genes...\n")
    go_up <- enrichGO(
        gene = upregulated_genes,
        OrgDb = org.Hs.eg.db,
        ont = "BP",
        pvalueCutoff = 0.05,
        readable = TRUE
    )
    
    write.csv(
        as.data.frame(go_up),
        file = file.path(enrichment_dir, "GO_enrichment_upregulated.csv")
    )
    cat("  Found", nrow(go_up), "GO terms for upregulated genes\n")
    
    if (nrow(go_up) > 0) {
        pdf(file.path(enrichment_dir, "GO_upregulated_dotplot.pdf"), width = 10, height = 8)
        print(dotplot(go_up, showCategory = 15, title = "GO Enrichment - Upregulated Genes"))
        dev.off()
    }
}

if (length(downregulated_genes) > 0) {
    cat("GO enrichment for downregulated genes...\n")
    go_down <- enrichGO(
        gene = downregulated_genes,
        OrgDb = org.Hs.eg.db,
        ont = "BP",
        pvalueCutoff = 0.05,
        readable = TRUE
    )
    
    write.csv(
        as.data.frame(go_down),
        file = file.path(enrichment_dir, "GO_enrichment_downregulated.csv")
    )
    cat("  Found", nrow(go_down), "GO terms for downregulated genes\n")
    
    if (nrow(go_down) > 0) {
        pdf(file.path(enrichment_dir, "GO_downregulated_dotplot.pdf"), width = 10, height = 8)
        print(dotplot(go_down, showCategory = 15, title = "GO Enrichment - Downregulated Genes"))
        dev.off()
    }
}

cat("\n")

# ============================
# PART 8: Create Summary Report
# ============================

cat("Creating summary report...\n\n")

summary_text <- sprintf(
    "FUNCTIONAL ENRICHMENT ANALYSIS SUMMARY\n%s\n\nGene Statistics:\n  Total genes analyzed: %d\n  Upregulated: %d\n  Downregulated: %d\n\nGene Ontology Results:\n  Significant BP terms: %d\n\nKEGG Pathway Results:\n  Significant pathways: %d\n\nOutput Files:\n  - GO_enrichment_BP.csv\n  - KEGG_enrichment.csv\n  - GO_enrichment_upregulated.csv\n  - GO_enrichment_downregulated.csv\n  - Visualization plots (PDF)\n",
    strrep("=", 50),
    length(all_genes),
    length(upregulated_genes),
    length(downregulated_genes),
    nrow(go_results),
    nrow(kegg_results)
)

cat(summary_text)

# Save summary report
writeLines(
    summary_text,
    file.path(enrichment_dir, "enrichment_summary.txt")
)

# ============================
# PART 9: Display Top Results
# ============================

cat("\nTop 10 GO Terms (by adjusted p-value):\n")
cat(strrep("-", 80), "\n")
if (nrow(go_results) > 0) {
    top_go <- head(go_results, 10)[, c("Description", "p.adjust", "Count")]
    print(as.data.frame(top_go))
} else {
    cat("No GO terms found\n")
}

cat("\n")
cat("Top 10 KEGG Pathways (by adjusted p-value):\n")
cat(strrep("-", 80), "\n")
if (nrow(kegg_results) > 0) {
    top_kegg <- head(kegg_results, 10)[, c("Description", "p.adjust", "Count")]
    print(as.data.frame(top_kegg))
} else {
    cat("No KEGG pathways found\n")
}

cat("\n")
cat("========================================\n")
cat("Enrichment analysis completed!\n")
cat("Results saved to:", enrichment_dir, "\n")
cat("========================================\n")
