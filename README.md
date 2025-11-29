# Bulk-RNA-seq
RNA sequencing (RNA-Seq) is a powerful and widely used technique for studying the
transcriptome, the complete set of RNA transcripts produced by the genome under specific
conditions or in a specific cell type. Unlike traditional microarrays, RNA-Seq provides a high
resolution, quantitative, and unbiased view of gene expression, enabling researchers to
identify differentially expressed genes, discover novel transcripts, detect alternative splicing
events, and quantify gene expression across the entire genome
<img width="885" height="940" alt="image" src="https://github.com/user-attachments/assets/8a5dd817-70b7-403f-a895-7088be3a17b8" />

The RNA-Seq data analysis process can be broadly divided into two main stages: the
transformation of raw sequencing reads into a count matrix, and the downstream statistical
and functional interpretation of gene expression data.

It begins with quality control of raw sequencing reads in FASTQ format. Tools like FastQC
are commonly used to evaluate sequencing quality, identify adapter contamination, and
detect technical biases. If necessary, low-quality reads and adapter sequences are trimmed
using programs such as Trimmomatic or Cutadapt. The cleaned reads are then aligned to a
reference genome or transcriptome using aligners like STAR or HISAT2, or alternatively
pseudoaligned using tools such as Salmon or Kallisto. Following alignment, gene- or
transcript-level quantification is performed with tools like featureCounts or HTSeq, resulting
in a count matrix that represents the number of reads mapped to each gene across all
samples.

The downstream analysis of this count data to derive meaningful biological insights.
Differential gene expression (DGE) analysis is typically performed using statistical
frameworks such as DESeq2, edgeR, or limma to identify genes that show significant
changes in expression between conditions or experimental groups. The results are often
visualized using principal component analysis (PCA), heatmaps, and volcano plots to
explore sample relationships and highlight key regulatory genes. To interpret the biological
relevance of differentially expressed genes, functional enrichment analysis is carried out
using gene ontology (GO) and pathway analysis tools such as clusterProfiler or GSEA. This
integrative workflow allows researchers to move from raw sequencing data to a deeper
understanding of cellular states, molecular functions, and regulatory mechanisms underlying
the biological system under study
