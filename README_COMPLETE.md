# Bulk RNA-Seq Analysis Pipeline

A comprehensive workflow for analyzing bulk RNA-sequencing data from raw FASTQ files to differential gene expression and functional enrichment analysis.

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Installation](#installation)
- [Project Structure](#project-structure)
- [Quick Start](#quick-start)
- [Workflow Steps](#workflow-steps)
- [Configuration](#configuration)
- [Output Files](#output-files)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [License](#license)

## Overview

This pipeline performs a complete bulk RNA-sequencing analysis workflow including:

1. **Quality Control**: Assessment of raw sequencing data using FastQC and MultiQC
2. **Read Processing**: Adapter and quality trimming with Trimmomatic
3. **Alignment**: Mapping reads to the reference genome using HISAT2
4. **Quantification**: Gene-level read counting with featureCounts
5. **Differential Expression**: Statistical analysis using DESeq2
6. **Functional Analysis**: Gene Ontology and KEGG pathway enrichment with clusterProfiler

## Features

✅ Automated environment setup with conda  
✅ Parallel processing support for faster analysis  
✅ Comprehensive quality control reports  
✅ Publication-ready visualizations  
✅ Functional enrichment analysis  
✅ Detailed configuration options  
✅ Error handling and logging  
✅ Sample metadata tracking  

## Installation

### Prerequisites

- Linux/macOS system
- Conda (Anaconda or Miniconda)
- Bash shell
- ~100-200 GB disk space (depending on sample count and reference genome)
- 16+ GB RAM recommended

### Quick Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/Bulk-RNA-seq.git
cd Bulk-RNA-seq

# Run setup script (creates conda environment and directory structure)
bash setup.sh

# Activate the conda environment
conda activate bulk_env

# Load environment variables
source BulkTranscriptomics/.env
```

### Manual Installation

```bash
# Create conda environment from file
conda env create -f environment.yml

# Activate environment
conda activate bulk_env
```

## Project Structure

```
BulkTranscriptomics/
├── FASTQ/                 # Raw sequencing data (FASTQ files)
├── FASTQC_Results/        # Quality control reports
├── TRIMMED/               # Trimmed FASTQ files
├── ALIGN/                 # Alignment intermediate files
├── BAM/                   # Sorted BAM files and indices
├── REFERENCE/             # Reference genome and annotations
├── COUNTS/                # Read count matrices
├── RESULTS/               # Analysis results
│   └── functional_enrichment/
├── SCRIPTS/               # Analysis scripts
├── DATA/                  # Metadata and sample information
└── .env                   # Environment variables
```

## Quick Start

### 1. Prepare Input Data

**Download SRA data:**
```bash
# Download individual samples
fasterq-dump SRR11262284 --threads 10 --progress --split-files -O $FASTQ

# Or use GNU parallel for faster downloads
parallel -j 2 'fasterq-dump {} --threads 4 --split-files --progress -O $FASTQ' \
    ::: SRR11262284 SRR11262285 SRR11262286
```

**Or copy local FASTQ files:**
```bash
cp /path/to/local/fastq/*.fastq $FASTQ/
```

### 2. Download Reference Genome

```bash
# Download HISAT2 index
wget -P $REFERENCE https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
tar -xzf $REFERENCE/grch38_genome.tar.gz -C $REFERENCE

# Download GTF annotation
wget -P $REFERENCE https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz
gunzip $REFERENCE/Homo_sapiens.GRCh38.109.gtf.gz
```

### 3. Prepare Sample Metadata

Edit `DATA/sample_metadata.csv` with your sample information:

```csv
Sample_ID,SRR_Accession,Condition,Replicate,Tissue,Disease_Status
SAMPLE1,SRR11262284,Control,1,Tissue,Normal
SAMPLE2,SRR11262285,Treatment,1,Tissue,Treated
...
```

### 4. Run the Pipeline

```bash
# Step 1: Quality Control
bash $SCRIPTS/fastqc_analysis.sh

# Step 2: Trimming
bash $SCRIPTS/trimming.sh

# Step 3: Alignment
bash $SCRIPTS/alignment.sh

# Step 4: Quantification
bash $SCRIPTS/quantification.sh

# Step 5: Differential Expression Analysis
Rscript $SCRIPTS/dge_analysis.R

# Step 6: Functional Enrichment
Rscript $SCRIPTS/enrichment_analysis.R
```

## Workflow Steps

### Step 1: Quality Control (FastQC)

```bash
bash $SCRIPTS/fastqc_analysis.sh
```

**Outputs:**
- Individual FastQC HTML reports
- MultiQC aggregated report
- Quality metrics

**What to look for:**
- Per-base quality scores > Q20
- GC content within expected range
- No adapter contamination
- Acceptable duplication levels

### Step 2: Read Trimming (Trimmomatic)

```bash
bash $SCRIPTS/trimming.sh
```

**Parameters (adjustable in config.yaml):**
- Leading/trailing quality: 3
- Sliding window: 4:15
- Minimum length: 36 bp

**Outputs:**
- Paired-end trimmed FASTQ files
- Trimming statistics and logs

### Step 3: Read Alignment (HISAT2)

```bash
bash $SCRIPTS/alignment.sh
```

**Parameters:**
- Max intron length: 500,000 bp
- Threads: 16

**Outputs:**
- Sorted BAM files
- BAM indices (.bai)
- Alignment statistics (flagstat)

### Step 4: Gene Quantification (featureCounts)

```bash
bash $SCRIPTS/quantification.sh
```

**Parameters:**
- Feature type: exon
- ID type: gene_id
- Strand specific: No

**Outputs:**
- Count matrix (counts × samples)
- Assignment statistics
- Summary file

### Step 5: Differential Expression Analysis (DESeq2)

```bash
Rscript $SCRIPTS/dge_analysis.R
```

**Outputs:**
- DEG results table (all genes)
- Filtered results (padj < 0.05)
- Visualizations:
  - PCA plot
  - MA plot
  - Volcano plot
  - Heatmap
  - Expression boxplots

### Step 6: Functional Enrichment (clusterProfiler)

```bash
Rscript $SCRIPTS/enrichment_analysis.R
```

**Outputs:**
- GO enrichment results (BP, MF, CC)
- KEGG pathway analysis
- Enrichment visualizations:
  - Dotplots
  - Barplots
  - Comparative analysis

## Configuration

Edit `config.yaml` to customize pipeline parameters:

```yaml
# Reference genome settings
reference:
  genome_version: "GRCh38"
  gtf_version: "Ensembl v109"

# Quality trimming thresholds
trimming:
  leading_quality: 3
  trailing_quality: 3
  min_length_after_trim: 36

# Alignment parameters
alignment:
  max_intron_length: 500000
  threads: 16

# Statistical analysis
deseq2:
  alpha: 0.05
  lfc_threshold: 0.5

# Enrichment analysis
enrichment:
  padj_cutoff: 0.05
  ontology: ["BP", "MF", "CC"]
```

## Output Files

### Quality Control
- `FASTQC_Results/multiqc_report.html` - Interactive QC summary

### Alignment
- `BAM/*.bam` - Sorted BAM files
- `BAM/*.bam.bai` - BAM indices
- `BAM/*.flagstat.txt` - Alignment statistics

### Quantification
- `COUNTS/count_matrix.txt` - Full count matrix
- `COUNTS/count_matrix_clean.txt` - Cleaned matrix (genes × samples)
- `COUNTS/featureCounts.log` - Counting log

### Differential Expression
- `RESULTS/DEG_results_all.csv` - All genes with statistics
- `RESULTS/DEG_results_significant.csv` - Significant genes (padj < 0.05)
- `RESULTS/*.pdf` - Plots (PCA, MA, Volcano, Heatmap)

### Functional Enrichment
- `RESULTS/functional_enrichment/GO_enrichment_BP.csv` - GO terms
- `RESULTS/functional_enrichment/KEGG_enrichment.csv` - KEGG pathways
- `RESULTS/functional_enrichment/*.pdf` - Enrichment plots

## Troubleshooting

### Issue: "fasterq-dump not found"
**Solution:** Install sra-tools in your conda environment
```bash
conda install -c bioconda sra-tools=3.0.7
```

### Issue: "HISAT2 index not found"
**Solution:** Download and extract the reference genome index
```bash
wget -P $REFERENCE https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
tar -xzf $REFERENCE/grch38_genome.tar.gz -C $REFERENCE
```

### Issue: "Low mapping percentage"
**Possible causes:**
- Wrong reference genome version
- Adapter contamination (check FastQC report)
- Poor sample quality
- Biological contamination

**Solutions:**
1. Check FastQC report for adapter content
2. Increase Trimmomatic stringency
3. Verify reference genome matches library
4. Consider removing low-quality samples

### Issue: "featureCounts assignment too low"
**Possible causes:**
- Strand-specific sequencing not accounted for
- Mismatched GTF version
- Incorrect SAM flag handling

**Solutions:**
1. Check `-s` parameter in config.yaml (0=unstranded, 1=stranded)
2. Verify GTF matches reference genome version
3. Check BAM file quality and indexes

### Issue: "R package installation fails"
**Solution:** Install from BiocManager
```bash
# In R console
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2")
```

## Tips and Best Practices

1. **Always review QC reports** before proceeding to downstream analysis
2. **Keep raw data** - Don't delete original FASTQ files
3. **Document metadata** - Accurate sample information is critical
4. **Use consistent naming** - Sample IDs should be informative
5. **Monitor disk space** - Pipeline generates large intermediate files
6. **Backup results** - Important data should be backed up
7. **Test with small dataset** - Validate pipeline before full analysis

## References

- HISAT2: Kim D, Paggi JM, Park C, et al. Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. Nat Biotechnol. 2019.
- DESeq2: Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol. 2014.
- clusterProfiler: Yu G, Wang LG, Han Y, He QY. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS. 2012.

## License

This project is licensed under the MIT License - see LICENSE file for details.

## Contributing

Contributions are welcome! Please see CONTRIBUTING.md for guidelines.

## Citation

If you use this pipeline in your research, please cite:

```
[Your citation here]
```

## Support

For issues, questions, or suggestions:
- Open an issue on GitHub
- Check the FAQ in the wiki
- Contact: your.email@institution.org

## Acknowledgments

- Built with best practices from Bioconda, Bioconductor, and RNA-seq community
- Reference genome from Ensembl
- Tools: FastQC, HISAT2, Trimmomatic, DESeq2, clusterProfiler
