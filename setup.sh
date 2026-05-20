#!/bin/bash

################################################################################
# RNA-Seq Bulk Transcriptomics Pipeline Setup Script
# This script sets up the conda environment and creates the required directory 
# structure for bulk RNA-seq analysis
################################################################################

set -e  # Exit on error

echo "=========================================="
echo "Bulk RNA-Seq Pipeline Setup"
echo "=========================================="
echo ""

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Step 1: Create Conda Environment
echo -e "${BLUE}[STEP 1]${NC} Creating Conda environment 'bulk_env'..."
conda create -n bulk_env -c conda-forge -c bioconda r-base=4.3.3 r-essentials -y

if [ $? -eq 0 ]; then
    echo -e "${GREEN}✓ Conda environment created successfully${NC}"
else
    echo -e "${YELLOW}✗ Conda environment creation failed or already exists${NC}"
fi

echo ""

# Step 2: Install bioinformatics tools
echo -e "${BLUE}[STEP 2]${NC} Installing bioinformatics tools..."
conda install -n bulk_env -c bioconda fastqc multiqc trimmomatic hisat2 samtools subread sra-tools=3.0.7 -y

if [ $? -eq 0 ]; then
    echo -e "${GREEN}✓ Bioinformatics tools installed successfully${NC}"
else
    echo -e "${YELLOW}✗ Tool installation failed${NC}"
fi

echo ""

# Step 3: Create project directory structure
echo -e "${BLUE}[STEP 3]${NC} Creating project directory structure..."

PROJECT_DIR="BulkTranscriptomics"

if [ ! -d "$PROJECT_DIR" ]; then
    mkdir -p "$PROJECT_DIR"/{FASTQ,FASTQC_Results,TRIMMED,ALIGN,BAM,REFERENCE,COUNTS,RESULTS,SCRIPTS,DATA}
    echo -e "${GREEN}✓ Project directories created:${NC}"
    echo "   - $PROJECT_DIR/FASTQ (raw sequencing data)"
    echo "   - $PROJECT_DIR/FASTQC_Results (quality control outputs)"
    echo "   - $PROJECT_DIR/TRIMMED (trimmed reads)"
    echo "   - $PROJECT_DIR/ALIGN (alignment intermediate files)"
    echo "   - $PROJECT_DIR/BAM (sorted BAM files)"
    echo "   - $PROJECT_DIR/REFERENCE (reference genome and annotations)"
    echo "   - $PROJECT_DIR/COUNTS (read count matrices)"
    echo "   - $PROJECT_DIR/RESULTS (analysis results)"
    echo "   - $PROJECT_DIR/SCRIPTS (analysis scripts)"
    echo "   - $PROJECT_DIR/DATA (metadata and supplementary data)"
else
    echo -e "${YELLOW}✓ Project directory already exists${NC}"
fi

echo ""

# Step 4: Export environment variables
echo -e "${BLUE}[STEP 4]${NC} Setting up environment variables..."

# Get absolute path of project directory
PROJECT_PATH="$(cd "$PROJECT_DIR" 2>/dev/null && pwd)" || PROJECT_PATH="$(pwd)/$PROJECT_DIR"

# Create .env file with environment variables
cat > "$PROJECT_DIR/.env" << EOF
#!/bin/bash
# Environment variables for Bulk RNA-Seq Pipeline
# Source this file with: source .env

export BULK_PROJECT="$PROJECT_PATH"
export FASTQ="\$BULK_PROJECT/FASTQ"
export FASTQC_Results="\$BULK_PROJECT/FASTQC_Results"
export TRIMMED="\$BULK_PROJECT/TRIMMED"
export ALIGN="\$BULK_PROJECT/ALIGN"
export BAM="\$BULK_PROJECT/BAM"
export REFERENCE="\$BULK_PROJECT/REFERENCE"
export COUNTS="\$BULK_PROJECT/COUNTS"
export RESULTS="\$BULK_PROJECT/RESULTS"
export SCRIPTS="\$BULK_PROJECT/SCRIPTS"
export DATA="\$BULK_PROJECT/DATA"

echo "Environment variables loaded from $PROJECT_PATH"
EOF

echo -e "${GREEN}✓ Environment variables file created (.env)${NC}"
echo "   Source with: source $PROJECT_DIR/.env"

echo ""

# Step 5: Copy scripts to SCRIPTS folder
echo -e "${BLUE}[STEP 5]${NC} Note on scripts..."
echo "Copy the provided analysis scripts to $PROJECT_DIR/SCRIPTS:"
echo "  - fastqc_analysis.sh"
echo "  - trimming.sh"
echo "  - alignment.sh"
echo "  - quantification.sh"
echo "  - dge_analysis.R"
echo "  - enrichment_analysis.R"

echo ""

# Step 6: Instructions
echo -e "${BLUE}[STEP 6]${NC} Next steps:"
echo "1. Activate the conda environment:"
echo "   ${YELLOW}conda activate bulk_env${NC}"
echo ""
echo "2. Load environment variables:"
echo "   ${YELLOW}source $PROJECT_DIR/.env${NC}"
echo ""
echo "3. Download reference genome and annotations:"
echo "   ${YELLOW}wget -P \$REFERENCE https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz${NC}"
echo "   ${YELLOW}tar -xzf \$REFERENCE/grch38_genome.tar.gz -C \$REFERENCE${NC}"
echo "   ${YELLOW}wget -P \$REFERENCE https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz${NC}"
echo "   ${YELLOW}gunzip \$REFERENCE/Homo_sapiens.GRCh38.109.gtf.gz${NC}"
echo ""
echo "4. Download SRA data:"
echo "   ${YELLOW}fasterq-dump <SRR_ID> --threads 10 --progress --split-files -O \$FASTQ${NC}"
echo ""
echo "5. Prepare sample metadata in $PROJECT_DIR/DATA/sample_metadata.csv"
echo "6. Run analysis scripts in order:"
echo "   ${YELLOW}bash \$SCRIPTS/fastqc_analysis.sh${NC}"
echo "   ${YELLOW}bash \$SCRIPTS/trimming.sh${NC}"
echo "   ${YELLOW}bash \$SCRIPTS/alignment.sh${NC}"
echo "   ${YELLOW}bash \$SCRIPTS/quantification.sh${NC}"

echo ""
echo -e "${GREEN}=========================================="
echo "Setup completed successfully!"
echo "==========================================${NC}"
