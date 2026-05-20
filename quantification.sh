#!/bin/bash

################################################################################
# Gene Quantification Script
# Uses featureCounts to count reads per gene
# Prerequisites: Load environment variables with: source .env
################################################################################

set -e

# Source environment variables
if [ ! -f ".env" ]; then
    echo "Error: .env file not found. Please run setup.sh first."
    exit 1
fi
source .env

echo "=========================================="
echo "STEP 04: Gene Quantification with featureCounts"
echo "=========================================="
echo ""

# Color output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

# Configuration
GTF_FILE="$REFERENCE/Homo_sapiens.GRCh38.109.gtf"
STRAND_SPECIFIC=0  # 0=unstranded, 1=forward, 2=reverse
THREADS=8
MIN_QUALITY=10

echo -e "${BLUE}Input directory: ${NC}$BAM"
echo -e "${BLUE}Output directory: ${NC}$COUNTS"
echo -e "${BLUE}GTF file: ${NC}$GTF_FILE"
echo ""
echo "Quantification parameters:"
echo "  Strand specific: $STRAND_SPECIFIC (0=unstranded)"
echo "  Threads: $THREADS"
echo "  Minimum mapping quality: $MIN_QUALITY"
echo ""

# Check if GTF file exists
if [ ! -f "$GTF_FILE" ]; then
    echo -e "${RED}Error: GTF annotation file not found at $GTF_FILE${NC}"
    echo -e "${YELLOW}Please download the annotation file first:${NC}"
    echo "  wget -P \$REFERENCE https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz"
    echo "  gunzip \$REFERENCE/Homo_sapiens.GRCh38.109.gtf.gz"
    exit 1
fi

# Check if BAM files exist
if [ -z "$(ls -A $BAM/*.bam 2>/dev/null | grep -v ".bam.bai")" ]; then
    echo -e "${RED}Error: No BAM files found in $BAM${NC}"
    exit 1
fi

# Get list of BAM files
BAM_FILES=$(ls $BAM/*.bam | grep -v ".bam.bai" | sort)
SAMPLE_COUNT=$(echo "$BAM_FILES" | wc -l)

echo -e "${BLUE}Found ${YELLOW}${SAMPLE_COUNT}${BLUE} BAM files${NC}"
echo ""

# Create a list of BAM files for featureCounts
BAM_LIST=$(echo "$BAM_FILES" | tr '\n' ' ')

# Run featureCounts
echo -e "${BLUE}Running featureCounts...${NC}"
echo "Command: featureCounts -p -B -C -s $STRAND_SPECIFIC -T $THREADS -a $GTF_FILE -o $COUNTS/count_matrix.txt $BAM_LIST"
echo ""

featureCounts \
    -p \
    -B \
    -C \
    -s $STRAND_SPECIFIC \
    -T $THREADS \
    -a "$GTF_FILE" \
    -o "$COUNTS/count_matrix.txt" \
    $BAM_FILES \
    2>&1 | tee "$COUNTS/featureCounts.log"

if [ $? -eq 0 ]; then
    echo -e "${GREEN}✓ featureCounts completed successfully${NC}"
else
    echo -e "${RED}✗ featureCounts failed${NC}"
    exit 1
fi

echo ""

# Post-processing: Create a summary table
echo -e "${BLUE}[POSTPROCESSING]${NC} Creating count summary..."

# Extract summary information
if [ -f "$COUNTS/count_matrix.txt.summary" ]; then
    echo "featureCounts Summary Statistics:"
    echo "=================================="
    cat "$COUNTS/count_matrix.txt.summary"
    echo ""
fi

# Create a cleaner count matrix by removing metadata columns
echo "Creating cleaned count matrix..."
cat "$COUNTS/count_matrix.txt" | cut -f1,7- > "$COUNTS/count_matrix_clean.txt"

if [ $? -eq 0 ]; then
    echo -e "${GREEN}✓ Cleaned count matrix created${NC}"
    echo "  File: $COUNTS/count_matrix_clean.txt"
else
    echo -e "${RED}✗ Failed to create cleaned count matrix${NC}"
fi

echo ""

# Generate statistics
echo -e "${BLUE}[STATISTICS]${NC} Count Matrix Statistics"
echo "=========================================="
echo ""

# Count total genes and samples
N_GENES=$(wc -l < "$COUNTS/count_matrix.txt" | awk '{print $1-1}')
N_SAMPLES=$(head -1 "$COUNTS/count_matrix.txt" | awk '{print NF-6}')

echo "Number of genes: $N_GENES"
echo "Number of samples: $N_SAMPLES"
echo ""

# Show first few rows
echo "First 5 rows of count matrix:"
head -6 "$COUNTS/count_matrix_clean.txt" | column -t
echo ""

# Calculate mapping statistics per sample
echo "Read counts per sample:"
echo "Sample | Assigned Reads | Percentage"
echo "-------|----------------|----------"
tail -n +2 "$COUNTS/count_matrix.txt.summary" | head -$SAMPLE_COUNT | while read line; do
    SAMPLE=$(echo "$line" | awk '{print $1}' | sed 's|.*\/||' | sed 's|\.bam||')
    ASSIGNED=$(echo "$line" | awk '{print $2}')
    TOTAL=$(echo "$line" | awk '{for(i=2;i<=NF;i++) sum+=$i} END {print sum}')
    PCT=$(echo "scale=2; ($ASSIGNED / $TOTAL) * 100" | bc)
    printf "%-6s | %14s | %6.2f%%\n" "$SAMPLE" "$ASSIGNED" "$PCT"
done

echo ""
echo -e "${YELLOW}Next step: Perform differential expression analysis${NC}"
echo "  Rscript \$SCRIPTS/dge_analysis.R"
echo ""

echo -e "${GREEN}=========================================="
echo "Gene Quantification Complete!"
echo "==========================================${NC}"
