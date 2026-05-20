#!/bin/bash

################################################################################
# Read Alignment Script
# Uses HISAT2 to align trimmed reads to the reference genome
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
echo "STEP 03: Read Alignment with HISAT2"
echo "=========================================="
echo ""

# Color output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

# Configuration
HISAT2_INDEX="$REFERENCE/grch38_genome/genome"
THREADS=16
MAX_INTRON=500000

echo -e "${BLUE}Input directory: ${NC}$TRIMMED"
echo -e "${BLUE}Output directory: ${NC}$BAM"
echo -e "${BLUE}HISAT2 index: ${NC}$HISAT2_INDEX"
echo ""
echo "Alignment parameters:"
echo "  Threads: $THREADS"
echo "  Max intron length: $MAX_INTRON"
echo ""

# Check if HISAT2 index exists
if [ ! -f "${HISAT2_INDEX}.1.ht2" ]; then
    echo -e "${RED}Error: HISAT2 index not found at $HISAT2_INDEX${NC}"
    echo -e "${YELLOW}Please download the reference genome first:${NC}"
    echo "  wget -P \$REFERENCE https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz"
    echo "  tar -xzf \$REFERENCE/grch38_genome.tar.gz -C \$REFERENCE"
    exit 1
fi

# Check if trimmed files exist
if [ -z "$(ls -A $TRIMMED 2>/dev/null)" ]; then
    echo -e "${RED}Error: No trimmed FASTQ files found in $TRIMMED${NC}"
    exit 1
fi

# Get list of samples
SAMPLES=$(ls $TRIMMED/*_1_paired.fastq 2>/dev/null | sed 's/_1_paired\.fastq$//' | xargs -n1 basename)

if [ -z "$SAMPLES" ]; then
    echo -e "${RED}Error: No paired trimmed FASTQ files found${NC}"
    exit 1
fi

SAMPLE_COUNT=$(echo "$SAMPLES" | wc -l)
echo -e "${BLUE}Found ${YELLOW}${SAMPLE_COUNT}${BLUE} samples to align${NC}"
echo ""

# Counter for progress
COUNTER=0

# Process each sample
for SAMPLE in $SAMPLES; do
    COUNTER=$((COUNTER + 1))
    
    R1="$TRIMMED/${SAMPLE}_1_paired.fastq"
    R2="$TRIMMED/${SAMPLE}_2_paired.fastq"
    SAM="$ALIGN/${SAMPLE}.sam"
    BAM_UNSORTED="$ALIGN/${SAMPLE}.unsorted.bam"
    BAM_SORTED="$BAM/${SAMPLE}.bam"
    BAM_INDEX="$BAM/${SAMPLE}.bam.bai"
    
    echo -e "${BLUE}[$COUNTER/$SAMPLE_COUNT]${NC} Aligning: $SAMPLE"
    
    # Check if input files exist
    if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
        echo -e "${YELLOW}Warning: Trimmed files not found for $SAMPLE, skipping...${NC}"
        continue
    fi
    
    # Step 1: Alignment with HISAT2
    echo "  Step 1: HISAT2 alignment..."
    hisat2 \
        -x "$HISAT2_INDEX" \
        -1 "$R1" \
        -2 "$R2" \
        --max-intronlen $MAX_INTRON \
        --threads $THREADS \
        --met-file "$ALIGN/${SAMPLE}.metrics.txt" \
        -S "$SAM" \
        2>&1 | tee "$ALIGN/${SAMPLE}_alignment.log"
    
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}✓ HISAT2 alignment completed${NC}"
    else
        echo -e "${RED}✗ HISAT2 alignment failed${NC}"
        rm -f "$SAM"
        exit 1
    fi
    
    # Step 2: Convert SAM to BAM
    echo "  Step 2: Converting SAM to BAM..."
    samtools view -@ $THREADS -b -S "$SAM" > "$BAM_UNSORTED"
    
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}✓ SAM to BAM conversion completed${NC}"
    else
        echo -e "${RED}✗ SAM to BAM conversion failed${NC}"
        exit 1
    fi
    
    # Step 3: Sort BAM file
    echo "  Step 3: Sorting BAM file..."
    samtools sort -@ $THREADS "$BAM_UNSORTED" -o "$BAM_SORTED"
    
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}✓ BAM sorting completed${NC}"
    else
        echo -e "${RED}✗ BAM sorting failed${NC}"
        exit 1
    fi
    
    # Step 4: Index BAM file
    echo "  Step 4: Indexing BAM file..."
    samtools index -@ $THREADS "$BAM_SORTED"
    
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}✓ BAM indexing completed${NC}"
    else
        echo -e "${RED}✗ BAM indexing failed${NC}"
        exit 1
    fi
    
    # Step 5: Get alignment statistics
    echo "  Step 5: Calculating alignment statistics..."
    samtools flagstat "$BAM_SORTED" > "$BAM/${SAMPLE}.flagstat.txt"
    
    # Clean up intermediate files
    rm -f "$SAM" "$BAM_UNSORTED"
    
    echo -e "${GREEN}✓ Completed for $SAMPLE${NC}"
    echo ""
done

echo ""

# Summary statistics
echo -e "${BLUE}[SUMMARY]${NC} Alignment Statistics"
echo "=========================================="
echo ""
echo -e "${YELLOW}Flagstat results:${NC}"
echo ""
for FLAGSTAT in $BAM/*.flagstat.txt; do
    SAMPLE=$(basename "$FLAGSTAT" .flagstat.txt)
    echo "Sample: $SAMPLE"
    head -3 "$FLAGSTAT"
    echo ""
done

echo ""
echo -e "${YELLOW}Next step: Quantify gene expression${NC}"
echo "  bash \$SCRIPTS/quantification.sh"
echo ""

echo -e "${GREEN}=========================================="
echo "Alignment Complete!"
echo "==========================================${NC}"
