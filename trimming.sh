#!/bin/bash

################################################################################
# Read Trimming Script
# Uses Trimmomatic to remove adapters and low-quality sequences
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
echo "STEP 02: Read Trimming with Trimmomatic"
echo "=========================================="
echo ""

# Color output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

# Check input files
if [ -z "$(ls -A $FASTQ 2>/dev/null)" ]; then
    echo -e "${RED}Error: No FASTQ files found in $FASTQ${NC}"
    exit 1
fi

# Configuration
ADAPTERS="TruSeq3-PE.fa"
LEADING=3
TRAILING=3
SLIDINGWINDOW="4:15"
MINLEN=36
THREADS=8

echo -e "${BLUE}Input directory: ${NC}$FASTQ"
echo -e "${BLUE}Output directory: ${NC}$TRIMMED"
echo ""
echo "Trimming parameters:"
echo "  Leading quality: $LEADING"
echo "  Trailing quality: $TRAILING"
echo "  Sliding window: $SLIDINGWINDOW"
echo "  Minimum length: $MINLEN"
echo "  Threads: $THREADS"
echo ""

# Get list of paired-end FASTQ files
# Assumes files are named as: sample_1.fastq and sample_2.fastq
SAMPLES=$(ls $FASTQ/*_1.fastq 2>/dev/null | sed 's/_1\.fastq$//' | xargs -n1 basename)

if [ -z "$SAMPLES" ]; then
    echo -e "${RED}Error: No paired-end FASTQ files found (expected format: *_1.fastq)${NC}"
    echo "Files found: $(ls $FASTQ/*.fastq | head -5)"
    exit 1
fi

SAMPLE_COUNT=$(echo "$SAMPLES" | wc -l)
echo -e "${BLUE}Found ${YELLOW}${SAMPLE_COUNT}${BLUE} sample pairs${NC}"
echo ""

# Counter for progress
COUNTER=0

# Process each sample
for SAMPLE in $SAMPLES; do
    COUNTER=$((COUNTER + 1))
    
    R1_INPUT="$FASTQ/${SAMPLE}_1.fastq"
    R2_INPUT="$FASTQ/${SAMPLE}_2.fastq"
    
    R1_PAIRED="$TRIMMED/${SAMPLE}_1_paired.fastq"
    R1_UNPAIRED="$TRIMMED/${SAMPLE}_1_unpaired.fastq"
    R2_PAIRED="$TRIMMED/${SAMPLE}_2_paired.fastq"
    R2_UNPAIRED="$TRIMMED/${SAMPLE}_2_unpaired.fastq"
    
    echo -e "${BLUE}[$COUNTER/$SAMPLE_COUNT]${NC} Processing: $SAMPLE"
    
    # Check if input files exist
    if [ ! -f "$R1_INPUT" ] || [ ! -f "$R2_INPUT" ]; then
        echo -e "${YELLOW}Warning: Paired files not found for $SAMPLE, skipping...${NC}"
        continue
    fi
    
    # Run Trimmomatic
    trimmomatic PE \
        -threads $THREADS \
        -phred33 \
        "$R1_INPUT" "$R2_INPUT" \
        "$R1_PAIRED" "$R1_UNPAIRED" \
        "$R2_PAIRED" "$R2_UNPAIRED" \
        LEADING:$LEADING \
        TRAILING:$TRAILING \
        SLIDINGWINDOW:$SLIDINGWINDOW \
        MINLEN:$MINLEN \
        2>&1 | tee "$TRIMMED/${SAMPLE}_trimming.log"
    
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}✓ Trimming completed for $SAMPLE${NC}"
    else
        echo -e "${RED}✗ Trimming failed for $SAMPLE${NC}"
        exit 1
    fi
    
    echo ""
done

echo ""

# Summary statistics
echo -e "${BLUE}[SUMMARY]${NC} Trimming Statistics"
echo "=========================================="
echo ""
echo "Original read files:"
ls -lh $FASTQ/*_1.fastq | awk '{print $9, $5}' | head -3
echo ""
echo "Trimmed paired files:"
ls -lh $TRIMMED/*_1_paired.fastq | awk '{print $9, $5}' | head -3
echo ""

# Count reads before and after
echo -e "${YELLOW}Reads per sample:${NC}"
echo "Sample | Original R1 | Trimmed R1 | Reduction (%)"
echo "-------|-------------|-----------|---------------"

for SAMPLE in $SAMPLES; do
    R1_INPUT="$FASTQ/${SAMPLE}_1.fastq"
    R1_PAIRED="$TRIMMED/${SAMPLE}_1_paired.fastq"
    
    if [ -f "$R1_INPUT" ] && [ -f "$R1_PAIRED" ]; then
        ORIG=$(wc -l < "$R1_INPUT" | awk '{print $1/4}' | bc)
        TRIM=$(wc -l < "$R1_PAIRED" | awk '{print $1/4}' | bc)
        REDUCTION=$(echo "scale=2; (($ORIG - $TRIM) / $ORIG) * 100" | bc)
        printf "%-6s | %11d | %9d | %6.2f%%\n" "$SAMPLE" "$ORIG" "$TRIM" "$REDUCTION"
    fi
done | head -5

echo ""
echo -e "${YELLOW}Next step: Align trimmed reads to reference genome${NC}"
echo "  bash \$SCRIPTS/alignment.sh"
echo ""

echo -e "${GREEN}=========================================="
echo "Trimming Complete!"
echo "==========================================${NC}"
