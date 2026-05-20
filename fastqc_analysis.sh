#!/bin/bash

################################################################################
# FastQC Quality Control Analysis Script
# Performs quality control checks on raw FASTQ files
# Prerequisites: Load environment variables with: source .env
################################################################################

set -e  # Exit on error

# Source environment variables
if [ ! -f ".env" ]; then
    echo "Error: .env file not found. Please run setup.sh first."
    exit 1
fi
source .env

echo "=========================================="
echo "STEP 01: FastQC Quality Control Analysis"
echo "=========================================="
echo ""

# Color output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

# Check if FASTQ files exist
if [ -z "$(ls -A $FASTQ 2>/dev/null)" ]; then
    echo -e "${RED}Error: No FASTQ files found in $FASTQ${NC}"
    echo "Please download FASTQ files first using: fasterq-dump <SRR_ID> -O \$FASTQ"
    exit 1
fi

echo -e "${BLUE}Input directory: ${NC}$FASTQ"
echo -e "${BLUE}Output directory: ${NC}$FASTQC_Results"
echo ""

# Count FASTQ files
FASTQ_COUNT=$(ls $FASTQ/*.fastq 2>/dev/null | wc -l)
echo -e "${BLUE}Found ${YELLOW}${FASTQ_COUNT}${BLUE} FASTQ files${NC}"
echo ""

# Step 1: Run FastQC on all FASTQ files
echo -e "${BLUE}[SUBSTEP 1]${NC} Running FastQC on all samples..."
echo "Command: fastqc $FASTQ/*.fastq -o $FASTQC_Results -t 10"
echo ""

fastqc $FASTQ/*.fastq -o $FASTQC_Results -t 10

if [ $? -eq 0 ]; then
    echo -e "${GREEN}✓ FastQC analysis completed successfully${NC}"
    
    # Count generated reports
    HTML_COUNT=$(ls $FASTQC_Results/*.html 2>/dev/null | wc -l)
    echo "  Generated $HTML_COUNT FastQC reports"
else
    echo -e "${RED}✗ FastQC analysis failed${NC}"
    exit 1
fi

echo ""

# Step 2: Run MultiQC to aggregate results
echo -e "${BLUE}[SUBSTEP 2]${NC} Aggregating results with MultiQC..."
echo "Command: multiqc $FASTQC_Results -o $FASTQC_Results"
echo ""

multiqc $FASTQC_Results -o $FASTQC_Results --force

if [ $? -eq 0 ]; then
    echo -e "${GREEN}✓ MultiQC report generated successfully${NC}"
    echo "  Report: $FASTQC_Results/multiqc_report.html"
else
    echo -e "${RED}✗ MultiQC report generation failed${NC}"
    exit 1
fi

echo ""

# Step 3: Summary statistics
echo -e "${BLUE}[SUBSTEP 3]${NC} Quality Control Summary"
echo "=========================================="

# Extract key statistics from FastQC output
echo ""
echo "Files analyzed:"
ls -1 $FASTQ/*.fastq | while read f; do
    basename "$f"
done | head -5
if [ $(ls $FASTQ/*.fastq | wc -l) -gt 5 ]; then
    echo "... and $(($(ls $FASTQ/*.fastq | wc -l) - 5)) more files"
fi

echo ""
echo -e "${YELLOW}Next step: Review the MultiQC report${NC}"
echo "  Open in browser: $FASTQC_Results/multiqc_report.html"
echo ""
echo -e "${YELLOW}If quality is acceptable, proceed with read trimming:${NC}"
echo "  bash \$SCRIPTS/trimming.sh"
echo ""

echo -e "${GREEN}=========================================="
echo "FastQC Analysis Complete!"
echo "==========================================${NC}"
