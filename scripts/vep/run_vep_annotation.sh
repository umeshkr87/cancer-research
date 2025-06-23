#!/bin/bash
#SBATCH --job-name=vep_annotation
#SBATCH --partition=IllinoisComputes
#SBATCH --account=aa107-ic
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=32G
#SBATCH --output=/u/aa107/uiuc-cancer-research/scripts/vep/vep_%j.out
#SBATCH --error=/u/aa107/uiuc-cancer-research/scripts/vep/vep_%j.err

# VEP Annotation Script for Prostate Cancer Variants
# Author: TabNet Prostate Cancer Research Team
# Date: $(date '+%Y-%m-%d')

set -e  # Exit on any error

# === CONFIGURATION ===
PROJECT_DIR="/u/aa107/uiuc-cancer-research"
INPUT_VCF="${PROJECT_DIR}/data/processed/merged_vcf/merged_prostate_variants.vcf"
OUTPUT_DIR="${PROJECT_DIR}/data/processed/vep"
LOG_FILE="${OUTPUT_DIR}/vep_annotation.log"

echo "=== VEP Annotation Pipeline Started: $(date) ===" | tee $LOG_FILE

# === SETUP ENVIRONMENT ===
module load anaconda3
eval "$(conda shell.bash hook)"

# === FIX CORRUPTED ENVIRONMENT ===
ENV_PATH="$HOME/.conda/envs/tabnet-prostate"
if [ -d "$ENV_PATH" ] && [ ! -f "$ENV_PATH/conda-meta/history" ]; then
    echo "Fixing corrupted conda environment..." | tee -a $LOG_FILE
    rm -rf $ENV_PATH
fi

# === CREATE/ACTIVATE ENVIRONMENT ===
if [ ! -d "$ENV_PATH" ]; then
    echo "Creating tabnet-prostate environment..." | tee -a $LOG_FILE
    conda create -n tabnet-prostate python=3.11 -y
fi

conda activate tabnet-prostate

# === INSTALL VEP IF NEEDED ===
if ! command -v vep &> /dev/null; then
    echo "Installing VEP via conda..." | tee -a $LOG_FILE
    conda install -c bioconda ensembl-vep -y
    
    # Install VEP cache
    echo "Installing VEP cache..." | tee -a $LOG_FILE
    vep_install -a cf -s homo_sapiens -y GRCh38
fi

# === VALIDATE INPUT ===
echo "Validating input file..." | tee -a $LOG_FILE
if [ ! -f "$INPUT_VCF" ]; then
    echo "ERROR: Input VCF not found: $INPUT_VCF" | tee -a $LOG_FILE
    exit 1
fi

INPUT_COUNT=$(grep -v "^#" $INPUT_VCF | wc -l)
echo "Input variants: $INPUT_COUNT" | tee -a $LOG_FILE

# === RUN VEP ANNOTATION ===
echo "Running VEP annotation..." | tee -a $LOG_FILE

vep \
    --input_file $INPUT_VCF \
    --output_file ${OUTPUT_DIR}/vep_annotated.vcf \
    --format vcf --vcf --force_overwrite \
    --assembly GRCh38 --offline --cache \
    --everything \
    --af_gnomad --af_1kg --af_esp \
    --sift b --polyphen b --ccds --uniprot --hgvs --symbol \
    --numbers --domains --regulatory --canonical --protein \
    --biotype --tsl --appris \
    --check_existing --pubmed \
    --fork 8 \
    --stats_file ${OUTPUT_DIR}/vep_summary.html 2>&1 | tee -a $LOG_FILE

# === VALIDATION REPORT ===
echo "=== VALIDATION REPORT ===" | tee ${OUTPUT_DIR}/validation_report.txt
echo "Generated: $(date)" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "" | tee -a ${OUTPUT_DIR}/validation_report.txt

# Check output file exists
if [ -f "${OUTPUT_DIR}/vep_annotated.vcf" ]; then
    OUTPUT_COUNT=$(grep -v "^#" ${OUTPUT_DIR}/vep_annotated.vcf | wc -l)
    echo "✅ VEP annotation completed successfully" | tee -a ${OUTPUT_DIR}/validation_report.txt
    echo "Input variants: $INPUT_COUNT" | tee -a ${OUTPUT_DIR}/validation_report.txt
    echo "Output variants: $OUTPUT_COUNT" | tee -a ${OUTPUT_DIR}/validation_report.txt
    
    # Calculate success rate
    if [ $INPUT_COUNT -eq $OUTPUT_COUNT ]; then
        echo "✅ Variant count preserved: 100%" | tee -a ${OUTPUT_DIR}/validation_report.txt
    else
        RATE=$(echo "scale=2; $OUTPUT_COUNT * 100 / $INPUT_COUNT" | bc -l)
        echo "⚠️  Variant preservation rate: ${RATE}%" | tee -a ${OUTPUT_DIR}/validation_report.txt
    fi
    
    # Check file size
    OUTPUT_SIZE=$(du -h ${OUTPUT_DIR}/vep_annotated.vcf | cut -f1)
    echo "Output file size: $OUTPUT_SIZE" | tee -a ${OUTPUT_DIR}/validation_report.txt
    
    # Check for CSQ annotations
    CSQ_COUNT=$(grep "CSQ=" ${OUTPUT_DIR}/vep_annotated.vcf | wc -l)
    echo "Variants with VEP annotations: $CSQ_COUNT" | tee -a ${OUTPUT_DIR}/validation_report.txt
    
    if [ $CSQ_COUNT -gt 0 ]; then
        ANNOTATION_RATE=$(echo "scale=2; $CSQ_COUNT * 100 / $OUTPUT_COUNT" | bc -l)
        echo "✅ Annotation coverage: ${ANNOTATION_RATE}%" | tee -a ${OUTPUT_DIR}/validation_report.txt
    else
        echo "❌ No VEP annotations found" | tee -a ${OUTPUT_DIR}/validation_report.txt
    fi
    
else
    echo "❌ VEP annotation failed - output file not found" | tee -a ${OUTPUT_DIR}/validation_report.txt
    exit 1
fi

echo "" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "=== FILES GENERATED ===" | tee -a ${OUTPUT_DIR}/validation_report.txt
ls -lh ${OUTPUT_DIR}/ | tee -a ${OUTPUT_DIR}/validation_report.txt

echo "" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "=== NEXT STEPS ===" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "1. Review validation report: ${OUTPUT_DIR}/validation_report.txt" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "2. Check VEP summary: ${OUTPUT_DIR}/vep_summary.html" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "3. Proceed to TabNet training with: ${OUTPUT_DIR}/vep_annotated.vcf" | tee -a ${OUTPUT_DIR}/validation_report.txt

echo "=== VEP Pipeline Completed: $(date) ===" | tee -a $LOG_FILE
echo "Validation report: ${OUTPUT_DIR}/validation_report.txt"