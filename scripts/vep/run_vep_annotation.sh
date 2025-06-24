#!/bin/bash
#SBATCH --job-name=vep_annotation
#SBATCH --partition=IllinoisComputes
#SBATCH --account=aa107-ic
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=64G
#SBATCH --output=/u/aa107/uiuc-cancer-research/scripts/vep/vep_%j.out
#SBATCH --error=/u/aa107/uiuc-cancer-research/scripts/vep/vep_%j.err

# Robust VEP Annotation Script with Cache Mode
# Updated with correct container version and improved reliability

set -e  # Exit on any error

# === CONFIGURATION ===
PROJECT_DIR="/u/aa107/uiuc-cancer-research"
INPUT_VCF="${PROJECT_DIR}/data/processed/merged_vcf/merged_prostate_variants.vcf"
OUTPUT_DIR="${PROJECT_DIR}/data/processed/vep"
CACHE_DIR="${OUTPUT_DIR}/vep_cache"
LOG_FILE="${OUTPUT_DIR}/vep_annotation.log"
VEP_CONTAINER="${OUTPUT_DIR}/vep-biocontainer.sif"

echo "=== VEP Cache-Based Pipeline Started: $(date) ===" | tee $LOG_FILE

# === CREATE OUTPUT DIRECTORY ===
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${CACHE_DIR}"

# === PULL CORRECT BIOCONTAINER VEP ===
if [ ! -f "$VEP_CONTAINER" ]; then
    echo "Downloading correct BioContainers VEP (114.1)..." | tee -a $LOG_FILE
    cd "${OUTPUT_DIR}"
    apptainer pull vep-biocontainer.sif docker://quay.io/biocontainers/ensembl-vep:114.1--pl5321h2a3209d_0
    echo "✅ BioContainers VEP downloaded" | tee -a $LOG_FILE
else
    echo "✅ BioContainers VEP container already exists" | tee -a $LOG_FILE
fi

# === VALIDATE INPUT ===
echo "Validating input file..." | tee -a $LOG_FILE
if [ ! -f "$INPUT_VCF" ]; then
    echo "❌ ERROR: Input VCF not found: $INPUT_VCF" | tee -a $LOG_FILE
    exit 1
fi

INPUT_COUNT=$(grep -v "^#" $INPUT_VCF | wc -l)
echo "✅ Input variants: $INPUT_COUNT" | tee -a $LOG_FILE

if [ $INPUT_COUNT -gt 300000 ]; then
    echo "⚠️  Large dataset detected ($INPUT_COUNT variants) - using optimized settings" | tee -a $LOG_FILE
fi

# === INSTALL VEP CACHE (OFFLINE MODE) ===
echo "Setting up VEP cache for offline operation..." | tee -a $LOG_FILE

if [ ! -d "${CACHE_DIR}/homo_sapiens" ]; then
    echo "Installing VEP cache for human GRCh38..." | tee -a $LOG_FILE
    
    apptainer exec \
        --bind ${PROJECT_DIR}:${PROJECT_DIR} \
        $VEP_CONTAINER \
        vep_install \
        -a cf \
        -s homo_sapiens \
        -y GRCh38 \
        -c ${CACHE_DIR} \
        --CONVERT 2>&1 | tee -a $LOG_FILE
    
    if [ $? -eq 0 ]; then
        echo "✅ VEP cache installation completed" | tee -a $LOG_FILE
    else
        echo "❌ VEP cache installation failed, trying alternative method..." | tee -a $LOG_FILE
        
        # Fallback: Download cache manually
        echo "Downloading cache manually..." | tee -a $LOG_FILE
        wget -q ftp://ftp.ensembl.org/pub/release-112/variation/indexed_vep_cache/homo_sapiens_vep_112_GRCh38.tar.gz -P ${CACHE_DIR}/
        cd ${CACHE_DIR}
        tar -xzf homo_sapiens_vep_112_GRCh38.tar.gz
        rm homo_sapiens_vep_112_GRCh38.tar.gz
        echo "✅ Manual cache download completed" | tee -a $LOG_FILE
    fi
else
    echo "✅ VEP cache already exists" | tee -a $LOG_FILE
fi

# === TEST VEP WITH SMALL SAMPLE ===
echo "Testing VEP with small sample..." | tee -a $LOG_FILE
head -1000 $INPUT_VCF > ${OUTPUT_DIR}/test_sample.vcf

apptainer exec \
    --bind ${PROJECT_DIR}:${PROJECT_DIR} \
    $VEP_CONTAINER \
    vep \
    --input_file ${OUTPUT_DIR}/test_sample.vcf \
    --output_file ${OUTPUT_DIR}/test_output.vcf \
    --format vcf --vcf --force_overwrite \
    --species homo_sapiens --assembly GRCh38 \
    --cache --dir_cache ${CACHE_DIR} \
    --offline \
    --sift b --polyphen b --symbol --numbers \
    --canonical --protein --biotype \
    --fork 4 2>&1 | tee -a $LOG_FILE

# Check if test worked
if [ -f "${OUTPUT_DIR}/test_output.vcf" ]; then
    TEST_VARIANTS=$(grep -v "^#" ${OUTPUT_DIR}/test_output.vcf | wc -l)
    echo "✅ Test successful: $TEST_VARIANTS variants processed" | tee -a $LOG_FILE
    rm -f ${OUTPUT_DIR}/test_sample.vcf ${OUTPUT_DIR}/test_output.vcf
else
    echo "❌ Test failed - check VEP configuration" | tee -a $LOG_FILE
    exit 1
fi

# === RUN FULL VEP ANNOTATION (CACHE MODE) ===
echo "Running full VEP annotation with cache..." | tee -a $LOG_FILE

# Determine optimal fork count based on dataset size
if [ $INPUT_COUNT -gt 100000 ]; then
    FORK_COUNT=16
    echo "Using 16 forks for large dataset" | tee -a $LOG_FILE
else
    FORK_COUNT=8
    echo "Using 8 forks for standard dataset" | tee -a $LOG_FILE
fi

apptainer exec \
    --bind ${PROJECT_DIR}:${PROJECT_DIR} \
    $VEP_CONTAINER \
    vep \
    --input_file $INPUT_VCF \
    --output_file ${OUTPUT_DIR}/vep_annotated.vcf \
    --format vcf --vcf --force_overwrite \
    --species homo_sapiens --assembly GRCh38 \
    --cache --dir_cache ${CACHE_DIR} \
    --offline \
    --sift b --polyphen b --symbol --numbers --biotype \
    --canonical --protein --ccds --uniprot --domains \
    --regulatory --variant_class \
    --af --af_1kg --af_gnomad \
    --pubmed --var_synonyms \
    --fork $FORK_COUNT \
    --buffer_size 5000 \
    --stats_file ${OUTPUT_DIR}/vep_summary.html 2>&1 | tee -a $LOG_FILE

# === VALIDATION ===
echo "=== VALIDATION ===" | tee ${OUTPUT_DIR}/validation_report.txt
echo "Generated: $(date)" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "" | tee -a ${OUTPUT_DIR}/validation_report.txt

if [ -f "${OUTPUT_DIR}/vep_annotated.vcf" ]; then
    OUTPUT_COUNT=$(grep -v "^#" ${OUTPUT_DIR}/vep_annotated.vcf | wc -l)
    echo "✅ VEP annotation completed" | tee -a ${OUTPUT_DIR}/validation_report.txt
    echo "Input variants: $INPUT_COUNT" | tee -a ${OUTPUT_DIR}/validation_report.txt
    echo "Output variants: $OUTPUT_COUNT" | tee -a ${OUTPUT_DIR}/validation_report.txt
    
    # Check annotation quality
    if [ $OUTPUT_COUNT -gt 0 ]; then
        # Check for CSQ annotations
        CSQ_COUNT=$(grep -c "CSQ=" ${OUTPUT_DIR}/vep_annotated.vcf || echo "0")
        if [ $CSQ_COUNT -gt 0 ]; then
            echo "✅ VEP annotations found: $CSQ_COUNT variants with CSQ data" | tee -a ${OUTPUT_DIR}/validation_report.txt
            
            # Sample annotation check
            SAMPLE_CSQ=$(grep -v "^#" ${OUTPUT_DIR}/vep_annotated.vcf | head -1 | grep -o "CSQ=[^;]*" | head -1)
            echo "Sample annotation: $SAMPLE_CSQ" | tee -a ${OUTPUT_DIR}/validation_report.txt
            
            # Count functional consequences
            MISSENSE_COUNT=$(grep -c "missense_variant" ${OUTPUT_DIR}/vep_annotated.vcf || echo "0")
            STOP_GAINED_COUNT=$(grep -c "stop_gained" ${OUTPUT_DIR}/vep_annotated.vcf || echo "0")
            SPLICE_COUNT=$(grep -c "splice" ${OUTPUT_DIR}/vep_annotated.vcf || echo "0")
            
            echo "Functional consequences found:" | tee -a ${OUTPUT_DIR}/validation_report.txt
            echo "  Missense variants: $MISSENSE_COUNT" | tee -a ${OUTPUT_DIR}/validation_report.txt
            echo "  Stop gained: $STOP_GAINED_COUNT" | tee -a ${OUTPUT_DIR}/validation_report.txt
            echo "  Splice variants: $SPLICE_COUNT" | tee -a ${OUTPUT_DIR}/validation_report.txt
            
        else
            echo "⚠️ No CSQ annotations found - checking for other VEP fields" | tee -a ${OUTPUT_DIR}/validation_report.txt
        fi
    else
        echo "❌ No variants in output file" | tee -a ${OUTPUT_DIR}/validation_report.txt
        exit 1
    fi
    
    OUTPUT_SIZE=$(du -h ${OUTPUT_DIR}/vep_annotated.vcf | cut -f1)
    echo "Output file size: $OUTPUT_SIZE" | tee -a ${OUTPUT_DIR}/validation_report.txt
    
    # Success rate calculation
    if [ $INPUT_COUNT -gt 0 ]; then
        SUCCESS_RATE=$(echo "scale=1; $OUTPUT_COUNT * 100 / $INPUT_COUNT" | bc -l)
        echo "Success rate: ${SUCCESS_RATE}%" | tee -a ${OUTPUT_DIR}/validation_report.txt
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
echo "3. Convert VCF to CSV for TabNet: python scripts/vep/vcf_to_tabnet.py" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "4. Proceed with TabNet training" | tee -a ${OUTPUT_DIR}/validation_report.txt

echo "=== VEP Cache-Based Pipeline Completed: $(date) ===" | tee -a $LOG_FILE

# === FINAL STATUS CHECK ===
if [ -f "${OUTPUT_DIR}/vep_annotated.vcf" ] && [ $OUTPUT_COUNT -gt 0 ]; then
    echo "✅ SUCCESS! VEP annotation completed with $OUTPUT_COUNT variants"
    echo "📁 Main output: ${OUTPUT_DIR}/vep_annotated.vcf"
    echo "📊 Summary: ${OUTPUT_DIR}/vep_summary.html"
    echo "📋 Report: ${OUTPUT_DIR}/validation_report.txt"
else
    echo "❌ FAILED! Check logs and validation report"
    exit 1
fi