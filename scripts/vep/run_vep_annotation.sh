#!/bin/bash
#SBATCH --job-name=vep_annotation
#SBATCH --partition=IllinoisComputes
#SBATCH --account=aa107-ic
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=32G
#SBATCH --output=/u/aa107/uiuc-cancer-research/scripts/vep/vep_%j.out
#SBATCH --error=/u/aa107/uiuc-cancer-research/scripts/vep/vep_%j.err

# Simple VEP Annotation Script using BioContainers
# Uses quay.io/biocontainers for better dependency resolution

set -e  # Exit on any error

# === CONFIGURATION ===
PROJECT_DIR="/u/aa107/uiuc-cancer-research"
INPUT_VCF="${PROJECT_DIR}/data/processed/merged_vcf/merged_prostate_variants.vcf"
OUTPUT_DIR="${PROJECT_DIR}/data/processed/vep"
LOG_FILE="${OUTPUT_DIR}/vep_annotation.log"
VEP_CONTAINER="${OUTPUT_DIR}/vep-biocontainer.sif"

echo "=== VEP BioContainer Pipeline Started: $(date) ===" | tee $LOG_FILE

# === CREATE OUTPUT DIRECTORY ===
mkdir -p "${OUTPUT_DIR}"

# === PULL BIOCONTAINER VEP IF NOT EXISTS ===
if [ ! -f "$VEP_CONTAINER" ]; then
    echo "Downloading BioContainers VEP..." | tee -a $LOG_FILE
    cd "${OUTPUT_DIR}"
    apptainer pull vep-biocontainer.sif docker://quay.io/biocontainers/ensembl-vep:112.0_cv1
    echo "✅ BioContainers VEP downloaded" | tee -a $LOG_FILE
else
    echo "✅ BioContainers VEP container already exists" | tee -a $LOG_FILE
fi

# === TEST DATABASE CONNECTIVITY ===
echo "Testing database connectivity..." | tee -a $LOG_FILE
DB_TEST=$(apptainer exec --bind ${PROJECT_DIR}:${PROJECT_DIR} $VEP_CONTAINER \
    perl -e "use DBI; my \$dbh = DBI->connect('DBI:mysql:host=ensembldb.ensembl.org;port=3306;database=homo_sapiens_core_112_38', 'anonymous', ''); print 'SUCCESS' if \$dbh; \$dbh->disconnect() if \$dbh;" 2>/dev/null || echo "FAILED")

if [ "$DB_TEST" = "SUCCESS" ]; then
    echo "✅ Database connectivity test passed" | tee -a $LOG_FILE
else
    echo "❌ Database connectivity test failed - trying alternative approach" | tee -a $LOG_FILE
fi

# === VALIDATE INPUT ===
echo "Validating input file..." | tee -a $LOG_FILE
if [ ! -f "$INPUT_VCF" ]; then
    echo "❌ ERROR: Input VCF not found: $INPUT_VCF" | tee -a $LOG_FILE
    exit 1
fi

INPUT_COUNT=$(grep -v "^#" $INPUT_VCF | wc -l)
echo "✅ Input variants: $INPUT_COUNT" | tee -a $LOG_FILE

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
    --database --host ensembldb.ensembl.org --port 3306 \
    --sift b --polyphen b --symbol --numbers \
    --canonical --protein --biotype \
    --fork 2 2>&1 | tee -a $LOG_FILE

# Check if test worked
if [ -f "${OUTPUT_DIR}/test_output.vcf" ]; then
    TEST_VARIANTS=$(grep -v "^#" ${OUTPUT_DIR}/test_output.vcf | wc -l)
    echo "✅ Test successful: $TEST_VARIANTS variants processed" | tee -a $LOG_FILE
    rm -f ${OUTPUT_DIR}/test_sample.vcf ${OUTPUT_DIR}/test_output.vcf
else
    echo "❌ Test failed - switching to offline mode would require cache download" | tee -a $LOG_FILE
    echo "For this test, we'll proceed with online mode and troubleshoot if needed" | tee -a $LOG_FILE
fi

# === RUN FULL VEP ANNOTATION ===
echo "Running full VEP annotation..." | tee -a $LOG_FILE

apptainer exec \
    --bind ${PROJECT_DIR}:${PROJECT_DIR} \
    $VEP_CONTAINER \
    vep \
    --input_file $INPUT_VCF \
    --output_file ${OUTPUT_DIR}/vep_annotated.vcf \
    --format vcf --vcf --force_overwrite \
    --species homo_sapiens --assembly GRCh38 \
    --database --host ensembldb.ensembl.org --port 3306 \
    --sift b --polyphen b --symbol --numbers --biotype \
    --canonical --protein --ccds --uniprot --domains \
    --regulatory --variant_class \
    --fork 8 \
    --stats_file ${OUTPUT_DIR}/vep_summary.html 2>&1 | tee -a $LOG_FILE

# === SIMPLE VALIDATION ===
echo "=== SIMPLE VALIDATION ===" | tee ${OUTPUT_DIR}/validation_report.txt
echo "Generated: $(date)" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "" | tee -a ${OUTPUT_DIR}/validation_report.txt

if [ -f "${OUTPUT_DIR}/vep_annotated.vcf" ]; then
    OUTPUT_COUNT=$(grep -v "^#" ${OUTPUT_DIR}/vep_annotated.vcf | wc -l)
    echo "✅ VEP annotation completed" | tee -a ${OUTPUT_DIR}/validation_report.txt
    echo "Input variants: $INPUT_COUNT" | tee -a ${OUTPUT_DIR}/validation_report.txt
    echo "Output variants: $OUTPUT_COUNT" | tee -a ${OUTPUT_DIR}/validation_report.txt
    
    # Simple annotation check
    if [ $OUTPUT_COUNT -gt 0 ]; then
        # Check for CSQ annotations (VEP consequence annotations)
        CSQ_COUNT=$(grep -c "CSQ=" ${OUTPUT_DIR}/vep_annotated.vcf || echo "0")
        if [ $CSQ_COUNT -gt 0 ]; then
            echo "✅ VEP annotations found: $CSQ_COUNT variants with CSQ data" | tee -a ${OUTPUT_DIR}/validation_report.txt
        else
            echo "⚠️ No CSQ annotations found - checking for other VEP annotations" | tee -a ${OUTPUT_DIR}/validation_report.txt
            # Check for any VEP-added INFO fields
            VEP_INFO=$(grep -v "^#" ${OUTPUT_DIR}/vep_annotated.vcf | head -1 | cut -f8)
            echo "Sample INFO field: $VEP_INFO" | tee -a ${OUTPUT_DIR}/validation_report.txt
        fi
    else
        echo "❌ No variants in output file" | tee -a ${OUTPUT_DIR}/validation_report.txt
    fi
    
    OUTPUT_SIZE=$(du -h ${OUTPUT_DIR}/vep_annotated.vcf | cut -f1)
    echo "Output file size: $OUTPUT_SIZE" | tee -a ${OUTPUT_DIR}/validation_report.txt
    
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
echo "3. Inspect annotated VCF: head -50 ${OUTPUT_DIR}/vep_annotated.vcf" | tee -a ${OUTPUT_DIR}/validation_report.txt

echo "=== VEP BioContainer Pipeline Completed: $(date) ===" | tee -a $LOG_FILE

# === FINAL STATUS CHECK ===
if [ -f "${OUTPUT_DIR}/vep_annotated.vcf" ] && [ $OUTPUT_COUNT -gt 0 ]; then
    echo "✅ SUCCESS! VEP annotation completed with $OUTPUT_COUNT variants"
else
    echo "❌ FAILED! Check logs and validation report"
    exit 1
fi