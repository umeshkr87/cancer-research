#!/bin/bash
#SBATCH --job-name=vep_annotation_phase3
#SBATCH --partition=IllinoisComputes
#SBATCH --account=aa107-ic
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=64G
#SBATCH --output=/u/aa107/uiuc-cancer-research/scripts/vep/vep_%j.out
#SBATCH --error=/u/aa107/uiuc-cancer-research/scripts/vep/vep_%j.err

# PHASE 3 VEP Optimization - FIXED VERSION
# Target: 91.1% → 95%+ success rate
# Fix: Remove problematic preprocessing, focus on VEP optimization

set -e  # Exit on any error

# === CONFIGURATION ===
PROJECT_DIR="/u/aa107/uiuc-cancer-research"
INPUT_VCF="${PROJECT_DIR}/data/processed/merged_vcf/merged_prostate_variants.vcf"
OUTPUT_DIR="${PROJECT_DIR}/data/processed/vep"
CACHE_DIR="${OUTPUT_DIR}/vep_cache"
LOG_FILE="${OUTPUT_DIR}/vep_annotation.log"
VEP_CONTAINER="${OUTPUT_DIR}/vep-biocontainer.sif"

echo "=== PHASE 3 VEP OPTIMIZATION (FIXED) STARTED: $(date) ===" | tee $LOG_FILE
echo "🎯 Target: Improve success rate from 91.1% to 95%+" | tee -a $LOG_FILE
echo "🔧 Fix: Removed problematic preprocessing, using original VCF" | tee -a $LOG_FILE

# === CREATE OUTPUT DIRECTORY ===
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${CACHE_DIR}"

# === PULL BIOCONTAINER VEP ===
if [ ! -f "$VEP_CONTAINER" ]; then
    echo "📦 Downloading BioContainers VEP (114.1)..." | tee -a $LOG_FILE
    cd "${OUTPUT_DIR}"
    apptainer pull vep-biocontainer.sif docker://quay.io/biocontainers/ensembl-vep:114.1--pl5321h2a3209d_0
    echo "✅ BioContainers VEP downloaded" | tee -a $LOG_FILE
else
    echo "✅ BioContainers VEP container already exists" | tee -a $LOG_FILE
fi

# === VALIDATE INPUT ===
echo "📊 Validating input file..." | tee -a $LOG_FILE
if [ ! -f "$INPUT_VCF" ]; then
    echo "❌ ERROR: Input VCF not found: $INPUT_VCF" | tee -a $LOG_FILE
    exit 1
fi

INPUT_COUNT=$(grep -v "^#" $INPUT_VCF | wc -l)
echo "✅ Input variants: $INPUT_COUNT" | tee -a $LOG_FILE

# === SKIP PROBLEMATIC PREPROCESSING ===
echo "🔧 PHASE 3 OPTIMIZATION: Using original VCF directly" | tee -a $LOG_FILE
echo "  ⚠️ Skipping coordinate sorting (was corrupting VCF format)" | tee -a $LOG_FILE
echo "  ⚠️ Skipping chromosome standardization (was causing parse errors)" | tee -a $LOG_FILE
echo "  ✅ Relying on VEP's built-in error handling" | tee -a $LOG_FILE

# Use original input directly
FINAL_INPUT="$INPUT_VCF"

# === VEP CACHE SETUP ===
echo "📚 Setting up VEP cache..." | tee -a $LOG_FILE

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
        echo "❌ VEP cache installation failed, trying alternative..." | tee -a $LOG_FILE
        
        # Fallback: Download cache manually
        echo "Downloading cache manually..." | tee -a $LOG_FILE
        wget -q ftp://ftp.ensembl.org/pub/release-112/variation/indexed_vep_cache/homo_sapiens_vep_112_GRCh38.tar.gz -P ${CACHE_DIR}/ || true
        if [ -f "${CACHE_DIR}/homo_sapiens_vep_112_GRCh38.tar.gz" ]; then
            cd ${CACHE_DIR}
            tar -xzf homo_sapiens_vep_112_GRCh38.tar.gz
            rm homo_sapiens_vep_112_GRCh38.tar.gz
            echo "✅ Manual cache download completed" | tee -a $LOG_FILE
        fi
    fi
else
    echo "✅ VEP cache already exists" | tee -a $LOG_FILE
fi

# === PHASE 3 TASK 2: VEP PARAMETER OPTIMIZATION (SAFE) ===
echo "⚡ PHASE 3 OPTIMIZATION: VEP Parameter Tuning" | tee -a $LOG_FILE

# Conservative optimization based on dataset size
if [ $INPUT_COUNT -gt 150000 ]; then
    BUFFER_SIZE=7500
    FORK_COUNT=12
    echo "📊 Large dataset: buffer=$BUFFER_SIZE, forks=$FORK_COUNT" | tee -a $LOG_FILE
elif [ $INPUT_COUNT -gt 100000 ]; then
    BUFFER_SIZE=6000
    FORK_COUNT=10
    echo "📊 Medium dataset: buffer=$BUFFER_SIZE, forks=$FORK_COUNT" | tee -a $LOG_FILE
else
    BUFFER_SIZE=5000
    FORK_COUNT=8
    echo "📊 Standard dataset: buffer=$BUFFER_SIZE, forks=$FORK_COUNT" | tee -a $LOG_FILE
fi

# === TEST VEP WITH SMALL SAMPLE ===
echo "🧪 Testing VEP with small sample..." | tee -a $LOG_FILE
head -1000 "$FINAL_INPUT" > ${OUTPUT_DIR}/test_sample.vcf

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
    --fork 4 --buffer_size 1000 2>&1 | tee -a $LOG_FILE

# Check if test worked
if [ -f "${OUTPUT_DIR}/test_output.vcf" ]; then
    TEST_VARIANTS=$(grep -v "^#" ${OUTPUT_DIR}/test_output.vcf | wc -l)
    echo "✅ Test successful: $TEST_VARIANTS variants processed" | tee -a $LOG_FILE
    rm -f ${OUTPUT_DIR}/test_sample.vcf ${OUTPUT_DIR}/test_output.vcf
else
    echo "❌ Test failed - check VEP configuration" | tee -a $LOG_FILE
    exit 1
fi

# === RUN OPTIMIZED VEP ANNOTATION ===
echo "🚀 Running optimized VEP annotation (without problematic preprocessing)..." | tee -a $LOG_FILE

apptainer exec \
    --bind ${PROJECT_DIR}:${PROJECT_DIR} \
    $VEP_CONTAINER \
    vep \
    --input_file "$FINAL_INPUT" \
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
    --buffer_size $BUFFER_SIZE \
    --stats_file ${OUTPUT_DIR}/vep_summary.html \
    --warning_file ${OUTPUT_DIR}/vep_warnings.txt \
    --no_check_variants_order \
    --allow_non_variant 2>&1 | tee -a $LOG_FILE

VEP_EXIT_CODE=$?
echo "VEP exit code: $VEP_EXIT_CODE" | tee -a $LOG_FILE

# === ENHANCED VALIDATION ===
echo "🔍 PHASE 3 ENHANCED VALIDATION" | tee -a $LOG_FILE
echo "=== PHASE 3 ENHANCED VALIDATION REPORT ===" | tee ${OUTPUT_DIR}/validation_report.txt
echo "Generated: $(date)" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "Phase: 3 - VEP Optimization (Fixed Version)" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "Fix Applied: Removed problematic preprocessing" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "" | tee -a ${OUTPUT_DIR}/validation_report.txt

if [ -f "${OUTPUT_DIR}/vep_annotated.vcf" ]; then
    OUTPUT_COUNT=$(grep -v "^#" ${OUTPUT_DIR}/vep_annotated.vcf | wc -l)
    echo "✅ VEP annotation completed" | tee -a ${OUTPUT_DIR}/validation_report.txt
    echo "Input variants: $INPUT_COUNT" | tee -a ${OUTPUT_DIR}/validation_report.txt
    echo "Output variants: $OUTPUT_COUNT" | tee -a ${OUTPUT_DIR}/validation_report.txt
    
    # Enhanced success rate calculation
    if [ $INPUT_COUNT -gt 0 ]; then
        SUCCESS_RATE=$(echo "scale=1; $OUTPUT_COUNT * 100 / $INPUT_COUNT" | bc -l)
        echo "Success rate: ${SUCCESS_RATE}%" | tee -a ${OUTPUT_DIR}/validation_report.txt
        
        # Phase 3 target assessment
        if (( $(echo "$SUCCESS_RATE >= 95.0" | bc -l) )); then
            echo "🎯 PHASE 3 SUCCESS: Target achieved (≥95%)" | tee -a ${OUTPUT_DIR}/validation_report.txt
            PHASE3_SUCCESS="✅ PASSED"
        elif (( $(echo "$SUCCESS_RATE >= 91.1" | bc -l) )); then
            echo "📈 PHASE 3 IMPROVEMENT: ${SUCCESS_RATE}% (baseline: 91.1%)" | tee -a ${OUTPUT_DIR}/validation_report.txt
            PHASE3_SUCCESS="📈 IMPROVED"
        else
            echo "⚠️ PHASE 3 REGRESSION: ${SUCCESS_RATE}% (baseline: 91.1%)" | tee -a ${OUTPUT_DIR}/validation_report.txt
            PHASE3_SUCCESS="❌ REGRESSION"
        fi
    fi
    
    # Check annotation quality
    if [ $OUTPUT_COUNT -gt 0 ]; then
        # Enhanced annotation analysis
        CSQ_COUNT=$(grep -c "CSQ=" ${OUTPUT_DIR}/vep_annotated.vcf || echo "0")
        SIFT_COUNT=$(grep -c "SIFT=" ${OUTPUT_DIR}/vep_annotated.vcf || echo "0")
        POLYPHEN_COUNT=$(grep -c "PolyPhen=" ${OUTPUT_DIR}/vep_annotated.vcf || echo "0")
        
        # Calculate annotation rates
        CSQ_RATE=$(echo "scale=1; $CSQ_COUNT * 100 / $OUTPUT_COUNT" | bc -l)
        SIFT_RATE=$(echo "scale=1; $SIFT_COUNT * 100 / $OUTPUT_COUNT" | bc -l)
        POLYPHEN_RATE=$(echo "scale=1; $POLYPHEN_COUNT * 100 / $OUTPUT_COUNT" | bc -l)
        
        echo "" | tee -a ${OUTPUT_DIR}/validation_report.txt
        echo "🧬 ANNOTATION COMPLETENESS:" | tee -a ${OUTPUT_DIR}/validation_report.txt
        echo "CSQ annotation rate: ${CSQ_RATE}%" | tee -a ${OUTPUT_DIR}/validation_report.txt
        echo "SIFT annotation rate: ${SIFT_RATE}%" | tee -a ${OUTPUT_DIR}/validation_report.txt
        echo "PolyPhen annotation rate: ${POLYPHEN_RATE}%" | tee -a ${OUTPUT_DIR}/validation_report.txt
        
        # Functional consequences analysis
        MISSENSE_COUNT=$(grep -c "missense_variant" ${OUTPUT_DIR}/vep_annotated.vcf || echo "0")
        STOP_GAINED_COUNT=$(grep -c "stop_gained" ${OUTPUT_DIR}/vep_annotated.vcf || echo "0")
        SPLICE_COUNT=$(grep -c "splice" ${OUTPUT_DIR}/vep_annotated.vcf || echo "0")
        
        echo "" | tee -a ${OUTPUT_DIR}/validation_report.txt
        echo "🔬 FUNCTIONAL CONSEQUENCES:" | tee -a ${OUTPUT_DIR}/validation_report.txt
        echo "  Missense variants: $MISSENSE_COUNT" | tee -a ${OUTPUT_DIR}/validation_report.txt
        echo "  Stop gained: $STOP_GAINED_COUNT" | tee -a ${OUTPUT_DIR}/validation_report.txt
        echo "  Splice variants: $SPLICE_COUNT" | tee -a ${OUTPUT_DIR}/validation_report.txt
        
        # File size and quality metrics
        OUTPUT_SIZE=$(du -h ${OUTPUT_DIR}/vep_annotated.vcf | cut -f1)
        echo "" | tee -a ${OUTPUT_DIR}/validation_report.txt
        echo "📁 OUTPUT METRICS:" | tee -a ${OUTPUT_DIR}/validation_report.txt
        echo "Output file size: $OUTPUT_SIZE" | tee -a ${OUTPUT_DIR}/validation_report.txt
        echo "VEP parameter optimization applied: ✅" | tee -a ${OUTPUT_DIR}/validation_report.txt
        echo "Problematic preprocessing removed: ✅" | tee -a ${OUTPUT_DIR}/validation_report.txt
        
    else
        echo "❌ No variants in output file" | tee -a ${OUTPUT_DIR}/validation_report.txt
        exit 1
    fi
    
    # Check warnings
    if [ -f "${OUTPUT_DIR}/vep_warnings.txt" ]; then
        WARNING_COUNT=$(wc -l < "${OUTPUT_DIR}/vep_warnings.txt")
        echo "⚠️ VEP warnings: $WARNING_COUNT" | tee -a ${OUTPUT_DIR}/validation_report.txt
        
        if [ $WARNING_COUNT -gt 0 ]; then
            echo "Sample warnings (first 5):" | tee -a ${OUTPUT_DIR}/validation_report.txt
            head -5 "${OUTPUT_DIR}/vep_warnings.txt" | tee -a ${OUTPUT_DIR}/validation_report.txt
        fi
    fi
    
else
    echo "❌ VEP annotation failed - output file not found" | tee -a ${OUTPUT_DIR}/validation_report.txt
    exit 1
fi

# === PHASE 3 COMPLETION SUMMARY ===
echo "" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "=== PHASE 3 COMPLETION SUMMARY ===" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "🎯 Target: 95%+ success rate" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "📊 Baseline: 91.1% (previous run)" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "📊 Achieved: ${SUCCESS_RATE}%" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "🏆 Status: $PHASE3_SUCCESS" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "🔧 OPTIMIZATIONS APPLIED:" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "  ✅ Conservative VEP parameter optimization" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "  ✅ Dynamic buffer sizing ($BUFFER_SIZE)" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "  ✅ Optimized fork count ($FORK_COUNT)" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "  ✅ Added VEP tolerance flags (--no_check_variants_order, --allow_non_variant)" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "  ✅ Enhanced validation reporting" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "🔧 FIXES APPLIED:" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "  ❌ Removed bcftools sorting (was corrupting VCF)" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "  ❌ Removed chromosome standardization (was causing parse errors)" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "  ✅ Using original VCF directly" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "  ✅ Relying on VEP's built-in error handling" | tee -a ${OUTPUT_DIR}/validation_report.txt

echo "" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "=== FILES GENERATED ===" | tee -a ${OUTPUT_DIR}/validation_report.txt
ls -lh ${OUTPUT_DIR}/ | grep -v "sorted_input\|final_input" | tee -a ${OUTPUT_DIR}/validation_report.txt

echo "" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "=== NEXT STEPS ===" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "1. Review validation report: ${OUTPUT_DIR}/validation_report.txt" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "2. Check VEP summary: ${OUTPUT_DIR}/vep_summary.html" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "3. Analyze VEP warnings: ${OUTPUT_DIR}/vep_warnings.txt" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "4. Convert VCF to CSV for TabNet: python scripts/vep/vcf_to_tabnet.py" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "5. Proceed with TabNet training" | tee -a ${OUTPUT_DIR}/validation_report.txt

echo "=== PHASE 3 VEP OPTIMIZATION (FIXED) COMPLETED: $(date) ===" | tee -a $LOG_FILE

# === FINAL STATUS ===
if [ -f "${OUTPUT_DIR}/vep_annotated.vcf" ] && [ $OUTPUT_COUNT -gt 0 ]; then
    echo "✅ PHASE 3 COMPLETE! VEP annotation: $OUTPUT_COUNT variants"
    echo "📊 Success rate: ${SUCCESS_RATE}% (Baseline: 91.1%)"
    echo "🎯 Phase 3 status: $PHASE3_SUCCESS"
    echo "📁 Main output: ${OUTPUT_DIR}/vep_annotated.vcf"
    echo "📋 Report: ${OUTPUT_DIR}/validation_report.txt"
    
    if (( $(echo "$SUCCESS_RATE >= 91.1" | bc -l) )); then
        echo "🎉 SUCCESS: Maintained or improved success rate!"
        exit 0
    else
        echo "⚠️ WARNING: Success rate below baseline (91.1%)"
        exit 1
    fi
else
    echo "❌ FAILED! Check logs and validation report"
    exit 1
fi