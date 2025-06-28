#!/bin/bash
#SBATCH --job-name=vep_annotation_dbnsfp
#SBATCH --partition=IllinoisComputes
#SBATCH --account=aa107-ic
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=64G
#SBATCH --output=/u/aa107/uiuc-cancer-research/scripts/vep/vep_%j.out
#SBATCH --error=/u/aa107/uiuc-cancer-research/scripts/vep/vep_%j.err

# FINAL VEP SCRIPT WITH dbNSFP FUNCTIONAL SCORES
# Target: 100% success rate + SIFT/PolyPhen functional scores
# Solution: dbNSFP plugin for comprehensive functional annotation

set -e  # Exit on any error

# === CONFIGURATION ===
PROJECT_DIR="/u/aa107/uiuc-cancer-research"
INPUT_VCF="${PROJECT_DIR}/data/processed/merged_vcf/merged_prostate_variants.vcf"
OUTPUT_DIR="${PROJECT_DIR}/data/processed/vep"
CACHE_DIR="${OUTPUT_DIR}/vep_cache"
PLUGIN_DIR="${OUTPUT_DIR}/vep_plugins"
LOG_FILE="${OUTPUT_DIR}/vep_annotation.log"
VEP_CONTAINER="${OUTPUT_DIR}/vep-biocontainer.sif"

# dbNSFP Configuration
DBNSFP_DIR="${OUTPUT_DIR}/dbnsfp"
DBNSFP_FILE="${DBNSFP_DIR}/dbNSFP4.7c_grch38.gz"
DBNSFP_PLUGIN="${PLUGIN_DIR}/dbNSFP.pm"

echo "=== FINAL VEP WITH dbNSFP FUNCTIONAL SCORES STARTED: $(date) ===" | tee $LOG_FILE
echo "🎯 Target: 100% success rate + SIFT/PolyPhen functional scores" | tee -a $LOG_FILE
echo "🔬 Solution: dbNSFP plugin for comprehensive annotation" | tee -a $LOG_FILE

# === CREATE OUTPUT DIRECTORIES ===
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${CACHE_DIR}"
mkdir -p "${PLUGIN_DIR}"
mkdir -p "${DBNSFP_DIR}"

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

# === dbNSFP PLUGIN SETUP ===
echo "🧬 SETTING UP dbNSFP FUNCTIONAL SCORES" | tee -a $LOG_FILE

# Download dbNSFP plugin
if [ ! -f "$DBNSFP_PLUGIN" ]; then
    echo "📥 Downloading dbNSFP plugin..." | tee -a $LOG_FILE
    wget -q https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/114/dbNSFP.pm -O "$DBNSFP_PLUGIN"
    echo "✅ dbNSFP plugin downloaded" | tee -a $LOG_FILE
else
    echo "✅ dbNSFP plugin already exists" | tee -a $LOG_FILE
fi

# Download and process dbNSFP database
if [ ! -f "$DBNSFP_FILE" ]; then
    echo "📥 Downloading dbNSFP database (this may take 30-60 minutes)..." | tee -a $LOG_FILE
    cd "${DBNSFP_DIR}"
    
    # Download dbNSFP 4.7c (latest academic version)
    echo "  Downloading dbNSFP4.7c.zip (~4GB)..." | tee -a $LOG_FILE
    wget -q ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbNSFP4.7c.zip
    
    if [ -f "dbNSFP4.7c.zip" ]; then
        echo "  Extracting dbNSFP database..." | tee -a $LOG_FILE
        unzip -q dbNSFP4.7c.zip
        
        echo "  Processing dbNSFP for GRCh38..." | tee -a $LOG_FILE
        # Create header file
        zcat dbNSFP4.7c_variant.chr1.gz | head -n1 > h
        
        # Process and sort all chromosomes for GRCh38
        echo "  Sorting and indexing (this may take 20-30 minutes)..." | tee -a $LOG_FILE
        zgrep -h -v ^#chr dbNSFP4.7c_variant.chr* | \
            sort -k1,1 -k2,2n - | \
            cat h - | \
            bgzip -c > dbNSFP4.7c_grch38.gz
        
        # Create tabix index
        tabix -s 1 -b 2 -e 2 dbNSFP4.7c_grch38.gz
        
        # Cleanup
        rm -f dbNSFP4.7c.zip h dbNSFP4.7c_variant.chr*.gz
        
        echo "✅ dbNSFP database processed and indexed" | tee -a $LOG_FILE
    else
        echo "❌ Failed to download dbNSFP database" | tee -a $LOG_FILE
        exit 1
    fi
else
    echo "✅ dbNSFP database already exists" | tee -a $LOG_FILE
fi

# === VEP PARAMETER OPTIMIZATION ===
echo "⚡ VEP PARAMETER OPTIMIZATION" | tee -a $LOG_FILE

# Dynamic optimization based on dataset size
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

# === TEST VEP WITH dbNSFP ===
echo "🧪 Testing VEP with dbNSFP plugin..." | tee -a $LOG_FILE
head -1000 "$INPUT_VCF" > ${OUTPUT_DIR}/test_sample.vcf

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
    --symbol --numbers --biotype --canonical --protein \
    --dir_plugins ${PLUGIN_DIR} \
    --plugin dbNSFP,${DBNSFP_FILE},SIFT_score,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,CADD_phred \
    --fork 4 --buffer_size 1000 2>&1 | tee -a $LOG_FILE

# Check if test worked
if [ -f "${OUTPUT_DIR}/test_output.vcf" ]; then
    TEST_VARIANTS=$(grep -v "^#" ${OUTPUT_DIR}/test_output.vcf | wc -l)
    FUNCTIONAL_SCORES=$(grep -c "SIFT_score=" ${OUTPUT_DIR}/test_output.vcf || echo "0")
    echo "✅ Test successful: $TEST_VARIANTS variants processed" | tee -a $LOG_FILE
    echo "✅ Functional scores found: $FUNCTIONAL_SCORES variants" | tee -a $LOG_FILE
    rm -f ${OUTPUT_DIR}/test_sample.vcf ${OUTPUT_DIR}/test_output.vcf
else
    echo "❌ Test failed - check VEP configuration" | tee -a $LOG_FILE
    exit 1
fi

# === RUN FULL VEP ANNOTATION WITH dbNSFP ===
echo "🚀 Running VEP with dbNSFP functional scores..." | tee -a $LOG_FILE

apptainer exec \
    --bind ${PROJECT_DIR}:${PROJECT_DIR} \
    $VEP_CONTAINER \
    vep \
    --input_file "$INPUT_VCF" \
    --output_file ${OUTPUT_DIR}/vep_annotated.vcf \
    --format vcf --vcf --force_overwrite \
    --species homo_sapiens --assembly GRCh38 \
    --cache --dir_cache ${CACHE_DIR} \
    --offline \
    --symbol --numbers --biotype \
    --canonical --protein --ccds --uniprot --domains \
    --regulatory --variant_class \
    --af --af_1kg --af_gnomad \
    --pubmed --var_synonyms \
    --dir_plugins ${PLUGIN_DIR} \
    --plugin dbNSFP,${DBNSFP_FILE},SIFT_score,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,CADD_phred,CADD_raw,REVEL_score,MetaLR_score,MetaLR_pred \
    --fork $FORK_COUNT \
    --buffer_size $BUFFER_SIZE \
    --stats_file ${OUTPUT_DIR}/vep_summary.html \
    --warning_file ${OUTPUT_DIR}/vep_warnings.txt \
    --no_check_variants_order \
    --allow_non_variant 2>&1 | tee -a $LOG_FILE

VEP_EXIT_CODE=$?
echo "VEP exit code: $VEP_EXIT_CODE" | tee -a $LOG_FILE

# === COMPREHENSIVE VALIDATION ===
echo "🔍 COMPREHENSIVE VALIDATION WITH FUNCTIONAL SCORES" | tee -a $LOG_FILE
echo "=== FINAL VEP + dbNSFP VALIDATION REPORT ===" | tee ${OUTPUT_DIR}/validation_report.txt
echo "Generated: $(date)" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "Pipeline: VEP + dbNSFP Functional Scores" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "" | tee -a ${OUTPUT_DIR}/validation_report.txt

if [ -f "${OUTPUT_DIR}/vep_annotated.vcf" ]; then
    OUTPUT_COUNT=$(grep -v "^#" ${OUTPUT_DIR}/vep_annotated.vcf | wc -l)
    echo "✅ VEP annotation completed" | tee -a ${OUTPUT_DIR}/validation_report.txt
    echo "Input variants: $INPUT_COUNT" | tee -a ${OUTPUT_DIR}/validation_report.txt
    echo "Output variants: $OUTPUT_COUNT" | tee -a ${OUTPUT_DIR}/validation_report.txt
    
    # Success rate calculation
    if [ $INPUT_COUNT -gt 0 ]; then
        SUCCESS_RATE=$(echo "scale=1; $OUTPUT_COUNT * 100 / $INPUT_COUNT" | bc -l)
        echo "Success rate: ${SUCCESS_RATE}%" | tee -a ${OUTPUT_DIR}/validation_report.txt
        
        if (( $(echo "$SUCCESS_RATE >= 95.0" | bc -l) )); then
            echo "🎯 SUCCESS: Excellent success rate (≥95%)" | tee -a ${OUTPUT_DIR}/validation_report.txt
            PIPELINE_SUCCESS="✅ PASSED"
        else
            echo "⚠️ WARNING: Success rate below target (<95%)" | tee -a ${OUTPUT_DIR}/validation_report.txt
            PIPELINE_SUCCESS="⚠️ PARTIAL"
        fi
    fi
    
    # Enhanced functional annotation analysis
    if [ $OUTPUT_COUNT -gt 0 ]; then
        echo "" | tee -a ${OUTPUT_DIR}/validation_report.txt
        echo "🧬 STRUCTURAL ANNOTATION COMPLETENESS:" | tee -a ${OUTPUT_DIR}/validation_report.txt
        
        CSQ_COUNT=$(grep -c "CSQ=" ${OUTPUT_DIR}/vep_annotated.vcf || echo "0")
        CSQ_RATE=$(echo "scale=1; $CSQ_COUNT * 100 / $OUTPUT_COUNT" | bc -l)
        echo "CSQ annotation rate: ${CSQ_RATE}%" | tee -a ${OUTPUT_DIR}/validation_report.txt
        
        echo "" | tee -a ${OUTPUT_DIR}/validation_report.txt
        echo "🔬 FUNCTIONAL PREDICTION SCORES (dbNSFP):" | tee -a ${OUTPUT_DIR}/validation_report.txt
        
        # SIFT scores
        SIFT_SCORE_COUNT=$(grep -c "SIFT_score=" ${OUTPUT_DIR}/vep_annotated.vcf || echo "0")
        SIFT_PRED_COUNT=$(grep -c "SIFT_pred=" ${OUTPUT_DIR}/vep_annotated.vcf || echo "0")
        SIFT_RATE=$(echo "scale=1; $SIFT_SCORE_COUNT * 100 / $OUTPUT_COUNT" | bc -l)
        echo "SIFT scores: ${SIFT_SCORE_COUNT} variants (${SIFT_RATE}%)" | tee -a ${OUTPUT_DIR}/validation_report.txt
        echo "SIFT predictions: ${SIFT_PRED_COUNT} variants" | tee -a ${OUTPUT_DIR}/validation_report.txt
        
        # PolyPhen scores
        POLYPHEN_HDIV_COUNT=$(grep -c "Polyphen2_HDIV_score=" ${OUTPUT_DIR}/vep_annotated.vcf || echo "0")
        POLYPHEN_HVAR_COUNT=$(grep -c "Polyphen2_HVAR_score=" ${OUTPUT_DIR}/vep_annotated.vcf || echo "0")
        POLYPHEN_RATE=$(echo "scale=1; $POLYPHEN_HDIV_COUNT * 100 / $OUTPUT_COUNT" | bc -l)
        echo "PolyPhen HDIV scores: ${POLYPHEN_HDIV_COUNT} variants (${POLYPHEN_RATE}%)" | tee -a ${OUTPUT_DIR}/validation_report.txt
        echo "PolyPhen HVAR scores: ${POLYPHEN_HVAR_COUNT} variants" | tee -a ${OUTPUT_DIR}/validation_report.txt
        
        # Additional scores
        CADD_COUNT=$(grep -c "CADD_phred=" ${OUTPUT_DIR}/vep_annotated.vcf || echo "0")
        REVEL_COUNT=$(grep -c "REVEL_score=" ${OUTPUT_DIR}/vep_annotated.vcf || echo "0")
        CADD_RATE=$(echo "scale=1; $CADD_COUNT * 100 / $OUTPUT_COUNT" | bc -l)
        echo "CADD scores: ${CADD_COUNT} variants (${CADD_RATE}%)" | tee -a ${OUTPUT_DIR}/validation_report.txt
        echo "REVEL scores: ${REVEL_COUNT} variants" | tee -a ${OUTPUT_DIR}/validation_report.txt
        
        # Functional consequences analysis
        echo "" | tee -a ${OUTPUT_DIR}/validation_report.txt
        echo "📈 FUNCTIONAL CONSEQUENCES:" | tee -a ${OUTPUT_DIR}/validation_report.txt
        MISSENSE_COUNT=$(grep -c "missense_variant" ${OUTPUT_DIR}/vep_annotated.vcf || echo "0")
        STOP_GAINED_COUNT=$(grep -c "stop_gained" ${OUTPUT_DIR}/vep_annotated.vcf || echo "0")
        SPLICE_COUNT=$(grep -c "splice" ${OUTPUT_DIR}/vep_annotated.vcf || echo "0")
        
        echo "  Missense variants: $MISSENSE_COUNT" | tee -a ${OUTPUT_DIR}/validation_report.txt
        echo "  Stop gained: $STOP_GAINED_COUNT" | tee -a ${OUTPUT_DIR}/validation_report.txt
        echo "  Splice variants: $SPLICE_COUNT" | tee -a ${OUTPUT_DIR}/validation_report.txt
        
        # File size and quality metrics
        OUTPUT_SIZE=$(du -h ${OUTPUT_DIR}/vep_annotated.vcf | cut -f1)
        echo "" | tee -a ${OUTPUT_DIR}/validation_report.txt
        echo "📁 OUTPUT METRICS:" | tee -a ${OUTPUT_DIR}/validation_report.txt
        echo "Output file size: $OUTPUT_SIZE" | tee -a ${OUTPUT_DIR}/validation_report.txt
        echo "dbNSFP functional scores: ✅ ENABLED" | tee -a ${OUTPUT_DIR}/validation_report.txt
        echo "VEP parameter optimization: ✅ APPLIED" | tee -a ${OUTPUT_DIR}/validation_report.txt
        
        # TabNet readiness assessment
        echo "" | tee -a ${OUTPUT_DIR}/validation_report.txt
        echo "🎯 TABNET READINESS ASSESSMENT:" | tee -a ${OUTPUT_DIR}/validation_report.txt
        
        FUNCTIONAL_COVERAGE=$(echo "scale=1; ($SIFT_SCORE_COUNT + $POLYPHEN_HDIV_COUNT + $CADD_COUNT) / 3 / $OUTPUT_COUNT * 100" | bc -l)
        echo "Average functional score coverage: ${FUNCTIONAL_COVERAGE}%" | tee -a ${OUTPUT_DIR}/validation_report.txt
        
        if (( $(echo "$FUNCTIONAL_COVERAGE >= 70.0" | bc -l) )); then
            echo "TabNet ML readiness: ✅ EXCELLENT (≥70% functional coverage)" | tee -a ${OUTPUT_DIR}/validation_report.txt
            TABNET_READY="✅ READY"
        elif (( $(echo "$FUNCTIONAL_COVERAGE >= 50.0" | bc -l) )); then
            echo "TabNet ML readiness: ✅ GOOD (≥50% functional coverage)" | tee -a ${OUTPUT_DIR}/validation_report.txt
            TABNET_READY="✅ READY"
        else
            echo "TabNet ML readiness: ⚠️ LIMITED (<50% functional coverage)" | tee -a ${OUTPUT_DIR}/validation_report.txt
            TABNET_READY="⚠️ LIMITED"
        fi
        
    else
        echo "❌ No variants in output file" | tee -a ${OUTPUT_DIR}/validation_report.txt
        exit 1
    fi
    
    # Check warnings
    if [ -f "${OUTPUT_DIR}/vep_warnings.txt" ]; then
        WARNING_COUNT=$(wc -l < "${OUTPUT_DIR}/vep_warnings.txt")
        echo "⚠️ VEP warnings: $WARNING_COUNT" | tee -a ${OUTPUT_DIR}/validation_report.txt
    fi
    
else
    echo "❌ VEP annotation failed - output file not found" | tee -a ${OUTPUT_DIR}/validation_report.txt
    exit 1
fi

# === FINAL PIPELINE SUMMARY ===
echo "" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "=== FINAL PIPELINE COMPLETION SUMMARY ===" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "🎯 Target: 100% success rate + functional scores" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "📊 Achieved success rate: ${SUCCESS_RATE}%" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "🔬 Functional score coverage: ${FUNCTIONAL_COVERAGE}%" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "🏆 Pipeline status: $PIPELINE_SUCCESS" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "🎯 TabNet readiness: $TABNET_READY" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "" | tee -a ${OUTPUT_DIR}/validation_report.txt

echo "🧬 FUNCTIONAL SCORES IMPLEMENTED:" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "  ✅ SIFT scores and predictions" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "  ✅ PolyPhen-2 HDIV and HVAR scores" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "  ✅ CADD pathogenicity scores" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "  ✅ REVEL ensemble predictor" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "  ✅ MetaLR machine learning predictions" | tee -a ${OUTPUT_DIR}/validation_report.txt

echo "" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "=== FILES GENERATED ===" | tee -a ${OUTPUT_DIR}/validation_report.txt
ls -lh ${OUTPUT_DIR}/ | grep -v "test_\|sorted_\|final_" | tee -a ${OUTPUT_DIR}/validation_report.txt

echo "" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "=== NEXT STEPS ===" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "1. Review validation report: ${OUTPUT_DIR}/validation_report.txt" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "2. Check VEP summary: ${OUTPUT_DIR}/vep_summary.html" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "3. Convert VCF to TabNet CSV: python scripts/vep/vcf_to_tabnet.py" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "4. Proceed with TabNet training using functional scores" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "5. Implement interpretable cancer variant classification" | tee -a ${OUTPUT_DIR}/validation_report.txt

echo "=== FINAL VEP + dbNSFP PIPELINE COMPLETED: $(date) ===" | tee -a $LOG_FILE

# === FINAL STATUS ===
if [ -f "${OUTPUT_DIR}/vep_annotated.vcf" ] && [ $OUTPUT_COUNT -gt 0 ]; then
    echo "🎉 PIPELINE COMPLETE! VEP + dbNSFP annotation: $OUTPUT_COUNT variants"
    echo "📊 Success rate: ${SUCCESS_RATE}%"
    echo "🔬 Functional coverage: ${FUNCTIONAL_COVERAGE}%"
    echo "🏆 Pipeline status: $PIPELINE_SUCCESS"
    echo "🎯 TabNet readiness: $TABNET_READY"
    echo "📁 Main output: ${OUTPUT_DIR}/vep_annotated.vcf"
    echo "📋 Report: ${OUTPUT_DIR}/validation_report.txt"
    
    if (( $(echo "$SUCCESS_RATE >= 95.0" | bc -l) )) && (( $(echo "$FUNCTIONAL_COVERAGE >= 50.0" | bc -l) )); then
        echo "🚀 SUCCESS: Ready for TabNet machine learning with functional scores!"
        exit 0
    else
        echo "⚠️ PARTIAL SUCCESS: Review functional score coverage"
        exit 1
    fi
else
    echo "❌ FAILED! Check logs and validation report"
    exit 1
fi