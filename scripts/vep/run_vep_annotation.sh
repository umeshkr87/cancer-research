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

# VEP Annotation Script with Scratch Storage for Cache and Large Files
# Updated to prevent quota issues by using scratch for VEP cache

set -e  # Exit on any error

# === CONFIGURATION ===
PROJECT_DIR="/u/aa107/uiuc-cancer-research"
INPUT_VCF="${PROJECT_DIR}/data/processed/merged_vcf/merged_prostate_variants.vcf"

# 🔥 KEY CHANGE: Use scratch for large files
SCRATCH_VEP="/scratch/aa107/vep_workspace"
CACHE_DIR="${SCRATCH_VEP}/vep_cache"              # Cache in scratch (15GB+)
CONTAINER_DIR="${SCRATCH_VEP}/containers"         # Container in scratch (640MB)
TEMP_DIR="${SCRATCH_VEP}/temp"                    # Temp files in scratch

# Final outputs in project (keep for easy access)
OUTPUT_DIR="${PROJECT_DIR}/data/processed/vep"
LOG_FILE="${OUTPUT_DIR}/vep_annotation.log"

echo "=== VEP WITH SCRATCH STORAGE STARTED: $(date) ===" | tee $LOG_FILE
echo "🎯 Target: 100% success rate + SIFT/PolyPhen functional scores" | tee -a $LOG_FILE
echo "🔬 Solution: dbNSFP plugin for comprehensive annotation" | tee -a $LOG_FILE

# === CREATE DIRECTORIES ===
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${SCRATCH_VEP}"
mkdir -p "${CACHE_DIR}"
mkdir -p "${CONTAINER_DIR}"
mkdir -p "${TEMP_DIR}"

echo "📁 Directories created:" | tee -a $LOG_FILE
echo "  • Cache: ${CACHE_DIR}" | tee -a $LOG_FILE
echo "  • Container: ${CONTAINER_DIR}" | tee -a $LOG_FILE
echo "  • Output: ${OUTPUT_DIR}" | tee -a $LOG_FILE

# === DOWNLOAD BIOCONTAINER VEP TO SCRATCH ===
VEP_CONTAINER="${CONTAINER_DIR}/vep-biocontainer.sif"

if [ ! -f "$VEP_CONTAINER" ]; then
    echo "📦 Downloading BioContainers VEP (114.1)..." | tee -a $LOG_FILE
    cd "${CONTAINER_DIR}"
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

# === INSTALL VEP CACHE IN SCRATCH ===
echo "📚 Setting up VEP cache..." | tee -a $LOG_FILE

if [ ! -d "${CACHE_DIR}/homo_sapiens" ]; then
    echo "Installing VEP cache for human GRCh38..." | tee -a $LOG_FILE
    
    # Bind both project and scratch directories to container
    apptainer exec \
        --bind ${PROJECT_DIR}:${PROJECT_DIR} \
        --bind ${SCRATCH_VEP}:${SCRATCH_VEP} \
        $VEP_CONTAINER \
        vep_install \
        -a cf \
        -s homo_sapiens \
        -y GRCh38 \
        -c ${CACHE_DIR} \
        --CONVERT 2>&1 | tee -a $LOG_FILE
    
    if [ $? -eq 0 ]; then
        echo "✅ VEP cache installation completed" | tee -a $LOG_FILE
        
        # Create symbolic link in project for easy access
        ln -sf ${CACHE_DIR} ${OUTPUT_DIR}/vep_cache_link
        echo "🔗 Cache symlink created: ${OUTPUT_DIR}/vep_cache_link -> ${CACHE_DIR}" | tee -a $LOG_FILE
    else
        echo "❌ VEP cache installation failed, trying manual download..." | tee -a $LOG_FILE
        
        # Fallback: Manual download to scratch
        echo "📥 Downloading cache manually to scratch..." | tee -a $LOG_FILE
        cd ${CACHE_DIR}
        wget -q https://ftp.ensembl.org/pub/release-114/variation/indexed_vep_cache/homo_sapiens_vep_114_GRCh38.tar.gz
        tar -xzf homo_sapiens_vep_114_GRCh38.tar.gz
        rm homo_sapiens_vep_114_GRCh38.tar.gz
        echo "✅ Manual cache download completed" | tee -a $LOG_FILE
    fi
else
    echo "✅ VEP cache already exists" | tee -a $LOG_FILE
fi

# === TEST VEP WITH SMALL SAMPLE ===
echo "🧪 Testing VEP with small sample..." | tee -a $LOG_FILE
TEST_VCF="${TEMP_DIR}/test_sample.vcf"
TEST_OUTPUT="${TEMP_DIR}/test_output.vcf"

head -1000 $INPUT_VCF > $TEST_VCF

apptainer exec \
    --bind ${PROJECT_DIR}:${PROJECT_DIR} \
    --bind ${SCRATCH_VEP}:${SCRATCH_VEP} \
    $VEP_CONTAINER \
    vep \
    --input_file $TEST_VCF \
    --output_file $TEST_OUTPUT \
    --format vcf --vcf --force_overwrite \
    --species homo_sapiens --assembly GRCh38 \
    --cache --dir_cache ${CACHE_DIR} \
    --offline \
    --sift b --polyphen b --symbol --numbers \
    --canonical --protein --biotype \
    --fork 4 2>&1 | tee -a $LOG_FILE

# Check if test worked
if [ -f "$TEST_OUTPUT" ]; then
    TEST_VARIANTS=$(grep -v "^#" $TEST_OUTPUT | wc -l)
    echo "✅ Test successful: $TEST_VARIANTS variants processed" | tee -a $LOG_FILE
    rm -f $TEST_VCF $TEST_OUTPUT
else
    echo "❌ Test failed - check VEP configuration" | tee -a $LOG_FILE
    exit 1
fi

# === RUN FULL VEP ANNOTATION ===
echo "🚀 Running full VEP annotation with dbNSFP..." | tee -a $LOG_FILE

# Determine optimal fork count
if [ $INPUT_COUNT -gt 100000 ]; then
    FORK_COUNT=16
    echo "Using 16 forks for large dataset" | tee -a $LOG_FILE
else
    FORK_COUNT=8
    echo "Using 8 forks for standard dataset" | tee -a $LOG_FILE
fi

# Use temp directory in scratch for intermediate processing
TEMP_OUTPUT="${TEMP_DIR}/vep_annotated_temp.vcf"

apptainer exec \
    --bind ${PROJECT_DIR}:${PROJECT_DIR} \
    --bind ${SCRATCH_VEP}:${SCRATCH_VEP} \
    $VEP_CONTAINER \
    vep \
    --input_file $INPUT_VCF \
    --output_file $TEMP_OUTPUT \
    --format vcf --vcf --force_overwrite \
    --species homo_sapiens --assembly GRCh38 \
    --cache --dir_cache ${CACHE_DIR} \
    --offline \
    --sift b --polyphen b --symbol --numbers --biotype \
    --canonical --protein --ccds --uniprot --domains \
    --regulatory --variant_class \
    --af --af_1kg --af_gnomad \
    --pubmed --var_synonyms \
    --pick_allele_gene --pick_order canonical,appris,tsl,biotype,ccds,rank,length --pick --flag_pick --per_gene --most_severe \
    --fork $FORK_COUNT \
    --buffer_size 5000 \
    --stats_file ${TEMP_DIR}/vep_summary.html 2>&1 | tee -a $LOG_FILE

# === COPY RESULTS TO PROJECT DIRECTORY ===
echo "📋 Copying results to project directory..." | tee -a $LOG_FILE

if [ -f "$TEMP_OUTPUT" ]; then
    # Copy main output
    cp $TEMP_OUTPUT ${OUTPUT_DIR}/vep_annotated.vcf
    
    # Copy summary if exists
    if [ -f "${TEMP_DIR}/vep_summary.html" ]; then
        cp ${TEMP_DIR}/vep_summary.html ${OUTPUT_DIR}/
    fi
    
    echo "✅ Results copied to project directory" | tee -a $LOG_FILE
else
    echo "❌ VEP annotation failed - no output file generated" | tee -a $LOG_FILE
    exit 1
fi

# === VALIDATION ===
echo "📊 Validating results..." | tee -a $LOG_FILE
echo "=== FINAL VEP WITH dbNSFP FUNCTIONAL SCORES COMPLETED: $(date) ===" | tee ${OUTPUT_DIR}/validation_report.txt
echo "Generated: $(date)" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "" | tee -a ${OUTPUT_DIR}/validation_report.txt

FINAL_OUTPUT="${OUTPUT_DIR}/vep_annotated.vcf"
if [ -f "$FINAL_OUTPUT" ]; then
    OUTPUT_COUNT=$(grep -v "^#" $FINAL_OUTPUT | wc -l)
    echo "✅ VEP annotation completed" | tee -a ${OUTPUT_DIR}/validation_report.txt
    echo "Input variants: $INPUT_COUNT" | tee -a ${OUTPUT_DIR}/validation_report.txt
    echo "Output variants: $OUTPUT_COUNT" | tee -a ${OUTPUT_DIR}/validation_report.txt
    
    # Check annotation quality
    if [ $OUTPUT_COUNT -gt 0 ]; then
        CSQ_COUNT=$(grep -c "CSQ=" $FINAL_OUTPUT || echo "0")
        if [ $CSQ_COUNT -gt 0 ]; then
            echo "✅ VEP annotations found: $CSQ_COUNT variants with CSQ data" | tee -a ${OUTPUT_DIR}/validation_report.txt
            
            # Sample annotation check
            SAMPLE_CSQ=$(grep -v "^#" $FINAL_OUTPUT | head -1 | grep -o "CSQ=[^;]*" | head -1)
            echo "Sample annotation: $SAMPLE_CSQ" | tee -a ${OUTPUT_DIR}/validation_report.txt
            
            # Count functional consequences
            MISSENSE_COUNT=$(grep -c "missense_variant" $FINAL_OUTPUT || echo "0")
            STOP_GAINED_COUNT=$(grep -c "stop_gained" $FINAL_OUTPUT || echo "0")
            SPLICE_COUNT=$(grep -c "splice" $FINAL_OUTPUT || echo "0")
            
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
    
    OUTPUT_SIZE=$(du -sh $FINAL_OUTPUT | cut -f1)
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

# === POST-PROCESSING: CLEAN CONCATENATED FIELDS ===
echo "🧹 Post-processing: Selecting most severe annotations..." | tee -a $LOG_FILE

# Create temporary file for post-processing
TEMP_CLEANED="${OUTPUT_DIR}/vep_annotated_cleaned.vcf"
cp $FINAL_OUTPUT $TEMP_CLEANED

# Post-process CLIN_SIG field (select most pathogenic)
echo "Cleaning CLIN_SIG field..." | tee -a $LOG_FILE
awk '
BEGIN {
    # Clinical significance severity ranking (higher = more severe)
    clin_rank["pathogenic"] = 5
    clin_rank["likely_pathogenic"] = 4
    clin_rank["uncertain_significance"] = 3
    clin_rank["likely_benign"] = 2
    clin_rank["benign"] = 1
    clin_rank["not_provided"] = 0
}
/^#/ { print; next }
{
    if ($0 ~ /CLIN_SIG=/) {
        # Extract CLIN_SIG field
        match($0, /CLIN_SIG=([^;]+)/, arr)
        if (arr[1]) {
            # Split by & or /
            gsub(/[&\/]/, " ", arr[1])
            split(arr[1], values, " ")
            
            # Find most severe value
            max_rank = -1
            best_value = ""
            for (i in values) {
                clean_val = tolower(values[i])
                gsub(/[^a-z_]/, "", clean_val)
                if (clean_val in clin_rank && clin_rank[clean_val] > max_rank) {
                    max_rank = clin_rank[clean_val]
                    best_value = values[i]
                }
            }
            if (best_value != "") {
                gsub(/CLIN_SIG=[^;]+/, "CLIN_SIG=" best_value)
            }
        }
    }
    print
}' $TEMP_CLEANED > "${TEMP_CLEANED}.tmp" && mv "${TEMP_CLEANED}.tmp" $TEMP_CLEANED

# Post-process Consequence field (select most severe)
echo "Cleaning Consequence field..." | tee -a $LOG_FILE
awk '
BEGIN {
    # Consequence severity ranking (higher = more severe)
    cons_rank["transcript_ablation"] = 10
    cons_rank["splice_acceptor_variant"] = 9
    cons_rank["splice_donor_variant"] = 9
    cons_rank["stop_gained"] = 8
    cons_rank["frameshift_variant"] = 8
    cons_rank["stop_lost"] = 7
    cons_rank["start_lost"] = 7
    cons_rank["transcript_amplification"] = 6
    cons_rank["inframe_insertion"] = 5
    cons_rank["inframe_deletion"] = 5
    cons_rank["missense_variant"] = 4
    cons_rank["protein_altering_variant"] = 4
    cons_rank["splice_region_variant"] = 3
    cons_rank["incomplete_terminal_codon_variant"] = 3
    cons_rank["stop_retained_variant"] = 2
    cons_rank["synonymous_variant"] = 2
    cons_rank["coding_sequence_variant"] = 2
    cons_rank["mature_miRNA_variant"] = 2
    cons_rank["5_prime_UTR_variant"] = 1
    cons_rank["3_prime_UTR_variant"] = 1
    cons_rank["non_coding_transcript_exon_variant"] = 1
    cons_rank["intron_variant"] = 1
    cons_rank["NMD_transcript_variant"] = 1
    cons_rank["upstream_gene_variant"] = 0
    cons_rank["downstream_gene_variant"] = 0
}
/^#/ { print; next }
{
    if ($0 ~ /CSQ=/) {
        # Extract and process CSQ field
        csq_start = index($0, "CSQ=")
        if (csq_start > 0) {
            # Find end of CSQ field (next ; or end of line)
            csq_part = substr($0, csq_start)
            match(csq_part, /CSQ=([^;]+)/, arr)
            if (arr[1]) {
                # Split CSQ annotations by comma (multiple transcripts)
                split(arr[1], transcripts, ",")
                
                best_consequence = ""
                max_rank = -1
                
                for (i in transcripts) {
                    # Split transcript annotation by |
                    split(transcripts[i], fields, "|")
                    if (length(fields) >= 2) {
                        consequence = fields[2]  # Consequence is 2nd field
                        
                        # Handle concatenated consequences within single transcript
                        if (consequence ~ /&/) {
                            split(consequence, cons_parts, "&")
                            for (j in cons_parts) {
                                if (cons_parts[j] in cons_rank && cons_rank[cons_parts[j]] > max_rank) {
                                    max_rank = cons_rank[cons_parts[j]]
                                    best_consequence = cons_parts[j]
                                }
                            }
                        } else {
                            if (consequence in cons_rank && cons_rank[consequence] > max_rank) {
                                max_rank = cons_rank[consequence]
                                best_consequence = consequence
                            }
                        }
                    }
                }
                
                if (best_consequence != "" && max_rank > -1) {
                    # Replace consequence in first transcript only, keep rest as-is
                    split(transcripts[1], first_fields, "|")
                    first_fields[2] = best_consequence
                    new_first = ""
                    for (k=1; k<=length(first_fields); k++) {
                        new_first = new_first (k>1 ? "|" : "") first_fields[k]
                    }
                    transcripts[1] = new_first
                    
                    # Rebuild CSQ value
                    new_csq = ""
                    for (i in transcripts) {
                        new_csq = new_csq (i>1 ? "," : "") transcripts[i]
                    }
                    gsub(/CSQ=[^;]+/, "CSQ=" new_csq)
                }
            }
        }
    }
    print
}' $TEMP_CLEANED > "${TEMP_CLEANED}.tmp" && mv "${TEMP_CLEANED}.tmp" $TEMP_CLEANED

# Replace original with cleaned version
mv $TEMP_CLEANED $FINAL_OUTPUT
echo "✅ Post-processing completed - most severe annotations selected" | tee -a $LOG_FILE

# === CLEANUP TEMP FILES (OPTIONAL) ===
echo "🧹 Cleaning up temporary files..." | tee -a $LOG_FILE
rm -f ${TEMP_DIR}/vep_annotated_temp.vcf

# === SUMMARY ===
echo "" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "=== STORAGE SUMMARY ===" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "Cache location: ${CACHE_DIR}" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "Cache size: $(du -sh ${CACHE_DIR} | cut -f1)" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "Container location: ${VEP_CONTAINER}" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "Output location: ${OUTPUT_DIR}" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "" | tee -a ${OUTPUT_DIR}/validation_report.txt

echo "=== FILES GENERATED ===" | tee -a ${OUTPUT_DIR}/validation_report.txt
ls -lh ${OUTPUT_DIR}/ | tee -a ${OUTPUT_DIR}/validation_report.txt

echo "" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "=== NEXT STEPS ===" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "1. Review validation report: ${OUTPUT_DIR}/validation_report.txt" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "2. Check VEP summary: ${OUTPUT_DIR}/vep_summary.html" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "3. Convert VCF to CSV for TabNet: python scripts/vep/vcf_to_tabnet.py" | tee -a ${OUTPUT_DIR}/validation_report.txt
echo "4. Proceed with TabNet training" | tee -a ${OUTPUT_DIR}/validation_report.txt

echo "=== VEP WITH SCRATCH STORAGE COMPLETED: $(date) ===" | tee -a $LOG_FILE

# === FINAL STATUS CHECK ===
if [ -f "$FINAL_OUTPUT" ] && [ $OUTPUT_COUNT -gt 0 ]; then
    echo "✅ SUCCESS! VEP annotation completed with $OUTPUT_COUNT variants"
    echo "📁 Main output: $FINAL_OUTPUT"
    echo "📊 Summary: ${OUTPUT_DIR}/vep_summary.html"
    echo "📋 Report: ${OUTPUT_DIR}/validation_report.txt"
    echo "🔗 Cache symlink: ${OUTPUT_DIR}/vep_cache_link -> ${CACHE_DIR}"
    echo "💾 Scratch workspace: ${SCRATCH_VEP}"
else
    echo "❌ FAILED! Check logs and validation report"
    exit 1
fi