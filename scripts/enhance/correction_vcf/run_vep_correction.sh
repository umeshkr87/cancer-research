#!/bin/bash
#SBATCH --job-name=vep_correction
#SBATCH --partition=IllinoisComputes
#SBATCH --account=aa107-ic
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:30:00
#SBATCH --output=/u/aa107/uiuc-cancer-research/scripts/enhance/correction_vcf/correction_%j.out
#SBATCH --error=/u/aa107/uiuc-cancer-research/scripts/enhance/correction_vcf/correction_%j.err

# =============================================================================
# VEP CONCATENATION CORRECTION RUNNER
# Fixes concatenated fields in VEP annotated VCF files for TabNet ML training
# =============================================================================

set -euo pipefail  # Exit on error, undefined vars, pipe failures

echo "🧬 VEP CONCATENATION CORRECTION"
echo "=============================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "Start time: $(date)"
echo ""

# =============================================================================
# ENVIRONMENT SETUP
# =============================================================================

echo "🔧 SETTING UP ENVIRONMENT"
echo "========================="

# Set project root and navigate
export PROJECT_ROOT="/u/aa107/uiuc-cancer-research"
cd $PROJECT_ROOT

# Script locations
SCRIPT_DIR="$PROJECT_ROOT/scripts/enhance/correction_vcf"
CORRECTION_SCRIPT="$SCRIPT_DIR/post_process_vep_concatenation.py"

# Data file paths
INPUT_VCF="/u/aa107/uiuc-cancer-research/data/processed/vep/vep_annotated.vcf"
OUTPUT_VCF="/u/aa107/uiuc-cancer-research/data/processed/vep/vep_annotated_clean.vcf"
BACKUP_VCF="/u/aa107/uiuc-cancer-research/data/processed/vep/vep_annotated_original_backup.vcf"

echo "✅ Project root: $PROJECT_ROOT"
echo "✅ Script directory: $SCRIPT_DIR"
echo "✅ Input VCF: $INPUT_VCF"
echo "✅ Output VCF: $OUTPUT_VCF"

# =============================================================================
# PRE-PROCESSING VALIDATION
# =============================================================================

echo ""
echo "🔍 PRE-PROCESSING VALIDATION"
echo "==========================="

# Check if input VCF exists
if [ ! -f "$INPUT_VCF" ]; then
    echo "❌ ERROR: Input VCF not found: $INPUT_VCF"
    echo "💡 Run VEP annotation pipeline first"
    exit 1
fi

# Get file information
FILE_SIZE=$(du -sh "$INPUT_VCF" | cut -f1)
LINE_COUNT=$(wc -l < "$INPUT_VCF")
VARIANT_COUNT=$(grep -v "^#" "$INPUT_VCF" | wc -l)

echo "✅ Input VCF found"
echo "   📏 Size: $FILE_SIZE"
echo "   📊 Total lines: $((LINE_COUNT))"
echo "   🧬 Variants: $((VARIANT_COUNT))"

# Validate VCF structure
if grep -q "##INFO=<ID=CSQ" "$INPUT_VCF"; then
    echo "✅ CSQ annotations detected"
    CSQ_FIELDS=$(grep "##INFO=<ID=CSQ" "$INPUT_VCF" | grep -o "Format: [^\"]*" | cut -d: -f2- | tr '|' '\n' | wc -l)
    echo "   🏷️  CSQ fields: $CSQ_FIELDS"
else
    echo "❌ ERROR: No CSQ annotations found in VCF"
    echo "💡 Ensure VEP annotation was successful"
    exit 1
fi

# Check for concatenation patterns in key fields
echo ""
echo "🔍 DETECTING CONCATENATION PATTERNS"
echo "==================================="

CONCAT_SAMPLES=$(grep -v "^#" "$INPUT_VCF" | head -1000 | grep -o "CSQ=[^;]*" | grep -c "&" || echo "0")
echo "📊 Concatenation samples found in first 1000 variants: $CONCAT_SAMPLES"

if [ $CONCAT_SAMPLES -eq 0 ]; then
    echo "⚠️  WARNING: No concatenation detected in sample. Script may not be needed."
    echo "💡 Continuing anyway for thorough processing..."
fi

# =============================================================================
# BACKUP ORIGINAL FILE
# =============================================================================

echo ""
echo "💾 CREATING BACKUP"
echo "=================="

if [ ! -f "$BACKUP_VCF" ]; then
    echo "Creating backup of original VCF..."
    cp "$INPUT_VCF" "$BACKUP_VCF"
    echo "✅ Backup created: $BACKUP_VCF"
else
    echo "✅ Backup already exists: $BACKUP_VCF"
fi

# =============================================================================
# SYSTEM RESOURCE CHECK
# =============================================================================

echo ""
echo "💻 SYSTEM RESOURCE CHECK"
echo "========================"

echo "Available CPU cores: $(nproc)"
echo "Available memory: $(free -h | grep '^Mem:' | awk '{print $2}')"
echo "Available disk space: $(df -h $PROJECT_ROOT | tail -1 | awk '{print $4}')"

# Estimate processing requirements
VCF_SIZE_MB=$(du -m "$INPUT_VCF" | cut -f1)
ESTIMATED_TIME_MIN=$((VCF_SIZE_MB / 10))  # Rough estimate: 10MB per minute

echo "VCF file size: ${VCF_SIZE_MB}MB"
echo "Estimated processing time: ~${ESTIMATED_TIME_MIN} minutes"

# =============================================================================
# EXECUTE CORRECTION SCRIPT
# =============================================================================

echo ""
echo "🚀 RUNNING VEP CONCATENATION CORRECTION"
echo "======================================"

# Verify correction script exists
if [ ! -f "$CORRECTION_SCRIPT" ]; then
    echo "❌ ERROR: Correction script not found: $CORRECTION_SCRIPT"
    exit 1
fi

# Make script executable
chmod +x "$CORRECTION_SCRIPT"

# Run the correction with timing
echo "Starting VEP correction..."
START_TIME=$(date +%s)

if python3 "$CORRECTION_SCRIPT"; then
    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    echo ""
    echo "✅ VEP CORRECTION COMPLETED SUCCESSFULLY!"
    echo "⏱️  Execution time: $((DURATION / 60)) minutes $((DURATION % 60)) seconds"
    CORRECTION_STATUS="SUCCESS"
else
    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    echo ""
    echo "❌ VEP CORRECTION FAILED!"
    echo "⏱️  Execution time: $((DURATION / 60)) minutes $((DURATION % 60)) seconds"
    echo "Check error logs for details."
    CORRECTION_STATUS="FAILED"
    exit 1
fi

# =============================================================================
# POST-PROCESSING VALIDATION
# =============================================================================

echo ""
echo "🔍 POST-PROCESSING VALIDATION"
echo "============================="

if [ -f "$OUTPUT_VCF" ]; then
    OUTPUT_SIZE=$(du -sh "$OUTPUT_VCF" | cut -f1)
    OUTPUT_LINES=$(wc -l < "$OUTPUT_VCF")
    OUTPUT_VARIANTS=$(grep -v "^#" "$OUTPUT_VCF" | wc -l)
    
    echo "✅ Output VCF created successfully"
    echo "   📏 Size: $OUTPUT_SIZE"
    echo "   📊 Total lines: $OUTPUT_LINES"
    echo "   🧬 Variants: $OUTPUT_VARIANTS"
    
    # Compare variant counts
    if [ $VARIANT_COUNT -eq $OUTPUT_VARIANTS ]; then
        echo "✅ Variant count preserved: $VARIANT_COUNT variants"
    else
        echo "⚠️  WARNING: Variant count changed from $VARIANT_COUNT to $OUTPUT_VARIANTS"
    fi
    
    # Quick concatenation check
    POST_CONCAT=$(grep -v "^#" "$OUTPUT_VCF" | head -1000 | grep -o "CSQ=[^;]*" | grep -c "&" || echo "0")
    echo "📊 Concatenation samples remaining: $POST_CONCAT (was $CONCAT_SAMPLES)"
    
    if [ $POST_CONCAT -lt $CONCAT_SAMPLES ]; then
        echo "✅ Concatenation reduction successful"
    else
        echo "⚠️  WARNING: Concatenation not reduced as expected"
    fi
    
else
    echo "❌ ERROR: Output VCF not created!"
    exit 1
fi

# =============================================================================
# SUMMARY AND NEXT STEPS
# =============================================================================

echo ""
echo "📈 CORRECTION SUMMARY"
echo "===================="
echo "Status: $CORRECTION_STATUS"
echo "Input file: $INPUT_VCF ($FILE_SIZE)"
echo "Output file: $OUTPUT_VCF ($OUTPUT_SIZE)"
echo "Processing time: $((DURATION / 60)) min $((DURATION % 60)) sec"
echo "Concatenation reduction: $CONCAT_SAMPLES → $POST_CONCAT"

echo ""
echo "🎯 NEXT STEPS"
echo "============"
echo "1. Run validation audit on cleaned VCF:"
echo "   cd $PROJECT_ROOT"
echo "   python scripts/validation/comprehensive_column_audit.py"
echo ""
echo "2. Compare before/after concatenation rates"
echo ""
echo "3. If validation passes, proceed with TabNet training:"
echo "   python scripts/enhance/tabnet_vcf_to_csv/vcf_to_tabnet_csv.py"
echo ""
echo "4. Files generated:"
echo "   ✅ Original backup: $BACKUP_VCF"
echo "   ✅ Cleaned VCF: $OUTPUT_VCF"

echo ""
echo "🎉 VEP CONCATENATION CORRECTION PIPELINE COMPLETED!"
echo "Job finished: $(date)"