#!/bin/bash
#SBATCH --job-name=cancer_data_audit
#SBATCH --partition=IllinoisComputes
#SBATCH --account=aa107-ic
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/u/aa107/uiuc-cancer-research/scripts/validation/audit_%j.out
#SBATCH --error=/u/aa107/uiuc-cancer-research/scripts/validation/audit_%j.err

# =============================================================================
# PROSTATE CANCER DATA QUALITY AUDIT RUNNER
# Phase 1: Comprehensive Column Analysis
# =============================================================================

set -euo pipefail  # Exit on error, undefined vars, pipe failures

echo "🧬 PROSTATE CANCER DATA QUALITY AUDIT"
echo "Phase 1: Comprehensive Column Analysis"
echo "======================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "Start time: $(date)"
echo ""

# =============================================================================
# ENVIRONMENT SETUP
# =============================================================================

echo "🔧 SETTING UP ENVIRONMENT"
echo "========================="

# Set project root
export PROJECT_ROOT="/u/aa107/uiuc-cancer-research"
cd $PROJECT_ROOT

# Create validation scripts directory
VALIDATION_DIR="$PROJECT_ROOT/scripts/validation"
mkdir -p $VALIDATION_DIR

# Create results directory
RESULTS_DIR="$PROJECT_ROOT/results/validation/phase1_audit"
mkdir -p $RESULTS_DIR

echo "✅ Project root: $PROJECT_ROOT"
echo "✅ Validation directory: $VALIDATION_DIR"
echo "✅ Results directory: $RESULTS_DIR"

# =============================================================================
# DEPENDENCY CHECK
# =============================================================================

echo ""
echo "📋 CHECKING DEPENDENCIES"
echo "========================"

# Load required modules
module load anaconda3

# Initialize conda for bash (required for compute nodes)
echo "Initializing conda..."
eval "$(conda shell.bash hook)"

# Activate conda environment (create if doesn't exist)
CONDA_ENV="tabnet-prostate"
if conda info --envs | grep -q $CONDA_ENV; then
    echo "✅ Activating existing environment: $CONDA_ENV"
    conda activate $CONDA_ENV
else
    echo "🔨 Creating new conda environment: $CONDA_ENV"
    conda create -y -n $CONDA_ENV python=3.9
    conda activate $CONDA_ENV
fi

# Install required packages
echo "📦 Installing/updating Python packages..."
pip install --upgrade pip
pip install pandas numpy matplotlib seaborn scikit-learn
pip install jupyter notebook ipywidgets
pip install openpyxl xlsxwriter  # For Excel file handling

# Verify environment activation
echo "✅ Environment verification:"
echo "  Python location: $(which python)"
echo "  Python version: $(python --version)"
echo "  Conda environment: $CONDA_DEFAULT_ENV"

# =============================================================================
# FILE VALIDATION
# =============================================================================

echo ""
echo "📁 VALIDATING INPUT FILES"
echo "========================="

# Define file paths
ENHANCED_FILE="$PROJECT_ROOT/data/processed/tabnet_csv/prostate_variants_tabnet_enhanced.csv"
MERGED_FILE="$PROJECT_ROOT/data/processed/merged/merged_prostate_variants.csv"

# Check enhanced file (required)
if [ -f "$ENHANCED_FILE" ]; then
    FILE_SIZE=$(du -sh "$ENHANCED_FILE" | cut -f1)
    ROW_COUNT=$(wc -l < "$ENHANCED_FILE")
    echo "✅ Enhanced dataset found: $FILE_SIZE, $((ROW_COUNT-1)) variants"
else
    echo "❌ ERROR: Enhanced dataset not found: $ENHANCED_FILE"
    echo "💡 This file should be generated after VEP annotation and functional enhancement"
    exit 1
fi

# Check merged file (optional, for comparison)
if [ -f "$MERGED_FILE" ]; then
    FILE_SIZE=$(du -sh "$MERGED_FILE" | cut -f1)
    ROW_COUNT=$(wc -l < "$MERGED_FILE")
    echo "✅ Merged dataset found: $FILE_SIZE, $((ROW_COUNT-1)) variants"
else
    echo "⚠️  Merged dataset not found (comparison analysis will be limited)"
fi

# =============================================================================
# CREATE AUDIT SCRIPT
# =============================================================================

echo ""
echo "📝 CREATING AUDIT SCRIPT"
echo "========================"

# Create the comprehensive audit script in the validation directory
AUDIT_SCRIPT="$VALIDATION_DIR/comprehensive_column_audit.py"

if [ ! -f "$AUDIT_SCRIPT" ]; then
    echo "⚠️  Audit script not found. Please ensure comprehensive_column_audit.py is in:"
    echo "   $VALIDATION_DIR/"
    echo ""
    echo "💡 You can create it by copying the script from the project artifacts or run:"
    echo "   cp /path/to/comprehensive_column_audit.py $VALIDATION_DIR/"
    exit 1
else
    echo "✅ Audit script found: $AUDIT_SCRIPT"
fi

# Make script executable
chmod +x "$AUDIT_SCRIPT"

# =============================================================================
# SYSTEM RESOURCE CHECK
# =============================================================================

echo ""
echo "💻 SYSTEM RESOURCE CHECK"
echo "========================"

echo "Available CPU cores: $(nproc)"
echo "Available memory: $(free -h | grep '^Mem:' | awk '{print $2}')"
echo "Available disk space: $(df -h $PROJECT_ROOT | tail -1 | awk '{print $4}')"

# Check memory requirements for large dataset
ENHANCED_SIZE_MB=$(du -m "$ENHANCED_FILE" | cut -f1)
ESTIMATED_MEMORY_MB=$((ENHANCED_SIZE_MB * 5))  # Estimate 5x file size for processing

echo "Dataset size: ${ENHANCED_SIZE_MB}MB"
echo "Estimated memory needed: ${ESTIMATED_MEMORY_MB}MB"

if [ $ESTIMATED_MEMORY_MB -gt 30000 ]; then
    echo "⚠️  Large dataset detected. Processing may take significant time."
fi

# =============================================================================
# RUN COMPREHENSIVE AUDIT
# =============================================================================

echo ""
echo "🚀 RUNNING COMPREHENSIVE DATA QUALITY AUDIT"
echo "==========================================="

# Set Python path
export PYTHONPATH="$PROJECT_ROOT:${PYTHONPATH:-}"

# Run the audit with error handling
echo "Starting audit execution..."
START_TIME=$(date +%s)

if python "$AUDIT_SCRIPT"; then
    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    echo ""
    echo "✅ AUDIT COMPLETED SUCCESSFULLY!"
    echo "Execution time: ${DURATION} seconds"
    AUDIT_STATUS="SUCCESS"
else
    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    echo ""
    echo "❌ AUDIT FAILED!"
    echo "Execution time: ${DURATION} seconds"
    echo "Check error logs for details."
    AUDIT_STATUS="FAILED"
fi

# =============================================================================
# RESULTS VALIDATION
# =============================================================================

echo ""
echo "📊 VALIDATING RESULTS"
echo "===================="

# Check if output files were generated
JSON_RESULT="$RESULTS_DIR/comprehensive_audit_results.json"
REPORT_RESULT="$RESULTS_DIR/data_quality_audit_report.txt"
CSV_RESULT="$RESULTS_DIR/column_summary.csv"

RESULTS_GENERATED=0

if [ -f "$JSON_RESULT" ]; then
    FILE_SIZE=$(du -sh "$JSON_RESULT" | cut -f1)
    echo "✅ JSON results: $FILE_SIZE"
    RESULTS_GENERATED=$((RESULTS_GENERATED + 1))
else
    echo "❌ JSON results not generated"
fi

if [ -f "$REPORT_RESULT" ]; then
    FILE_SIZE=$(du -sh "$REPORT_RESULT" | cut -f1)
    echo "✅ Text report: $FILE_SIZE"
    RESULTS_GENERATED=$((RESULTS_GENERATED + 1))
else
    echo "❌ Text report not generated"
fi

if [ -f "$CSV_RESULT" ]; then
    FILE_SIZE=$(du -sh "$CSV_RESULT" | cut -f1)
    echo "✅ CSV summary: $FILE_SIZE"
    RESULTS_GENERATED=$((RESULTS_GENERATED + 1))
else
    echo "❌ CSV summary not generated"
fi

# =============================================================================
# GENERATE EXECUTION SUMMARY
# =============================================================================

echo ""
echo "📋 GENERATING EXECUTION SUMMARY"
echo "==============================="

SUMMARY_FILE="$RESULTS_DIR/audit_execution_summary.txt"

cat > "$SUMMARY_FILE" << EOF
PROSTATE CANCER DATA QUALITY AUDIT - EXECUTION SUMMARY
=====================================================

Job Information:
- Job ID: $SLURM_JOB_ID
- Node: $SLURMD_NODENAME
- Start Time: $(date)
- Execution Duration: ${DURATION} seconds
- Status: $AUDIT_STATUS

Input Files:
- Enhanced Dataset: $ENHANCED_FILE
- Enhanced Dataset Size: $(du -sh "$ENHANCED_FILE" | cut -f1)
- Enhanced Dataset Variants: $(echo $(($(wc -l < "$ENHANCED_FILE")-1)))

Output Files Generated: $RESULTS_GENERATED/3
- JSON Results: $([ -f "$JSON_RESULT" ] && echo "✅" || echo "❌")
- Text Report: $([ -f "$REPORT_RESULT" ] && echo "✅" || echo "❌")
- CSV Summary: $([ -f "$CSV_RESULT" ] && echo "✅" || echo "❌")

Environment:
- Python Version: $(python --version)
- Conda Environment: $CONDA_ENV
- Project Root: $PROJECT_ROOT

Next Steps:
1. Review generated reports in: $RESULTS_DIR
2. Analyze priority recommendations
3. Plan Phase 2: Legitimate vs. Problematic Classification
4. Consider VEP parameter optimization based on findings

Contact: aa107@illinois.edu
EOF

echo "📄 Execution summary saved: $SUMMARY_FILE"

# =============================================================================
# FINAL STATUS AND RECOMMENDATIONS
# =============================================================================

echo ""
echo "🎯 FINAL STATUS"
echo "==============="

if [ "$AUDIT_STATUS" = "SUCCESS" ] && [ $RESULTS_GENERATED -eq 3 ]; then
    echo "✅ SUCCESS! PHASE 1 AUDIT COMPLETED!"
    echo ""
    echo "📊 Key deliverables generated:"
    echo "   - Comprehensive column analysis (89 columns)"
    echo "   - Concatenation pattern identification"
    echo "   - Clinical relevance scoring"
    echo "   - Pipeline source attribution"
    echo "   - Priority recommendations"
    echo ""
    echo "🔍 Next actions:"
    echo "   1. Review report: $REPORT_RESULT"
    echo "   2. Analyze CSV summary: $CSV_RESULT"
    echo "   3. Implement priority recommendations"
    echo "   4. Proceed to Phase 2 analysis"
    echo ""
    echo "📁 All results in: $RESULTS_DIR"
    
elif [ "$AUDIT_STATUS" = "SUCCESS" ]; then
    echo "⚠️  AUDIT COMPLETED WITH ISSUES"
    echo "Some output files were not generated. Check logs for details."
    echo "📁 Partial results in: $RESULTS_DIR"
    
else
    echo "❌ AUDIT FAILED"
    echo "Check error logs and fix issues before rerunning."
    echo ""
    echo "💡 Common issues:"
    echo "   - Missing input files"
    echo "   - Insufficient memory"
    echo "   - Python package conflicts"
    echo "   - File permission issues"
    echo ""
    echo "📞 Support: aa107@illinois.edu"
fi

echo ""
echo "🕐 End time: $(date)"
echo ""

# Set exit code based on audit status
if [ "$AUDIT_STATUS" = "SUCCESS" ] && [ $RESULTS_GENERATED -eq 3 ]; then
    exit 0
else
    exit 1
fi