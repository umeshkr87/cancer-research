#!/bin/bash
#SBATCH --job-name=tabnet_validation_enhanced
#SBATCH --partition=IllinoisComputes-GPU
#SBATCH --account=aa107-ic
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --gres=gpu:1
#SBATCH --time=01:00:00
#SBATCH --output=/u/aa107/uiuc-cancer-research/logs/tabnet_validation_enhanced_%j.out
#SBATCH --error=/u/aa107/uiuc-cancer-research/logs/tabnet_validation_enhanced_%j.err

# ================================================================
# Enhanced TabNet Validation Script
# Validates AlphaMissense integration and data leakage elimination
# ================================================================

set -euo pipefail

echo "🧬 ENHANCED TABNET VALIDATION ON CAMPUS CLUSTER"
echo "================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "GPU Info: $(nvidia-smi --query-gpu=name --format=csv,noheader,nounits)"
echo ""

# === CONFIGURATION ===
PROJECT_DIR="/u/aa107/uiuc-cancer-research"
LOG_DIR="${PROJECT_DIR}/logs"
RESULTS_DIR="${PROJECT_DIR}/results/validation"
ENHANCED_DATASET="${PROJECT_DIR}/data/processed/tabnet_csv/prostate_variants_tabnet_enhanced.csv"

# Create necessary directories
mkdir -p "${LOG_DIR}"
mkdir -p "${RESULTS_DIR}"

# === ENVIRONMENT SETUP ===
echo "🔧 ENVIRONMENT SETUP"
echo "-----------------"

# Navigate to project root
cd "${PROJECT_DIR}"
echo "Working directory: $(pwd)"

# Load anaconda module
echo "Loading anaconda module..."
module load anaconda3

# Initialize conda for batch environment
echo "Initializing conda..."
eval "$(conda shell.bash hook)"

# Activate environment
echo "Activating tabnet-prostate environment..."
conda activate tabnet-prostate || {
    echo "❌ Failed to activate conda environment"
    echo "💡 Create environment first: conda create -n tabnet-prostate python=3.11"
    exit 1
}

# Verify environment
echo "✅ Environment verification:"
echo "  Python: $(which python)"
echo "  Python version: $(python --version)"
echo "  Conda env: $CONDA_DEFAULT_ENV"
echo ""

# === PREREQUISITE CHECKS ===
echo "📋 PREREQUISITE CHECKS"
echo "--------------------"

# Check enhanced dataset exists
if [ ! -f "${ENHANCED_DATASET}" ]; then
    echo "❌ Enhanced dataset not found: ${ENHANCED_DATASET}"
    echo "💡 Run AlphaMissense enhancement first:"
    echo "   sbatch scripts/enhance/functional_enhancement/run_functional_imputation.sh"
    exit 1
else
    DATASET_SIZE=$(du -h "${ENHANCED_DATASET}" | cut -f1)
    DATASET_LINES=$(wc -l < "${ENHANCED_DATASET}")
    echo "✅ Enhanced dataset found: ${ENHANCED_DATASET}"
    echo "  Size: ${DATASET_SIZE}"
    echo "  Lines: ${DATASET_LINES}"
fi

# Check Python dependencies
echo "🐍 Checking Python dependencies..."
python -c "import torch; print(f'PyTorch: {torch.__version__}')" || {
    echo "❌ PyTorch not available"
    exit 1
}

python -c "import pytorch_tabnet; print('TabNet: OK')" || {
    echo "❌ PyTorch TabNet not available"
    echo "💡 Install: pip install pytorch-tabnet"
    exit 1
}

python -c "import pandas, numpy, sklearn; print('Core packages: OK')" || {
    echo "❌ Core packages not available"
    exit 1
}

echo "✅ All dependencies available"
echo ""

# === STEP 1: ENVIRONMENT TEST ===
echo "🧪 STEP 1: ENVIRONMENT TESTING"
echo "=============================="

echo "Running environment tests..."
python src/model/tests/test_environment.py

if [ $? -eq 0 ]; then
    echo "✅ Environment tests passed"
else
    echo "❌ Environment tests failed"
    echo "Check logs for details"
    exit 1
fi
echo ""

# === STEP 2: ENHANCED DATASET VALIDATION ===
echo "🔍 STEP 2: ENHANCED DATASET VALIDATION"
echo "===================================="

echo "Validating enhanced dataset structure..."
python -c "
import pandas as pd
import sys

try:
    df = pd.read_csv('${ENHANCED_DATASET}')
    print(f'Dataset loaded: {len(df):,} variants × {len(df.columns)} features')
    
    # Check for data leakage features (should be ABSENT)
    leakage_features = ['functional_pathogenicity', 'sift_confidence', 'polyphen_confidence']
    leakage_found = [f for f in leakage_features if f in df.columns]
    
    if leakage_found:
        print(f'❌ CRITICAL: Data leakage features found: {leakage_found}')
        sys.exit(1)
    else:
        print('✅ No data leakage features detected')
    
    # Check for AlphaMissense features (should be PRESENT)
    am_features = ['alphamissense_pathogenicity', 'alphamissense_class']
    missing_am = [f for f in am_features if f not in df.columns]
    
    if missing_am:
        print(f'❌ CRITICAL: AlphaMissense features missing: {missing_am}')
        sys.exit(1)
    else:
        print('✅ AlphaMissense features present')
    
    # Check coverage
    am_coverage = df['alphamissense_pathogenicity'].notna().sum()
    coverage_rate = am_coverage / len(df) * 100
    print(f'📊 AlphaMissense coverage: {am_coverage:,} variants ({coverage_rate:.1f}%)')
    
    if coverage_rate < 30:
        print('⚠️  Low AlphaMissense coverage')
    else:
        print('✅ Good AlphaMissense coverage')
        
    print('✅ Enhanced dataset validation passed')
    
except Exception as e:
    print(f'❌ Dataset validation failed: {e}')
    sys.exit(1)
"

if [ $? -eq 0 ]; then
    echo "✅ Enhanced dataset validation passed"
else
    echo "❌ Enhanced dataset validation failed"
    exit 1
fi
echo ""

# === STEP 3: COMPREHENSIVE TABNET VALIDATION ===
echo "🔥 STEP 3: COMPREHENSIVE TABNET VALIDATION"
echo "========================================"

echo "Running comprehensive TabNet validation..."
python src/model/validation/validate_tabnet.py

VALIDATION_EXIT_CODE=$?

if [ $VALIDATION_EXIT_CODE -eq 0 ]; then
    echo "✅ TabNet validation completed successfully"
else
    echo "❌ TabNet validation failed (exit code: $VALIDATION_EXIT_CODE)"
    echo "Check validation logs for details"
fi
echo ""

# === STEP 4: RESULTS SUMMARY ===
echo "📊 STEP 4: VALIDATION RESULTS SUMMARY"
echo "==================================="

# Find the latest validation report
LATEST_REPORT=$(find "${RESULTS_DIR}" -name "enhanced_tabnet_validation_*.json" -type f -printf '%T@ %p\n' | sort -n | tail -1 | cut -d' ' -f2-)

if [ -f "$LATEST_REPORT" ]; then
    echo "📋 Latest validation report: $LATEST_REPORT"
    
    # Extract key results
    echo "🔍 Key Results:"
    python -c "
import json
import sys

try:
    with open('$LATEST_REPORT', 'r') as f:
        results = json.load(f)
    
    # Summary
    summary = results.get('summary', {})
    print(f\"  Data leakage detected: {'❌ YES' if summary.get('data_leakage_detected') else '✅ NO'}\")
    print(f\"  AlphaMissense integrated: {'✅ YES' if summary.get('alphamissense_integrated') else '❌ NO'}\")
    
    # Performance
    if 'tabnet_validation' in results:
        tabnet_acc = results['tabnet_validation']['mean_accuracy']
        print(f'  TabNet accuracy: {tabnet_acc:.3f}')
        
        if tabnet_acc > 0.95:
            print('  Status: ❌ SUSPICIOUS - Investigate data leakage')
        elif tabnet_acc > 0.75:
            print('  Status: ✅ EXCELLENT - Ready for production')
        else:
            print('  Status: ✅ GOOD - Realistic performance')
    
except Exception as e:
    print(f'Error reading report: {e}')
    sys.exit(1)
"
else
    echo "⚠️  No validation report found"
fi

# Log file locations
echo ""
echo "📁 OUTPUT LOCATIONS:"
echo "  Validation report: ${RESULTS_DIR}/"
echo "  Job output: ${LOG_DIR}/tabnet_validation_enhanced_${SLURM_JOB_ID}.out"
echo "  Job errors: ${LOG_DIR}/tabnet_validation_enhanced_${SLURM_JOB_ID}.err"

# === FINAL STATUS ===
echo ""
echo "=== VALIDATION COMPLETED: $(date) ==="

if [ $VALIDATION_EXIT_CODE -eq 0 ]; then
    echo "✅ SUCCESS: Enhanced TabNet validation completed"
    echo ""
    echo "🎯 NEXT STEPS:"
    echo "  1. Review validation report in ${RESULTS_DIR}/"
    echo "  2. If no data leakage detected, proceed with full training"
    echo "  3. Use enhanced dataset for production TabNet model"
    echo "  4. Submit training job: sbatch scripts/train_tabnet_enhanced.sbatch"
    
    exit 0
else
    echo "❌ FAILED: Enhanced TabNet validation failed"
    echo ""
    echo "🔧 TROUBLESHOOTING:"
    echo "  1. Check job output logs for detailed errors"
    echo "  2. Verify enhanced dataset was created correctly"
    echo "  3. Ensure all dependencies are installed"
    echo "  4. Contact support if issues persist"
    
    exit 1
fi