#!/bin/bash
# ================================================================
# TabNet Validation Runner - Enhanced Version
# Properly sets up conda environment and runs validation
# ================================================================

set -euo pipefail

echo "🔥 TABNET VALIDATION RUNNER - ENHANCED VERSION"
echo "==============================================="
echo "Setting up environment and running validation..."
echo ""

# === CONFIGURATION ===
PROJECT_DIR="/u/aa107/uiuc-cancer-research"
VALIDATION_SCRIPT="${PROJECT_DIR}/src/model/validation/validate_tabnet.py"
RESULTS_DIR="${PROJECT_DIR}/results/validation"

# === ENVIRONMENT SETUP ===
echo "🔧 ENVIRONMENT SETUP"
echo "-------------------"

# Navigate to project directory
cd "${PROJECT_DIR}"
echo "Working directory: $(pwd)"

# Create results directory
mkdir -p "${RESULTS_DIR}"
echo "Results directory: ${RESULTS_DIR}"

# Load anaconda module (required on campus cluster)
echo "Loading anaconda module..."
module load anaconda3

# Initialize conda for shell script
echo "Initializing conda..."
eval "$(conda shell.bash hook)"

# Check if environment exists
if conda env list | grep -q "tabnet-prostate"; then
    echo "✅ Found tabnet-prostate environment"
else
    echo "❌ tabnet-prostate environment not found"
    echo "💡 Create it first by running:"
    echo "   ./src/model/tests/run_tabnet_tests.sh"
    exit 1
fi

# Activate the conda environment
echo "Activating tabnet-prostate environment..."
conda activate tabnet-prostate

# Verify activation
echo "✅ Environment activated:"
echo "  Python: $(which python)"
echo "  Python version: $(python --version)"
echo "  Conda env: $CONDA_DEFAULT_ENV"
echo ""

# === DEPENDENCY CHECK ===
echo "📦 CHECKING VALIDATION DEPENDENCIES"
echo "----------------------------------"

echo "Checking core packages..."
python -c "
import sys
try:
    import torch
    import pytorch_tabnet
    import sklearn
    import pandas as pd
    import numpy as np
    print('  ✅ All core packages available')
    print(f'  📦 PyTorch: {torch.__version__}')
    print(f'  📦 Sklearn: {sklearn.__version__}')
    print(f'  📦 Pandas: {pd.__version__}')
except ImportError as e:
    print(f'  ❌ Missing package: {e}')
    sys.exit(1)
"

if [ $? -ne 0 ]; then
    echo "❌ Missing dependencies - run test runner first"
    exit 1
fi

echo "✅ All validation dependencies available"
echo ""

# === ENHANCED DATASET CHECK ===
echo "🧬 CHECKING ENHANCED DATASET"
echo "---------------------------"

ENHANCED_DATASET="${PROJECT_DIR}/data/processed/tabnet_csv/prostate_variants_tabnet_enhanced.csv"

if [ -f "$ENHANCED_DATASET" ]; then
    DATASET_SIZE=$(du -h "$ENHANCED_DATASET" | cut -f1)
    DATASET_LINES=$(wc -l < "$ENHANCED_DATASET")
    echo "✅ Enhanced dataset found:"
    echo "  📁 File: $ENHANCED_DATASET"
    echo "  📏 Size: $DATASET_SIZE"
    echo "  📊 Lines: $DATASET_LINES"
    
    # Quick validation of dataset structure
    echo "🔍 Quick dataset structure check..."
    python -c "
import pandas as pd
import sys

try:
    df = pd.read_csv('$ENHANCED_DATASET', low_memory=False)
    print(f'  📊 Loaded: {len(df):,} variants × {len(df.columns)} features')
    
    # Check for critical features
    leakage_features = ['functional_pathogenicity', 'sift_confidence', 'polyphen_confidence']
    leakage_found = [f for f in leakage_features if f in df.columns]
    
    am_features = ['alphamissense_pathogenicity', 'alphamissense_class']
    am_missing = [f for f in am_features if f not in df.columns]
    
    if leakage_found:
        print(f'  ❌ Data leakage features found: {leakage_found}')
        sys.exit(1)
    elif am_missing:
        print(f'  ❌ AlphaMissense features missing: {am_missing}')
        sys.exit(1)
    else:
        print('  ✅ Dataset structure looks good')
        
except Exception as e:
    print(f'  ❌ Dataset check failed: {e}')
    sys.exit(1)
"
    
    if [ $? -ne 0 ]; then
        echo "❌ Dataset structure check failed"
        exit 1
    fi
    
else
    echo "❌ Enhanced dataset not found: $ENHANCED_DATASET"
    echo "💡 Run AlphaMissense enhancement first:"
    echo "   sbatch scripts/enhance/functional_enhancement/run_functional_imputation.sh"
    exit 1
fi

echo ""

# === RUN VALIDATION ===
echo "🧬 RUNNING ENHANCED TABNET VALIDATION"
echo "===================================="

if [ -f "$VALIDATION_SCRIPT" ]; then
    echo "Running validation script: $VALIDATION_SCRIPT"
    echo ""
    
    # Run the validation script
    python "$VALIDATION_SCRIPT"
    
    VALIDATION_EXIT_CODE=$?
    
    echo ""
    echo "=== VALIDATION RESULTS ==="
    
    if [ $VALIDATION_EXIT_CODE -eq 0 ]; then
        echo "✅ VALIDATION COMPLETED SUCCESSFULLY!"
        echo ""
        
        # Check for latest report
        LATEST_REPORT=$(find "${RESULTS_DIR}" -name "enhanced_tabnet_validation_*.json" -type f -printf '%T@ %p\n' 2>/dev/null | sort -n | tail -1 | cut -d' ' -f2- || echo "")
        
        if [ -n "$LATEST_REPORT" ] && [ -f "$LATEST_REPORT" ]; then
            echo "📋 Latest validation report: $LATEST_REPORT"
            
            # Extract key results
            echo "🔍 Key Validation Results:"
            python -c "
import json
try:
    with open('$LATEST_REPORT', 'r') as f:
        results = json.load(f)
    
    # Summary
    summary = results.get('summary', {})
    data_leakage = summary.get('data_leakage_detected', False)
    am_integrated = summary.get('alphamissense_integrated', False)
    
    print(f'  Data leakage detected: {\"❌ YES\" if data_leakage else \"✅ NO\"}')
    print(f'  AlphaMissense integrated: {\"✅ YES\" if am_integrated else \"❌ NO\"}')
    
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
    print(f'  ⚠️  Could not parse report: {e}')
"
        else
            echo "⚠️  No validation report found in ${RESULTS_DIR}"
        fi
        
        echo ""
        echo "🎯 NEXT STEPS:"
        echo "  1. Review validation report (if generated)"
        echo "  2. If no data leakage detected, proceed with training:"
        echo "     python src/model/tabnet_prostate_variant_classifier.py"
        echo "  3. Or submit cluster training job"
        
        exit 0
        
    else
        echo "❌ VALIDATION FAILED (exit code: $VALIDATION_EXIT_CODE)"
        echo ""
        echo "🔧 TROUBLESHOOTING:"
        echo "  1. Check the validation output above for specific errors"
        echo "  2. Ensure enhanced dataset is properly formatted"
        echo "  3. Verify all dependencies are installed correctly"
        echo "  4. Try running tests first: ./src/model/tests/run_tabnet_tests.sh"
        
        exit 1
    fi
    
else
    echo "❌ Validation script not found: $VALIDATION_SCRIPT"
    echo "💡 Make sure you're in the project root directory"
    exit 1
fi