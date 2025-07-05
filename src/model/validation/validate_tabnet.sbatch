#!/bin/bash
#SBATCH --job-name=enhanced_tabnet_validation
#SBATCH --account=aces
#SBATCH --partition=gpu-a100
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --gres=gpu:1
#SBATCH --time=1:00:00
#SBATCH --output=enhanced_tabnet_validation_%j.out
#SBATCH --error=enhanced_tabnet_validation_%j.err

echo "🧬 ENHANCED TABNET VALIDATION ON CAMPUS CLUSTER"
echo "================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "Started: $(date)"
echo "GPU Info: $(nvidia-smi --query-gpu=name --format=csv,noheader,nounits)"
echo ""

# === CONFIGURATION ===
PROJECT_DIR="/u/aa107/uiuc-cancer-research"
VALIDATION_SCRIPT="${PROJECT_DIR}/src/model/validation/validate_tabnet.py"
CLEAN_DATASET="${PROJECT_DIR}/data/processed/tabnet_csv/prostate_variants_tabnet_clean.csv"

# === ENVIRONMENT SETUP ===
echo "🔧 ENVIRONMENT SETUP"
echo "-----------------"

# Navigate to project directory
cd "${PROJECT_DIR}"
echo "Working directory: $(pwd)"

# Load anaconda module
echo "Loading anaconda module..."
module load anaconda3

# Initialize conda
echo "Initializing conda..."
eval "$(conda shell.bash hook)"

# Activate environment
echo "Activating tabnet-prostate environment..."
conda activate tabnet-prostate

# Verify environment
echo "✅ Environment verification:"
echo "  Python: $(which python)"
echo "  Python version: $(python --version)"
echo "  Conda env: $CONDA_DEFAULT_ENV"
echo ""

# === PREREQUISITE CHECKS ===
echo "📋 PREREQUISITE CHECKS"
echo "--------------------"

# Check clean dataset
if [ ! -f "${CLEAN_DATASET}" ]; then
    echo "❌ Clean dataset not found: ${CLEAN_DATASET}"
    echo "💡 Create clean dataset first:"
    echo "   python3 -c \"import pandas as pd; df = pd.read_csv('/u/aa107/uiuc-cancer-research/data/processed/tabnet_csv/prostate_variants_tabnet_enhanced.csv', low_memory=False); df_clean = df.drop(columns=['sift_prediction', 'polyphen_prediction'] if 'sift_prediction' in df.columns else []); df_clean.to_csv('${CLEAN_DATASET}', index=False); print('Clean dataset created')\""
    exit 1
else
    DATASET_SIZE=$(du -h "${CLEAN_DATASET}" | cut -f1)
    DATASET_LINES=$(wc -l < "${CLEAN_DATASET}")
    echo "✅ Clean dataset found: ${CLEAN_DATASET}"
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

if [ $? -ne 0 ]; then
    echo "❌ Environment tests failed"
    exit 1
fi

echo "✅ Environment tests passed"

# === STEP 2: ENHANCED DATASET VALIDATION ===
echo ""
echo "🔍 STEP 2: ENHANCED DATASET VALIDATION"
echo "======================================"

echo "Validating clean dataset structure..."
python -c "
import pandas as pd

# Load clean dataset
df = pd.read_csv('${CLEAN_DATASET}', low_memory=False)
print(f'Dataset loaded: {df.shape[0]:,} variants × {df.shape[1]} features')

# Check for data leakage
leakage_features = ['sift_prediction', 'polyphen_prediction', 'functional_pathogenicity']
leakage_found = [f for f in leakage_features if f in df.columns]

if leakage_found:
    print(f'❌ Data leakage features found: {leakage_found}')
    exit(1)
else:
    print('✅ No data leakage features detected')

# Check AlphaMissense
am_features = ['alphamissense_pathogenicity', 'alphamissense_class']
missing_am = [f for f in am_features if f not in df.columns]

if missing_am:
    print(f'❌ AlphaMissense features missing: {missing_am}')
else:
    print('✅ AlphaMissense features present')
    am_coverage = df['alphamissense_pathogenicity'].notna().sum()
    coverage_pct = am_coverage / len(df) * 100
    print(f'📊 AlphaMissense coverage: {am_coverage:,} variants ({coverage_pct:.1f}%)')
    
    if coverage_pct >= 30:
        print('✅ Good AlphaMissense coverage')
    else:
        print('⚠️  Low AlphaMissense coverage')

print('✅ Enhanced dataset validation passed')
"

if [ $? -ne 0 ]; then
    echo "❌ Enhanced dataset validation failed"
    exit 1
fi

echo "✅ Enhanced dataset validation passed"

# === STEP 3: COMPREHENSIVE TABNET VALIDATION ===
echo ""
echo "🔥 STEP 3: COMPREHENSIVE TABNET VALIDATION"
echo "=========================================="

echo "Running comprehensive TabNet validation..."

# Run the fixed validation script
python "${VALIDATION_SCRIPT}"

VALIDATION_EXIT_CODE=$?

echo ""
echo "=== VALIDATION RESULTS ==="

if [ $VALIDATION_EXIT_CODE -eq 0 ]; then
    echo "✅ VALIDATION COMPLETED SUCCESSFULLY!"
    echo ""
    
    # Check for latest report
    RESULTS_DIR="${PROJECT_DIR}/results/validation"
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
    if 'baseline_validation' in results:
        baseline_acc = results['baseline_validation']['mean_accuracy']
        print(f'  Baseline accuracy: {baseline_acc:.3f}')
        
        if baseline_acc > 0.95:
            print('  Baseline status: ❌ SUSPICIOUS - Check for leakage')
        elif baseline_acc > 0.75:
            print('  Baseline status: ✅ EXCELLENT')
        else:
            print('  Baseline status: ✅ REALISTIC')
    
    if 'tabnet_validation' in results:
        tabnet_acc = results['tabnet_validation']['mean_accuracy']
        print(f'  TabNet accuracy: {tabnet_acc:.3f}')
        
        if tabnet_acc > 0.95:
            print('  TabNet status: ❌ SUSPICIOUS - Check for leakage')
        elif tabnet_acc > 0.75:
            print('  TabNet status: ✅ EXCELLENT')
        else:
            print('  TabNet status: ✅ REALISTIC')
    
except Exception as e:
    print(f'  ⚠️  Could not parse report: {e}')
"
    else
        echo "⚠️  No validation report found in ${RESULTS_DIR}"
    fi
    
    echo ""
    echo "🎯 NEXT STEPS:"
    echo "  1. Review validation report for detailed results"
    echo "  2. If no data leakage detected, proceed with full training:"
    echo "     python src/model/tabnet_prostate_variant_classifier.py"
    echo "  3. Or submit full training cluster job"
    
    exit 0
    
else
    echo "❌ VALIDATION FAILED (exit code: $VALIDATION_EXIT_CODE)"
    echo ""
    echo "🔧 TROUBLESHOOTING:"
    echo "  1. Check the validation output above for specific errors"
    echo "  2. Ensure clean dataset is properly formatted"
    echo "  3. Verify all dependencies are installed correctly"
    echo "  4. Check that custom TabNet model can be imported"
    
    exit 1
fi