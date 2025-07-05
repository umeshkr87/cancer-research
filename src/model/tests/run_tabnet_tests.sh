#!/bin/bash
# ================================================================
# TabNet Test Runner - Enhanced Version
# Properly sets up conda environment and runs tests
# ================================================================

set -euo pipefail

echo "🧪 TABNET TEST RUNNER - ENHANCED VERSION"
echo "=================================================="
echo "Setting up environment and running tests..."
echo ""

# === CONFIGURATION ===
PROJECT_DIR="/u/aa107/uiuc-cancer-research"
TEST_SCRIPT="${PROJECT_DIR}/src/model/tests/test_environment.py"

# === ENVIRONMENT SETUP ===
echo "🔧 ENVIRONMENT SETUP"
echo "-------------------"

# Navigate to project directory
cd "${PROJECT_DIR}"
echo "Working directory: $(pwd)"

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
    echo "💡 Create it first with:"
    echo "   module load anaconda3"
    echo "   conda create -n tabnet-prostate python=3.11 -y"
    echo "   conda activate tabnet-prostate"
    echo "   pip install torch torchvision pytorch-tabnet scikit-learn pandas numpy matplotlib seaborn"
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
echo "📦 CHECKING DEPENDENCIES"
echo "-----------------------"

# Check core dependencies
echo "Checking PyTorch..."
python -c "import torch; print(f'  ✅ PyTorch: {torch.__version__}')" || {
    echo "  ❌ PyTorch missing - installing..."
    conda install pytorch torchvision torchaudio pytorch-cuda=11.8 -c pytorch -c nvidia -y
}

echo "Checking TabNet..."
python -c "import pytorch_tabnet; print('  ✅ TabNet: OK')" || {
    echo "  ❌ TabNet missing - installing..."
    pip install pytorch-tabnet
}

echo "Checking sklearn..."
python -c "import sklearn; print(f'  ✅ Sklearn: {sklearn.__version__}')" || {
    echo "  ❌ Sklearn missing - installing..."
    pip install scikit-learn
}

echo "Checking pandas/numpy..."
python -c "import pandas as pd, numpy as np; print(f'  ✅ Pandas: {pd.__version__}, NumPy: {np.__version__}')" || {
    echo "  ❌ Pandas/NumPy missing - installing..."
    pip install pandas numpy
}

echo "✅ All dependencies available"
echo ""

# === RUN TESTS ===
echo "🧪 RUNNING ENHANCED TABNET TESTS"
echo "==============================="

if [ -f "$TEST_SCRIPT" ]; then
    echo "Running: $TEST_SCRIPT"
    echo ""
    
    # Run the test script
    python "$TEST_SCRIPT"
    
    TEST_EXIT_CODE=$?
    
    echo ""
    echo "=== TEST RESULTS ==="
    
    if [ $TEST_EXIT_CODE -eq 0 ]; then
        echo "✅ ALL TESTS PASSED!"
        echo "🎯 Environment is ready for TabNet training"
        echo ""
        echo "Next steps:"
        echo "  1. Run validation: python src/model/validation/validate_tabnet.py"
        echo "  2. Or submit cluster job: sbatch src/model/validation/validate_tabnet.sbatch"
        exit 0
    else
        echo "❌ SOME TESTS FAILED (exit code: $TEST_EXIT_CODE)"
        echo "🔧 Check the test output above for specific issues"
        echo ""
        echo "Common solutions:"
        echo "  - Ensure enhanced dataset exists"
        echo "  - Check conda environment has all packages"
        echo "  - Verify file paths are correct"
        exit 1
    fi
    
else
    echo "❌ Test script not found: $TEST_SCRIPT"
    echo "💡 Make sure you're in the project root directory"
    exit 1
fi