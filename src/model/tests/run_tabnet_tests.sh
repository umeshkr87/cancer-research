#!/bin/bash
set -e

echo "🧬 TabNet Environment Validation - Campus Cluster"
echo "==============================================="
echo "Starting comprehensive TabNet environment validation..."

# Configuration
PROJECT_ROOT="/u/aa107/uiuc-cancer-research"
cd "$PROJECT_ROOT"

echo "📁 Project directory: $(pwd)"
echo "🕐 Started at: $(date)"

# Step 1: Environment Setup
echo ""
echo "🔧 STEP 1: ENVIRONMENT SETUP"
echo "============================"

# Load anaconda module
echo "Loading anaconda module..."
module load anaconda3

# Initialize conda for this session (key fix)
echo "Initializing conda for bash..."
eval "$(conda shell.bash hook)"

# Check if environment exists
echo "Checking conda environments..."
if conda env list | grep -q "tabnet-prostate"; then
    echo "✅ tabnet-prostate environment found"
else
    echo "❌ tabnet-prostate environment not found"
    echo "🔧 Creating environment..."
    conda create -n tabnet-prostate python=3.11 -y
    echo "✅ Environment created"
fi

# Activate environment
echo "Activating tabnet-prostate environment..."
conda activate tabnet-prostate

# Verify activation
echo "✅ Active environment: $CONDA_DEFAULT_ENV"
echo "✅ Python location: $(which python)"
echo "✅ Python version: $(python --version)"

# Install required packages if missing
echo "Checking/installing required packages..."
python -c "
import subprocess
import sys

packages = ['torch', 'pandas', 'numpy', 'scikit-learn', 'pytorch-tabnet', 'matplotlib', 'seaborn']
missing = []

for package in packages:
    try:
        __import__(package.replace('-', '_'))
        print(f'✅ {package}')
    except ImportError:
        missing.append(package)
        print(f'❌ {package} (missing)')

if missing:
    print(f'📦 Installing missing packages: {missing}')
    for package in missing:
        try:
            if package == 'torch':
                # Install PyTorch with CUDA support
                subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'torch', 'torchvision', 'torchaudio', '--index-url', 'https://download.pytorch.org/whl/cu118'], timeout=300)
            else:
                subprocess.check_call([sys.executable, '-m', 'pip', 'install', package], timeout=120)
            print(f'✅ {package} installed')
        except Exception as e:
            print(f'❌ Failed to install {package}: {e}')
else:
    print('✅ All packages available')
"

echo ""
echo "🧪 STEP 2: LOCAL ENVIRONMENT TESTS"
echo "=================================="
echo "Running local tests with activated environment..."

# Run local tests in activated environment
python src/model/tests/test_environment.py

LOCAL_EXIT_CODE=$?

echo ""
echo "🚀 STEP 3: GPU ENVIRONMENT TESTS"
echo "==============================="

if [ $LOCAL_EXIT_CODE -eq 0 ]; then
    echo "✅ Local tests passed - submitting GPU tests..."
else
    echo "⚠️  Local tests had issues - submitting GPU tests anyway..."
fi

# Submit GPU tests
cd src/model/tests
echo "Submitting GPU test job..."

JOB_ID=$(sbatch test_gpu_environment.sbatch | awk '{print $4}')
echo "✅ GPU test job submitted: $JOB_ID"

# Wait for completion with timeout
echo "Waiting for GPU test completion (max 10 minutes)..."
TIMEOUT=600  # 10 minutes
ELAPSED=0
INTERVAL=15

while [ $ELAPSED -lt $TIMEOUT ]; do
    if [ $(squeue -j $JOB_ID -h 2>/dev/null | wc -l) -eq 0 ]; then
        echo "✅ GPU test completed"
        break
    fi
    
    echo "  ⏳ Still running... (${ELAPSED}s/${TIMEOUT}s)"
    sleep $INTERVAL
    ELAPSED=$((ELAPSED + INTERVAL))
done

if [ $ELAPSED -ge $TIMEOUT ]; then
    echo "⏰ GPU test timeout - check manually with: squeue -j $JOB_ID"
    echo "📄 Check results later with: cat env_test_${JOB_ID}.out"
    exit 1
fi

echo ""
echo "📊 STEP 4: RESULTS SUMMARY"
echo "========================="

# Show GPU test results
if [ -f "env_test_${JOB_ID}.out" ]; then
    echo "GPU Test Results:"
    echo "----------------"
    cat "env_test_${JOB_ID}.out"
    
    # Clean up
    rm -f "env_test_${JOB_ID}.out" "env_test_${JOB_ID}.err" 2>/dev/null
    echo ""
    echo "🧹 Cleaned up output files"
else
    echo "❌ GPU test output file not found"
fi

echo ""
echo "🎉 VALIDATION COMPLETED"
echo "======================"
echo "Completed at: $(date)"

# Final assessment
if [ $LOCAL_EXIT_CODE -eq 0 ]; then
    echo "✅ Local tests: PASSED"
else
    echo "⚠️  Local tests: ISSUES DETECTED"
fi

echo ""
echo "🚀 NEXT STEPS FOR WEEK 5:"
echo "========================"
echo "1. Run validation: ./src/model/tests/run_tabnet_tests.sh"
echo "2. Quick validation: python src/model/validate_tabnet.py"  
echo "3. Start optimization: python scripts/optimization/hyperparameter_optimization.py"
echo ""
echo "📚 Documentation: src/model/tests/README.md"