#!/bin/bash
set -e

echo "🚀 TabNet H100 Optimization Runner"
echo "=================================="

cd /u/aa107/uiuc-cancer-research

# Setup environment
echo "🔧 Setting up environment..."
module load anaconda3
eval "$(conda shell.bash hook)"
conda activate tabnet-prostate

echo "✅ Environment: $CONDA_DEFAULT_ENV"

# Check GPU availability
echo "🔍 Checking GPU access..."
if python -c "import torch; print('✅ CUDA available' if torch.cuda.is_available() else '❌ No GPU')"; then
    echo "🚀 Starting 36-hour optimization..."
    python scripts/optimization/hyperparameter_optimization.py
else
    echo "⚠️  No GPU detected - submitting SLURM job instead..."
    echo "📋 Use: sbatch scripts/optimization/optimize_tabnet.sbatch"
fi

echo "✅ Process completed!"