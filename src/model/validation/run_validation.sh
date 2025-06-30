#!/bin/bash
set -e

echo "🧬 TabNet Validation Runner"
echo "=========================="

cd /u/aa107/uiuc-cancer-research

# Setup environment
echo "🔧 Setting up environment..."
module load anaconda3
eval "$(conda shell.bash hook)"
conda activate tabnet-prostate

echo "✅ Environment: $CONDA_DEFAULT_ENV"

# Run validation
echo "🔄 Running TabNet validation..."
python src/model/validation/validate_tabnet.py

echo "✅ Validation completed!"