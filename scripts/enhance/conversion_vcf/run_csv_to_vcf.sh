#!/bin/bash
set -e

echo "🧬 Simple CSV to VCF Converter"
echo "=============================="

# Setup environment
module load anaconda3 2>/dev/null || echo "Using system Python"

if command -v conda >/dev/null 2>&1; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate tabnet-prostate 2>/dev/null || echo "Using base environment"
fi

# Check dependencies
python3 -c "import pandas" 2>/dev/null || {
    echo "Installing pandas..."
    pip install pandas --quiet
}

# Run conversion
echo "Running conversion..."
cd /u/aa107/uiuc-cancer-research

python3 scripts/enhance/conversion_vcf/csv_to_vcf.py

echo "✅ Conversion complete!"
echo "📁 Output: data/processed/merged_vcf/merged_prostate_variants.vcf"
echo "🔄 Ready for VarStack annotation"