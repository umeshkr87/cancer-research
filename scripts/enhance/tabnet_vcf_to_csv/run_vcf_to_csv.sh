#!/bin/bash
# VCF to TabNet CSV Converter Runner
# Location: /u/aa107/uiuc-cancer-research/scripts/enhance/run_vcf_to_csv.sh
# Usage: bash run_vcf_to_csv.sh [test|full]

set -e

echo "🧬 VEP VCF to TabNet CSV Converter"
echo "=================================="

# Setup environment
cd /u/aa107/uiuc-cancer-research

# Check if Python is available
if command -v python3 >/dev/null 2>&1; then
    PYTHON=python3
elif command -v python >/dev/null 2>&1; then
    PYTHON=python
else
    echo "❌ Python not found! Please load Python module"
    exit 1
fi

echo "🐍 Using Python: $(which $PYTHON)"

# Check for required Python packages
echo "📦 Checking Python dependencies..."
$PYTHON -c "import pandas; print('✅ pandas available')" || {
    echo "❌ pandas not found. Installing..."
    pip install pandas --user
}

$PYTHON -c "import numpy; print('✅ numpy available')" || {
    echo "❌ numpy not found. Installing..."
    pip install numpy --user
}

# Input and output paths
INPUT_VCF="/u/aa107/uiuc-cancer-research/data/processed/vep/vep_annotated.vcf"
OUTPUT_CSV="/u/aa107/uiuc-cancer-research/data/processed/tabnet_csv/prostate_variants_tabnet.csv"

# Check if input exists
if [ ! -f "$INPUT_VCF" ]; then
    echo "❌ Input VCF not found: $INPUT_VCF"
    echo "Please run VEP annotation first"
    exit 1
fi

echo "✅ Input VCF found: $INPUT_VCF"
echo "📏 Input size: $(du -sh $INPUT_VCF | cut -f1)"

# Create output directory
mkdir -p "$(dirname $OUTPUT_CSV)"
echo "📁 Output directory: $(dirname $OUTPUT_CSV)"
echo

# Determine processing mode
MODE=${1:-full}

case $MODE in
    "test")
        echo "🧪 Running in TEST mode (first 10,000 variants)"
        echo "⚡ Faster processing for validation"
        $PYTHON scripts/enhance/tabnet_vcf_to_csv/vcf_to_tabnet_csv.py \
            --input "$INPUT_VCF" \
            --output "${OUTPUT_CSV%.csv}_test.csv" \
            --max-variants 10000
        ;;
    "full")
        echo "🚀 Running in FULL mode (all variants)"
        echo "⏱️  This may take 10-15 minutes for 193K variants"
        $PYTHON scripts/enhance/tabnet_vcf_to_csv/vcf_to_tabnet_csv.py \
            --input "$INPUT_VCF" \
            --output "$OUTPUT_CSV"
        ;;
    *)
        echo "❌ Unknown mode: $MODE"
        echo "Usage: bash run_vcf_to_csv.sh [test|full]"
        echo "  test: Process first 10K variants for validation"
        echo "  full: Process all variants for training"
        exit 1
        ;;
esac

echo
echo "✅ Conversion completed!"

# Show output info
if [ -f "$OUTPUT_CSV" ]; then
    echo "📁 Output file: $OUTPUT_CSV"
    echo "📏 Output size: $(du -sh $OUTPUT_CSV | cut -f1)"
    echo "📊 Row count: $(tail -n +2 $OUTPUT_CSV | wc -l) variants"
elif [ -f "${OUTPUT_CSV%.csv}_test.csv" ]; then
    TEST_OUTPUT="${OUTPUT_CSV%.csv}_test.csv"
    echo "📁 Test output: $TEST_OUTPUT"
    echo "📏 Output size: $(du -sh $TEST_OUTPUT | cut -f1)"
    echo "📊 Row count: $(tail -n +2 $TEST_OUTPUT | wc -l) variants"
fi

echo
echo "🔄 Next Steps:"
echo "1. Review CSV file structure and features"
echo "2. Load into TabNet training pipeline"
echo "3. Configure feature selection for prostate cancer classification"
echo "4. Train interpretable model with attention mechanisms"

echo
echo "💡 Quick validation:"
echo "head -5 $OUTPUT_CSV"
echo "# This will show first 5 rows of your TabNet training data"