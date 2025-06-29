#!/bin/bash
# VEP VCF Analysis Runner
# Location: /u/aa107/uiuc-cancer-research/scripts/vep/validation/run_vep_analysis.sh
# Usage: bash run_vep_analysis.sh

set -e

echo "🧬 VEP VCF Analysis for TabNet Training"
echo "======================================"

# Setup environment
cd /u/aa107/uiuc-cancer-research

# Check if Python is available
if command -v python3 >/dev/null 2>&1; then
    PYTHON=python3
elif command -v python >/dev/null 2>&1; then
    PYTHON=python
else
    echo "❌ Python not found! Please load Python module or install Python"
    exit 1
fi

echo "🐍 Using Python: $(which $PYTHON)"
echo "📁 Project directory: $(pwd)"
echo

# Check if VCF file exists
VCF_FILE="/u/aa107/uiuc-cancer-research/data/processed/vep/vep_annotated.vcf"
if [ ! -f "$VCF_FILE" ]; then
    echo "❌ VCF file not found: $VCF_FILE"
    echo "Please run VEP annotation first."
    exit 1
fi

echo "✅ VCF file found: $VCF_FILE"
echo "📏 File size: $(du -sh $VCF_FILE | cut -f1)"
echo

# Run analysis
echo "🚀 Starting VEP analysis..."
echo "⏱️  This may take a few minutes for large files..."
echo

$PYTHON scripts/vep/validation/analyze_vep_vcf.py

echo
echo "✅ Analysis complete!"
echo "📄 Report saved to: /u/aa107/uiuc-cancer-research/data/processed/vep/vep_analysis_report.txt"
echo
echo "🔄 Next Steps:"
echo "1. Review functional scores availability in the report"
echo "2. Proceed with VCF to CSV conversion if ready"
echo "3. Start TabNet feature engineering"