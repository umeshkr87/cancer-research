#!/bin/bash
set -e

echo "🧬 Advanced ClinVar Prostate Cancer Variant Analyzer"
echo "=================================================="

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Paths
VCF_PATH_GRCH38="$PROJECT_ROOT/data/raw/variants/vcf_GRCh38/clinvar_vcf_GRCh38.vcf"
VCF_PATH_GRCH38_GZ="$PROJECT_ROOT/data/raw/variants/vcf_GRCh38/clinvar_vcf_GRCh38.vcf.gz"
OUTPUT_DIR="$PROJECT_ROOT/data/processed/clinvar_prostate"
SCRIPT_PATH="$SCRIPT_DIR/filter_clinvar_prostate.py"

echo "Project root: $PROJECT_ROOT"
echo "Script location: $SCRIPT_PATH"
echo "Output directory: $OUTPUT_DIR"

# Check for VCF file (try both compressed and uncompressed)
VCF_FILE=""
if [ -f "$VCF_PATH_GRCH38" ]; then
    VCF_FILE="$VCF_PATH_GRCH38"
    echo "Found uncompressed VCF: $VCF_FILE"
elif [ -f "$VCF_PATH_GRCH38_GZ" ]; then
    VCF_FILE="$VCF_PATH_GRCH38_GZ"
    echo "Found compressed VCF: $VCF_FILE"
else
    echo "❌ ClinVar VCF file not found!"
    echo ""
    echo "Expected locations:"
    echo "  $VCF_PATH_GRCH38"
    echo "  $VCF_PATH_GRCH38_GZ"
    echo ""
    echo "Please upload your ClinVar VCF file:"
    echo "  scp vcf_GRCh38/clinvar.vcf* aa107@cli-dtn.researchdata.illinois.edu:$PROJECT_ROOT/data/raw/variants/vcf_GRCh38/"
    exit 1
fi

# Check if script exists
if [ ! -f "$SCRIPT_PATH" ]; then
    echo "❌ Analysis script not found: $SCRIPT_PATH"
    echo "Please ensure analyze_and_prepare_clinvar.py is in the scripts/variants/ directory"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Load environment
echo "Loading environment..."
module load anaconda3 2>/dev/null || echo "Note: anaconda3 module not available, using system Python"

# Activate conda environment if it exists
if command -v conda >/dev/null 2>&1; then
    echo "Activating tabnet-prostate environment..."
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate tabnet-prostate 2>/dev/null || echo "Note: tabnet-prostate environment not found, using base"
fi

# Install required packages if needed
echo "Checking Python dependencies..."
python -c "import pandas, numpy, tqdm" 2>/dev/null || {
    echo "Installing required packages..."
    pip install pandas numpy tqdm --quiet
}

# Determine if this is a test run or full run
TEST_MODE=false
if [ "$1" = "--test" ]; then
    TEST_MODE=true
    echo "🧪 Running in TEST MODE (limited variants)"
fi

# Run the analysis
echo "Starting ClinVar analysis..."
echo "Input VCF: $VCF_FILE"
echo "Output directory: $OUTPUT_DIR"

if [ "$TEST_MODE" = true ]; then
    # Test run with limited variants
    python "$SCRIPT_PATH" \
        --vcf-path "$VCF_FILE" \
        --output-dir "$OUTPUT_DIR" \
        --max-variants 10000 \
        --log-file "$OUTPUT_DIR/clinvar_analysis_test.log"
else
    # Full run
    python "$SCRIPT_PATH" \
        --vcf-path "$VCF_FILE" \
        --output-dir "$OUTPUT_DIR" \
        --log-file "$OUTPUT_DIR/clinvar_analysis.log"
fi

# Check results
echo ""
echo "✅ Analysis completed!"
echo ""

# Display results summary
if [ -f "$OUTPUT_DIR/clinvar_prostate.csv" ]; then
    VARIANT_COUNT=$(tail -n +2 "$OUTPUT_DIR/clinvar_prostate.csv" | wc -l)
    echo "📊 Results Summary:"
    echo "  Main dataset: $OUTPUT_DIR/clinvar_prostate.csv"
    echo "  Variant count: $VARIANT_COUNT"
    
    # Quick quality checks
    echo ""
    echo "🔍 Quick Quality Checks:"
    
    # Check pathogenic distribution
    echo "  Pathogenic distribution:"
    tail -n +2 "$OUTPUT_DIR/clinvar_prostate.csv" | cut -d',' -f7 | sort | uniq -c | while read count value; do
        if [ "$value" = "1" ]; then
            echo "    Pathogenic: $count"
        elif [ "$value" = "0" ]; then
            echo "    Non-pathogenic: $count"
        fi
    done
    
    # Check for key genes
    echo "  Key prostate genes found:"
    KEY_GENES=("BRCA1" "BRCA2" "ATM" "AR" "PTEN" "TP53")
    for gene in "${KEY_GENES[@]}"; do
        count=$(grep -c "$gene" "$OUTPUT_DIR/clinvar_prostate.csv" || echo "0")
        echo "    $gene: $count variants"
    done
    
    echo ""
    echo "📋 Next Steps:"
    echo "  1. Review comprehensive results: $OUTPUT_DIR/clinvar_prostate_comprehensive.csv"
    echo "  2. Check analysis summary: $OUTPUT_DIR/clinvar_prostate_summary.txt"
    echo "  3. Integrate with TCGA-PRAD and COSMIC datasets"
    echo "  4. Proceed with TabNet feature engineering"
    
else
    echo "❌ Main output file not found. Check log for errors:"
    echo "   $OUTPUT_DIR/clinvar_analysis.log"
fi

echo ""
echo "🎯 Ready for TabNet Integration!"
echo "Main dataset: $OUTPUT_DIR/clinvar_prostate.csv"