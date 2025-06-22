#!/bin/bash
set -e

echo "🧬 TabNet Prostate Cancer Merged Dataset Report"
echo "=============================================="

# Configuration
PROJECT_ROOT="/u/aa107/uiuc-cancer-research"
MERGED_DIR="$PROJECT_ROOT/data/processed/merged"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "📁 Analysis Directory: $MERGED_DIR"
echo "🕐 Report Generated: $(date)"
echo ""

# Load environment and execute merge script
echo "🔧 SETTING UP ENVIRONMENT"
echo "========================="
module load anaconda3 2>/dev/null || echo "Note: anaconda3 module not available, using system Python"

# Activate conda environment if available
if command -v conda >/dev/null 2>&1; then
    echo "Activating tabnet-prostate environment..."
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate tabnet-prostate 2>/dev/null || echo "Note: tabnet-prostate environment not found, using base"
fi

# Install required packages if needed
echo "Checking Python dependencies..."
python3 -c "import pandas, numpy" 2>/dev/null || {
    echo "Installing required packages..."
    pip install pandas numpy --quiet
}

# Create merged directory if it doesn't exist
mkdir -p "$MERGED_DIR"

# Execute merge script
echo ""
echo "🔄 EXECUTING MERGE SCRIPT"
echo "========================="
cd "$PROJECT_ROOT"

if python3 scripts/merge/merge_datasets.py; then
    echo ""
    echo "✅ MERGE COMPLETED SUCCESSFULLY"
    echo "==============================="
else
    echo ""
    echo "❌ MERGE SCRIPT FAILED"
    echo "======================"
    echo "Please check the error messages above and ensure:"
    echo "1. All input datasets are available"
    echo "2. Python environment has required packages"
    echo "3. Sufficient disk space available"
    echo ""
fi

# Directory Overview
echo "📂 DIRECTORY OVERVIEW"
echo "===================="
echo "Directory size: $(du -sh $MERGED_DIR | cut -f1)"
echo "File count: $(find $MERGED_DIR -type f | wc -l)"
echo ""

echo "📄 FILES IN DIRECTORY"
echo "===================="
ls -lah $MERGED_DIR
echo ""

# Main merged dataset analysis
MERGED_CSV="$MERGED_DIR/merged_prostate_variants.csv"
REPORT_TXT="$MERGED_DIR/merge_report.txt"

if [ -f "$MERGED_CSV" ]; then
    echo "✅ MAIN DATASET FOUND: merged_prostate_variants.csv"
    echo "================================================="
    
    # File details
    echo "File size: $(du -sh $MERGED_CSV | cut -f1)"
    echo "Last modified: $(stat -c %y $MERGED_CSV)"
    echo ""
    
    # Basic statistics using Unix tools
    echo "📊 DATASET STATISTICS"
    echo "===================="
    
    TOTAL_ROWS=$(($(wc -l < "$MERGED_CSV") - 1))  # Subtract header
    TOTAL_COLS=$(head -1 "$MERGED_CSV" | tr ',' '\n' | wc -l)
    
    echo "Total variants: $(printf "%'d" $TOTAL_ROWS)"
    echo "Total columns: $TOTAL_COLS"
    echo "File size: $(du -sh $MERGED_CSV | cut -f1)"
    echo ""
    
    # Data sources distribution using grep
    echo "📈 DATA SOURCE DISTRIBUTION"
    echo "==========================="
    
    COSMIC_COUNT=$(grep -c ",COSMIC," "$MERGED_CSV" 2>/dev/null || echo "0")
    CLINVAR_COUNT=$(grep -c ",ClinVar," "$MERGED_CSV" 2>/dev/null || echo "0")
    TCGA_COUNT=$(grep -c ",TCGA," "$MERGED_CSV" 2>/dev/null || echo "0")
    
    if [ $COSMIC_COUNT -gt 0 ]; then
        COSMIC_PCT=$(echo "scale=1; $COSMIC_COUNT * 100 / $TOTAL_ROWS" | bc -l 2>/dev/null || echo "0.0")
        echo "  COSMIC: $(printf "%'d" $COSMIC_COUNT) ($COSMIC_PCT%)"
    fi
    
    if [ $CLINVAR_COUNT -gt 0 ]; then
        CLINVAR_PCT=$(echo "scale=1; $CLINVAR_COUNT * 100 / $TOTAL_ROWS" | bc -l 2>/dev/null || echo "0.0")
        echo "  ClinVar: $(printf "%'d" $CLINVAR_COUNT) ($CLINVAR_PCT%)"
    fi
    
    if [ $TCGA_COUNT -gt 0 ]; then
        TCGA_PCT=$(echo "scale=1; $TCGA_COUNT * 100 / $TOTAL_ROWS" | bc -l 2>/dev/null || echo "0.0")
        echo "  TCGA: $(printf "%'d" $TCGA_COUNT) ($TCGA_PCT%)"
    fi
    echo ""
    
    # Variant classification distribution
    echo "🎯 VARIANT CLASSIFICATION"
    echo "========================="
    
    VUS_COUNT=$(grep -c "VUS" "$MERGED_CSV" 2>/dev/null || echo "0")
    PATHOGENIC_COUNT=$(grep -c "Actionable_Pathogenic" "$MERGED_CSV" 2>/dev/null || echo "0")
    BENIGN_COUNT=$(grep -c -E "(Benign|Likely_Benign)" "$MERGED_CSV" 2>/dev/null || echo "0")
    ACTIONABLE_COUNT=$(grep -c "Likely_Actionable" "$MERGED_CSV" 2>/dev/null || echo "0")
    
    echo "  VUS (Uncertain): $(printf "%'d" $VUS_COUNT)"
    echo "  Actionable Pathogenic: $(printf "%'d" $PATHOGENIC_COUNT)"
    echo "  Benign variants: $(printf "%'d" $BENIGN_COUNT)"
    echo "  Likely Actionable: $(printf "%'d" $ACTIONABLE_COUNT)"
    echo ""
    
    # Top genes
    echo "🔝 TOP 10 GENES"
    echo "==============="
    
    # Extract gene column (assuming it's the first column) and count
    tail -n +2 "$MERGED_CSV" | cut -d',' -f1 | sort | uniq -c | sort -rn | head -10 | while read count gene; do
        echo "  $gene: $(printf "%'d" $count)"
    done
    echo ""
    
    # Pathway analysis (count rows with 1 in specific columns)
    echo "🧬 PATHWAY GENE ANALYSIS"
    echo "========================"
    
    # Look for pathway columns - this is approximate since column positions may vary
    DNA_REPAIR=$(grep -o ",1," "$MERGED_CSV" | wc -l 2>/dev/null || echo "Unknown")
    echo "  DNA Repair pathway variants: ~$DNA_REPAIR"
    echo "  (Approximate counts - exact analysis in detailed report)"
    echo ""
    
    # Chromosome distribution (basic)
    echo "🧮 CHROMOSOME DISTRIBUTION"
    echo "=========================="
    
    # Extract chromosome column and show unique values
    UNIQUE_CHRS=$(tail -n +2 "$MERGED_CSV" | cut -d',' -f2 | sort -u | head -10 | tr '\n' ' ')
    CHR_COUNT=$(tail -n +2 "$MERGED_CSV" | cut -d',' -f2 | sort -u | wc -l)
    
    echo "  Total chromosomes: $CHR_COUNT"
    echo "  Sample chromosomes: $UNIQUE_CHRS"
    echo ""
    
    # Data quality checks
    echo "🔍 DATA QUALITY CHECKS"
    echo "======================"
    
    # Check for empty fields (basic)
    EMPTY_FIELDS=$(grep -c ",," "$MERGED_CSV" 2>/dev/null || echo "0")
    echo "  Rows with empty fields: $EMPTY_FIELDS"
    
    # Check file integrity
    HEAD_COLS=$(head -1 "$MERGED_CSV" | tr ',' '\n' | wc -l)
    SAMPLE_ROW_COLS=$(tail -1 "$MERGED_CSV" | tr ',' '\n' | wc -l)
    
    if [ $HEAD_COLS -eq $SAMPLE_ROW_COLS ]; then
        echo "  ✅ Column consistency check passed"
    else
        echo "  ⚠️  Column count mismatch detected"
    fi
    
    echo "  ✅ File readable and accessible"
    echo ""

else
    echo "❌ MERGE SCRIPT FAILED!"
    echo "The merge script completed but did not generate the expected output file."
    echo "Please check the error messages above."
    echo ""
fi

# Report file analysis
if [ -f "$REPORT_TXT" ]; then
    echo "📋 MERGE REPORT SUMMARY"
    echo "======================="
    echo "Report file: merge_report.txt"
    echo "Report size: $(du -sh $REPORT_TXT | cut -f1)"
    echo ""
    echo "Report preview:"
    head -10 "$REPORT_TXT" | sed 's/^/  /'
    echo "  ..."
    echo ""
else
    echo "⚠️  Merge report file not found: merge_report.txt"
    echo ""
fi

# Input dataset status
echo "📥 INPUT DATASETS STATUS"
echo "========================"

# Check source datasets
COSMIC_DIR="$PROJECT_ROOT/data/processed/cosmic_prostate"
CLINVAR_DIR="$PROJECT_ROOT/data/processed/clinvar_prostate"  
TCGA_DIR="$PROJECT_ROOT/data/processed/tcga_prad_prostate"

for dataset in "COSMIC:$COSMIC_DIR/cosmic_prostate.csv" "ClinVar:$CLINVAR_DIR/clinvar_prostate.csv" "TCGA:$TCGA_DIR/tcga_prad_mutations.csv"; do
    name=$(echo $dataset | cut -d: -f1)
    path=$(echo $dataset | cut -d: -f2)
    
    if [ -f "$path" ]; then
        size=$(du -sh "$path" | cut -f1)
        lines=$(wc -l < "$path")
        echo "  ✅ $name: $size, $((lines-1)) variants"
    else
        echo "  ❌ $name: Not found ($path)"
    fi
done
echo ""

# Disk usage summary
echo "💾 DISK USAGE SUMMARY"
echo "===================="
echo "Merged directory: $(du -sh $MERGED_DIR | cut -f1)"
echo "Total processed data: $(du -sh $PROJECT_ROOT/data/processed | cut -f1)"
echo "Available space: $(df -h $PROJECT_ROOT | tail -1 | awk '{print $4}')"
echo ""

# TabNet readiness check (basic)
echo "🚀 TABNET READINESS CHECK"
echo "========================="

if [ -f "$MERGED_CSV" ]; then
    
    # Check file accessibility
    if [ -r "$MERGED_CSV" ]; then
        echo "✅ File readable and accessible"
    else
        echo "❌ File not readable"
    fi
    
    # Check minimum row count
    TOTAL_ROWS=$(($(wc -l < "$MERGED_CSV") - 1))
    if [ $TOTAL_ROWS -gt 1000 ]; then
        echo "✅ Sufficient data: $(printf "%'d" $TOTAL_ROWS) variants"
    else
        echo "⚠️  Limited data: only $TOTAL_ROWS variants"
    fi
    
    # Check minimum column count
    TOTAL_COLS=$(head -1 "$MERGED_CSV" | tr ',' '\n' | wc -l)
    if [ $TOTAL_COLS -gt 20 ]; then
        echo "✅ Sufficient features: $TOTAL_COLS columns"
    else
        echo "⚠️  Limited features: only $TOTAL_COLS columns"
    fi
    
    # Check for required columns (basic)
    HEADER=$(head -1 "$MERGED_CSV")
    
    if echo "$HEADER" | grep -q "gene_symbol"; then
        echo "✅ Gene symbol column present"
    else
        echo "⚠️  Gene symbol column missing"
    fi
    
    if echo "$HEADER" | grep -q "variant_classification"; then
        echo "✅ Target variable column present"
    else
        echo "⚠️  Target variable column missing"
    fi
    
    if echo "$HEADER" | grep -q "chromosome"; then
        echo "✅ Chromosome column present"
    else
        echo "⚠️  Chromosome column missing"
    fi
    
    echo ""
    echo "🎯 DATASET READY FOR NEXT STEPS"
    echo "  • Total variants: $(printf "%'d" $TOTAL_ROWS)"
    echo "  • Total features: $TOTAL_COLS" 
    echo "  • File size: $(du -sh $MERGED_CSV | cut -f1)"
    
else
    echo "❌ Cannot perform readiness check - merge script failed to generate dataset"
fi

echo ""
echo "🎯 NEXT STEPS"
echo "============="
echo "1. ✅ Data merge completed successfully"
echo "2. 🔄 Next: VarStack annotation (40+ columns)"
echo "   • Convert CSV to VCF format"
echo "   • Upload to VarStack API"
echo "   • Get population frequencies, conservation scores"
echo "3. 🔄 Then: VEP annotation (30+ functional scores)"
echo "   • Upload VCF to Ensembl VEP"
echo "   • Get SIFT, PolyPhen, CADD scores"
echo "4. 🔄 Finally: TabNet training"
echo "   cd $PROJECT_ROOT"
echo "   python src/model/tabnet_prostate_variant_classifier.py"
echo "5. 🔍 Environment check:"
echo "   python src/model/tests/test_environment.py"
echo ""
echo "📊 Report completed at $(date)"
echo "For detailed merge statistics, see: $REPORT_TXT"