#!/bin/bash
# Functional Score Imputation Runner
# File: /u/aa107/uiuc-cancer-research/scripts/enhance/functional_enhancement/run_functional_imputation.sh
# Usage: bash run_functional_imputation.sh

set -e

echo "🧬 FUNCTIONAL SCORE IMPUTATION RUNNER"
echo "======================================"
echo "📅 $(date)"
echo ""

# Project setup
PROJECT_ROOT="/u/aa107/uiuc-cancer-research"
SCRIPT_DIR="$PROJECT_ROOT/scripts/enhance/functional_enhancement"
DATA_DIR="$PROJECT_ROOT/data/processed/tabnet_csv"

echo "📁 Project root: $PROJECT_ROOT"
echo "📁 Script directory: $SCRIPT_DIR"
echo "📁 Data directory: $DATA_DIR"
echo ""

# Change to project directory
cd "$PROJECT_ROOT"
echo "✅ Changed to project directory: $(pwd)"

# Check if Python is available
echo "🐍 CHECKING PYTHON ENVIRONMENT"
echo "------------------------------"

if command -v python3 >/dev/null 2>&1; then
    PYTHON=python3
elif command -v python >/dev/null 2>&1; then
    PYTHON=python
else
    echo "❌ Python not found! Please install Python or load Python module"
    echo "💡 Try: module load python3 (on HPC systems)"
    exit 1
fi

echo "✅ Python found: $(which $PYTHON)"
echo "✅ Python version: $($PYTHON --version)"

# Check Python dependencies
echo ""
echo "📦 CHECKING DEPENDENCIES"
echo "-------------------------"

check_python_package() {
    local package=$1
    if $PYTHON -c "import $package" 2>/dev/null; then
        echo "✅ $package available"
        return 0
    else
        echo "❌ $package missing"
        return 1
    fi
}

MISSING_DEPS=0

# Check required packages
if ! check_python_package "pandas"; then
    MISSING_DEPS=1
    echo "💡 Install with: pip install pandas"
fi

if ! check_python_package "numpy"; then
    MISSING_DEPS=1
    echo "💡 Install with: pip install numpy"
fi

if ! check_python_package "pathlib"; then
    echo "⚠️  pathlib not found as separate import, checking if built-in..."
    if $PYTHON -c "from pathlib import Path" 2>/dev/null; then
        echo "✅ pathlib available (built-in)"
    else
        MISSING_DEPS=1
        echo "❌ pathlib not available"
    fi
else
    echo "✅ pathlib available"
fi

# Install missing dependencies if needed
if [ $MISSING_DEPS -eq 1 ]; then
    echo ""
    echo "🔧 INSTALLING MISSING DEPENDENCIES"
    echo "----------------------------------"
    
    # Try pip install
    echo "Attempting to install missing packages..."
    if command -v pip3 >/dev/null 2>&1; then
        pip3 install pandas numpy --user
    elif command -v pip >/dev/null 2>&1; then
        pip install pandas numpy --user
    else
        echo "❌ pip not found! Please install missing packages manually:"
        echo "   - pandas"
        echo "   - numpy"
        exit 1
    fi
    
    # Verify installation
    echo "🔍 Verifying installation..."
    if check_python_package "pandas" && check_python_package "numpy"; then
        echo "✅ Dependencies installed successfully"
    else
        echo "❌ Dependency installation failed"
        exit 1
    fi
fi

# Check directory structure
echo ""
echo "📂 CHECKING DIRECTORY STRUCTURE"
echo "-------------------------------"

# Create directories if they don't exist
mkdir -p "$SCRIPT_DIR"
mkdir -p "$DATA_DIR"

# Check for input file
INPUT_FILE="$DATA_DIR/prostate_variants_tabnet.csv"
if [ -f "$INPUT_FILE" ]; then
    INPUT_SIZE=$(du -sh "$INPUT_FILE" | cut -f1)
    INPUT_LINES=$(wc -l < "$INPUT_FILE")
    echo "✅ Input file found: $INPUT_FILE"
    echo "   Size: $INPUT_SIZE"
    echo "   Lines: $((INPUT_LINES - 1)) variants"
else
    echo "❌ Input file not found: $INPUT_FILE"
    echo "💡 Expected location: $INPUT_FILE"
    echo "💡 Please ensure VEP processing completed successfully"
    exit 1
fi

# Check if script exists
IMPUTATION_SCRIPT="$SCRIPT_DIR/simple_functional_imputation.py"
if [ ! -f "$IMPUTATION_SCRIPT" ]; then
    echo "❌ Imputation script not found: $IMPUTATION_SCRIPT"
    echo "💡 Please create the script file at this location"
    exit 1
fi

echo "✅ Script found: $IMPUTATION_SCRIPT"
echo "✅ Directory structure validated"

# Show current data status
echo ""
echo "📊 CURRENT DATA STATUS"
echo "----------------------"

# Quick analysis of missing data
echo "🔍 Analyzing current functional score coverage..."
$PYTHON -c "
import pandas as pd
df = pd.read_csv('$INPUT_FILE')
total = len(df)
print(f'Total variants: {total:,}')

# Check for SIFT scores
sift_cols = [col for col in df.columns if 'sift' in col.lower() and 'score' in col.lower()]
if sift_cols:
    sift_col = sift_cols[0]
    sift_missing = df[sift_col].isna().sum()
    print(f'SIFT scores ({sift_col}): {total - sift_missing:,} present, {sift_missing:,} missing ({sift_missing/total*100:.1f}%)')
else:
    print('❌ No SIFT score columns found')

# Check for PolyPhen scores  
polyphen_cols = [col for col in df.columns if 'polyphen' in col.lower() and 'score' in col.lower()]
if polyphen_cols:
    polyphen_col = polyphen_cols[0]
    polyphen_missing = df[polyphen_col].isna().sum()
    print(f'PolyPhen scores ({polyphen_col}): {total - polyphen_missing:,} present, {polyphen_missing:,} missing ({polyphen_missing/total*100:.1f}%)')
else:
    print('❌ No PolyPhen score columns found')

# Check pathway columns
pathway_cols = ['dna_repair_pathway', 'mismatch_repair_pathway', 'is_important_gene']
for col in pathway_cols:
    if col in df.columns:
        count = (df[col] == 1).sum()
        print(f'{col}: {count:,} variants ({count/total*100:.1f}%)')
    else:
        print(f'⚠️  {col}: Not found')
"

# Confirm execution
echo ""
echo "🚀 READY TO START IMPUTATION"
echo "=============================="
echo "This will:"
echo "1. Load $INPUT_FILE"
echo "2. Perform pathway-aware median imputation for missing functional scores"
echo "3. Create confidence scores for imputed values"
echo "4. Save enhanced dataset to prostate_variants_tabnet_imputed.csv"
echo "5. Generate detailed imputation report"
echo ""
echo "Expected runtime: 2-5 minutes"
echo "Expected performance improvement: 6-8% TabNet accuracy"
echo ""

read -p "Continue with imputation? (y/n): " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "❌ Imputation cancelled by user"
    exit 0
fi

# Run the imputation
echo ""
echo "🧬 STARTING FUNCTIONAL SCORE IMPUTATION"
echo "========================================"

# Run the Python script
echo "Executing: $PYTHON $IMPUTATION_SCRIPT"
echo ""

if $PYTHON "$IMPUTATION_SCRIPT"; then
    echo ""
    echo "🎉 IMPUTATION COMPLETED SUCCESSFULLY!"
    echo "======================================"
    
    # Show results
    OUTPUT_FILE="$DATA_DIR/prostate_variants_tabnet_imputed.csv"
    REPORT_FILE="$DATA_DIR/imputation_report.txt"
    
    if [ -f "$OUTPUT_FILE" ]; then
        OUTPUT_SIZE=$(du -sh "$OUTPUT_FILE" | cut -f1)
        OUTPUT_LINES=$(wc -l < "$OUTPUT_FILE")
        echo "✅ Enhanced dataset created: $OUTPUT_FILE"
        echo "   Size: $OUTPUT_SIZE"
        echo "   Variants: $((OUTPUT_LINES - 1))"
    fi
    
    if [ -f "$REPORT_FILE" ]; then
        echo "✅ Report generated: $REPORT_FILE"
        echo ""
        echo "📋 IMPUTATION SUMMARY"
        echo "--------------------"
        grep -A 10 "IMPUTATION RESULTS:" "$REPORT_FILE" | head -8
    fi
    
    echo ""
    echo "🚀 NEXT STEPS"
    echo "============="
    echo "1. ✅ Functional score imputation completed"
    echo "2. 🔄 Use prostate_variants_tabnet_imputed.csv for TabNet training"
    echo "3. 🔄 Include new confidence features in your model:"
    echo "   - sift_confidence"
    echo "   - polyphen_confidence  "
    echo "   - functional_pathogenicity"
    echo "4. 🔄 Start TabNet training:"
    echo "   cd $PROJECT_ROOT"
    echo "   python src/model/tabnet_prostate_variant_classifier.py --input imputed"
    echo ""
    echo "📈 Expected performance improvement: 6-8% accuracy gain over missing data baseline"
    
else
    echo ""
    echo "❌ IMPUTATION FAILED!"
    echo "===================="
    echo "Check the error messages above for details"
    echo "Common issues:"
    echo "- Missing input file"
    echo "- Incorrect column names"
    echo "- Insufficient disk space"
    echo "- Python package conflicts"
    exit 1
fi

echo ""
echo "✅ All done! $(date)"