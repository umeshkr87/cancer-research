#!/bin/bash
# Corrected Validation Commands for Source CSV Files
# Run these to confirm REF/ALT allele data exists before merger

echo "=== COSMIC PROSTATE CSV VALIDATION ==="
COSMIC_FILE="/u/aa107/uiuc-cancer-research/data/processed/cosmic_prostate/cosmic_prostate.csv"

echo "1. Check COSMIC column headers:"
head -1 $COSMIC_FILE

echo -e "\n2. Look for allele-related columns:"
head -1 $COSMIC_FILE | tr ',' '\n' | grep -i -E "(allele|ref|alt|wt|mut)"

echo -e "\n3. Sample COSMIC allele data (first 5 rows):"
# GENOMIC_WT_ALLELE is column 18, GENOMIC_MUT_ALLELE is column 19
cut -d',' -f18,19 $COSMIC_FILE | head -5

echo -e "\n4. Count non-empty allele entries:"
tail -n +2 $COSMIC_FILE | cut -d',' -f18 | grep -v "^$" | wc -l

echo -e "\n=== CLINVAR PROSTATE CSV VALIDATION ==="
CLINVAR_FILE="/u/aa107/uiuc-cancer-research/data/processed/clinvar_prostate/clinvar_prostate.csv"

echo "1. Check ClinVar column headers:"
head -1 $CLINVAR_FILE

echo -e "\n2. Look for allele-related columns:"
head -1 $CLINVAR_FILE | tr ',' '\n' | grep -i -E "(reference|alternate|ref|alt|allele)"

echo -e "\n3. Sample ClinVar allele data (first 5 rows):"
# reference is column 4, alternate is column 5
cut -d',' -f4,5 $CLINVAR_FILE | head -5

echo -e "\n4. Count non-empty allele entries:"
tail -n +2 $CLINVAR_FILE | cut -d',' -f4 | grep -v "^$" | wc -l

echo -e "\n=== TCGA-PRAD CSV VALIDATION ==="
TCGA_FILE="/u/aa107/uiuc-cancer-research/data/processed/tcga_prad_prostate/tcga_prad_mutations.csv"

echo "1. Check TCGA column headers:"
head -1 $TCGA_FILE

echo -e "\n2. Look for allele-related columns:"
head -1 $TCGA_FILE | tr ',' '\n' | grep -i -E "(reference|tumor|allele|ref|alt)"

echo -e "\n3. Sample TCGA allele data (first 5 rows):"
# Reference_Allele is column 7, Tumor_Seq_Allele2 is column 8
cut -d',' -f7,8 $TCGA_FILE | head -5

echo -e "\n4. Count non-empty allele entries:"
tail -n +2 $TCGA_FILE | cut -d',' -f7 | grep -v "^$" | wc -l

echo -e "\n=== MERGED DATASET VALIDATION ==="
MERGED_FILE="/u/aa107/uiuc-cancer-research/data/processed/merged/merged_prostate_variants.csv"

if [ -f "$MERGED_FILE" ]; then
    echo "1. Check merged dataset has REF/ALT columns:"
    head -1 $MERGED_FILE | tr ',' '\n' | grep -E "(reference|alternate)"
    
    echo -e "\n2. Sample merged allele data (first 5 rows):"
    # Find column numbers for reference and alternate
    REF_COL=$(head -1 $MERGED_FILE | tr ',' '\n' | grep -n "reference" | cut -d':' -f1)
    ALT_COL=$(head -1 $MERGED_FILE | tr ',' '\n' | grep -n "alternate" | cut -d':' -f1)
    
    if [ ! -z "$REF_COL" ] && [ ! -z "$ALT_COL" ]; then
        cut -d',' -f${REF_COL},${ALT_COL} $MERGED_FILE | head -5
        
        echo -e "\n3. Count merged variants with both REF/ALT:"
        tail -n +2 $MERGED_FILE | cut -d',' -f${REF_COL},${ALT_COL} | grep -v "^,$\|^,\|,$" | wc -l
        
        echo -e "\n4. VEP compatibility check (sample alleles):"
        tail -n +2 $MERGED_FILE | cut -d',' -f${REF_COL},${ALT_COL} | grep -E "^[ATCG]+,[ATCG]+$" | head -5
    else
        echo "❌ Reference/alternate columns not found in merged dataset!"
    fi
else
    echo "❌ Merged dataset not found: $MERGED_FILE"
fi

echo -e "\n=== VALIDATION SUMMARY ==="
echo "✅ SUCCESS CRITERIA:"
echo "- COSMIC: Should show GENOMIC_WT_ALLELE, GENOMIC_MUT_ALLELE columns with A/T/C/G values"
echo "- ClinVar: Should show reference, alternate columns with nucleotide sequences"  
echo "- TCGA: Should show Reference_Allele, Tumor_Seq_Allele2 columns with A/T/C/G values"
echo "- MERGED: Should show reference, alternate columns with valid nucleotide data"
echo ""
echo "🎯 VEP READINESS:"
echo "- Merged dataset should have >200K variants with valid REF/ALT alleles"
echo "- Alleles should be nucleotide sequences (A/T/C/G), not placeholder values"
echo "- Ready for VCF conversion and VEP annotation pipeline"