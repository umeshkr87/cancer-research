#!/bin/bash
# Validation Commands for Source CSV Files
# Run these to confirm REF/ALT allele data exists before merger

echo "=== COSMIC PROSTATE CSV VALIDATION ==="
COSMIC_FILE="/u/aa107/uiuc-cancer-research/data/processed/cosmic_prostate/cosmic_prostate.csv"

echo "1. Check COSMIC column headers:"
head -1 $COSMIC_FILE

echo -e "\n2. Look for allele-related columns:"
head -1 $COSMIC_FILE | tr ',' '\n' | grep -i -E "(allele|ref|alt|wt|mut)"

echo -e "\n3. Sample COSMIC allele data (first 5 rows):"
csvcut -c "GENOMIC_WT_ALLELE,GENOMIC_MUT_ALLELE" $COSMIC_FILE | head -5

echo -e "\n4. Count non-empty allele entries:"
tail -n +2 $COSMIC_FILE | cut -d',' -f$(head -1 $COSMIC_FILE | tr ',' '\n' | grep -n "GENOMIC_WT_ALLELE" | cut -d':' -f1) | grep -v "^$" | wc -l

echo -e "\n=== CLINVAR PROSTATE CSV VALIDATION ==="
CLINVAR_FILE="/u/aa107/uiuc-cancer-research/data/processed/clinvar_prostate/clinvar_prostate.csv"

echo "1. Check ClinVar column headers:"
head -1 $CLINVAR_FILE

echo -e "\n2. Look for allele-related columns:"
head -1 $CLINVAR_FILE | tr ',' '\n' | grep -i -E "(reference|alternate|ref|alt|allele)"

echo -e "\n3. Sample ClinVar allele data (first 5 rows):"
csvcut -c "reference,alternate" $CLINVAR_FILE | head -5

echo -e "\n4. Count non-empty allele entries:"
tail -n +2 $CLINVAR_FILE | cut -d',' -f$(head -1 $CLINVAR_FILE | tr ',' '\n' | grep -n "reference" | cut -d':' -f1) | grep -v "^$" | wc -l

echo -e "\n=== TCGA-PRAD CSV VALIDATION ==="
TCGA_FILE="/u/aa107/uiuc-cancer-research/data/processed/tcga_prad_prostate/tcga_prad_mutations.csv"

echo "1. Check TCGA column headers:"
head -1 $TCGA_FILE

echo -e "\n2. Look for allele-related columns:"
head -1 $TCGA_FILE | tr ',' '\n' | grep -i -E "(reference|tumor|allele|ref|alt)"

echo -e "\n3. Sample TCGA allele data (first 5 rows):"
csvcut -c "Reference_Allele,Tumor_Seq_Allele2" $TCGA_FILE | head -5

echo -e "\n4. Count non-empty allele entries:"
tail -n +2 $TCGA_FILE | cut -d',' -f$(head -1 $TCGA_FILE | tr ',' '\n' | grep -n "Reference_Allele" | cut -d':' -f1) | grep -v "^$" | wc -l

echo -e "\n=== VALIDATION SUMMARY ==="
echo "Expected results for successful validation:"
echo "- COSMIC: Should show GENOMIC_WT_ALLELE, GENOMIC_MUT_ALLELE columns with A/T/C/G values"
echo "- ClinVar: Should show reference, alternate columns with nucleotide sequences"  
echo "- TCGA: Should show Reference_Allele, Tumor_Seq_Allele2 columns with A/T/C/G values"
echo ""
echo "If any file shows missing columns or empty values, that source needs fixing"
echo "If all show proper allele data, then merger script is the issue"

# Alternative quick check commands (if csvcut not available)
echo -e "\n=== ALTERNATIVE COMMANDS (if csvcut not installed) ==="
echo "# Quick column check for each file:"
echo "head -1 $COSMIC_FILE | tr ',' '\n' | nl"
echo "head -1 $CLINVAR_FILE | tr ',' '\n' | nl" 
echo "head -1 $TCGA_FILE | tr ',' '\n' | nl"

echo -e "\n# Sample data check (replace column numbers based on above):"
echo "cut -d',' -f[COL_NUM] $COSMIC_FILE | head -5"
echo "cut -d',' -f[COL_NUM] $CLINVAR_FILE | head -5"
echo "cut -d',' -f[COL_NUM] $TCGA_FILE | head -5"