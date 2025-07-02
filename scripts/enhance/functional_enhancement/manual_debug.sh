# ================================================================
# MANUAL DEBUGGING COMMANDS FOR ALPHAMISSENSE INTEGRATION
# Run these on the UIUC Campus Cluster to diagnose the 0% coverage issue
# ================================================================

# === 1. ANALYZE VEP DATA FORMAT ===
echo "=== VEP DATA ANALYSIS ==="
VEP_FILE="/u/aa107/uiuc-cancer-research/data/processed/tabnet_csv/prostate_variants_tabnet.csv"

echo "1. Check VEP file exists and size:"
ls -lh "$VEP_FILE"

echo -e "\n2. Check VEP column headers:"
head -1 "$VEP_FILE" | tr ',' '\n' | grep -E "(chromosome|position|reference|alternate)" | nl

echo -e "\n3. Sample VEP chromosome formats:"
tail -n +2 "$VEP_FILE" | cut -d',' -f1 | head -10 | sort | uniq

echo -e "\n4. Sample VEP positions:"
tail -n +2 "$VEP_FILE" | cut -d',' -f2 | head -5

echo -e "\n5. Sample VEP alleles (REF,ALT):"
# Assuming columns 3,4 are reference_allele,alternate_allele
tail -n +2 "$VEP_FILE" | cut -d',' -f3,4 | head -5

echo -e "\n6. Create sample VEP lookup keys:"
tail -n +2 "$VEP_FILE" | head -5 | awk -F',' '{print $1"_"$2"_"$3"_"$4}'

# === 2. ANALYZE ALPHAMISSENSE DATA FORMAT ===
echo -e "\n=== ALPHAMISSENSE DATA ANALYSIS ==="
AM_FILE="/u/aa107/scratch/alphamissense/AlphaMissense_hg38.tsv"

echo "1. Check AlphaMissense file exists and size:"
ls -lh "$AM_FILE"

echo -e "\n2. Check AlphaMissense header:"
head -5 "$AM_FILE" | grep "^#CHROM"

echo -e "\n3. Sample AlphaMissense chromosome formats:"
tail -n +4 "$AM_FILE" | cut -f1 | head -10 | sort | uniq

echo -e "\n4. Sample AlphaMissense positions:"
tail -n +4 "$AM_FILE" | cut -f2 | head -5

echo -e "\n5. Sample AlphaMissense alleles (REF,ALT):"
tail -n +4 "$AM_FILE" | cut -f3,4 | head -5

echo -e "\n6. Create sample AlphaMissense lookup keys:"
tail -n +4 "$AM_FILE" | head -5 | awk -F'\t' '{print $1"_"$2"_"$3"_"$4}'

# === 3. DIRECT COMPARISON ===
echo -e "\n=== DIRECT FORMAT COMPARISON ==="

echo "VEP sample lookup keys:"
tail -n +2 "$VEP_FILE" | head -3 | awk -F',' '{print $1"_"$2"_"$3"_"$4}'

echo -e "\nAlphaMissense sample lookup keys:"
tail -n +4 "$AM_FILE" | head -3 | awk -F'\t' '{print $1"_"$2"_"$3"_"$4}'

# === 4. TEST CHROMOSOME NORMALIZATION ===
echo -e "\n=== CHROMOSOME NORMALIZATION TEST ==="

echo "VEP chromosomes (original):"
tail -n +2 "$VEP_FILE" | cut -d',' -f1 | head -5

echo -e "\nVEP chromosomes (remove chr prefix if present):"
tail -n +2 "$VEP_FILE" | cut -d',' -f1 | head -5 | sed 's/^chr//i'

echo -e "\nAlphaMissense chromosomes (original):"
tail -n +4 "$AM_FILE" | cut -f1 | head -5

echo -e "\nAlphaMissense chromosomes (remove chr prefix if present):"
tail -n +4 "$AM_FILE" | cut -f1 | head -5 | sed 's/^chr//i'

# === 5. SPECIFIC VARIANT LOOKUP TEST ===
echo -e "\n=== SPECIFIC VARIANT LOOKUP TEST ==="

# Get first variant from VEP
echo "First VEP variant:"
FIRST_VEP=$(tail -n +2 "$VEP_FILE" | head -1)
echo "$FIRST_VEP"

# Extract components
CHR=$(echo "$FIRST_VEP" | cut -d',' -f1)
POS=$(echo "$FIRST_VEP" | cut -d',' -f2)
REF=$(echo "$FIRST_VEP" | cut -d',' -f3)
ALT=$(echo "$FIRST_VEP" | cut -d',' -f4)

echo "Components: CHR=$CHR, POS=$POS, REF=$REF, ALT=$ALT"

# Search for this variant in AlphaMissense (exact match)
echo -e "\nSearching for exact match in AlphaMissense:"
grep -F "$CHR	$POS	$REF	$ALT" "$AM_FILE" | head -1

# Search with chromosome normalization
CHR_NORM=$(echo "$CHR" | sed 's/^chr//i')
echo -e "\nSearching with normalized chromosome ($CHR_NORM):"
grep -F "$CHR_NORM	$POS	$REF	$ALT" "$AM_FILE" | head -1

# Search with chr prefix added
CHR_CHR="chr$CHR_NORM"
echo -e "\nSearching with chr prefix ($CHR_CHR):"
grep -F "$CHR_CHR	$POS	$REF	$ALT" "$AM_FILE" | head -1

# === 6. QUICK STATISTICS ===
echo -e "\n=== QUICK STATISTICS ==="

echo "VEP total variants:"
tail -n +2 "$VEP_FILE" | wc -l

echo "AlphaMissense total entries:"
tail -n +4 "$AM_FILE" | wc -l

echo "VEP unique chromosomes:"
tail -n +2 "$VEP_FILE" | cut -d',' -f1 | sort | uniq | tr '\n' ' '

echo -e "\nAlphaMissense unique chromosomes:"
tail -n +4 "$AM_FILE" | cut -f1 | sort | uniq | head -10 | tr '\n' ' '

# === 7. PYTHON QUICK TEST ===
echo -e "\n=== PYTHON QUICK TEST ==="

python3 << 'EOF'
import pandas as pd
import sys

print("🐍 Python Quick Diagnostic")

# Load sample data
vep_file = "/u/aa107/uiuc-cancer-research/data/processed/tabnet_csv/prostate_variants_tabnet.csv"
am_file = "/u/aa107/scratch/alphamissense/AlphaMissense_hg38.tsv"

try:
    # VEP data
    vep_sample = pd.read_csv(vep_file, nrows=5)
    print(f"✅ VEP sample loaded: {len(vep_sample)} rows")
    print(f"VEP columns: {list(vep_sample.columns[:5])}")
    
    # AlphaMissense data
    am_sample = pd.read_csv(am_file, sep='\t', skiprows=3, nrows=5)
    print(f"✅ AlphaMissense sample loaded: {len(am_sample)} rows")
    print(f"AM columns: {list(am_sample.columns[:5])}")
    
    # Show sample data types
    print(f"\nVEP data types:")
    print(f"  chromosome: {vep_sample['chromosome'].dtype}")
    print(f"  position: {vep_sample['position'].dtype}")
    
    print(f"\nAM data types:")
    print(f"  #CHROM: {am_sample['#CHROM'].dtype}")
    print(f"  POS: {am_sample['POS'].dtype}")
    
except Exception as e:
    print(f"❌ Error: {e}")
    sys.exit(1)
EOF

echo -e "\n=== DEBUGGING COMPLETE ==="
echo "Run these commands to identify the format mismatch causing 0% coverage."
echo "Then use the comprehensive debug script to fix the issue."