#!/usr/bin/env python3
"""
Corrected CSV to VCF Converter for TabNet Prostate Cancer Dataset
Reads actual REF/ALT alleles from merged CSV instead of using placeholders
Location: /u/aa107/uiuc-cancer-research/scripts/enhance/conversion_vcf/csv_to_vcf_corrected.py
"""

import pandas as pd
import sys
from pathlib import Path
from datetime import datetime
import re

def validate_allele(allele):
    """Validate nucleotide allele sequence"""
    if pd.isna(allele) or str(allele).strip() == '':
        return None
    
    # Clean and validate
    allele_clean = str(allele).strip().upper()
    
    # Check for valid nucleotides only
    if re.match(r'^[ATCG-]+$', allele_clean) and len(allele_clean) > 0:
        return allele_clean
    
    return None

def csv_to_vcf(input_csv, output_vcf):
    """Convert CSV with real REF/ALT alleles to VCF format"""
    
    print(f"🔄 Converting {input_csv} to VCF format with real alleles...")
    
    # Load data
    df = pd.read_csv(input_csv, low_memory=False)
    print(f"✅ Loaded {len(df):,} variants")
    
    # Validate required columns
    required_cols = ['chromosome', 'position', 'gene_symbol', 'reference', 'alternate']
    missing_cols = [col for col in required_cols if col not in df.columns]
    
    if missing_cols:
        print(f"❌ Missing required columns: {missing_cols}")
        return False
    
    # Clean and validate data
    print("🧹 Cleaning and validating variant data...")
    
    # Clean chromosomes and positions
    df = df.dropna(subset=['chromosome', 'position', 'gene_symbol'])
    df['chromosome'] = df['chromosome'].astype(str).str.replace('chr', '', case=False)
    df['position'] = df['position'].astype(int)
    
    # Validate and clean alleles
    df['ref_clean'] = df['reference'].apply(validate_allele)
    df['alt_clean'] = df['alternate'].apply(validate_allele)
    
    # Keep only variants with valid REF/ALT alleles
    valid_variants = df.dropna(subset=['ref_clean', 'alt_clean'])
    
    print(f"✅ Valid variants with REF/ALT alleles: {len(valid_variants):,}/{len(df):,} ({len(valid_variants)/len(df)*100:.1f}%)")
    
    if len(valid_variants) == 0:
        print("❌ No valid variants found with proper REF/ALT alleles")
        return False
    
    # Create output directory
    Path(output_vcf).parent.mkdir(parents=True, exist_ok=True)
    
    # Write VCF
    print("📝 Writing VCF file...")
    with open(output_vcf, 'w') as f:
        # Write VCF header
        f.write("##fileformat=VCFv4.2\n")
        f.write(f"##fileDate={datetime.now().strftime('%Y%m%d')}\n")
        f.write("##reference=GRCh38\n")
        f.write("##source=TabNet_Prostate_Cancer_Dataset\n")
        f.write("##INFO=<ID=GENE,Number=1,Type=String,Description=\"Gene symbol\">\n")
        f.write("##INFO=<ID=SOURCE,Number=1,Type=String,Description=\"Data source\">\n")
        f.write("##INFO=<ID=CLASS,Number=1,Type=String,Description=\"Variant classification\">\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        
        # Write variants
        written_count = 0
        for _, row in valid_variants.iterrows():
            chrom = str(row['chromosome'])
            pos = int(row['position'])
            ref = row['ref_clean']
            alt = row['alt_clean']
            
            # Create variant ID
            variant_id = f"{chrom}:{pos}:{ref}:{alt}"
            
            # Standard VCF fields
            qual = "."
            filter_val = "PASS"
            
            # Build INFO field
            info_parts = [f"GENE={row['gene_symbol']}"]
            
            if 'data_source' in row and pd.notna(row['data_source']):
                info_parts.append(f"SOURCE={row['data_source']}")
            
            if 'variant_classification' in row and pd.notna(row['variant_classification']):
                info_parts.append(f"CLASS={row['variant_classification']}")
            
            info = ";".join(info_parts)
            
            # Write VCF line
            f.write(f"{chrom}\t{pos}\t{variant_id}\t{ref}\t{alt}\t{qual}\t{filter_val}\t{info}\n")
            written_count += 1
    
    print(f"✅ VCF created: {output_vcf}")
    print(f"📊 Written {written_count:,} variants with valid REF/ALT alleles")
    print(f"📁 File size: {Path(output_vcf).stat().st_size / 1024 / 1024:.1f} MB")
    
    return True

def main():
    """Main function"""
    
    # Define paths
    project_root = Path(__file__).parent.parent.parent.parent
    input_csv = project_root / "data/processed/merged/merged_prostate_variants.csv"
    output_vcf = project_root / "data/processed/merged_vcf/merged_prostate_variants.vcf"
    
    print("🧬 Corrected CSV to VCF Converter")
    print("=================================")
    print("Now reads actual REF/ALT alleles from merged dataset!")
    print(f"Input: {input_csv}")
    print(f"Output: {output_vcf}")
    
    # Check input
    if not input_csv.exists():
        print(f"❌ Input file not found: {input_csv}")
        return False
    
    # Convert
    try:
        success = csv_to_vcf(input_csv, output_vcf)
        if success:
            print(f"\n🎯 SUCCESS! VCF ready for VEP annotation")
            print(f"📁 VCF file: {output_vcf}")
            print(f"🚀 Expected VEP result: >200K variants processed (vs previous 0)")
        return success
    except Exception as e:
        print(f"❌ Conversion failed: {e}")
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)