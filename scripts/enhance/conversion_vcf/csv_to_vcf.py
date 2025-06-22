#!/usr/bin/env python3
"""
Simple CSV to VCF Converter for TabNet Prostate Cancer Dataset
Location: /u/aa107/uiuc-cancer-research/scripts/enhance/conversion_vcf/csv_to_vcf.py
"""

import pandas as pd
import sys
from pathlib import Path
from datetime import datetime

def csv_to_vcf(input_csv, output_vcf):
    """Simple, reliable CSV to VCF conversion"""
    
    print(f"🔄 Converting {input_csv} to VCF format...")
    
    # Load data
    df = pd.read_csv(input_csv, low_memory=False)
    print(f"✅ Loaded {len(df):,} variants")
    
    # Clean and validate essential columns
    df = df.dropna(subset=['chromosome', 'position', 'gene_symbol'])
    df['chromosome'] = df['chromosome'].astype(str).str.replace('chr', '', case=False)
    df['position'] = df['position'].astype(int)
    
    print(f"✅ {len(df):,} variants with valid coordinates")
    
    # Create output directory
    Path(output_vcf).parent.mkdir(parents=True, exist_ok=True)
    
    # Write VCF
    with open(output_vcf, 'w') as f:
        # Write header
        f.write("##fileformat=VCFv4.2\n")
        f.write(f"##fileDate={datetime.now().strftime('%Y%m%d')}\n")
        f.write("##reference=GRCh38\n")
        f.write("##INFO=<ID=GENE,Number=1,Type=String,Description=\"Gene symbol\">\n")
        f.write("##INFO=<ID=SOURCE,Number=1,Type=String,Description=\"Data source\">\n")
        f.write("##INFO=<ID=CLASS,Number=1,Type=String,Description=\"Variant classification\">\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        
        # Write variants
        for _, row in df.iterrows():
            chrom = row['chromosome']
            pos = int(row['position'])
            variant_id = f"{chrom}:{pos}:{row['gene_symbol']}"
            ref = "N"  # Placeholder - annotation tools will add real alleles
            alt = "."  # Placeholder
            qual = "."
            filter_val = "PASS"
            
            # Build INFO
            info_parts = [f"GENE={row['gene_symbol']}"]
            if 'data_source' in row and pd.notna(row['data_source']):
                info_parts.append(f"SOURCE={row['data_source']}")
            if 'variant_classification' in row and pd.notna(row['variant_classification']):
                info_parts.append(f"CLASS={row['variant_classification']}")
            
            info = ";".join(info_parts)
            
            f.write(f"{chrom}\t{pos}\t{variant_id}\t{ref}\t{alt}\t{qual}\t{filter_val}\t{info}\n")
    
    print(f"✅ VCF created: {output_vcf}")
    print(f"📊 File size: {Path(output_vcf).stat().st_size / 1024 / 1024:.1f} MB")
    
    return True

def main():
    """Main function"""
    
    # Define paths
    project_root = Path(__file__).parent.parent.parent.parent
    input_csv = project_root / "data/processed/merged/merged_prostate_variants.csv"
    output_vcf = project_root / "data/processed/merged_vcf/merged_prostate_variants.vcf"
    
    print("🧬 Simple CSV to VCF Converter")
    print("==============================")
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
            print(f"\n🎯 SUCCESS! Ready for VarStack annotation")
            print(f"📁 VCF file: {output_vcf}")
        return success
    except Exception as e:
        print(f"❌ Conversion failed: {e}")
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)