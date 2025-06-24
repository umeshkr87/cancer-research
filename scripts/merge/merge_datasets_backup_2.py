#!/usr/bin/env python3
"""
Enhanced merger for COSMIC, ClinVar, and TCGA-PRAD datasets.
Now includes REF/ALT allele column mapping for VEP compatibility.
Location: /u/aa107/uiuc-cancer-research/scripts/merge/merge_datasets.py
"""

import pandas as pd
import numpy as np
import re
from pathlib import Path

def clean_chromosome_column(df):
    """Clean chromosome column to handle mixed types"""
    # Convert to string and handle NaN values
    df['chromosome'] = df['chromosome'].astype(str).replace('nan', 'Unknown')
    # Remove any 'chr' prefix if present
    df['chromosome'] = df['chromosome'].str.replace('chr', '', case=False)
    # Replace 'Unknown' back to NaN for proper handling
    df['chromosome'] = df['chromosome'].replace('Unknown', np.nan)
    return df

def validate_allele(allele_value):
    """Validate and clean allele sequences"""
    if pd.isna(allele_value) or str(allele_value).strip() == '':
        return None
    
    # Convert to string and clean
    allele = str(allele_value).strip().upper()
    
    # Handle common placeholder values
    if allele in ['', '.', '-', 'NULL', 'NA', 'NAN']:
        return None
    
    # Remove any non-ATCG characters but preserve valid indel notation
    # Allow A, T, C, G and - for deletions
    cleaned_allele = re.sub(r'[^ATCG-]', '', allele)
    
    # Return None if completely invalid
    if not cleaned_allele or cleaned_allele == '-':
        return None
        
    return cleaned_allele

def process_allele_columns(df, ref_col, alt_col):
    """Process and validate reference and alternate allele columns"""
    # Create standardized columns
    df['reference'] = df[ref_col].apply(validate_allele) if ref_col in df.columns else None
    df['alternate'] = df[alt_col].apply(validate_allele) if alt_col in df.columns else None
    
    # Count valid alleles for reporting
    valid_ref = df['reference'].notna().sum()
    valid_alt = df['alternate'].notna().sum()
    both_valid = (df['reference'].notna() & df['alternate'].notna()).sum()
    
    return df, valid_ref, valid_alt, both_valid

def main():
    """Main function to merge the three processed datasets with allele support"""
    
    print("🧬 Enhanced TabNet Prostate Cancer Dataset Merger")
    print("===============================================")
    print("Now with REF/ALT allele support for VEP compatibility!")
    
    # Define file paths
    project_root = Path(__file__).parent.parent.parent
    
    cosmic_path = project_root / "data/processed/cosmic_prostate/cosmic_prostate.csv"
    clinvar_path = project_root / "data/processed/clinvar_prostate/clinvar_prostate.csv"
    tcga_path = project_root / "data/processed/tcga_prad_prostate/tcga_prad_mutations.csv"
    
    # Check if files exist
    missing_files = []
    for name, path in [("COSMIC", cosmic_path), ("ClinVar", clinvar_path), ("TCGA", tcga_path)]:
        if not path.exists():
            missing_files.append(f"{name}: {path}")
        else:
            print(f"✅ Found {name}: {path}")
    
    if missing_files:
        print("\n❌ Missing required files:")
        for file in missing_files:
            print(f"   {file}")
        print("\nPlease ensure all processed datasets are available.")
        return False
    
    # Load datasets with proper dtype handling
    print("\n📊 Loading datasets...")
    
    # Load COSMIC
    cosmic_df = pd.read_csv(cosmic_path, low_memory=False)
    cosmic_df['data_source'] = 'COSMIC'
    print(f"   COSMIC: {len(cosmic_df):,} variants")
    
    # Load ClinVar  
    clinvar_df = pd.read_csv(clinvar_path, low_memory=False)
    clinvar_df['data_source'] = 'ClinVar'
    print(f"   ClinVar: {len(clinvar_df):,} variants")
    
    # Load TCGA
    tcga_df = pd.read_csv(tcga_path, low_memory=False)
    tcga_df['data_source'] = 'TCGA'
    print(f"   TCGA: {len(tcga_df):,} variants")
    
    # Standardize column names and process alleles
    print("\n🔧 Standardizing columns and processing alleles...")
    
    # COSMIC columns + alleles
    cosmic_clean = cosmic_df.rename(columns={
        'gene': 'gene_symbol',
        'chr': 'chromosome',
        'pos': 'position'
    }).copy()
    
    # Process COSMIC alleles
    cosmic_clean, cosmic_ref, cosmic_alt, cosmic_both = process_allele_columns(
        cosmic_clean, 'GENOMIC_WT_ALLELE', 'GENOMIC_MUT_ALLELE'
    )
    print(f"   COSMIC alleles: {cosmic_both}/{len(cosmic_clean)} variants with both REF/ALT")
    
    # ClinVar columns + alleles (already correctly named)
    clinvar_clean = clinvar_df.rename(columns={
        'gene': 'gene_symbol',
        'chr': 'chromosome',
        'pos': 'position'
    }).copy()
    
    # Process ClinVar alleles (already named reference/alternate)
    clinvar_clean, clinvar_ref, clinvar_alt, clinvar_both = process_allele_columns(
        clinvar_clean, 'reference', 'alternate'
    )
    print(f"   ClinVar alleles: {clinvar_both}/{len(clinvar_clean)} variants with both REF/ALT")
    
    # TCGA columns + alleles
    tcga_clean = tcga_df.rename(columns={
        'Hugo_Symbol': 'gene_symbol',
        'Chromosome': 'chromosome', 
        'Start_Position': 'position'
    }).copy()
    
    # Process TCGA alleles
    tcga_clean, tcga_ref, tcga_alt, tcga_both = process_allele_columns(
        tcga_clean, 'Reference_Allele', 'Tumor_Seq_Allele2'
    )
    print(f"   TCGA alleles: {tcga_both}/{len(tcga_clean)} variants with both REF/ALT")
    
    # Clean chromosome columns for all datasets
    cosmic_clean = clean_chromosome_column(cosmic_clean)
    clinvar_clean = clean_chromosome_column(clinvar_clean)
    tcga_clean = clean_chromosome_column(tcga_clean)
    
    # Combine datasets
    print("\n🔀 Merging datasets...")
    all_variants = pd.concat([cosmic_clean, clinvar_clean, tcga_clean], ignore_index=True)
    print(f"   Total variants: {len(all_variants):,}")
    
    # Validate merged allele data
    total_with_alleles = (all_variants['reference'].notna() & all_variants['alternate'].notna()).sum()
    allele_coverage = (total_with_alleles / len(all_variants)) * 100
    print(f"   Variants with both REF/ALT: {total_with_alleles:,} ({allele_coverage:.1f}%)")
    
    # Add prostate pathway features
    print("\n🧬 Adding prostate pathway features...")
    
    # Define prostate cancer gene sets
    dna_repair_genes = {'BRCA1', 'BRCA2', 'ATM', 'CHEK2', 'PALB2', 'RAD51D', 'BRIP1', 'FANCA', 'MLH1', 'MSH2', 'MSH6', 'PMS2', 'NBN'}
    hormone_genes = {'AR', 'CYP17A1', 'SRD5A2', 'CYP19A1', 'ESR1', 'ESR2'}
    pi3k_genes = {'PTEN', 'PIK3CA', 'PIK3R1', 'AKT1', 'AKT2', 'AKT3', 'MTOR', 'TSC1', 'TSC2'}
    core_prostate_genes = {'AR', 'PTEN', 'TP53', 'MYC', 'ERG', 'SPOP', 'FOXA1', 'CHD1'}
    
    # Add pathway binary features
    all_variants['dna_repair_gene'] = all_variants['gene_symbol'].isin(dna_repair_genes).astype(int)
    all_variants['hormone_pathway_gene'] = all_variants['gene_symbol'].isin(hormone_genes).astype(int)
    all_variants['pi3k_pathway_gene'] = all_variants['gene_symbol'].isin(pi3k_genes).astype(int)
    all_variants['core_prostate_gene'] = all_variants['gene_symbol'].isin(core_prostate_genes).astype(int)
    
    # Add therapeutic targets
    parp_targets = dna_repair_genes
    all_variants['parp_inhibitor_target'] = all_variants['gene_symbol'].isin(parp_targets).astype(int)
    all_variants['hormone_therapy_target'] = all_variants['gene_symbol'].isin(hormone_genes).astype(int)
    
    # Create variant classification for TabNet
    print("\n🎯 Creating variant classification...")
    
    def classify_variant(row):
        # Use ClinVar clinical significance if available
        if 'clinical_significance' in row and pd.notna(row['clinical_significance']):
            sig = str(row['clinical_significance']).lower()
            if 'pathogenic' in sig and 'likely' not in sig:
                return 'Actionable_Pathogenic'
            elif 'likely pathogenic' in sig:
                return 'Likely_Actionable'
            elif 'benign' in sig and 'likely' not in sig:
                return 'Benign'
            elif 'likely benign' in sig:
                return 'Likely_Benign'
            else:
                return 'VUS'
        
        # For variants without ClinVar annotation
        if row['data_source'] == 'COSMIC':
            return 'Likely_Actionable'
        elif row['data_source'] == 'TCGA':
            return 'VUS'
        else:
            return 'VUS'
    
    all_variants['variant_classification'] = all_variants.apply(classify_variant, axis=1)
    
    # Add basic genomic features
    all_variants['chromosome_num'] = pd.to_numeric(all_variants['chromosome'], errors='coerce')
    all_variants['is_sex_chromosome'] = all_variants['chromosome'].isin(['X', 'Y']).astype(int)
    
    # Save merged dataset
    print("\n💾 Saving enhanced merged dataset...")
    output_dir = project_root / "data/processed/merged"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    output_file = output_dir / "merged_prostate_variants.csv"
    all_variants.to_csv(output_file, index=False)
    
    # Generate comprehensive report with allele statistics
    print("\n📋 Generating enhanced report...")
    report_file = output_dir / "merge_report.txt"
    
    # Clean chromosome data for reporting
    valid_chromosomes = all_variants['chromosome'].dropna().unique()
    # Convert to string and sort properly
    chr_list = []
    for chr_val in valid_chromosomes:
        chr_str = str(chr_val)
        if chr_str.isdigit():
            chr_list.append(int(chr_str))
        else:
            chr_list.append(chr_str)
    
    # Sort chromosomes properly (numbers first, then letters)
    numeric_chrs = sorted([c for c in chr_list if isinstance(c, int)])
    alpha_chrs = sorted([c for c in chr_list if isinstance(c, str)])
    sorted_chromosomes = [str(c) for c in numeric_chrs] + alpha_chrs
    
    with open(report_file, 'w') as f:
        f.write("Enhanced TabNet Prostate Cancer Dataset Merge Report\n")
        f.write("===================================================\n\n")
        
        f.write("📊 DATASET SUMMARY\n")
        f.write("-" * 50 + "\n")
        f.write(f"Total variants: {len(all_variants):,}\n")
        f.write(f"Unique genes: {all_variants['gene_symbol'].nunique()}\n")
        f.write(f"Unique chromosomes: {sorted_chromosomes}\n")
        f.write(f"Data sources: {', '.join(all_variants['data_source'].unique())}\n\n")
        
        f.write("🧬 ALLELE DATA ANALYSIS (NEW)\n")
        f.write("-" * 50 + "\n")
        f.write(f"Variants with REF allele: {all_variants['reference'].notna().sum():,}\n")
        f.write(f"Variants with ALT allele: {all_variants['alternate'].notna().sum():,}\n")
        f.write(f"Variants with both REF/ALT: {total_with_alleles:,} ({allele_coverage:.1f}%)\n")
        f.write(f"VEP compatibility: {total_with_alleles:,} variants ready for annotation\n\n")
        
        f.write("📊 ALLELE DATA BY SOURCE\n")
        f.write("-" * 50 + "\n")
        f.write(f"COSMIC: {cosmic_both}/{len(cosmic_clean)} variants ({(cosmic_both/len(cosmic_clean)*100):.1f}%)\n")
        f.write(f"ClinVar: {clinvar_both}/{len(clinvar_clean)} variants ({(clinvar_both/len(clinvar_clean)*100):.1f}%)\n")
        f.write(f"TCGA: {tcga_both}/{len(tcga_clean)} variants ({(tcga_both/len(tcga_clean)*100):.1f}%)\n\n")
        
        f.write("📈 SOURCE DISTRIBUTION\n")
        f.write("-" * 50 + "\n")
        source_counts = all_variants['data_source'].value_counts()
        for source, count in source_counts.items():
            percentage = (count / len(all_variants)) * 100
            f.write(f"{source}: {count:,} ({percentage:.1f}%)\n")
        f.write("\n")
        
        f.write("🎯 VARIANT CLASSIFICATION DISTRIBUTION\n")
        f.write("-" * 50 + "\n")
        class_counts = all_variants['variant_classification'].value_counts()
        for class_name, count in class_counts.items():
            percentage = (count / len(all_variants)) * 100
            f.write(f"{class_name}: {count:,} ({percentage:.1f}%)\n")
        f.write("\n")
        
        f.write("🧬 PATHWAY GENE ANALYSIS\n")
        f.write("-" * 50 + "\n")
        pathway_features = ['dna_repair_gene', 'hormone_pathway_gene', 'pi3k_pathway_gene', 'core_prostate_gene']
        for feature in pathway_features:
            count = all_variants[feature].sum()
            percentage = (count / len(all_variants)) * 100
            pathway_name = feature.replace('_gene', '').replace('_', ' ').title()
            f.write(f"{pathway_name}: {count:,} variants ({percentage:.1f}%)\n")
        f.write("\n")
        
        f.write("💊 THERAPEUTIC TARGETS\n")
        f.write("-" * 50 + "\n")
        parp_count = all_variants['parp_inhibitor_target'].sum()
        hormone_count = all_variants['hormone_therapy_target'].sum()
        f.write(f"PARP Inhibitor Targets: {parp_count:,}\n")
        f.write(f"Hormone Therapy Targets: {hormone_count:,}\n\n")
        
        f.write("🔝 TOP 20 GENES\n")
        f.write("-" * 50 + "\n")
        top_genes = all_variants['gene_symbol'].value_counts().head(20)
        for gene, count in top_genes.items():
            f.write(f"{gene}: {count:,}\n")
        f.write("\n")
        
        f.write("📂 OUTPUT FILES\n")
        f.write("-" * 50 + "\n")
        f.write(f"Enhanced merged dataset: {output_file}\n")
        f.write(f"This report: {report_file}\n")
        f.write(f"Features ready for TabNet training: {len(all_variants.columns)}\n")
        f.write(f"Target variable: variant_classification (5 classes)\n")
        f.write(f"NEW: REF/ALT alleles included for VEP compatibility\n")
    
    print(f"   Report saved: {report_file}")
    
    # Print summary
    print(f"\n✅ SUCCESS! Enhanced merged dataset saved to:")
    print(f"   {output_file}")
    print(f"\n📊 Dataset Summary:")
    print(f"   Total variants: {len(all_variants):,}")
    print(f"   Unique genes: {all_variants['gene_symbol'].nunique()}")
    print(f"   Data sources: {', '.join(all_variants['data_source'].unique())}")
    print(f"   🧬 NEW: Variants with REF/ALT alleles: {total_with_alleles:,} ({allele_coverage:.1f}%)")
    
    print(f"\n🎯 Variant Classification Distribution:")
    for class_name, count in all_variants['variant_classification'].value_counts().items():
        print(f"   {class_name}: {count:,}")
    
    print(f"\n🧬 Pathway Gene Distribution:")
    pathway_features = ['dna_repair_gene', 'hormone_pathway_gene', 'pi3k_pathway_gene', 'core_prostate_gene']
    for feature in pathway_features:
        count = all_variants[feature].sum()
        print(f"   {feature}: {count:,}")
    
    print(f"\n🚀 VEP READY! Dataset now contains REF/ALT alleles for annotation!")
    print(f"   Previous VEP result: 0/213,533 variants processed")
    print(f"   Expected VEP result: {total_with_alleles:,} variants ready for processing")
    
    return True

if __name__ == "__main__":
    success = main()
    if not success:
        exit(1)