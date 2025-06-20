#!/usr/bin/env python3

import pandas as pd
import numpy as np

def filter_prostate_samples(cosmic_df):
    """Filter COSMIC data for prostate cancer samples"""
    
    print(f"Original COSMIC data: {len(cosmic_df)} mutations")
    
    # Method 1: Filter by TCGA-PRAD samples (most reliable)
    tcga_prad_samples = cosmic_df[
        cosmic_df['SAMPLE_NAME'].str.contains('TCGA-', na=False) &
        cosmic_df['SAMPLE_NAME'].str.contains('PRAD|prostate', case=False, na=False)
    ]
    
    # Method 2: Filter by known prostate cancer cell lines
    prostate_cell_lines = [
        'PC3', 'LNCaP', 'DU145', '22Rv1', 'VCaP', 'LAPC4', 
        'C4-2', 'MDA-PCa', 'RWPE-1', 'PWR-1E'
    ]
    
    cell_line_samples = cosmic_df[
        cosmic_df['SAMPLE_NAME'].str.contains('|'.join(prostate_cell_lines), case=False, na=False)
    ]
    
    # Method 3: You'll need to check COSMIC_PHENOTYPE_ID mapping
    # For now, we'll use a placeholder - you may need to download the phenotype mapping file
    
    # Combine all prostate samples
    prostate_samples = pd.concat([tcga_prad_samples, cell_line_samples]).drop_duplicates()
    
    print(f"Prostate cancer mutations found: {len(prostate_samples)}")
    print(f"Unique prostate samples: {prostate_samples['SAMPLE_NAME'].nunique()}")
    
    return prostate_samples

def standardize_cosmic_columns(cosmic_df):
    """Standardize COSMIC columns for merging with TCGA/ClinVar"""
    
    cosmic_clean = cosmic_df.copy()
    
    # Standardize column names for merging
    cosmic_clean['gene'] = cosmic_clean['GENE_SYMBOL']
    cosmic_clean['chr'] = cosmic_clean['CHROMOSOME'].astype(str)
    cosmic_clean['pos'] = cosmic_clean['GENOME_START']
    
    # Select relevant columns for analysis
    columns_to_keep = [
        'gene', 'chr', 'pos',
        'GENE_SYMBOL', 'CHROMOSOME', 'GENOME_START', 'GENOME_STOP',
        'MUTATION_CDS', 'MUTATION_AA', 'MUTATION_DESCRIPTION',
        'MUTATION_SOMATIC_STATUS', 'SAMPLE_NAME', 'COSMIC_PHENOTYPE_ID',
        'GENOMIC_MUTATION_ID', 'HGVSP', 'HGVSC', 'HGVSG',
        'GENOMIC_WT_ALLELE', 'GENOMIC_MUT_ALLELE'
    ]
    
    # Keep only existing columns
    existing_columns = [col for col in columns_to_keep if col in cosmic_clean.columns]
    cosmic_filtered = cosmic_clean[existing_columns]
    
    return cosmic_filtered

def main():
    """Main function to create cosmic_prostate.csv"""
    
    print("Loading COSMIC mutations data...")
    
    # Load the full COSMIC file
    cosmic_file = "../../data/raw/variants/Cosmic_MutantCensus_v102_GRCh38.tsv"
    
    try:
        # Read COSMIC file (large file, may take time)
        cosmic_df = pd.read_csv(cosmic_file, sep='\t', low_memory=False)
        print(f"Loaded COSMIC data: {cosmic_df.shape}")
        
    except FileNotFoundError:
        print(f"COSMIC file not found at: {cosmic_file}")
        print("Please place your COSMIC file in: /u/aa107/uiuc-cancer-research/data/raw/variants/")
        return
    
    # Filter for prostate cancer
    prostate_df = filter_prostate_samples(cosmic_df)
    
    if len(prostate_df) == 0:
        print("No prostate cancer samples found!")
        print("Sample names preview:")
        print(cosmic_df['SAMPLE_NAME'].head(20).tolist())
        return
    
    # Standardize columns
    prostate_clean = standardize_cosmic_columns(prostate_df)
    
    # Save filtered prostate data
    output_file = "../../data/raw/variants/cosmic_prostate.csv"
    prostate_clean.to_csv(output_file, index=False)
    
    print(f"\nProstate COSMIC data saved: {output_file}")
    print(f"Shape: {prostate_clean.shape}")
    print(f"Genes: {prostate_clean['gene'].nunique()}")
    print(f"Samples: {prostate_clean['SAMPLE_NAME'].nunique()}")
    
    # Show sample preview
    print("\nSample data preview:")
    print(prostate_clean[['gene', 'chr', 'pos', 'MUTATION_DESCRIPTION', 'SAMPLE_NAME']].head())
    
    # Show mutation types
    print("\nMutation types:")
    print(prostate_clean['MUTATION_DESCRIPTION'].value_counts().head())

if __name__ == "__main__":
    main()