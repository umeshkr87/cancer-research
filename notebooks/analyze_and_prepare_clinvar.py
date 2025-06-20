import pandas as pd
import numpy as np
import os
import sys
from pathlib import Path
import subprocess
import re
from tqdm import tqdm
import logging

# Add the project root directory to the Python path
project_root = Path(__file__).parent.parent
sys.path.append(str(project_root))

def setup_logging():
    """Setup basic logging configuration"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    return logging.getLogger(__name__)

def extract_prostate_cancer_variants(vcf_path, output_dir):
    """
    Extract prostate cancer-related variants from ClinVar VCF file using bcftools.
    
    Args:
        vcf_path (str): Path to the ClinVar VCF file
        output_dir (str): Directory to save the extracted variants
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Define output file path
    output_file = os.path.join(output_dir, 'prostate_cancer_variants.vcf')
    
    # Use bcftools to filter variants
    cmd = f"bcftools view {vcf_path} | grep -i 'prostate.*cancer\\|prostate.*neoplasm\\|prostate.*tumor' > {output_file}"
    
    try:
        subprocess.run(cmd, shell=True, check=True)
        print(f"Successfully extracted prostate cancer variants to {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error extracting variants: {e}")
        return None
    
    return output_file

def analyze_prostate_cancer_variants(vcf_file):
    """
    Analyze prostate cancer variants from the filtered VCF file.
    
    Args:
        vcf_file (str): Path to the filtered VCF file containing prostate cancer variants
    """
    variants = []
    
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
                
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue
                
            chrom, pos, id_, ref, alt, qual, filter_, info = fields
            
            # Extract relevant information from INFO field
            info_dict = dict(item.split('=') for item in info.split(';') if '=' in item)
            
            # Extract clinical significance
            clinsig = info_dict.get('CLNSIG', 'Unknown')
            
            # Extract gene name if available
            gene = info_dict.get('GENEINFO', '').split(':')[0] if 'GENEINFO' in info_dict else 'Unknown'
            
            # Extract variant type
            variant_type = info_dict.get('CLNVC', 'Unknown')
            
            variants.append({
                'chromosome': chrom,
                'position': pos,
                'variant_id': id_,
                'reference': ref,
                'alternate': alt,
                'clinical_significance': clinsig,
                'gene': gene,
                'variant_type': variant_type
            })
    
    return pd.DataFrame(variants)

def prepare_for_merging(df):
    """
    Prepare ClinVar data for merging with TCGA-PRAD and COSMIC.
    Adds a unique variant_id and binary features for pathogenicity.
    """
    # Create a unique variant identifier
    df['variant_id'] = df.apply(
        lambda x: f"{x['chromosome']}:{x['position']}:{x['reference']}:{x['alternate']}", axis=1
    )

    # Create binary features for clinical significance
    pathogenic_terms = [
        'Pathogenic', 'Likely_pathogenic', 'Pathogenic/Likely_pathogenic'
    ]
    df['is_pathogenic'] = df['clinical_significance'].apply(
        lambda x: 1 if any(term in str(x) for term in pathogenic_terms) else 0
    )

    # Select and rename columns for merging
    merge_df = df[[
        'variant_id', 'chromosome', 'position', 'reference', 'alternate',
        'gene', 'is_pathogenic', 'clinical_significance', 'variant_type'
    ]].copy()
    return merge_df

def main():
    # Setup logging
    logger = setup_logging()
    logger.info("Starting ClinVar prostate cancer variant analysis")
    
    # Define paths
    vcf_path = os.path.join(project_root, 'data/raw/clinvar/vcf_GRCh38/clinvar.vcf.gz')
    output_dir = os.path.join(project_root, 'data/processed/clinvar/prostate_cancer')
    
    # Extract prostate cancer variants
    logger.info("Extracting prostate cancer variants from ClinVar VCF")
    filtered_vcf = extract_prostate_cancer_variants(vcf_path, output_dir)
    
    if filtered_vcf:
        # Analyze the variants
        logger.info("Analyzing prostate cancer variants")
        variants_df = analyze_prostate_cancer_variants(filtered_vcf)
        
        # Save the analysis results
        output_file = os.path.join(output_dir, 'prostate_cancer_variants_analysis.csv')
        variants_df.to_csv(output_file, index=False)
        logger.info(f"Analysis results saved to {output_file}")
        
        # Prepare for merging
        logger.info("Preparing ClinVar data for merging")
        merge_df = prepare_for_merging(variants_df)
        merge_file = os.path.join(output_dir, 'clinvar_prostate_ready_for_merge.csv')
        merge_df.to_csv(merge_file, index=False)
        logger.info(f"Ready-to-merge file saved to {merge_file}")
        
        # Print summary statistics
        print("\nProstate Cancer Variants Analysis Summary:")
        print(f"Total variants found: {len(variants_df)}")
        print("\nClinical Significance Distribution:")
        print(variants_df['clinical_significance'].value_counts())
        print("\nTop Genes with Variants:")
        print(variants_df['gene'].value_counts().head(10))
        print("\nVariant Type Distribution:")
        print(variants_df['variant_type'].value_counts())
        print(f"\nReady-to-merge file saved to: {merge_file}")
    else:
        logger.error("Failed to extract prostate cancer variants")

if __name__ == "__main__":
    main() 