#!/usr/bin/env python3
"""
Comprehensive ClinVar Prostate Cancer Variant Analyzer
Corrected version addressing filtering limitations and TabNet integration requirements
"""

import pandas as pd
import numpy as np
import os
import sys
import gzip
import re
from pathlib import Path
import subprocess
import logging
from typing import Dict, List, Tuple, Optional
from tqdm import tqdm
import argparse

def setup_logging(log_file: str = "clinvar_analysis.log") -> logging.Logger:
    """Setup comprehensive logging configuration"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger(__name__)

class ProstateCancerGenes:
    """Comprehensive prostate cancer gene collections for filtering"""
    
    # Core prostate cancer genes
    CORE_PROSTATE_GENES = {
        'AR', 'PTEN', 'TP53', 'MYC', 'ERG', 'ETV1', 'ETV4', 'ETV5',
        'TMPRSS2', 'SPINK1', 'CHD1', 'SPOP', 'FOXA1', 'IDH1'
    }
    
    # DNA repair pathway genes (PARP inhibitor targets)
    DNA_REPAIR_GENES = {
        'BRCA1', 'BRCA2', 'ATM', 'CHEK2', 'PALB2', 'RAD51D', 'BRIP1', 
        'FANCA', 'MLH1', 'MSH2', 'MSH6', 'PMS2', 'NBN', 'RAD51C', 
        'RAD51B', 'BARD1', 'FANCM', 'FANCI', 'FANCD2', 'RAD54L'
    }
    
    # Hormone pathway genes
    HORMONE_PATHWAY_GENES = {
        'CYP17A1', 'SRD5A2', 'CYP19A1', 'ESR1', 'ESR2', 'CYP11A1',
        'HSD17B3', 'HSD3B1', 'HSD3B2', 'STAR', 'CYP21A2'
    }
    
    # PI3K/AKT/mTOR pathway genes
    PI3K_PATHWAY_GENES = {
        'PIK3CA', 'PIK3R1', 'AKT1', 'AKT2', 'AKT3', 'MTOR', 'TSC1', 
        'TSC2', 'PTEN', 'STK11', 'RICTOR', 'RAPTOR'
    }
    
    # Additional prostate cancer associated genes
    ADDITIONAL_PROSTATE_GENES = {
        'HOXB13', 'KLK3', 'CDH1', 'RB1', 'NCOR1', 'NCOR2', 'NKX3-1',
        'GSTP1', 'CDKN1B', 'CDKN2A', 'MEN1', 'VHL', 'FLCN', 'FH'
    }
    
    # Cancer predisposition genes relevant to prostate cancer
    CANCER_PREDISPOSITION_GENES = {
        'BRCA1', 'BRCA2', 'TP53', 'CHEK2', 'ATM', 'PALB2', 'NBN',
        'MUTYH', 'MSH2', 'MLH1', 'MSH6', 'PMS2', 'EPCAM'
    }
    
    @classmethod
    def get_all_prostate_genes(cls) -> set:
        """Return comprehensive set of all prostate cancer related genes"""
        return (cls.CORE_PROSTATE_GENES | 
                cls.DNA_REPAIR_GENES | 
                cls.HORMONE_PATHWAY_GENES |
                cls.PI3K_PATHWAY_GENES |
                cls.ADDITIONAL_PROSTATE_GENES |
                cls.CANCER_PREDISPOSITION_GENES)
    
    @classmethod
    def get_gene_pathway_mapping(cls) -> Dict[str, str]:
        """Return gene to pathway mapping for feature engineering"""
        mapping = {}
        for gene in cls.CORE_PROSTATE_GENES:
            mapping[gene] = 'core_prostate'
        for gene in cls.DNA_REPAIR_GENES:
            mapping[gene] = 'dna_repair'
        for gene in cls.HORMONE_PATHWAY_GENES:
            mapping[gene] = 'hormone_pathway'
        for gene in cls.PI3K_PATHWAY_GENES:
            mapping[gene] = 'pi3k_pathway'
        for gene in cls.ADDITIONAL_PROSTATE_GENES:
            mapping[gene] = 'prostate_associated'
        for gene in cls.CANCER_PREDISPOSITION_GENES:
            mapping[gene] = 'cancer_predisposition'
        return mapping

class ClinVarVCFParser:
    """Robust ClinVar VCF parser with comprehensive INFO field extraction"""
    
    def __init__(self, logger: logging.Logger):
        self.logger = logger
        self.prostate_genes = ProstateCancerGenes.get_all_prostate_genes()
        
    def parse_info_field(self, info_str: str) -> Dict[str, str]:
        """Parse VCF INFO field into dictionary with robust error handling"""
        info_dict = {}
        
        # Handle edge cases in INFO parsing
        if not info_str or info_str == '.':
            return info_dict
            
        try:
            for item in info_str.split(';'):
                if '=' in item:
                    key, value = item.split('=', 1)
                    info_dict[key] = value
                else:
                    # Flag fields (no value)
                    info_dict[item] = 'True'
        except Exception as e:
            self.logger.warning(f"Error parsing INFO field: {info_str[:100]}... Error: {e}")
            
        return info_dict
    
    def extract_gene_symbols(self, geneinfo: str) -> List[str]:
        """Extract gene symbols from GENEINFO field with robust parsing"""
        if not geneinfo or geneinfo == '.' or geneinfo == '':
            return []
        
        genes = []
        try:
            # GENEINFO format: "GENE1:ID1|GENE2:ID2" or variations
            for gene_pair in geneinfo.split('|'):
                if ':' in gene_pair:
                    gene_symbol = gene_pair.split(':')[0].strip()
                    if gene_symbol and gene_symbol != '.':
                        genes.append(gene_symbol)
                elif gene_pair.strip() and gene_pair.strip() != '.':
                    # Handle cases where gene symbol is provided without ID
                    genes.append(gene_pair.strip())
        except Exception as e:
            self.logger.warning(f"Error parsing GENEINFO: {geneinfo}. Error: {e}")
            
        return genes
    
    def extract_molecular_consequences(self, mc_field: str) -> List[str]:
        """Extract molecular consequences from MC field"""
        if not mc_field or mc_field == '.' or mc_field == '':
            return []
        
        consequences = []
        try:
            # MC format: "SO:ID|consequence_term,SO:ID2|consequence_term2"
            for item in mc_field.split(','):
                if '|' in item:
                    consequence = item.split('|')[1].strip()
                    if consequence and consequence != '.':
                        consequences.append(consequence)
        except Exception as e:
            self.logger.warning(f"Error parsing MC field: {mc_field}. Error: {e}")
            
        return consequences
    
    def clean_clinical_significance(self, clnsig: str) -> str:
        """Clean and standardize clinical significance values"""
        if not clnsig or clnsig == '.':
            return 'Unknown'
        
        # Replace underscores and normalize
        clean_sig = clnsig.replace('_', ' ').strip()
        
        # Standardize common variants
        standardization_map = {
            'Pathogenic/Likely pathogenic': 'Pathogenic',
            'Likely pathogenic/Pathogenic': 'Pathogenic',
            'Benign/Likely benign': 'Benign',
            'Likely benign/Benign': 'Benign',
            'Uncertain significance': 'Uncertain_significance',
            'not provided': 'Not_provided',
            'drug response': 'Drug_response',
            'risk factor': 'Risk_factor'
        }
        
        return standardization_map.get(clean_sig, clean_sig)
    
    def is_prostate_relevant_by_gene(self, gene_list: List[str]) -> bool:
        """Check if variant is prostate-relevant based on gene membership"""
        return any(gene in self.prostate_genes for gene in gene_list)
    
    def is_prostate_relevant_by_disease(self, disease_name: str, clinical_sig: str) -> bool:
        """Check if variant is prostate-relevant based on disease annotations"""
        if not disease_name:
            return False
            
        disease_lower = disease_name.lower()
        
        # Primary prostate terms
        prostate_terms = [
            'prostate', 'prostatic', 'prostate cancer', 'prostate carcinoma',
            'prostate adenocarcinoma', 'prostate neoplasm', 'prostate tumor',
            'benign prostatic hyperplasia', 'bph', 'androgen', 'psa'
        ]
        
        # Check for prostate-specific disease terms
        if any(term in disease_lower for term in prostate_terms):
            return True
        
        # Check for male-specific cancers that may be relevant
        male_cancer_terms = ['male breast cancer', 'hereditary cancer syndrome']
        if any(term in disease_lower for term in male_cancer_terms):
            return True
            
        return False
    
    def determine_pathogenicity_score(self, clinical_sig: str) -> int:
        """Convert clinical significance to binary pathogenicity score"""
        if not clinical_sig:
            return 0
            
        sig_lower = clinical_sig.lower()
        
        # Clearly pathogenic
        if any(term in sig_lower for term in ['pathogenic', 'likely pathogenic']):
            return 1
        
        # Clearly benign
        if any(term in sig_lower for term in ['benign', 'likely benign']):
            return 0
        
        # Special cases that may be clinically relevant
        if any(term in sig_lower for term in ['drug response', 'risk factor']):
            return 1
        
        # Everything else (VUS, etc.)
        return 0

class ClinVarProstateCancerAnalyzer:
    """Main analyzer class for extracting prostate cancer variants from ClinVar"""
    
    def __init__(self, output_dir: str, logger: logging.Logger):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.logger = logger
        self.parser = ClinVarVCFParser(logger)
        
    def process_vcf_file(self, vcf_path: str, max_variants: Optional[int] = None) -> pd.DataFrame:
        """Process ClinVar VCF file with comprehensive filtering"""
        
        self.logger.info(f"Processing ClinVar VCF: {vcf_path}")
        
        variants = []
        total_processed = 0
        prostate_relevant_count = 0
        
        # Determine if file is gzipped
        open_func = gzip.open if vcf_path.endswith('.gz') else open
        
        try:
            with open_func(vcf_path, 'rt') as f:
                for line_num, line in enumerate(tqdm(f, desc="Processing VCF")):
                    # Skip header lines
                    if line.startswith('#'):
                        continue
                    
                    total_processed += 1
                    
                    # Optional limit for testing
                    if max_variants and total_processed > max_variants:
                        break
                    
                    # Parse VCF line
                    variant_data = self._parse_vcf_line(line, line_num)
                    if not variant_data:
                        continue
                    
                    # Check if prostate-relevant
                    if self._is_variant_prostate_relevant(variant_data):
                        variants.append(variant_data)
                        prostate_relevant_count += 1
                    
                    # Progress reporting
                    if total_processed % 100000 == 0:
                        self.logger.info(f"Processed {total_processed:,} variants, "
                                       f"found {prostate_relevant_count:,} prostate-relevant")
        
        except Exception as e:
            self.logger.error(f"Error processing VCF file: {e}")
            return pd.DataFrame()
        
        self.logger.info(f"Final results: {total_processed:,} total variants processed, "
                        f"{prostate_relevant_count:,} prostate-relevant variants found")
        
        return pd.DataFrame(variants)
    
    def _parse_vcf_line(self, line: str, line_num: int) -> Optional[Dict]:
        """Parse a single VCF line into structured data"""
        try:
            fields = line.strip().split('\t')
            if len(fields) < 8:
                return None
            
            chrom, pos, var_id, ref, alt, qual, filter_val, info_str = fields[:8]
            
            # Parse INFO field
            info_dict = self.parser.parse_info_field(info_str)
            
            # Extract gene information
            gene_list = self.parser.extract_gene_symbols(info_dict.get('GENEINFO', ''))
            
            # Extract clinical information
            clinical_sig = self.parser.clean_clinical_significance(info_dict.get('CLNSIG', ''))
            disease_name = info_dict.get('CLNDN', '').replace('_', ' ')
            
            # Extract additional fields
            molecular_consequences = self.parser.extract_molecular_consequences(info_dict.get('MC', ''))
            review_status = info_dict.get('CLNREVSTAT', '').replace('_', ' ')
            variant_type = info_dict.get('CLNVC', '')
            rs_id = info_dict.get('RS', '')
            allele_id = info_dict.get('ALLELEID', '')
            
            # Create structured variant record
            variant_record = {
                'variant_id': f"{chrom}:{pos}:{ref}:{alt}",
                'chromosome': chrom.replace('chr', ''),  # Standardize chromosome naming
                'position': int(pos),
                'reference': ref,
                'alternate': alt,
                'gene': '|'.join(gene_list) if gene_list else '',
                'gene_list': gene_list,  # Keep for internal processing
                'clinical_significance': clinical_sig,
                'disease_name': disease_name,
                'molecular_consequence': '|'.join(molecular_consequences),
                'review_status': review_status,
                'variant_type': variant_type,
                'rs_id': rs_id,
                'allele_id': allele_id,
                'is_pathogenic': self.parser.determine_pathogenicity_score(clinical_sig)
            }
            
            return variant_record
            
        except Exception as e:
            self.logger.warning(f"Error parsing line {line_num}: {e}")
            return None
    
    def _is_variant_prostate_relevant(self, variant_data: Dict) -> bool:
        """Determine if variant is relevant to prostate cancer research"""
        
        # Gene-based filtering (primary filter)
        if self.parser.is_prostate_relevant_by_gene(variant_data['gene_list']):
            return True
        
        # Disease-based filtering (secondary filter)
        if self.parser.is_prostate_relevant_by_disease(
            variant_data['disease_name'], 
            variant_data['clinical_significance']
        ):
            return True
        
        return False
    
    def add_pathway_annotations(self, df: pd.DataFrame) -> pd.DataFrame:
        """Add pathway annotations for TabNet feature engineering"""
        
        gene_pathway_mapping = ProstateCancerGenes.get_gene_pathway_mapping()
        
        # Initialize pathway columns
        pathway_columns = ['core_prostate', 'dna_repair', 'hormone_pathway', 
                          'pi3k_pathway', 'prostate_associated', 'cancer_predisposition']
        
        for pathway in pathway_columns:
            df[f'{pathway}_gene'] = 0
        
        # Annotate pathways for each variant
        for idx, row in df.iterrows():
            genes = row['gene'].split('|') if row['gene'] else []
            
            for gene in genes:
                if gene in gene_pathway_mapping:
                    pathway = gene_pathway_mapping[gene]
                    if f'{pathway}_gene' in df.columns:
                        df.at[idx, f'{pathway}_gene'] = 1
        
        return df
    
    def generate_summary_statistics(self, df: pd.DataFrame) -> Dict:
        """Generate comprehensive summary statistics"""
        
        stats = {
            'total_variants': len(df),
            'unique_genes': len(set([gene for gene_str in df['gene'] 
                                   for gene in gene_str.split('|') if gene])),
            'chromosomes': sorted(df['chromosome'].unique()),
            'pathogenic_variants': df['is_pathogenic'].sum(),
            'pathogenic_percentage': (df['is_pathogenic'].sum() / len(df)) * 100 if len(df) > 0 else 0
        }
        
        # Clinical significance distribution
        stats['clinical_significance_dist'] = df['clinical_significance'].value_counts().to_dict()
        
        # Top genes
        gene_counts = {}
        for gene_str in df['gene']:
            for gene in gene_str.split('|'):
                if gene and gene in ProstateCancerGenes.get_all_prostate_genes():
                    gene_counts[gene] = gene_counts.get(gene, 0) + 1
        
        stats['top_genes'] = dict(sorted(gene_counts.items(), 
                                       key=lambda x: x[1], reverse=True)[:20])
        
        # Pathway distribution
        pathway_columns = [col for col in df.columns if col.endswith('_gene')]
        stats['pathway_distribution'] = {col: df[col].sum() for col in pathway_columns}
        
        return stats
    
    def save_results(self, df: pd.DataFrame, base_filename: str = "clinvar_prostate") -> Dict[str, str]:
        """Save results in multiple formats for downstream analysis"""
        
        output_files = {}
        
        # Clean dataset for TabNet
        tabnet_columns = [
            'variant_id', 'chromosome', 'position', 'reference', 'alternate',
            'gene', 'is_pathogenic', 'clinical_significance', 'variant_type'
        ]
        
        # Add pathway columns
        pathway_columns = [col for col in df.columns if col.endswith('_gene')]
        tabnet_columns.extend(pathway_columns)
        
        tabnet_df = df[tabnet_columns].copy()
        
        # Save main CSV for TabNet integration
        main_csv = self.output_dir / f"{base_filename}.csv"
        tabnet_df.to_csv(main_csv, index=False)
        output_files['main_csv'] = str(main_csv)
        
        # Save comprehensive dataset with all annotations
        comprehensive_csv = self.output_dir / f"{base_filename}_comprehensive.csv"
        df.to_csv(comprehensive_csv, index=False)
        output_files['comprehensive_csv'] = str(comprehensive_csv)
        
        # Save summary statistics
        stats = self.generate_summary_statistics(df)
        stats_file = self.output_dir / f"{base_filename}_summary.txt"
        
        with open(stats_file, 'w') as f:
            f.write("ClinVar Prostate Cancer Variant Analysis Summary\n")
            f.write("=" * 50 + "\n\n")
            
            f.write(f"Total variants: {stats['total_variants']:,}\n")
            f.write(f"Unique genes: {stats['unique_genes']}\n")
            f.write(f"Pathogenic variants: {stats['pathogenic_variants']:,} "
                   f"({stats['pathogenic_percentage']:.1f}%)\n")
            f.write(f"Chromosomes: {', '.join(stats['chromosomes'])}\n\n")
            
            f.write("Clinical Significance Distribution:\n")
            for sig, count in stats['clinical_significance_dist'].items():
                f.write(f"  {sig}: {count:,}\n")
            
            f.write("\nTop 20 Prostate Cancer Genes:\n")
            for gene, count in stats['top_genes'].items():
                f.write(f"  {gene}: {count:,}\n")
            
            f.write("\nPathway Distribution:\n")
            for pathway, count in stats['pathway_distribution'].items():
                f.write(f"  {pathway}: {count:,}\n")
        
        output_files['summary'] = str(stats_file)
        
        return output_files

def main():
    """Main function with command-line interface"""
    
    parser = argparse.ArgumentParser(
        description="Extract and analyze prostate cancer variants from ClinVar VCF"
    )
    parser.add_argument(
        "--vcf-path", 
        required=True,
        help="Path to ClinVar VCF file (can be gzipped)"
    )
    parser.add_argument(
        "--output-dir",
        default="./clinvar_prostate_analysis",
        help="Output directory for results"
    )
    parser.add_argument(
        "--max-variants",
        type=int,
        help="Maximum number of variants to process (for testing)"
    )
    parser.add_argument(
        "--log-file",
        default="clinvar_analysis.log",
        help="Log file path"
    )
    
    args = parser.parse_args()
    
    # Setup logging
    logger = setup_logging(args.log_file)
    logger.info("Starting ClinVar Prostate Cancer Variant Analysis")
    
    # Validate input file
    if not os.path.exists(args.vcf_path):
        logger.error(f"VCF file not found: {args.vcf_path}")
        sys.exit(1)
    
    # Initialize analyzer
    analyzer = ClinVarProstateCancerAnalyzer(args.output_dir, logger)
    
    try:
        # Process VCF file
        logger.info(f"Processing VCF: {args.vcf_path}")
        variants_df = analyzer.process_vcf_file(args.vcf_path, args.max_variants)
        
        if variants_df.empty:
            logger.error("No prostate-relevant variants found!")
            sys.exit(1)
        
        # Add pathway annotations
        logger.info("Adding pathway annotations")
        variants_df = analyzer.add_pathway_annotations(variants_df)
        
        # Save results
        logger.info("Saving results")
        output_files = analyzer.save_results(variants_df)
        
        # Print summary
        logger.info("Analysis completed successfully!")
        logger.info(f"Results saved to:")
        for file_type, filepath in output_files.items():
            logger.info(f"  {file_type}: {filepath}")
        
        # Print quick stats
        stats = analyzer.generate_summary_statistics(variants_df)
        print(f"\n🎯 Analysis Complete!")
        print(f"📊 {stats['total_variants']:,} prostate-relevant variants found")
        print(f"🧬 {stats['unique_genes']} unique genes")
        print(f"⚡ {stats['pathogenic_variants']:,} pathogenic variants ({stats['pathogenic_percentage']:.1f}%)")
        print(f"📁 Main output: {output_files['main_csv']}")
        print(f"📋 Ready for TabNet integration!")
        
    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()