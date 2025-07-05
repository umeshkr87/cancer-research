#!/usr/bin/env python3
"""
VEP VCF to TabNet CSV Converter
Converts VEP annotated VCF file to CSV format optimized for TabNet training

Location: /u/aa107/uiuc-cancer-research/scripts/enhance/vcf_to_tabnet_csv.py
Usage: python vcf_to_tabnet_csv.py
"""

import re
import gzip
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict
import argparse
from datetime import datetime

def open_vcf_file(vcf_path):
    """Open VCF file (handles both .vcf and .vcf.gz)"""
    if str(vcf_path).endswith('.gz'):
        return gzip.open(vcf_path, 'rt')
    else:
        return open(vcf_path, 'r')

def parse_csq_header(header_line):
    """Extract CSQ field names from VEP header"""
    match = re.search(r'Format: ([^"]+)', header_line)
    if match:
        return match.group(1).split('|')
    return []

def parse_variant_basic_info(fields):
    """Extract basic variant information from VCF fields"""
    chrom, pos, var_id, ref, alt, qual, filter_field, info = fields[:8]
    
    return {
        'chromosome': chrom,
        'position': int(pos),
        'variant_id': var_id,
        'reference_allele': ref,
        'alternate_allele': alt,
        'quality_score': float(qual) if qual != '.' else np.nan,
        'filter_status': filter_field,
        'ref_length': len(ref),
        'alt_length': len(alt),
        'variant_type': classify_variant_type(ref, alt)
    }

def classify_variant_type(ref, alt):
    """Classify variant type based on REF and ALT alleles"""
    if len(ref) == 1 and len(alt) == 1:
        return 'SNV'
    elif len(ref) > len(alt):
        return 'deletion'
    elif len(ref) < len(alt):
        return 'insertion'
    else:
        return 'complex'

def parse_csq_annotation(csq_string, csq_fields):
    """Parse CSQ annotation and return structured data"""
    if not csq_string or not csq_fields:
        return {}
    
    # Split by comma for multiple transcripts, take first (usually canonical)
    transcripts = csq_string.split(',')
    if not transcripts:
        return {}
    
    # Parse first transcript
    values = transcripts[0].split('|')
    annotation = {}
    
    for i, field in enumerate(csq_fields):
        if i < len(values):
            value = values[i].strip()
            annotation[field] = value if value and value != '.' else np.nan
        else:
            annotation[field] = np.nan
    
    return annotation

def engineer_features(basic_info, csq_annotation):
    """Engineer additional features for TabNet training"""
    features = {}
    
    # Variant size features
    features['variant_size'] = abs(basic_info['alt_length'] - basic_info['ref_length'])
    features['is_indel'] = 1 if basic_info['variant_type'] in ['insertion', 'deletion'] else 0
    features['is_snv'] = 1 if basic_info['variant_type'] == 'SNV' else 0
    
    # Impact encoding
    impact_map = {'HIGH': 4, 'MODERATE': 3, 'LOW': 2, 'MODIFIER': 1}
    features['impact_score'] = impact_map.get(csq_annotation.get('IMPACT'), 0)
    
    # Consequence severity (simplified)
    consequence = csq_annotation.get('Consequence', '')
    features['is_lof'] = 1 if any(term in consequence for term in ['stop_gained', 'frameshift', 'splice_donor', 'splice_acceptor']) else 0
    features['is_missense'] = 1 if 'missense_variant' in consequence else 0
    features['is_synonymous'] = 1 if 'synonymous_variant' in consequence else 0
    
    # Functional prediction parsing
    features = parse_functional_scores(csq_annotation, features)
    
    # Gene importance (based on our analysis)
    important_genes = {'BRCA1', 'BRCA2', 'ATM', 'TP53', 'AR', 'MSH2', 'MSH6', 'MLH1', 'PALB2', 'BRIP1'}
    features['is_important_gene'] = 1 if csq_annotation.get('SYMBOL') in important_genes else 0
    
    # Therapeutic pathway indicators
    dna_repair_genes = {'BRCA1', 'BRCA2', 'ATM', 'PALB2', 'BRIP1', 'RAD51C', 'RAD51D'}
    mismatch_repair_genes = {'MSH2', 'MSH6', 'MLH1', 'PMS2'}
    hormone_genes = {'AR', 'ESR1', 'ESR2'}
    
    features['dna_repair_pathway'] = 1 if csq_annotation.get('SYMBOL') in dna_repair_genes else 0
    features['mismatch_repair_pathway'] = 1 if csq_annotation.get('SYMBOL') in mismatch_repair_genes else 0
    features['hormone_pathway'] = 1 if csq_annotation.get('SYMBOL') in hormone_genes else 0
    
    return features

def parse_functional_scores(csq_annotation, features):
    """Parse and encode functional prediction scores"""
    
    # SIFT score processing
    sift_val = csq_annotation.get('SIFT', '')
    if pd.notna(sift_val) and sift_val:
        # SIFT can be "deleterious(0.01)" or just "0.01" or "deleterious"
        if '(' in sift_val and ')' in sift_val:
            # Extract score from parentheses
            score_match = re.search(r'\(([\d.]+)\)', sift_val)
            if score_match:
                features['sift_score'] = float(score_match.group(1))
            features['sift_prediction'] = 1 if 'deleterious' in sift_val.lower() else 0
        elif sift_val.replace('.', '').isdigit():
            # Just a number
            score = float(sift_val)
            features['sift_score'] = score
            features['sift_prediction'] = 1 if score < 0.05 else 0
        else:
            # Just prediction text
            features['sift_prediction'] = 1 if 'deleterious' in sift_val.lower() else 0
            features['sift_score'] = np.nan
    else:
        features['sift_score'] = np.nan
        features['sift_prediction'] = np.nan
    
    # PolyPhen score processing
    polyphen_val = csq_annotation.get('PolyPhen', '')
    if pd.notna(polyphen_val) and polyphen_val:
        # PolyPhen can be "probably_damaging(0.99)" or similar
        if '(' in polyphen_val and ')' in polyphen_val:
            score_match = re.search(r'\(([\d.]+)\)', polyphen_val)
            if score_match:
                features['polyphen_score'] = float(score_match.group(1))
            
            if 'probably_damaging' in polyphen_val.lower():
                features['polyphen_prediction'] = 2
            elif 'possibly_damaging' in polyphen_val.lower():
                features['polyphen_prediction'] = 1
            else:
                features['polyphen_prediction'] = 0
        else:
            features['polyphen_score'] = np.nan
            if 'probably_damaging' in polyphen_val.lower():
                features['polyphen_prediction'] = 2
            elif 'possibly_damaging' in polyphen_val.lower():
                features['polyphen_prediction'] = 1
            else:
                features['polyphen_prediction'] = 0
    else:
        features['polyphen_score'] = np.nan
        features['polyphen_prediction'] = np.nan
    
    # Allele frequency processing
    af_1kg = csq_annotation.get('AF', '')
    if pd.notna(af_1kg) and af_1kg:
        try:
            features['af_1kg'] = float(af_1kg)
            features['is_rare'] = 1 if float(af_1kg) < 0.01 else 0
        except:
            features['af_1kg'] = np.nan
            features['is_rare'] = np.nan
    else:
        features['af_1kg'] = np.nan
        features['is_rare'] = np.nan
    
    return features

def process_vcf_to_csv(vcf_path, output_path, max_variants=None):
    """Main function to convert VCF to CSV"""
    
    print(f"🔄 Converting VCF to TabNet CSV format...")
    print(f"📁 Input: {vcf_path}")
    print(f"📁 Output: {output_path}")
    print()
    
    # Initialize data storage
    all_data = []
    csq_fields = []
    total_variants = 0
    processed_variants = 0
    
    try:
        with open_vcf_file(vcf_path) as vcf:
            for line_num, line in enumerate(vcf, 1):
                line = line.strip()
                
                # Parse header for CSQ fields
                if line.startswith('##'):
                    if 'ID=CSQ' in line:
                        csq_fields = parse_csq_header(line)
                        print(f"📋 Found {len(csq_fields)} CSQ fields")
                
                elif line.startswith('#CHROM'):
                    print(f"🚀 Starting variant processing...")
                    print()
                
                # Process variants
                elif line and not line.startswith('#'):
                    total_variants += 1
                    
                    # Progress indicator
                    if total_variants % 25000 == 0:
                        print(f"📈 Processed {total_variants:,} variants...")
                    
                    # Parse VCF line
                    fields = line.split('\t')
                    if len(fields) >= 8:
                        
                        # Extract basic variant info
                        basic_info = parse_variant_basic_info(fields)
                        
                        # Extract CSQ annotation
                        info_field = fields[7]
                        csq_match = re.search(r'CSQ=([^;]+)', info_field)
                        
                        if csq_match:
                            csq_string = csq_match.group(1)
                            csq_annotation = parse_csq_annotation(csq_string, csq_fields)
                            
                            # Engineer features
                            engineered_features = engineer_features(basic_info, csq_annotation)
                            
                            # Combine all data
                            variant_data = {**basic_info, **csq_annotation, **engineered_features}
                            all_data.append(variant_data)
                            processed_variants += 1
                        
                        # Limit processing for testing
                        if max_variants and total_variants >= max_variants:
                            print(f"⚠️  Processing limited to {max_variants:,} variants for testing")
                            break
    
    except Exception as e:
        print(f"❌ Error processing VCF: {e}")
        return False
    
    if not all_data:
        print("❌ No variant data extracted!")
        return False
    
    # Convert to DataFrame
    print(f"\n📊 Converting {len(all_data):,} variants to DataFrame...")
    df = pd.DataFrame(all_data)
    
    # Data quality summary
    print(f"✅ Successfully processed {processed_variants:,}/{total_variants:,} variants")
    print(f"📊 DataFrame shape: {df.shape}")
    print(f"📋 Total features: {df.shape[1]}")
    
    # Show feature types
    print(f"\n🔍 Feature summary:")
    print(f"   Numeric features: {df.select_dtypes(include=[np.number]).shape[1]}")
    print(f"   Text features: {df.select_dtypes(include=[object]).shape[1]}")
    print(f"   Missing data: {df.isnull().sum().sum():,} cells")
    
    # Create output directory
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Save to CSV
    print(f"\n💾 Saving to CSV...")
    df.to_csv(output_path, index=False)
    
    # File size
    file_size_mb = output_path.stat().st_size / (1024 * 1024)
    print(f"✅ CSV saved successfully!")
    print(f"📁 File: {output_path}")
    print(f"📏 Size: {file_size_mb:.1f} MB")
    
    # Show sample of key features for validation
    print(f"\n🔬 Sample of key features:")
    key_cols = ['chromosome', 'position', 'SYMBOL', 'Consequence', 'IMPACT', 
                'sift_score', 'polyphen_score', 'impact_score', 'is_important_gene']
    available_cols = [col for col in key_cols if col in df.columns]
    
    if available_cols:
        sample_df = df[available_cols].head()
        print(sample_df.to_string(index=False))
    
    return True

def main():
    """Main function with command line arguments"""
    parser = argparse.ArgumentParser(description='Convert VEP VCF to TabNet CSV')
    parser.add_argument('--input', '-i', 
                       default='/u/aa107/uiuc-cancer-research/data/processed/vep/vep_annotated.vcf',
                       help='Input VCF file path')
    parser.add_argument('--output', '-o',
                       default='/u/aa107/uiuc-cancer-research/data/processed/tabnet_csv/prostate_variants_tabnet.csv',
                       help='Output CSV file path')
    parser.add_argument('--max-variants', '-m', type=int,
                       help='Maximum number of variants to process (for testing)')
    
    args = parser.parse_args()
    
    print("🧬 VEP VCF to TabNet CSV Converter")
    print("=" * 50)
    print(f"🎯 Purpose: Prepare data for interpretable prostate cancer classification")
    print(f"📅 Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()
    
    # File paths
    vcf_path = Path(args.input)
    output_path = Path(args.output)
    
    # Validate input
    if not vcf_path.exists():
        print(f"❌ Input VCF file not found: {vcf_path}")
        return False
    
    # Run conversion
    success = process_vcf_to_csv(vcf_path, output_path, args.max_variants)
    
    if success:
        print(f"\n🎉 Conversion completed successfully!")
        print(f"📊 CSV ready for TabNet training")
        print(f"🔄 Next steps:")
        print(f"   1. Load CSV into TabNet training pipeline")
        print(f"   2. Configure feature selection and engineering")
        print(f"   3. Train interpretable model with attention mechanisms")
    else:
        print(f"\n❌ Conversion failed!")
        print(f"🔍 Check error messages above for troubleshooting")
    
    return success

if __name__ == "__main__":
    main()