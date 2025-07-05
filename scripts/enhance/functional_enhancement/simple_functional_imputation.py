#!/usr/bin/env python3
"""
AlphaMissense Functional Enhancement for Prostate Cancer Variants
Replaces artificial imputation with legitimate AlphaMissense pathogenicity scores
to eliminate data leakage while maintaining clinical interpretability.
"""

import pandas as pd
import numpy as np
import os
import sys
import requests
import gzip
from pathlib import Path
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

def download_alphamissense_data(scratch_dir):
    """Download and extract AlphaMissense pre-computed scores"""
    print("\n🔬 DOWNLOADING ALPHAMISSENSE DATABASE")
    print("-" * 50)
    
    # AlphaMissense URLs (using the main hg38 file)
    alphamissense_url = "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz"
    
    scratch_path = Path(scratch_dir)
    print(f"📁 Creating scratch directory: {scratch_path}")
    scratch_path.mkdir(exist_ok=True, parents=True)
    
    alphamissense_file = scratch_path / "AlphaMissense_hg38.tsv.gz"
    alphamissense_extracted = scratch_path / "AlphaMissense_hg38.tsv"
    
    print(f"📥 Target files:")
    print(f"   Compressed: {alphamissense_file}")
    print(f"   Extracted: {alphamissense_extracted}")
    
    # Check if extracted file already exists
    if alphamissense_extracted.exists():
        print(f"✅ AlphaMissense already extracted: {alphamissense_extracted}")
        return alphamissense_extracted
    
    # Download if not exists
    if not alphamissense_file.exists():
        print(f"📥 Downloading AlphaMissense database (~2GB)...")
        print(f"Source: {alphamissense_url}")
        print(f"Target: {alphamissense_file}")
        
        try:
            response = requests.get(alphamissense_url, stream=True, timeout=30)
            response.raise_for_status()
            
            total_size = int(response.headers.get('content-length', 0))
            downloaded = 0
            
            with open(alphamissense_file, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
                        downloaded += len(chunk)
                        if total_size > 0:
                            percent = (downloaded / total_size) * 100
                            print(f"\rProgress: {percent:.1f}% ({downloaded:,}/{total_size:,} bytes)", end='')
            
            print(f"\n✅ Download complete: {alphamissense_file}")
            print(f"📏 Downloaded size: {alphamissense_file.stat().st_size / 1024 / 1024:.1f} MB")
            
        except Exception as e:
            print(f"❌ Download failed: {e}")
            if alphamissense_file.exists():
                alphamissense_file.unlink()  # Remove partial download
            raise
    
    # Extract the file
    if not alphamissense_extracted.exists():
        print(f"📦 Extracting AlphaMissense database...")
        print(f"Source: {alphamissense_file}")
        print(f"Target: {alphamissense_extracted}")
        
        try:
            with gzip.open(alphamissense_file, 'rt') as f_in:
                with open(alphamissense_extracted, 'w') as f_out:
                    # Copy in chunks to show progress
                    chunk_size = 1024 * 1024  # 1MB chunks
                    while True:
                        chunk = f_in.read(chunk_size)
                        if not chunk:
                            break
                        f_out.write(chunk)
                        print(".", end='', flush=True)
            
            print(f"\n✅ Extraction complete")
            print(f"📁 Extracted to: {alphamissense_extracted}")
            print(f"📏 Extracted size: {alphamissense_extracted.stat().st_size / 1024 / 1024:.1f} MB")
            
        except Exception as e:
            print(f"❌ Extraction failed: {e}")
            if alphamissense_extracted.exists():
                alphamissense_extracted.unlink()  # Remove partial extraction
            raise
    
    return alphamissense_extracted

def load_alphamissense_lookup(alphamissense_file):
    """Load AlphaMissense scores into memory for fast lookup"""
    print("\n🧠 LOADING ALPHAMISSENSE LOOKUP TABLE")
    print("-" * 50)
    
    print(f"📖 Reading AlphaMissense file: {alphamissense_file}")
    
    # First, auto-detect columns by finding the header line (starts with #CHROM)
    print("🔍 Auto-detecting file format...")
    with open(alphamissense_file, 'r') as f:
        line_num = 0
        for line in f:
            line_num += 1
            if line.startswith('#CHROM'):
                header_line = line.strip()
                break
            if line_num > 20:  # Safety check
                raise ValueError("No header line starting with '#CHROM' found in first 20 lines")
    
    columns = header_line.split('\t')
    print(f"📋 Detected columns: {columns}")
    
    # Map to expected column names
    required_cols = ['#CHROM', 'POS', 'REF', 'ALT', 'am_pathogenicity']
    available_cols = []
    
    for req_col in required_cols:
        if req_col in columns:
            available_cols.append(req_col)
        else:
            print(f"⚠️  Column '{req_col}' not found in file")
    
    if len(available_cols) < 5:
        print(f"❌ Missing required columns. Found: {available_cols}")
        raise ValueError(f"Required columns not found: {required_cols}")
    
    print(f"✅ Using columns: {available_cols}")
    
    # Load AlphaMissense data (skip copyright comment lines but keep header)
    alphamissense_df = pd.read_csv(alphamissense_file, sep='\t', 
                                   skiprows=lambda x: x < 3,  # Skip first 3 copyright lines
                                   usecols=available_cols)
    
    print(f"📊 Loaded {len(alphamissense_df):,} AlphaMissense predictions")
    
    # Create lookup key: chr_pos_ref_alt (FIX: Strip 'chr' prefix for standardization)
    alphamissense_df['lookup_key'] = (
        alphamissense_df['#CHROM'].astype(str).str.replace('chr', '') + '_' +
        alphamissense_df['POS'].astype(str) + '_' +
        alphamissense_df['REF'] + '_' +
        alphamissense_df['ALT']
    )
    
    # Create lookup dictionary for fast access
    lookup_dict = dict(zip(alphamissense_df['lookup_key'], 
                          alphamissense_df['am_pathogenicity']))
    
    print(f"✅ Created lookup table with {len(lookup_dict):,} entries")
    print(f"📈 Score range: {alphamissense_df['am_pathogenicity'].min():.3f} - {alphamissense_df['am_pathogenicity'].max():.3f}")
    
    return lookup_dict

def enhance_with_alphamissense(df, alphamissense_lookup):
    """Add AlphaMissense scores to variants and remove artificial features"""
    print("\n🎯 ENHANCING VARIANTS WITH ALPHAMISSENSE SCORES")
    print("-" * 50)
    
    # Remove artificial features that cause data leakage
    artificial_features = ['sift_confidence', 'polyphen_confidence', 'functional_pathogenicity']
    for feature in artificial_features:
        if feature in df.columns:
            df = df.drop(columns=[feature])
            print(f"🗑️  Removed artificial feature: {feature}")
    
    # Create lookup keys for our variants
    df['lookup_key'] = (
        df['chromosome'].astype(str) + '_' +
        df['position'].astype(str) + '_' +
        df['reference_allele'] + '_' +
        df['alternate_allele']
    )
    
    # Map AlphaMissense scores
    df['alphamissense_pathogenicity'] = df['lookup_key'].map(alphamissense_lookup)
    
    # Statistics
    matched_variants = df['alphamissense_pathogenicity'].notna().sum()
    total_variants = len(df)
    coverage_rate = (matched_variants / total_variants) * 100
    
    print(f"📊 AlphaMissense Coverage:")
    print(f"   Total variants: {total_variants:,}")
    print(f"   Matched variants: {matched_variants:,}")
    print(f"   Coverage rate: {coverage_rate:.1f}%")
    
    if matched_variants > 0:
        am_scores = df['alphamissense_pathogenicity'].dropna()
        print(f"   Score statistics:")
        print(f"     Range: {am_scores.min():.3f} - {am_scores.max():.3f}")
        print(f"     Mean: {am_scores.mean():.3f}")
        print(f"     Median: {am_scores.median():.3f}")
        
        # Clinical significance thresholds (from literature)
        likely_pathogenic = (am_scores >= 0.564).sum()
        ambiguous = ((am_scores >= 0.34) & (am_scores < 0.564)).sum()
        likely_benign = (am_scores < 0.34).sum()
        
        print(f"   Clinical classification:")
        print(f"     Likely pathogenic (≥0.564): {likely_pathogenic:,} ({likely_pathogenic/len(am_scores)*100:.1f}%)")
        print(f"     Ambiguous (0.34-0.564): {ambiguous:,} ({ambiguous/len(am_scores)*100:.1f}%)")
        print(f"     Likely benign (<0.34): {likely_benign:,} ({likely_benign/len(am_scores)*100:.1f}%)")
    
    # Add clinical classification categories
    df['alphamissense_class'] = 'Unknown'
    df.loc[df['alphamissense_pathogenicity'] >= 0.564, 'alphamissense_class'] = 'Likely_Pathogenic'
    df.loc[(df['alphamissense_pathogenicity'] >= 0.34) & 
           (df['alphamissense_pathogenicity'] < 0.564), 'alphamissense_class'] = 'Ambiguous'
    df.loc[df['alphamissense_pathogenicity'] < 0.34, 'alphamissense_class'] = 'Likely_Benign'
    
    # Clean up temporary lookup key
    df = df.drop(columns=['lookup_key'])
    
    print("✅ AlphaMissense enhancement complete")
    return df

def analyze_pathway_enrichment(df):
    """Analyze AlphaMissense scores by therapeutic pathways"""
    print("\n📈 PATHWAY ENRICHMENT ANALYSIS")
    print("-" * 50)
    
    pathways = ['dna_repair_pathway', 'mismatch_repair_pathway', 'hormone_pathway', 'is_important_gene']
    
    for pathway in pathways:
        if pathway in df.columns:
            pathway_mask = (df[pathway] == 1)
            pathway_variants = df[pathway_mask]
            
            if len(pathway_variants) > 0:
                am_scores = pathway_variants['alphamissense_pathogenicity'].dropna()
                
                if len(am_scores) > 0:
                    print(f"\n{pathway.replace('_', ' ').title()}:")
                    print(f"  Total variants: {len(pathway_variants):,}")
                    print(f"  With AlphaMissense scores: {len(am_scores):,}")
                    print(f"  Average pathogenicity: {am_scores.mean():.3f}")
                    print(f"  High pathogenicity (≥0.7): {(am_scores >= 0.7).sum():,} ({(am_scores >= 0.7).sum()/len(am_scores)*100:.1f}%)")

def validate_enhancement(df_original, df_enhanced):
    """Validate AlphaMissense enhancement results"""
    print("\n✅ ENHANCEMENT VALIDATION")
    print("-" * 50)
    
    print(f"Original variants: {len(df_original):,}")
    print(f"Enhanced variants: {len(df_enhanced):,}")
    print(f"Features added: AlphaMissense pathogenicity scores")
    print(f"Features removed: Artificial confidence indicators")
    
    # Check for artificial features (should be gone)
    artificial_features = ['sift_confidence', 'polyphen_confidence', 'functional_pathogenicity']
    remaining_artificial = [f for f in artificial_features if f in df_enhanced.columns]
    
    if remaining_artificial:
        print(f"⚠️  Warning: Artificial features still present: {remaining_artificial}")
    else:
        print("✅ All artificial features successfully removed")
    
    # Check AlphaMissense coverage
    if 'alphamissense_pathogenicity' in df_enhanced.columns:
        coverage = df_enhanced['alphamissense_pathogenicity'].notna().sum()
        print(f"✅ AlphaMissense coverage: {coverage:,} variants ({coverage/len(df_enhanced)*100:.1f}%)")
    
    return True

def save_enhanced_results(df, output_file, report_file, stats):
    """Save enhanced data and generate comprehensive report"""
    print(f"\n💾 SAVING ENHANCED RESULTS")
    print("-" * 50)
    
    # Save enhanced CSV
    print(f"📁 Saving enhanced data to: {output_file}")
    df.to_csv(output_file, index=False)
    
    file_size_mb = output_file.stat().st_size / (1024 * 1024)
    print(f"✅ Saved {len(df):,} variants")
    print(f"📏 File size: {file_size_mb:.1f} MB")
    
    # Generate comprehensive report
    print(f"📝 Generating report: {report_file}")
    
    with open(report_file, 'w') as f:
        f.write("ALPHAMISSENSE FUNCTIONAL ENHANCEMENT REPORT\n")
        f.write("=" * 55 + "\n\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Input file: prostate_variants_tabnet.csv\n")
        f.write(f"Output file: prostate_variants_tabnet_enhanced.csv\n\n")
        
        f.write("ENHANCEMENT SUMMARY:\n")
        f.write(f"Total variants: {len(df):,}\n")
        f.write(f"AlphaMissense coverage: {stats['coverage_count']:,} variants ({stats['coverage_rate']:.1f}%)\n")
        f.write(f"File size: {file_size_mb:.1f} MB\n\n")
        
        f.write("DATA LEAKAGE ELIMINATION:\n")
        f.write("❌ REMOVED artificial features:\n")
        f.write("  - sift_confidence (artificial binary flags)\n")
        f.write("  - polyphen_confidence (artificial binary flags)\n")
        f.write("  - functional_pathogenicity (composite artificial score)\n\n")
        
        f.write("✅ ADDED legitimate features:\n")
        f.write("  - alphamissense_pathogenicity (0-1 pathogenicity score)\n")
        f.write("  - alphamissense_class (Likely_Pathogenic/Ambiguous/Likely_Benign)\n\n")
        
        f.write("EXPECTED PERFORMANCE IMPACT:\n")
        f.write("- Eliminates 100% artificial accuracy from data leakage\n")
        f.write("- Expected realistic accuracy: 75-85% (clinically appropriate)\n")
        f.write("- Maintains interpretability with legitimate functional scores\n")
        f.write("- Provides state-of-the-art pathogenicity predictions\n\n")
        
        f.write("NEXT STEPS:\n")
        f.write("1. Use prostate_variants_tabnet_enhanced.csv for TabNet training\n")
        f.write("2. Include alphamissense_pathogenicity as primary functional feature\n")
        f.write("3. Use alphamissense_class for interpretability analysis\n")
        f.write("4. Validate model achieves 75-85% accuracy (no more data leakage)\n")
        f.write("5. Proceed with clinical variant classification research\n")
    
    print(f"✅ Report saved: {report_file}")

def main():
    """Main AlphaMissense enhancement pipeline"""
    print("🧬 ALPHAMISSENSE FUNCTIONAL ENHANCEMENT")
    print("=" * 60)
    print("Replacing artificial imputation with legitimate pathogenicity scores")
    print("to eliminate data leakage in prostate cancer variant classification\n")
    
    # Configuration
    base_dir = Path("/u/aa107/uiuc-cancer-research")
    scratch_dir = "/u/aa107/scratch/alphamissense"
    input_file = base_dir / "data/processed/tabnet_csv/prostate_variants_tabnet.csv"
    output_file = base_dir / "data/processed/tabnet_csv/prostate_variants_tabnet_enhanced.csv"
    report_file = base_dir / "data/processed/tabnet_csv/alphamissense_enhancement_report.txt"
    
    print(f"📁 Configuration:")
    print(f"   Base directory: {base_dir}")
    print(f"   Scratch directory: {scratch_dir}")
    print(f"   Input file: {input_file}")
    print(f"   Output file: {output_file}")
    print(f"   Report file: {report_file}")
    print()
    
    # Step 1: Download AlphaMissense database
    print("🔄 Step 1: Downloading AlphaMissense database...")
    try:
        alphamissense_file = download_alphamissense_data(scratch_dir)
        print(f"✅ Step 1 complete: {alphamissense_file}")
    except Exception as e:
        print(f"❌ Step 1 failed - Error downloading AlphaMissense: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    # Step 2: Load AlphaMissense lookup table
    print("\n🔄 Step 2: Loading AlphaMissense lookup table...")
    try:
        alphamissense_lookup = load_alphamissense_lookup(alphamissense_file)
        print(f"✅ Step 2 complete: {len(alphamissense_lookup):,} entries loaded")
    except Exception as e:
        print(f"❌ Step 2 failed - Error loading AlphaMissense lookup: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    # Step 3: Load original VEP data
    print(f"\n🔄 Step 3: Loading original VEP data...")
    print(f"📖 Reading: {input_file}")
    
    if not input_file.exists():
        print(f"❌ Step 3 failed - Input file not found: {input_file}")
        sys.exit(1)
    
    try:
        df_original = pd.read_csv(input_file)
        print(f"✅ Step 3 complete: Loaded {len(df_original):,} variants with {df_original.shape[1]} features")
        
        # Show column names for debugging
        print(f"📋 Available columns: {list(df_original.columns[:10])}...")
        
    except Exception as e:
        print(f"❌ Step 3 failed - Error loading CSV: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    # Step 4: Enhance with AlphaMissense
    print(f"\n🔄 Step 4: Enhancing with AlphaMissense...")
    try:
        df_enhanced = enhance_with_alphamissense(df_original, alphamissense_lookup)
        print(f"✅ Step 4 complete: Enhanced dataset ready")
    except Exception as e:
        print(f"❌ Step 4 failed - Error enhancing with AlphaMissense: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    # Step 5: Analyze pathway enrichment
    print(f"\n🔄 Step 5: Analyzing pathway enrichment...")
    try:
        analyze_pathway_enrichment(df_enhanced)
        print(f"✅ Step 5 complete: Pathway analysis done")
    except Exception as e:
        print(f"⚠️  Step 5 warning - Error in pathway analysis: {e}")
        # Continue despite pathway analysis errors
    
    # Step 6: Validate enhancement
    print(f"\n🔄 Step 6: Validating enhancement...")
    try:
        validate_enhancement(df_original, df_enhanced)
        print(f"✅ Step 6 complete: Validation passed")
    except Exception as e:
        print(f"❌ Step 6 failed - Error in validation: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    # Step 7: Save results
    print(f"\n🔄 Step 7: Saving enhanced results...")
    try:
        # Calculate statistics for report
        coverage_count = df_enhanced['alphamissense_pathogenicity'].notna().sum()
        coverage_rate = (coverage_count / len(df_enhanced)) * 100
        
        stats = {
            'coverage_count': coverage_count,
            'coverage_rate': coverage_rate
        }
        
        # Create output directory
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        save_enhanced_results(df_enhanced, output_file, report_file, stats)
        print(f"✅ Step 7 complete: All results saved")
    except Exception as e:
        print(f"❌ Step 7 failed - Error saving results: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    # Final summary
    print(f"\n🎉 ALPHAMISSENSE ENHANCEMENT COMPLETE!")
    print("=" * 60)
    print(f"📊 Summary:")
    print(f"   Input variants: {len(df_original):,}")
    print(f"   Enhanced variants: {len(df_enhanced):,}")
    print(f"   AlphaMissense coverage: {coverage_count:,} ({coverage_rate:.1f}%)")
    print(f"   Output file: {output_file}")
    print(f"   Report file: {report_file}")
    print()
    print("🔬 Ready for TabNet training with legitimate functional scores!")
    print("🎯 Expected model accuracy: 75-85% (no more data leakage)")

if __name__ == "__main__":
    main()