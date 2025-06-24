#!/usr/bin/env python3
"""
Enhanced TabNet Prostate Cancer Dataset Merger
PHASE 2: Evidence-Based Classification & Class Balance Fix - FULLY CORRECTED

Target: Achieve 8-12% Likely_Actionable through evidence-based classification
"""

import pandas as pd
import numpy as np
from pathlib import Path
import os

def clean_chromosome_column(df):
    """Clean chromosome column for consistent formatting"""
    df = df.copy()
    if 'chromosome' in df.columns:
        # Remove 'chr' prefix and handle various formats
        df['chromosome'] = df['chromosome'].astype(str).str.replace('chr', '', case=False)
        df['chromosome'] = df['chromosome'].str.upper()
        
        # Handle special cases
        df['chromosome'] = df['chromosome'].replace({'23': 'X', '24': 'Y', 'MT': 'M'})
        
        # Remove any remaining non-standard values
        valid_chrs = [str(i) for i in range(1, 23)] + ['X', 'Y', 'M']
        df = df[df['chromosome'].isin(valid_chrs)]
    
    return df

def has_ref_alt_alleles(row):
    """Check if variant has both REF and ALT alleles for VEP compatibility"""
    ref_cols = ['Reference_Allele', 'GENOMIC_WT_ALLELE', 'ref', 'REF']
    alt_cols = ['Tumor_Seq_Allele2', 'GENOMIC_MUT_ALLELE', 'alt', 'ALT']
    
    has_ref = any(pd.notna(row.get(col)) and str(row.get(col)).strip() != '' for col in ref_cols)
    has_alt = any(pd.notna(row.get(col)) and str(row.get(col)).strip() != '' for col in alt_cols)
    
    return has_ref and has_alt

def classify_cosmic_variant(row):
    """CORRECTED: Evidence-based COSMIC variant classification"""
    print(f"🔍 Classifying COSMIC variant: {row.get('gene_symbol', 'UNKNOWN')}")
    
    if pd.isna(row.get('MUTATION_SOMATIC_STATUS')):
        return 'VUS'
    
    status = str(row['MUTATION_SOMATIC_STATUS']).lower()
    gene = str(row.get('gene_symbol', ''))
    
    # Get impact score if available
    impact_score = row.get('impact_score', 0)
    
    print(f"   Gene: {gene}, Status: {status}, Impact: {impact_score}")
    
    # Define therapeutic gene sets for enhanced classification
    parp_genes = {'BRCA1', 'BRCA2', 'ATM', 'CHEK2', 'PALB2', 'RAD51D', 'BRIP1', 'MLH1', 'MSH2', 'MSH6', 'PMS2', 'NBN', 'BARD1'}
    hormone_genes = {'AR', 'CYP17A1', 'SRD5A2', 'FOXA1', 'NCOR1', 'NCOR2'}
    high_priority_genes = parp_genes.union(hormone_genes).union({'TP53', 'PTEN', 'RB1'})
    
    # CRITICAL FIX: High-impact therapeutic variants should be Likely_Actionable
    if gene in high_priority_genes and impact_score >= 3:
        print(f"   → HIGH-IMPACT THERAPEUTIC: {gene} (score: {impact_score}) → Likely_Actionable")
        return 'Likely_Actionable'
    
    # High confidence pathogenic with therapeutic relevance
    if any(term in status for term in ['confirmed somatic', 'reported in another cancer sample as somatic']):
        if gene in high_priority_genes:
            print(f"   → CONFIRMED THERAPEUTIC: {gene} → Actionable_Pathogenic")
            return 'Actionable_Pathogenic'
        else:
            print(f"   → CONFIRMED NON-THERAPEUTIC: {gene} → Likely_Actionable")
            return 'Likely_Actionable'
    
    # Moderate confidence
    elif any(term in status for term in ['likely somatic', 'probable']):
        if gene in high_priority_genes:
            print(f"   → LIKELY THERAPEUTIC: {gene} → Likely_Actionable")
            return 'Likely_Actionable'
        else:
            print(f"   → LIKELY NON-THERAPEUTIC: {gene} → VUS")
            return 'VUS'
    
    # Unknown or low confidence - but check if therapeutic + high impact
    elif 'variant of unknown origin' in status:
        if gene in high_priority_genes and impact_score >= 2:  # Lower threshold for therapeutic genes
            print(f"   → UNKNOWN BUT THERAPEUTIC HIGH-IMPACT: {gene} → Likely_Actionable")
            return 'Likely_Actionable'
        else:
            print(f"   → UNKNOWN: {gene} → VUS")
            return 'VUS'
    
    # Default for other statuses
    else:
        print(f"   → DEFAULT: {gene} → VUS")
        return 'VUS'

def classify_clinvar_variant(row):
    """Evidence-based ClinVar variant classification using clinical_significance"""
    if pd.isna(row.get('clinical_significance')):
        return 'VUS'
    
    significance = str(row['clinical_significance']).lower()
    
    # High confidence pathogenic
    if 'pathogenic' in significance and 'conflicting' not in significance and 'likely' not in significance:
        return 'Actionable_Pathogenic'
    
    # Moderate confidence pathogenic
    elif 'likely pathogenic' in significance and 'conflicting' not in significance:
        return 'Likely_Actionable'
    
    # Benign variants
    elif 'benign' in significance and 'conflicting' not in significance:
        if 'likely' in significance:
            return 'Likely_Benign'
        else:
            return 'Benign'
    
    # Drug response (actionable)
    elif 'drug response' in significance:
        return 'Likely_Actionable'
    
    # Conflicting or uncertain
    elif any(term in significance for term in ['conflicting', 'uncertain', 'not provided']):
        return 'VUS'
    
    # Default
    else:
        return 'VUS'

def classify_tcga_variant(row):
    """Impact-based TCGA variant classification using Variant_Classification + therapeutic context"""
    if pd.isna(row.get('Variant_Classification')):
        return 'VUS'
    
    variant_class = str(row['Variant_Classification']).lower()
    gene = str(row.get('gene_symbol', ''))
    
    # Define therapeutic gene sets
    therapeutic_genes = {
        'BRCA1', 'BRCA2', 'ATM', 'CHEK2', 'PALB2', 'RAD51D', 'BRIP1', 'MLH1', 'MSH2', 'MSH6', 'PMS2', 'NBN',
        'AR', 'CYP17A1', 'SRD5A2', 'FOXA1', 'PTEN', 'PIK3CA', 'AKT1', 'MTOR', 'TP53', 'RB1'
    }
    
    # High impact variants - likely actionable if in therapeutic genes
    if any(term in variant_class for term in ['nonsense_mutation', 'frame_shift', 'splice_site']):
        if gene in therapeutic_genes:
            return 'Likely_Actionable'
        else:
            return 'VUS'  # High impact but not therapeutic target
    
    # Moderate impact variants - context dependent
    elif 'missense_mutation' in variant_class:
        if gene in therapeutic_genes:
            return 'Likely_Actionable'  # Therapeutic gene missense
        else:
            return 'VUS'  # Regular missense
    
    # Low impact variants
    elif any(term in variant_class for term in ['silent', 'synonymous']):
        return 'Likely_Benign'
    
    # In-frame indels in therapeutic genes
    elif 'in_frame' in variant_class and gene in therapeutic_genes:
        return 'Likely_Actionable'
    
    # Default for other variant types
    else:
        return 'VUS'

def apply_stratified_resampling(df, target_likely_actionable_pct=10.0):
    """CORRECTED: Apply stratified resampling to achieve target class balance"""
    
    # Check sklearn availability first
    try:
        from sklearn.utils import resample
        print("✅ sklearn.utils.resample available")
    except ImportError:
        print("❌ ERROR: sklearn not available - installing...")
        try:
            import subprocess
            subprocess.check_call(["pip", "install", "scikit-learn"])
            from sklearn.utils import resample
            print("✅ sklearn installed and imported successfully")
        except:
            print("❌ FATAL: Cannot install sklearn - skipping resampling")
            return df
    
    current_likely_actionable_pct = (df['variant_classification'] == 'Likely_Actionable').mean() * 100
    
    print(f"\n⚖️  CLASS BALANCE ANALYSIS:")
    print(f"   Current Likely_Actionable: {current_likely_actionable_pct:.1f}%")
    print(f"   Target Likely_Actionable: {target_likely_actionable_pct:.1f}%")
    print(f"   Trigger range: Outside 8.0-15.0%")
    
    # Enhanced trigger condition with explicit logic
    needs_resampling = current_likely_actionable_pct < 8.0 or current_likely_actionable_pct > 15.0
    print(f"   Needs resampling: {needs_resampling}")
    print(f"     - Too low: {current_likely_actionable_pct < 8.0} (< 8.0%)")
    print(f"     - Too high: {current_likely_actionable_pct > 15.0} (> 15.0%)")
    
    if needs_resampling:
        print("\n🚀 APPLYING STRATIFIED RESAMPLING...")
        
        # Show current distribution
        print("📊 Current distribution:")
        current_counts = df['variant_classification'].value_counts()
        for class_name, count in current_counts.items():
            pct = (count / len(df)) * 100
            print(f"   {class_name}: {count:,} ({pct:.1f}%)")
        
        # Target distribution for balanced TabNet training
        total_variants = len(df)
        target_distribution = {
            'Actionable_Pathogenic': int(0.15 * total_variants),  # 15%
            'Likely_Actionable': int(target_likely_actionable_pct/100 * total_variants),  # 10%
            'VUS': int(0.45 * total_variants),  # 45%
            'Likely_Benign': int(0.20 * total_variants),  # 20%
            'Benign': int(0.10 * total_variants)  # 10%
        }
        
        print(f"\n🎯 Target distribution:")
        for class_name, target_count in target_distribution.items():
            pct = (target_count / total_variants) * 100
            print(f"   {class_name}: {target_count:,} ({pct:.1f}%)")
        
        # Resample each class
        balanced_variants = []
        for class_name, target_count in target_distribution.items():
            class_variants = df[df['variant_classification'] == class_name]
            
            if len(class_variants) == 0:
                print(f"⚠️  WARNING: No variants found for class {class_name} - skipping")
                continue
                
            if len(class_variants) > target_count:
                # Downsample
                resampled = resample(class_variants, n_samples=target_count, random_state=42)
                print(f"  📉 {class_name}: {len(class_variants):,} → {target_count:,} (downsampled)")
            else:
                # Upsample
                resampled = resample(class_variants, n_samples=target_count, random_state=42, replace=True)
                print(f"  📈 {class_name}: {len(class_variants):,} → {target_count:,} (upsampled)")
            
            balanced_variants.append(resampled)
        
        if balanced_variants:
            result_df = pd.concat(balanced_variants, ignore_index=True)
            print(f"\n✅ RESAMPLING COMPLETED: {len(df):,} → {len(result_df):,} variants")
            
            # Verify final distribution
            print("\n📊 Final distribution:")
            final_counts = result_df['variant_classification'].value_counts()
            for class_name, count in final_counts.items():
                pct = (count / len(result_df)) * 100
                print(f"   {class_name}: {count:,} ({pct:.1f}%)")
            
            final_likely_actionable_pct = (result_df['variant_classification'] == 'Likely_Actionable').mean() * 100
            print(f"\n🎉 FINAL LIKELY_ACTIONABLE: {final_likely_actionable_pct:.1f}%")
            
            return result_df
        else:
            print("❌ No classes could be resampled - returning original dataset")
            return df
    else:
        print("\n✅ Class balance within acceptable range (8-15%) - no resampling needed")
        return df

def main():
    """Main merge function with enhanced classification logic"""
    
    print("🧬 Enhanced TabNet Prostate Cancer Dataset Merger")
    print("=" * 50)
    print("PHASE 2: Evidence-Based Classification & Class Balance Fix")
    print("=" * 50)
    
    # Get the project root directory
    script_dir = Path(__file__).parent
    project_root = script_dir.parent.parent
    
    # CORRECTED FILE PATHS - Based on actual repository structure
    file_paths = {
        'cosmic': project_root / "data" / "processed" / "cosmic_prostate" / "cosmic_prostate.csv",
        'clinvar': project_root / "data" / "processed" / "clinvar_prostate" / "clinvar_prostate.csv",
        'tcga': project_root / "data" / "processed" / "tcga_prad_prostate" / "tcga_prad_mutations.csv"
    }
    
    print(f"\n🔍 Checking for required input files...")
    
    # Check if all files exist
    missing_files = []
    for name, path in file_paths.items():
        if path.exists():
            size = path.stat().st_size
            print(f"   ✅ {name.upper()}: {path} ({size:,} bytes)")
        else:
            print(f"   ❌ {name.upper()}: {path} (NOT FOUND)")
            missing_files.append(f"{name.upper()}: {path}")
    
    if missing_files:
        print(f"\n❌ Missing required files:")
        for file in missing_files:
            print(f"   {file}")
        return False
    
    # Load datasets with proper dtype handling
    print("\n📊 Loading datasets...")
    
    # Load COSMIC
    cosmic_df = pd.read_csv(file_paths['cosmic'], low_memory=False)
    cosmic_df['data_source'] = 'COSMIC'
    print(f"   COSMIC: {len(cosmic_df):,} variants")
    
    # Load ClinVar  
    clinvar_df = pd.read_csv(file_paths['clinvar'], low_memory=False)
    clinvar_df['data_source'] = 'ClinVar'
    print(f"   ClinVar: {len(clinvar_df):,} variants")
    
    # Load TCGA
    tcga_df = pd.read_csv(file_paths['tcga'], low_memory=False)
    tcga_df['data_source'] = 'TCGA'
    print(f"   TCGA: {len(tcga_df):,} variants")
    
    # Standardize column names
    print("\n🔧 Standardizing columns...")
    
    # COSMIC columns - handle both 'gene' and 'GENE_SYMBOL'
    cosmic_clean = cosmic_df.copy()
    if 'gene' in cosmic_clean.columns:
        cosmic_clean['gene_symbol'] = cosmic_clean['gene']
    elif 'GENE_SYMBOL' in cosmic_clean.columns:
        cosmic_clean['gene_symbol'] = cosmic_clean['GENE_SYMBOL']
    
    if 'chr' in cosmic_clean.columns:
        cosmic_clean['chromosome'] = cosmic_clean['chr']
    elif 'CHROMOSOME' in cosmic_clean.columns:
        cosmic_clean['chromosome'] = cosmic_clean['CHROMOSOME']
        
    if 'pos' in cosmic_clean.columns:
        cosmic_clean['position'] = cosmic_clean['pos']
    elif 'GENOME_START' in cosmic_clean.columns:
        cosmic_clean['position'] = cosmic_clean['GENOME_START']
    
    # ClinVar columns - handle flexible naming
    clinvar_clean = clinvar_df.copy()
    if 'gene' in clinvar_clean.columns:
        clinvar_clean['gene_symbol'] = clinvar_clean['gene']
    elif 'Gene' in clinvar_clean.columns:
        clinvar_clean['gene_symbol'] = clinvar_clean['Gene']
    
    if 'chr' in clinvar_clean.columns:
        clinvar_clean['chromosome'] = clinvar_clean['chr']
    elif 'Chromosome' in clinvar_clean.columns:
        clinvar_clean['chromosome'] = clinvar_clean['Chromosome']
        
    if 'pos' in clinvar_clean.columns:
        clinvar_clean['position'] = clinvar_clean['pos']
    elif 'Position' in clinvar_clean.columns:
        clinvar_clean['position'] = clinvar_clean['Position']
    
    # TCGA columns - standard TCGA naming
    tcga_clean = tcga_df.copy()
    if 'Hugo_Symbol' in tcga_clean.columns:
        tcga_clean['gene_symbol'] = tcga_clean['Hugo_Symbol']
    elif 'Gene' in tcga_clean.columns:
        tcga_clean['gene_symbol'] = tcga_clean['Gene']
        
    if 'Chromosome' in tcga_clean.columns:
        tcga_clean['chromosome'] = tcga_clean['Chromosome']
    elif 'chr' in tcga_clean.columns:
        tcga_clean['chromosome'] = tcga_clean['chr']
        
    if 'Start_Position' in tcga_clean.columns:
        tcga_clean['position'] = tcga_clean['Start_Position']
    elif 'Position' in tcga_clean.columns:
        tcga_clean['position'] = tcga_clean['Position']
    
    # Clean chromosome columns for all datasets
    cosmic_clean = clean_chromosome_column(cosmic_clean)
    clinvar_clean = clean_chromosome_column(clinvar_clean)
    tcga_clean = clean_chromosome_column(tcga_clean)
    
    # PHASE 2: EVIDENCE-BASED CLASSIFICATION
    print("\n🎯 Applying evidence-based classification...")
    
    # COSMIC: Use MUTATION_SOMATIC_STATUS with enhanced logic
    print(f"\n🔬 COSMIC Classification (Debug Mode):")
    print(f"   Processing {len(cosmic_clean)} COSMIC variants...")
    
    # Apply classification with debug output for first few variants
    cosmic_clean['variant_classification'] = cosmic_clean.apply(classify_cosmic_variant, axis=1)
    cosmic_counts = cosmic_clean['variant_classification'].value_counts()
    print(f"   ✅ COSMIC classification complete: {dict(cosmic_counts)}")
    
    # ClinVar: Use clinical_significance
    print(f"\n🏥 ClinVar Classification:")
    clinvar_clean['variant_classification'] = clinvar_clean.apply(classify_clinvar_variant, axis=1)
    clinvar_counts = clinvar_clean['variant_classification'].value_counts()
    print(f"   ✅ ClinVar classification complete: {dict(clinvar_counts)}")
    
    # TCGA: Use Variant_Classification + therapeutic context
    print(f"\n🧬 TCGA Classification:")
    tcga_clean['variant_classification'] = tcga_clean.apply(classify_tcga_variant, axis=1)
    tcga_counts = tcga_clean['variant_classification'].value_counts()
    print(f"   ✅ TCGA classification complete: {dict(tcga_counts)}")
    
    # Combine datasets
    print("\n🔀 Merging datasets...")
    all_variants = pd.concat([cosmic_clean, clinvar_clean, tcga_clean], ignore_index=True)
    print(f"   Total variants: {len(all_variants):,}")
    
    # Check initial class distribution
    print("\n📊 Initial class distribution:")
    initial_counts = all_variants['variant_classification'].value_counts()
    for class_name, count in initial_counts.items():
        pct = (count / len(all_variants)) * 100
        print(f"   {class_name}: {count:,} ({pct:.1f}%)")
    
    # Add prostate pathway features
    print("\n🧬 Adding prostate pathway features...")
    
    # Define prostate cancer gene sets
    dna_repair_genes = {'BRCA1', 'BRCA2', 'ATM', 'CHEK2', 'PALB2', 'RAD51D', 'BRIP1', 'FANCA', 'MLH1', 'MSH2', 'MSH6', 'PMS2', 'NBN'}
    hormone_genes = {'AR', 'CYP17A1', 'SRD5A2', 'CYP19A1', 'ESR1', 'ESR2'}
    pi3k_genes = {'PTEN', 'PIK3CA', 'PIK3R1', 'AKT1', 'AKT2', 'AKT3', 'MTOR', 'TSC1', 'TSC2'}
    core_prostate_genes = {'AR', 'ERG', 'ETV1', 'TMPRSS2', 'NKX3-1', 'HOXB13', 'KLK3'}
    
    # Add pathway features
    all_variants['dna_repair_gene'] = all_variants['gene_symbol'].isin(dna_repair_genes)
    all_variants['hormone_pathway_gene'] = all_variants['gene_symbol'].isin(hormone_genes)
    all_variants['pi3k_pathway_gene'] = all_variants['gene_symbol'].isin(pi3k_genes)
    all_variants['core_prostate_gene'] = all_variants['gene_symbol'].isin(core_prostate_genes)
    
    # Add therapeutic target features
    all_variants['parp_inhibitor_target'] = all_variants['gene_symbol'].isin(dna_repair_genes)
    all_variants['hormone_therapy_target'] = all_variants['gene_symbol'].isin(hormone_genes)
    
    # Add allele information for VEP compatibility
    print("\n🧬 Adding allele data for VEP compatibility...")
    all_variants['has_ref_alt'] = all_variants.apply(has_ref_alt_alleles, axis=1)
    total_with_alleles = all_variants['has_ref_alt'].sum()
    allele_coverage = (total_with_alleles / len(all_variants)) * 100
    
    # PHASE 2: APPLY STRATIFIED RESAMPLING
    print("\n" + "="*60)
    print("PHASE 2 CRITICAL: STRATIFIED RESAMPLING")
    print("="*60)
    
    all_variants = apply_stratified_resampling(all_variants, target_likely_actionable_pct=10.0)
    
    # Final validation
    print("\n📊 FINAL CLASS DISTRIBUTION:")
    final_counts = all_variants['variant_classification'].value_counts()
    for class_name, count in final_counts.items():
        pct = (count / len(all_variants)) * 100
        status = "🎯" if class_name == 'Likely_Actionable' and 8 <= pct <= 12 else ""
        print(f"   {class_name}: {count:,} ({pct:.1f}%) {status}")
    
    # Save enhanced dataset
    output_dir = project_root / "data" / "processed" / "merged"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    output_file = output_dir / "merged_prostate_variants.csv"
    all_variants.to_csv(output_file, index=False)
    
    # Generate comprehensive report
    report_file = output_dir / "merge_report.txt"
    
    final_likely_actionable_pct = (all_variants['variant_classification'] == 'Likely_Actionable').mean() * 100
    
    with open(report_file, 'w') as f:
        f.write("Enhanced TabNet Prostate Cancer Dataset Merge Report\n")
        f.write("="*51 + "\n")
        
        f.write("📊 DATASET SUMMARY\n")
        f.write("-" * 50 + "\n")
        f.write(f"Total variants: {len(all_variants):,}\n")
        f.write(f"Unique genes: {all_variants['gene_symbol'].nunique()}\n")
        f.write(f"Unique chromosomes: {sorted(all_variants['chromosome'].unique())}\n")
        f.write(f"Data sources: {', '.join(all_variants['data_source'].unique())}\n\n")
        
        f.write("🧬 ALLELE DATA ANALYSIS (VEP READY)\n")
        f.write("-" * 50 + "\n")
        f.write(f"Variants with REF/ALT alleles: {total_with_alleles:,} ({allele_coverage:.1f}%)\n")
        f.write(f"VEP compatibility: {total_with_alleles:,} variants ready for annotation\n\n")
        
        f.write("📈 SOURCE DISTRIBUTION\n")
        f.write("-" * 50 + "\n")
        source_counts = all_variants['data_source'].value_counts()
        for source, count in source_counts.items():
            percentage = (count / len(all_variants)) * 100
            f.write(f"{source}: {count:,} ({percentage:.1f}%)\n")
        f.write("\n")
        
        f.write("🎯 VARIANT CLASSIFICATION DISTRIBUTION (PHASE 2 CORRECTED)\n")
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
        f.write(f"PHASE 2 CORRECTED: Evidence-based classification + resampling implemented\n")
        f.write(f"Class balance achieved: {final_likely_actionable_pct:.1f}% Likely_Actionable\n")
    
    print(f"\n   📄 Report saved: {report_file}")
    
    # Print summary
    print(f"\n" + "="*60)
    print(f"✅ PHASE 2 CORRECTED - FINAL RESULTS")
    print("="*60)
    print(f"Enhanced merged dataset saved to:")
    print(f"   {output_file}")
    print(f"\n📊 Dataset Summary:")
    print(f"   Total variants: {len(all_variants):,}")
    print(f"   Unique genes: {all_variants['gene_symbol'].nunique()}")
    print(f"   Data sources: {', '.join(all_variants['data_source'].unique())}")
    print(f"   🧬 VEP-ready variants: {total_with_alleles:,} ({allele_coverage:.1f}%)")
    
    print(f"\n🎯 PHASE 2 CLASS BALANCE RESULTS:")
    for class_name, count in final_counts.items():
        pct = (count / len(all_variants)) * 100
        if class_name == 'Likely_Actionable':
            if 8 <= pct <= 12:
                status = "✅ TARGET ACHIEVED"
            else:
                status = "❌ TARGET MISSED"
        else:
            status = ""
        print(f"   {class_name}: {count:,} ({pct:.1f}%) {status}")
    
    if 8 <= final_likely_actionable_pct <= 12:
        print(f"\n🎉 SUCCESS: {final_likely_actionable_pct:.1f}% Likely_Actionable (target: 8-12%)")
        print(f"🚀 READY FOR PHASE 3: VEP Optimization")
    else:
        print(f"\n⚠️  NEEDS REVIEW: {final_likely_actionable_pct:.1f}% Likely_Actionable (target: 8-12%)")
    
    return True

if __name__ == "__main__":
    success = main()
    if not success:
        exit(1)