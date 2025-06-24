#!/usr/bin/env python3

import pandas as pd
import numpy as np
import os

def define_therapeutic_gene_sets():
    """Define therapeutic target gene sets for prostate cancer"""
    
    # PARP Inhibitor Pathway Genes - Key targets
    parp_pathway_genes = {
        'BRCA1', 'BRCA2', 'ATM', 'CHEK2', 'PALB2', 'RAD51D', 'BRIP1', 
        'MLH1', 'MSH2', 'MSH6', 'PMS2', 'NBN', 'BARD1'
    }
    
    # Hormone Therapy Pathway Genes - Prostate specific
    hormone_pathway_genes = {
        'AR', 'CYP17A1', 'SRD5A2', 'FOXA1', 'NCOR1', 'NCOR2'
    }
    
    # PI3K/AKT/mTOR Pathway - Core genes
    pi3k_pathway_genes = {
        'PTEN', 'PIK3CA', 'AKT1', 'MTOR', 'TSC1', 'TSC2'
    }
    
    # Cell Cycle Control - Key genes
    cell_cycle_genes = {
        'TP53', 'RB1', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B'
    }
    
    # Emerging Targeted Therapy Genes - Prostate relevant
    emerging_targets = {
        'MYC', 'SPOP', 'ERG', 'ETV1', 'TMPRSS2', 'NKX3-1', 'HOXB13', 'ZFHX3'
    }
    
    return {
        'parp_pathway': parp_pathway_genes,
        'hormone_pathway': hormone_pathway_genes,
        'pi3k_pathway': pi3k_pathway_genes,
        'cell_cycle': cell_cycle_genes,
        'emerging_targets': emerging_targets
    }

def calculate_mutation_impact_score(mutation_desc):
    """Calculate impact score based on mutation description"""
    if pd.isna(mutation_desc):
        return 0
        
    mutation_str = str(mutation_desc).lower()
    
    # High Impact (Score 4)
    if any(term in mutation_str for term in ['nonsense', 'stop_gained', 'frameshift']):
        return 4
    
    # Moderate-High Impact (Score 3)
    elif any(term in mutation_str for term in ['splice_site', 'splice_region', 'inframe_deletion']):
        return 3
    
    # Moderate Impact (Score 2)
    elif 'missense' in mutation_str:
        return 2
    
    # Low Impact (Score 1)
    elif any(term in mutation_str for term in ['synonymous', 'silent']):
        return 1
    
    else:
        return 0

def has_prostate_context(row):
    """Check if a variant has prostate cancer context"""
    sample_name = str(row.get('SAMPLE_NAME', '')).lower()
    phenotype_id = str(row.get('COSMIC_PHENOTYPE_ID', ''))
    
    # TCGA-PRAD samples
    if 'tcga' in sample_name and ('prad' in sample_name or 'prostate' in sample_name):
        return True
    
    # Known prostate cell lines
    prostate_indicators = [
        'pc3', 'lncap', 'du145', '22rv1', 'vcap', 'lapc4', 
        'c4-2', 'mda-pca', 'prostate', 'npc3', 'red94', 'zhe06'
    ]
    
    for indicator in prostate_indicators:
        if indicator in sample_name:
            return True
    
    # Prostate phenotype IDs (if available)
    if phenotype_id in ['26', '592', '593', '594', '595', '596', '597', '598', '599', '600']:
        return True
    
    return False

def filter_prostate_samples(cosmic_df):
    """FIXED: Ensure ALL variants have prostate context"""
    
    print(f"Original COSMIC data: {len(cosmic_df)} mutations")
    
    # Get therapeutic gene sets
    gene_sets = define_therapeutic_gene_sets()
    all_therapeutic_genes = set()
    for gene_set in gene_sets.values():
        all_therapeutic_genes.update(gene_set)
    
    print(f"Targeting {len(all_therapeutic_genes)} therapeutic genes")
    
    # STEP 1: Original prostate-specific filtering (baseline)
    tcga_prad_samples = cosmic_df[
        cosmic_df['SAMPLE_NAME'].str.contains('TCGA-', na=False) &
        cosmic_df['SAMPLE_NAME'].str.contains('PRAD|prostate', case=False, na=False)
    ]
    
    prostate_cell_lines = [
        'PC3', 'LNCaP', 'DU145', '22Rv1', 'VCaP', 'LAPC4', 
        'C4-2', 'MDA-PCa', 'RWPE-1', 'PWR-1E', 'CWR22Rv1',
        'LNCaP104S', 'NPC3D', 'MDA-Pca-2b', 'PC31T', 'RED94',
        'ZHE06', 'NPC3F', 'PC35', 'LPC3p', 'BxPc3'
    ]
    
    cell_line_samples = cosmic_df[
        cosmic_df['SAMPLE_NAME'].str.contains('|'.join(prostate_cell_lines), case=False, na=False)
    ]
    
    # Baseline prostate samples
    baseline_prostate = pd.concat([tcga_prad_samples, cell_line_samples]).drop_duplicates()
    print(f"Baseline prostate samples: {len(baseline_prostate)}")
    
    # STEP 2: Add prostate context column to full dataset for filtering
    cosmic_df_with_context = cosmic_df.copy()
    cosmic_df_with_context['has_prostate_context'] = cosmic_df_with_context.apply(has_prostate_context, axis=1)
    cosmic_df_with_context['impact_score'] = cosmic_df_with_context['MUTATION_DESCRIPTION'].apply(calculate_mutation_impact_score)
    
    # STEP 3: Therapeutic variants with prostate context ONLY
    therapeutic_with_context = cosmic_df_with_context[
        cosmic_df_with_context['GENE_SYMBOL'].isin(all_therapeutic_genes) &
        cosmic_df_with_context['has_prostate_context']
    ]
    
    print(f"Therapeutic variants with prostate context: {len(therapeutic_with_context)}")
    
    # STEP 4: FIXED - High-impact therapeutic variants MUST have prostate context
    high_priority_genes = gene_sets['parp_pathway'].union(gene_sets['hormone_pathway']).union(gene_sets['cell_cycle'])
    
    # CRITICAL FIX: Add prostate context requirement
    high_impact_therapeutic = cosmic_df_with_context[
        cosmic_df_with_context['GENE_SYMBOL'].isin(high_priority_genes) &
        (cosmic_df_with_context['impact_score'] >= 3) &
        cosmic_df_with_context['MUTATION_SOMATIC_STATUS'].str.contains('somatic|confirmed', case=False, na=False) &
        cosmic_df_with_context['has_prostate_context']  # ← ADDED PROSTATE REQUIREMENT
    ]
    
    print(f"High-impact therapeutic variants (prostate-only): {len(high_impact_therapeutic)}")
    
    # STEP 5: Combine all PROSTATE-SPECIFIC sources
    all_candidates = pd.concat([
        baseline_prostate,
        therapeutic_with_context,
        high_impact_therapeutic
    ]).drop_duplicates()
    
    print(f"Combined prostate-specific candidates: {len(all_candidates)}")
    
    # STEP 6: Apply prioritization if needed
    if len(all_candidates) > 1000:
        print("Applying prioritization to optimize dataset size...")
        
        # Add annotations for prioritization
        all_candidates = all_candidates.copy()
        all_candidates['impact_score'] = all_candidates['MUTATION_DESCRIPTION'].apply(calculate_mutation_impact_score)
        all_candidates['has_prostate_context'] = all_candidates.apply(has_prostate_context, axis=1)
        
        # Add therapeutic pathway annotations
        all_candidates['is_parp_target'] = all_candidates['GENE_SYMBOL'].isin(gene_sets['parp_pathway'])
        all_candidates['is_hormone_target'] = all_candidates['GENE_SYMBOL'].isin(gene_sets['hormone_pathway'])
        all_candidates['is_pi3k_target'] = all_candidates['GENE_SYMBOL'].isin(gene_sets['pi3k_pathway'])
        all_candidates['is_cell_cycle_target'] = all_candidates['GENE_SYMBOL'].isin(gene_sets['cell_cycle'])
        all_candidates['is_emerging_target'] = all_candidates['GENE_SYMBOL'].isin(gene_sets['emerging_targets'])
        all_candidates['is_therapeutic_target'] = all_candidates['GENE_SYMBOL'].isin(all_therapeutic_genes)
        
        # Create priority score
        all_candidates['priority_score'] = 0
        
        # Prostate context bonus (all should have this now)
        all_candidates.loc[all_candidates['has_prostate_context'], 'priority_score'] += 5
        
        # Impact score
        all_candidates['priority_score'] += all_candidates['impact_score']
        
        # Therapeutic pathway bonuses
        all_candidates.loc[all_candidates['is_parp_target'], 'priority_score'] += 3
        all_candidates.loc[all_candidates['is_hormone_target'], 'priority_score'] += 3
        all_candidates.loc[all_candidates['is_pi3k_target'], 'priority_score'] += 2
        all_candidates.loc[all_candidates['is_cell_cycle_target'], 'priority_score'] += 2
        all_candidates.loc[all_candidates['is_emerging_target'], 'priority_score'] += 1
        
        # Somatic status bonus
        all_candidates.loc[all_candidates['MUTATION_SOMATIC_STATUS'].str.contains('confirmed', case=False, na=False), 'priority_score'] += 2
        
        # Sort by priority and take top variants
        all_candidates = all_candidates.sort_values('priority_score', ascending=False)
        
        # Keep top 1000 variants to maintain quality while ensuring prostate specificity
        final_samples = all_candidates.head(1000).copy()
        
        # Convert priority_score to actionability_score for consistency
        final_samples['actionability_score'] = final_samples['priority_score']
        
    else:
        final_samples = all_candidates.copy()
        
        # Add annotations
        final_samples['impact_score'] = final_samples['MUTATION_DESCRIPTION'].apply(calculate_mutation_impact_score)
        final_samples['has_prostate_context'] = final_samples.apply(has_prostate_context, axis=1)
        
        # Add therapeutic pathway annotations
        final_samples['is_parp_target'] = final_samples['GENE_SYMBOL'].isin(gene_sets['parp_pathway'])
        final_samples['is_hormone_target'] = final_samples['GENE_SYMBOL'].isin(gene_sets['hormone_pathway'])
        final_samples['is_pi3k_target'] = final_samples['GENE_SYMBOL'].isin(gene_sets['pi3k_pathway'])
        final_samples['is_cell_cycle_target'] = final_samples['GENE_SYMBOL'].isin(gene_sets['cell_cycle'])
        final_samples['is_emerging_target'] = final_samples['GENE_SYMBOL'].isin(gene_sets['emerging_targets'])
        final_samples['is_therapeutic_target'] = final_samples['GENE_SYMBOL'].isin(all_therapeutic_genes)
        
        # Create actionability score
        final_samples['actionability_score'] = 0
        final_samples.loc[final_samples['is_parp_target'], 'actionability_score'] += 3
        final_samples.loc[final_samples['is_hormone_target'], 'actionability_score'] += 3
        final_samples.loc[final_samples['is_pi3k_target'], 'actionability_score'] += 2
        final_samples.loc[final_samples['is_cell_cycle_target'], 'actionability_score'] += 2
        final_samples.loc[final_samples['is_emerging_target'], 'actionability_score'] += 1
        final_samples['actionability_score'] += final_samples['impact_score']
        
        # Sort by actionability
        final_samples = final_samples.sort_values('actionability_score', ascending=False)
    
    print(f"Final prostate-specific cancer mutations: {len(final_samples)}")
    print(f"Unique prostate samples: {final_samples['SAMPLE_NAME'].nunique()}")
    print(f"Therapeutic targets: {final_samples['is_therapeutic_target'].sum()}")
    print(f"Variants with prostate context: {final_samples['has_prostate_context'].sum()}")
    print(f"High impact variants: {(final_samples['impact_score'] >= 3).sum()}")
    
    return final_samples

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
        'GENOMIC_WT_ALLELE', 'GENOMIC_MUT_ALLELE',
        # Enhanced columns
        'impact_score', 'actionability_score', 'has_prostate_context',
        'is_parp_target', 'is_hormone_target', 'is_pi3k_target',
        'is_cell_cycle_target', 'is_emerging_target', 'is_therapeutic_target'
    ]
    
    # Keep only existing columns
    existing_columns = [col for col in columns_to_keep if col in cosmic_clean.columns]
    cosmic_filtered = cosmic_clean[existing_columns]
    
    return cosmic_filtered

def generate_enhancement_summary(prostate_df):
    """Generate summary of enhancements"""
    
    print("\n" + "="*50)
    print("PROSTATE-SPECIFIC COSMIC ANALYSIS SUMMARY")
    print("="*50)
    
    print(f"Total prostate-specific mutations: {len(prostate_df)}")
    print(f"Unique genes: {prostate_df['gene'].nunique()}")
    
    # Impact score distribution
    print("\nImpact Score Distribution:")
    for score in sorted(prostate_df['impact_score'].unique()):
        count = (prostate_df['impact_score'] == score).sum()
        print(f"  Score {score}: {count} ({count/len(prostate_df)*100:.1f}%)")
    
    # Prostate context
    prostate_context_count = prostate_df['has_prostate_context'].sum()
    print(f"\nProstate Context: {prostate_context_count}/{len(prostate_df)} ({prostate_context_count/len(prostate_df)*100:.1f}%)")
    
    # Therapeutic pathway representation
    print("\nTherapeutic Pathway Representation:")
    print(f"  PARP Inhibitor Targets: {prostate_df['is_parp_target'].sum()}")
    print(f"  Hormone Therapy Targets: {prostate_df['is_hormone_target'].sum()}")
    print(f"  PI3K/AKT/mTOR Targets: {prostate_df['is_pi3k_target'].sum()}")
    print(f"  Cell Cycle Targets: {prostate_df['is_cell_cycle_target'].sum()}")
    print(f"  Emerging Targets: {prostate_df['is_emerging_target'].sum()}")
    print(f"  Any Therapeutic Target: {prostate_df['is_therapeutic_target'].sum()}")
    
    # Top genes by actionability
    print("\nTop 15 Genes by Actionability Score:")
    top_genes = prostate_df.groupby('gene')['actionability_score'].max().sort_values(ascending=False).head(15)
    for gene, score in top_genes.items():
        count = (prostate_df['gene'] == gene).sum()
        print(f"  {gene}: {count} variants (max score: {score})")

def main():
    """Main function to create prostate-specific cosmic_prostate.csv"""
    
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
    
    # PROSTATE-SPECIFIC filter for prostate cancer
    prostate_df = filter_prostate_samples(cosmic_df)
    
    if len(prostate_df) == 0:
        print("No prostate cancer samples found!")
        print("Sample names preview:")
        print(cosmic_df['SAMPLE_NAME'].head(20).tolist())
        return
    
    # Standardize columns
    prostate_clean = standardize_cosmic_columns(prostate_df)
    
    # Create output directory and save to correct path
    output_file = "../../data/raw/variants/cosmic_prostate.csv"
    prostate_clean.to_csv(output_file, index=False)
    
    print(f"\nProstate-specific COSMIC data saved: {output_file}")
    print(f"Shape: {prostate_clean.shape}")
    print(f"Genes: {prostate_clean['gene'].nunique()}")
    print(f"Samples: {prostate_clean['SAMPLE_NAME'].nunique()}")
    
    # Generate enhancement summary
    generate_enhancement_summary(prostate_clean)
    
    # Show sample preview
    print("\nSample data preview:")
    display_columns = ['gene', 'chr', 'pos', 'MUTATION_DESCRIPTION', 'actionability_score', 'has_prostate_context']
    available_columns = [col for col in display_columns if col in prostate_clean.columns]
    print(prostate_clean[available_columns].head(10))
    
    # Show mutation types
    print("\nMutation types:")
    print(prostate_clean['MUTATION_DESCRIPTION'].value_counts().head(10))

if __name__ == "__main__":
    main()