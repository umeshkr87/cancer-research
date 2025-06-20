#!/usr/bin/env python3
"""
COSMIC Prostate Cancer Summary Generator
Creates summary.txt file equivalent to ClinVar format
"""

import pandas as pd
import numpy as np
from pathlib import Path
from collections import Counter

def generate_cosmic_summary(cosmic_csv_path, output_dir):
    """
    Generate comprehensive summary of COSMIC prostate cancer data
    
    Args:
        cosmic_csv_path: Path to cosmic_prostate.csv
        output_dir: Output directory for summary file
    """
    
    # Load data
    df = pd.read_csv(cosmic_csv_path)
    
    # Create output directory
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Summary file path
    summary_path = output_dir / "cosmic_prostate_summary.txt"
    
    # Define prostate cancer gene pathways
    dna_repair_genes = {'BRCA1', 'BRCA2', 'ATM', 'CHEK2', 'PALB2', 'RAD51D', 'BRIP1', 'FANCA', 'MLH1', 'MSH2', 'MSH6', 'PMS2', 'NBN'}
    pi3k_pathway_genes = {'PTEN', 'PIK3CA', 'PIK3R1', 'AKT1', 'AKT2', 'AKT3', 'MTOR', 'TSC1', 'TSC2'}
    ar_pathway_genes = {'AR', 'FOXA1', 'HOXB13', 'NKX3-1'}
    cell_cycle_genes = {'TP53', 'RB1', 'CDKN1B', 'CDKN2A', 'CDK4', 'CCND1'}
    chromatin_genes = {'ARID1A', 'KMT2C', 'KMT2D', 'CHD1', 'SPOP'}
    
    # Calculate statistics
    total_mutations = len(df)
    unique_genes = df['gene'].nunique()
    
    # Somatic mutations
    somatic_mutations = df[df['MUTATION_SOMATIC_STATUS'].str.contains('somatic', case=False, na=False)]
    somatic_count = len(somatic_mutations)
    somatic_percentage = (somatic_count / total_mutations) * 100
    
    # Chromosomes
    chromosomes = sorted([str(x) for x in df['chr'].unique() if pd.notna(x)])
    
    # Count mutations by pathway
    pathway_counts = {}
    for gene in df['gene']:
        if gene in dna_repair_genes:
            pathway_counts['DNA Repair Pathway'] = pathway_counts.get('DNA Repair Pathway', 0) + 1
        if gene in pi3k_pathway_genes:
            pathway_counts['PI3K/AKT/mTOR Pathway'] = pathway_counts.get('PI3K/AKT/mTOR Pathway', 0) + 1
        if gene in ar_pathway_genes:
            pathway_counts['Androgen Receptor Pathway'] = pathway_counts.get('Androgen Receptor Pathway', 0) + 1
        if gene in cell_cycle_genes:
            pathway_counts['Cell Cycle Control'] = pathway_counts.get('Cell Cycle Control', 0) + 1
        if gene in chromatin_genes:
            pathway_counts['Chromatin Remodeling'] = pathway_counts.get('Chromatin Remodeling', 0) + 1
    
    # Write summary file
    with open(summary_path, 'w') as f:
        f.write("COSMIC Prostate Cancer Mutation Analysis Summary\n")
        f.write("===============================================\n\n")
        
        # Basic statistics
        f.write(f"Total mutations: {total_mutations:,}\n")
        f.write(f"Unique genes: {unique_genes}\n")
        f.write(f"Somatic mutations: {somatic_count:,} ({somatic_percentage:.1f}%)\n")
        f.write(f"Chromosomes: {', '.join(chromosomes)}\n\n")
        
        # Mutation type distribution
        f.write("Mutation Type Distribution:\n")
        mutation_types = df['MUTATION_DESCRIPTION'].value_counts()
        for mut_type, count in mutation_types.items():
            f.write(f"  {mut_type}: {count}\n")
        f.write("\n")
        
        # Somatic status distribution
        f.write("Somatic Status Distribution:\n")
        somatic_status = df['MUTATION_SOMATIC_STATUS'].value_counts()
        for status, count in somatic_status.items():
            f.write(f"  {status}: {count}\n")
        f.write("\n")
        
        # Sample distribution
        f.write("Sample Distribution:\n")
        samples = df['SAMPLE_NAME'].value_counts()
        for sample, count in samples.items():
            f.write(f"  {sample}: {count}\n")
        f.write("\n")
        
        # Top genes
        f.write("Top 20 Prostate Cancer Genes:\n")
        top_genes = df['gene'].value_counts().head(20)
        for gene, count in top_genes.items():
            f.write(f"  {gene}: {count}\n")
        f.write("\n")
        
        # Pathway distribution
        f.write("Pathway Distribution:\n")
        for pathway, count in sorted(pathway_counts.items()):
            f.write(f"  {pathway}: {count}\n")
        f.write("\n")
        
        # Chromosome distribution
        f.write("Chromosome Distribution:\n")
        chr_counts = df['chr'].value_counts().sort_index()
        for chr_name, count in chr_counts.items():
            f.write(f"  Chr {chr_name}: {count}\n")
        f.write("\n")
        
        # Key therapeutic targets
        f.write("Key Therapeutic Targets:\n")
        
        # PARP inhibitor targets
        parp_targets = df[df['gene'].isin(dna_repair_genes)]
        f.write(f"  PARP Inhibitor Targets (BRCA1/2, ATM, CHEK2): {len(parp_targets)}\n")
        
        # Hormone therapy targets
        hormone_targets = df[df['gene'].isin(ar_pathway_genes)]
        f.write(f"  Hormone Therapy Targets (AR pathway): {len(hormone_targets)}\n")
        
        # DNA repair deficient (immunotherapy relevant)
        dna_repair_mutations = df[df['gene'].isin(dna_repair_genes)]
        f.write(f"  Immunotherapy Biomarkers (DNA repair deficient): {len(dna_repair_mutations)}\n")
        
        # All actionable targets
        actionable_genes = dna_repair_genes | pi3k_pathway_genes | ar_pathway_genes
        actionable_mutations = df[df['gene'].isin(actionable_genes)]
        f.write(f"  Targeted Therapy Candidates: {len(actionable_mutations)}\n")
        f.write("\n")
        
        # Clinical actionability assessment
        f.write("Clinical Actionability:\n")
        
        # Highly actionable (known therapeutic targets)
        highly_actionable = df[df['gene'].isin({'BRCA1', 'BRCA2', 'AR', 'PTEN', 'PIK3CA'})]
        highly_actionable_pct = (len(highly_actionable) / total_mutations) * 100
        f.write(f"  Highly Actionable: {len(highly_actionable)} ({highly_actionable_pct:.1f}%)\n")
        
        # Potentially actionable
        potentially_actionable = df[df['gene'].isin(actionable_genes)] 
        potentially_actionable = potentially_actionable[~potentially_actionable['gene'].isin({'BRCA1', 'BRCA2', 'AR', 'PTEN', 'PIK3CA'})]
        potentially_actionable_pct = (len(potentially_actionable) / total_mutations) * 100
        f.write(f"  Potentially Actionable: {len(potentially_actionable)} ({potentially_actionable_pct:.1f}%)\n")
        
        # Research interest
        research_interest = total_mutations - len(highly_actionable) - len(potentially_actionable)
        research_interest_pct = (research_interest / total_mutations) * 100
        f.write(f"  Research Interest: {research_interest} ({research_interest_pct:.1f}%)\n")
        f.write("\n")
        
        # Quality metrics
        f.write("Quality Metrics:\n")
        
        # Complete coordinates
        complete_coords = df[df['pos'].notna() & df['chr'].notna()]
        complete_coords_pct = (len(complete_coords) / total_mutations) * 100
        f.write(f"  Complete genomic coordinates: {len(complete_coords)} ({complete_coords_pct:.1f}%)\n")
        
        # Protein annotation
        protein_annotation = df[df['MUTATION_AA'].notna() & (df['MUTATION_AA'] != 'p.?')]
        protein_annotation_pct = (len(protein_annotation) / total_mutations) * 100
        f.write(f"  Protein-level annotation: {len(protein_annotation)} ({protein_annotation_pct:.1f}%)\n")
        
        # Cell line validation
        validated_samples = df[df['SAMPLE_NAME'].notna()]
        validated_samples_pct = (len(validated_samples) / total_mutations) * 100
        f.write(f"  Cell line validation: {len(validated_samples)} ({validated_samples_pct:.1f}%)\n")
        
        # Somatic confirmation
        f.write(f"  Somatic status confirmed: {somatic_count} ({somatic_percentage:.1f}%)\n")
    
    print(f"✅ COSMIC summary generated: {summary_path}")
    print(f"📊 Summary includes {total_mutations:,} mutations across {unique_genes} genes")
    
    return summary_path

def main():
    """Main function to generate COSMIC summary"""
    
    # Define paths
    cosmic_csv = "data/processed/cosmic_prostate/cosmic_prostate.csv"
    output_dir = "data/processed/cosmic_prostate"
    
    # Alternative paths to try
    possible_paths = [
        cosmic_csv,
        "cosmic_prostate.csv",
        "../data/processed/cosmic_prostate/cosmic_prostate.csv",
        "../../data/processed/cosmic_prostate/cosmic_prostate.csv"
    ]
    
    # Find the file
    file_found = False
    for path in possible_paths:
        if Path(path).exists():
            print(f"📁 Found COSMIC file: {path}")
            
            # Determine output directory
            if "data/processed/cosmic_prostate" in path:
                output_dir = str(Path(path).parent)
            else:
                output_dir = "."
            
            # Generate summary
            summary_path = generate_cosmic_summary(path, output_dir)
            
            print(f"\n🎯 SUMMARY GENERATED:")
            print(f"📄 File: {summary_path}")
            print(f"📂 Location: {Path(summary_path).absolute()}")
            
            file_found = True
            break
    
    if not file_found:
        print("❌ cosmic_prostate.csv file not found!")
        print("\nTried these paths:")
        for path in possible_paths:
            print(f"  - {path}")
        print("\nPlease run from the project root or adjust the path.")

if __name__ == "__main__":
    main()