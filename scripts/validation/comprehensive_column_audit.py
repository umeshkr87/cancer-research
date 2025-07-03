#!/usr/bin/env python3
"""
Comprehensive Data Quality Investigation Script
Phase 1: Column-by-Column Analysis for Prostate Cancer Variant Classification Pipeline

This script performs systematic analysis of all 89 columns in the enhanced prostate cancer dataset
to identify data quality issues, concatenation patterns, and clinical relevance scoring.

Author: PhD Research Student, University of Illinois
Contact: aa107@illinois.edu
"""

import pandas as pd
import numpy as np
import re
import os
import sys
from datetime import datetime
from collections import defaultdict, Counter
import json

class ProstateCancerDataAuditor:
    """Comprehensive data quality auditor for prostate cancer variant classification pipeline."""
    
    def __init__(self, enhanced_file_path, merged_file_path=None):
        """
        Initialize the auditor with file paths.
        
        Args:
            enhanced_file_path: Path to final enhanced dataset (89 columns)
            merged_file_path: Path to pre-VEP merged dataset (62 columns) for comparison
        """
        self.enhanced_file = enhanced_file_path
        self.merged_file = merged_file_path
        self.enhanced_df = None
        self.merged_df = None
        
        # Clinical relevance scoring dictionary
        self.clinical_scores = {
            'ESSENTIAL': 5,  # Critical for ML training
            'HIGH': 4,       # Important for clinical interpretation
            'MEDIUM': 3,     # Useful for analysis
            'LOW': 2,        # Optional/supplementary
            'ARTIFACT': 1    # Pipeline artifacts, not clinically relevant
        }
        
        # Concatenation patterns to detect
        self.concat_patterns = {
            'ampersand': r'&',
            'pipe': r'\|', 
            'comma': r',',
            'semicolon': r';',
            'slash': r'/',
            'colon': r'::'
        }
        
        # Results storage
        self.audit_results = {
            'summary': {},
            'column_analysis': {},
            'concatenation_analysis': {},
            'clinical_scoring': {},
            'pipeline_attribution': {},
            'recommendations': []
        }
    
    def load_datasets(self):
        """Load the datasets for analysis."""
        print("📊 Loading datasets...")
        
        try:
            # Load enhanced dataset (final output)
            if os.path.exists(self.enhanced_file):
                self.enhanced_df = pd.read_csv(self.enhanced_file, low_memory=False)
                print(f"✅ Enhanced dataset loaded: {self.enhanced_df.shape[0]:,} variants × {self.enhanced_df.shape[1]} columns")
            else:
                raise FileNotFoundError(f"Enhanced dataset not found: {self.enhanced_file}")
            
            # Load merged dataset for comparison (optional)
            if self.merged_file and os.path.exists(self.merged_file):
                self.merged_df = pd.read_csv(self.merged_file, low_memory=False)
                print(f"✅ Merged dataset loaded: {self.merged_df.shape[0]:,} variants × {self.merged_df.shape[1]} columns")
            else:
                print("⚠️  Merged dataset not available for comparison")
                
        except Exception as e:
            print(f"❌ Error loading datasets: {e}")
            sys.exit(1)
    
    def analyze_column_data_types(self, df, column_name):
        """Analyze data types and patterns in a specific column."""
        series = df[column_name]
        
        analysis = {
            'total_values': len(series),
            'non_null_values': series.notna().sum(),
            'null_values': series.isna().sum(),
            'unique_values': series.nunique(),
            'data_type': str(series.dtype)
        }
        
        # Analyze non-null values only
        non_null_series = series.dropna().astype(str)
        
        if len(non_null_series) > 0:
            # Detect concatenation patterns
            concat_analysis = {}
            for pattern_name, pattern_regex in self.concat_patterns.items():
                matches = non_null_series.str.contains(pattern_regex, regex=True, na=False)
                concat_count = matches.sum()
                if concat_count > 0:
                    concat_analysis[pattern_name] = {
                        'count': int(concat_count),
                        'percentage': float(concat_count / len(non_null_series) * 100),
                        'examples': non_null_series[matches].head(3).tolist()
                    }
            
            analysis['concatenation'] = concat_analysis
            
            # Value length statistics
            lengths = non_null_series.str.len()
            analysis['value_lengths'] = {
                'min': int(lengths.min()),
                'max': int(lengths.max()),
                'mean': float(lengths.mean()),
                'median': float(lengths.median())
            }
            
            # Sample values
            analysis['sample_values'] = non_null_series.head(5).tolist()
            
        return analysis
    
    def assign_clinical_relevance_score(self, column_name):
        """Assign clinical relevance score based on column name and biological importance."""
        column_lower = column_name.lower()
        
        # Essential columns for ML training
        essential_keywords = [
            'clinical_significance', 'clin_sig', 'pathogenicity', 'benign', 'pathogenic',
            'consequence', 'mutation_type', 'variant_class', 'gene_symbol', 'gene',
            'chromosome', 'position', 'reference', 'alternate', 'ref', 'alt'
        ]
        
        # High importance columns
        high_keywords = [
            'alphamissense', 'sift', 'polyphen', 'cadd', 'gerp', 'phred',
            'cosmic', 'clinvar', 'tcga', 'frequency', 'allele_freq',
            'exon', 'transcript', 'protein', 'amino_acid'
        ]
        
        # Medium importance columns
        medium_keywords = [
            'sample', 'tissue', 'tumor', 'normal', 'stage', 'grade',
            'treatment', 'response', 'survival', 'age', 'race'
        ]
        
        # Artifact/technical columns
        artifact_keywords = [
            'id', 'uuid', 'index', 'file', 'source', 'batch', 'version',
            'timestamp', 'checksum', 'flag', 'status', 'url', 'link'
        ]
        
        # Check for matches
        if any(keyword in column_lower for keyword in essential_keywords):
            return 'ESSENTIAL'
        elif any(keyword in column_lower for keyword in high_keywords):
            return 'HIGH'
        elif any(keyword in column_lower for keyword in medium_keywords):
            return 'MEDIUM'
        elif any(keyword in column_lower for keyword in artifact_keywords):
            return 'ARTIFACT'
        else:
            return 'LOW'
    
    def identify_pipeline_source(self, column_name, has_concatenation):
        """Identify which pipeline stage likely introduced data quality issues."""
        column_lower = column_name.lower()
        
        # VEP-specific columns (added during annotation)
        vep_columns = [
            'consequence', 'impact', 'symbol', 'gene', 'feature_type',
            'feature', 'biotype', 'exon', 'intron', 'hgvsc', 'hgvsp',
            'cdna_position', 'cds_position', 'protein_position', 'amino_acids',
            'codons', 'existing_variation', 'distance', 'strand', 'flags',
            'symbol_source', 'hgnc_id', 'canonical', 'mane', 'tsl', 'appris',
            'ccds', 'ensp', 'swissprot', 'trembl', 'uniparc', 'uniprot_isoform',
            'source', 'given_ref', 'used_ref', 'bam_edit', 'sift', 'polyphen',
            'domains', 'hgvs_offset', 'gmaf', 'afr_maf', 'amr_maf', 'eas_maf',
            'eur_maf', 'sas_maf', 'gnomadg_af', 'max_af', 'clin_sig',
            'somatic', 'pubmed', 'var_synonyms'
        ]
        
        # Enhancement-specific columns (added during functional enhancement)
        enhancement_columns = [
            'alphamissense', 'pathogenicity_score', 'confidence'
        ]
        
        # Source dataset columns (from original COSMIC/ClinVar/TCGA)
        source_columns = [
            'cosmic_id', 'clinvar_id', 'tcga_barcode', 'sample_id',
            'genomic_wt_allele', 'genomic_mut_allele', 'mutation_description',
            'clinical_significance', 'review_status', 'last_evaluated',
            'reference_allele', 'tumor_seq_allele', 'variant_classification'
        ]
        
        if any(vep_col in column_lower for vep_col in vep_columns):
            if has_concatenation:
                return 'VEP_ANNOTATION'
            else:
                return 'VEP_CLEAN'
        elif any(enh_col in column_lower for enh_col in enhancement_columns):
            return 'FUNCTIONAL_ENHANCEMENT'
        elif any(src_col in column_lower for src_col in source_columns):
            return 'SOURCE_DATASET'
        else:
            return 'MERGE_STAGE'
    
    def run_comprehensive_audit(self):
        """Execute the complete data quality audit."""
        print("\n🔍 STARTING COMPREHENSIVE COLUMN AUDIT")
        print("="*50)
        
        # Analyze each column in enhanced dataset
        total_columns = len(self.enhanced_df.columns)
        total_concatenated_columns = 0
        total_concatenated_values = 0
        
        for i, column in enumerate(self.enhanced_df.columns, 1):
            print(f"Analyzing column {i}/{total_columns}: {column}")
            
            # Perform detailed analysis
            column_analysis = self.analyze_column_data_types(self.enhanced_df, column)
            clinical_score = self.assign_clinical_relevance_score(column)
            
            # Check for concatenation
            has_concatenation = bool(column_analysis.get('concatenation', {}))
            if has_concatenation:
                total_concatenated_columns += 1
                # Count total concatenated values across all patterns
                concat_counts = [pattern['count'] for pattern in column_analysis['concatenation'].values()]
                total_concatenated_values += sum(concat_counts)
            
            pipeline_source = self.identify_pipeline_source(column, has_concatenation)
            
            # Store results
            self.audit_results['column_analysis'][column] = {
                'analysis': column_analysis,
                'clinical_relevance': clinical_score,
                'clinical_score': self.clinical_scores[clinical_score],
                'pipeline_source': pipeline_source,
                'has_concatenation': has_concatenation,
                'ml_ready': not has_concatenation and clinical_score in ['ESSENTIAL', 'HIGH', 'MEDIUM']
            }
        
        # Generate summary statistics
        self.audit_results['summary'] = {
            'total_columns': total_columns,
            'concatenated_columns': total_concatenated_columns,
            'concatenated_percentage': (total_concatenated_columns / total_columns) * 100,
            'total_variants': len(self.enhanced_df),
            'total_concatenated_values': total_concatenated_values,
            'audit_timestamp': datetime.now().isoformat()
        }
        
        # Perform additional analyses
        self._analyze_concatenation_patterns()
        self._score_clinical_importance()
        self._attribute_pipeline_sources()
        self._generate_recommendations()
        
        print(f"\n✅ Audit completed: {total_concatenated_columns}/{total_columns} columns have concatenation issues")
    
    def _analyze_concatenation_patterns(self):
        """Analyze patterns of concatenation across the dataset."""
        pattern_summary = defaultdict(int)
        source_concatenation = defaultdict(int)
        
        for column, data in self.audit_results['column_analysis'].items():
            if data['has_concatenation']:
                source_concatenation[data['pipeline_source']] += 1
                
                for pattern_name in data['analysis']['concatenation'].keys():
                    pattern_summary[pattern_name] += 1
        
        self.audit_results['concatenation_analysis'] = {
            'patterns_by_frequency': dict(pattern_summary),
            'sources_with_concatenation': dict(source_concatenation),
            'most_problematic_patterns': sorted(pattern_summary.items(), key=lambda x: x[1], reverse=True)[:5]
        }
    
    def _score_clinical_importance(self):
        """Score and rank columns by clinical importance."""
        clinical_ranking = defaultdict(list)
        
        for column, data in self.audit_results['column_analysis'].items():
            relevance = data['clinical_relevance']
            clinical_ranking[relevance].append({
                'column': column,
                'has_concatenation': data['has_concatenation'],
                'ml_ready': data['ml_ready']
            })
        
        self.audit_results['clinical_scoring'] = {
            'ranking_by_importance': dict(clinical_ranking),
            'essential_columns_with_issues': len([c for c in clinical_ranking['ESSENTIAL'] if c['has_concatenation']]),
            'total_ml_ready_columns': sum(1 for data in self.audit_results['column_analysis'].values() if data['ml_ready'])
        }
    
    def _attribute_pipeline_sources(self):
        """Attribute concatenation issues to specific pipeline stages."""
        source_attribution = defaultdict(list)
        
        for column, data in self.audit_results['column_analysis'].items():
            if data['has_concatenation']:
                source_attribution[data['pipeline_source']].append(column)
        
        self.audit_results['pipeline_attribution'] = {
            'issues_by_source': dict(source_attribution),
            'vep_introduced_issues': len(source_attribution.get('VEP_ANNOTATION', [])),
            'source_dataset_issues': len(source_attribution.get('SOURCE_DATASET', [])),
            'merge_stage_issues': len(source_attribution.get('MERGE_STAGE', []))
        }
    
    def _generate_recommendations(self):
        """Generate prioritized recommendations for data quality improvement."""
        recommendations = []
        
        # High priority: Essential columns with concatenation
        essential_issues = [
            col for col, data in self.audit_results['column_analysis'].items()
            if data['clinical_relevance'] == 'ESSENTIAL' and data['has_concatenation']
        ]
        
        if essential_issues:
            recommendations.append({
                'priority': 'CRITICAL',
                'action': 'Fix concatenation in essential columns',
                'columns': essential_issues,
                'rationale': 'These columns are critical for ML training but have concatenation issues'
            })
        
        # Medium priority: VEP parameter optimization
        vep_issues = self.audit_results['pipeline_attribution'].get('issues_by_source', {}).get('VEP_ANNOTATION', [])
        if len(vep_issues) > 10:
            recommendations.append({
                'priority': 'HIGH',
                'action': 'Optimize VEP parameters to reduce concatenation',
                'columns': vep_issues[:10],  # Show top 10
                'rationale': f'VEP introduced concatenation in {len(vep_issues)} columns'
            })
        
        # Lower priority: Clean up artifact columns
        artifact_columns = [
            col for col, data in self.audit_results['column_analysis'].items()
            if data['clinical_relevance'] == 'ARTIFACT'
        ]
        
        if len(artifact_columns) > 5:
            recommendations.append({
                'priority': 'LOW',
                'action': 'Consider removing artifact columns',
                'columns': artifact_columns,
                'rationale': 'These columns are not clinically relevant and may confuse ML models'
            })
        
        self.audit_results['recommendations'] = recommendations
    
    def save_results(self, output_dir):
        """Save comprehensive audit results to files."""
        os.makedirs(output_dir, exist_ok=True)
        
        # Save JSON results
        json_file = os.path.join(output_dir, 'comprehensive_audit_results.json')
        with open(json_file, 'w') as f:
            json.dump(self.audit_results, f, indent=2, default=str)
        
        # Save human-readable report
        report_file = os.path.join(output_dir, 'data_quality_audit_report.txt')
        self._generate_text_report(report_file)
        
        # Save CSV summary
        csv_file = os.path.join(output_dir, 'column_summary.csv')
        self._generate_csv_summary(csv_file)
        
        print(f"\n📄 Results saved:")
        print(f"   JSON: {json_file}")
        print(f"   Report: {report_file}")
        print(f"   CSV: {csv_file}")
    
    def _generate_text_report(self, filename):
        """Generate human-readable text report."""
        with open(filename, 'w') as f:
            f.write("PROSTATE CANCER DATA QUALITY AUDIT REPORT\n")
            f.write("="*50 + "\n\n")
            
            # Summary section
            summary = self.audit_results['summary']
            f.write(f"AUDIT SUMMARY\n")
            f.write(f"Timestamp: {summary['audit_timestamp']}\n")
            f.write(f"Total columns analyzed: {summary['total_columns']}\n")
            f.write(f"Columns with concatenation: {summary['concatenated_columns']} ({summary['concatenated_percentage']:.1f}%)\n")
            f.write(f"Total variants: {summary['total_variants']:,}\n\n")
            
            # Recommendations section
            f.write("PRIORITY RECOMMENDATIONS\n")
            f.write("-"*25 + "\n")
            for i, rec in enumerate(self.audit_results['recommendations'], 1):
                f.write(f"{i}. [{rec['priority']}] {rec['action']}\n")
                f.write(f"   Rationale: {rec['rationale']}\n")
                f.write(f"   Affected columns: {len(rec['columns'])}\n\n")
            
            # Pipeline attribution
            f.write("PIPELINE SOURCE ATTRIBUTION\n")
            f.write("-"*30 + "\n")
            attribution = self.audit_results['pipeline_attribution']
            for source, columns in attribution['issues_by_source'].items():
                f.write(f"{source}: {len(columns)} columns\n")
            f.write("\n")
            
            # Top problematic columns
            f.write("TOP 10 MOST PROBLEMATIC COLUMNS\n")
            f.write("-"*35 + "\n")
            problematic_columns = [
                (col, data) for col, data in self.audit_results['column_analysis'].items()
                if data['has_concatenation'] and data['clinical_relevance'] in ['ESSENTIAL', 'HIGH']
            ]
            problematic_columns.sort(key=lambda x: x[1]['clinical_score'], reverse=True)
            
            for i, (col, data) in enumerate(problematic_columns[:10], 1):
                f.write(f"{i:2d}. {col} ({data['clinical_relevance']})\n")
                concat_info = data['analysis']['concatenation']
                for pattern, info in concat_info.items():
                    f.write(f"     {pattern}: {info['count']} values ({info['percentage']:.1f}%)\n")
                f.write("\n")
    
    def _generate_csv_summary(self, filename):
        """Generate CSV summary for easy analysis."""
        rows = []
        for column, data in self.audit_results['column_analysis'].items():
            row = {
                'column_name': column,
                'clinical_relevance': data['clinical_relevance'],
                'clinical_score': data['clinical_score'],
                'pipeline_source': data['pipeline_source'],
                'has_concatenation': data['has_concatenation'],
                'ml_ready': data['ml_ready'],
                'total_values': data['analysis']['total_values'],
                'non_null_values': data['analysis']['non_null_values'],
                'unique_values': data['analysis']['unique_values'],
                'concatenation_patterns': len(data['analysis'].get('concatenation', {}))
            }
            rows.append(row)
        
        df = pd.DataFrame(rows)
        df.to_csv(filename, index=False)


def main():
    """Main execution function."""
    # File paths
    enhanced_file = "/u/aa107/uiuc-cancer-research/data/processed/tabnet_csv/prostate_variants_tabnet_enhanced.csv"
    merged_file = "/u/aa107/uiuc-cancer-research/data/processed/merged/merged_prostate_variants.csv"
    output_dir = "/u/aa107/uiuc-cancer-research/results/validation/phase1_audit"
    
    print("🧬 PROSTATE CANCER DATA QUALITY AUDITOR")
    print("Phase 1: Comprehensive Column Analysis")
    print("="*50)
    
    # Initialize auditor
    auditor = ProstateCancerDataAuditor(enhanced_file, merged_file)
    
    # Load datasets
    auditor.load_datasets()
    
    # Run comprehensive audit
    auditor.run_comprehensive_audit()
    
    # Save results
    auditor.save_results(output_dir)
    
    print("\n🎯 AUDIT COMPLETED SUCCESSFULLY!")
    print(f"Results saved to: {output_dir}")
    print("\nNext steps:")
    print("1. Review generated reports")
    print("2. Prioritize fixing essential columns with concatenation")
    print("3. Consider VEP parameter optimization")
    print("4. Proceed to Phase 2: Legitimate vs. Problematic Classification")


if __name__ == "__main__":
    main()