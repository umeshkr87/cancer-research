#!/usr/bin/env python3
"""
Comprehensive Data Quality Investigation Script - VCF MODE
Phase 1: Column-by-Column Analysis for Prostate Cancer Variant Classification Pipeline

This script performs systematic analysis of VCF CSQ fields to identify data quality issues,
concatenation patterns, and clinical relevance scoring at the VEP annotation level.

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

class ProstateCancerVCFAuditor:
    """Comprehensive VCF data quality auditor for prostate cancer variant classification pipeline."""
    
    def __init__(self, vcf_file_path, enhanced_csv_path=None):
        """
        Initialize the auditor with VCF file path.
        
        Args:
            vcf_file_path: Path to VEP annotated VCF file
            enhanced_csv_path: Path to enhanced CSV for comparison (optional)
        """
        self.vcf_file = vcf_file_path
        self.enhanced_csv = enhanced_csv_path
        self.vcf_data = []
        self.csq_fields = []
        self.csq_field_mapping = {}
        
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
    
    def parse_vcf_header(self):
        """Extract CSQ field information from VCF header."""
        print("📋 Parsing VCF header for CSQ field definitions...")
        
        with open(self.vcf_file, 'r') as f:
            for line in f:
                if line.startswith('##INFO=<ID=CSQ'):
                    # Extract Format field from CSQ header
                    match = re.search(r'Format: ([^"]+)', line)
                    if match:
                        self.csq_fields = match.group(1).split('|')
                        print(f"✅ Found {len(self.csq_fields)} CSQ fields")
                        
                        # Create field position mapping
                        for i, field in enumerate(self.csq_fields):
                            self.csq_field_mapping[field] = i
                        
                        return True
                elif line.startswith('#CHROM'):
                    break
        
        raise ValueError("CSQ header not found in VCF file")
    
    def load_vcf_data(self, max_variants=None):
        """Load and parse VCF data."""
        print("📊 Loading VCF data...")
        
        variant_count = 0
        annotated_count = 0
        
        with open(self.vcf_file, 'r') as f:
            for line in f:
                # Skip header lines
                if line.startswith('#'):
                    continue
                
                # Parse variant line
                fields = line.strip().split('\t')
                if len(fields) >= 8:
                    chrom, pos, var_id, ref, alt, qual, filter_val, info = fields[:8]
                    
                    # Extract CSQ annotation from INFO field
                    csq_match = re.search(r'CSQ=([^;]+)', info)
                    if csq_match:
                        csq_string = csq_match.group(1)
                        
                        # Process multiple transcripts (comma-separated)
                        transcripts = csq_string.split(',')
                        
                        # Use first transcript (usually canonical with --pick)
                        if transcripts:
                            csq_values = transcripts[0].split('|')
                            
                            # Create variant record
                            variant_record = {
                                'chromosome': chrom,
                                'position': int(pos),
                                'variant_id': var_id,
                                'reference_allele': ref,
                                'alternate_allele': alt,
                                'quality_score': qual,
                                'filter_status': filter_val,
                                'csq_values': csq_values,
                                'raw_csq': csq_string
                            }
                            
                            self.vcf_data.append(variant_record)
                            annotated_count += 1
                    
                    variant_count += 1
                    
                    # Progress indicator
                    if variant_count % 50000 == 0:
                        print(f"📈 Processed {variant_count:,} variants, {annotated_count:,} annotated...")
                    
                    # Optional limit for testing
                    if max_variants and variant_count >= max_variants:
                        break
        
        print(f"✅ Loaded {variant_count:,} variants, {annotated_count:,} with CSQ annotations")
        return variant_count, annotated_count
    
    def get_clinical_relevance(self, field_name):
        """Determine clinical relevance of a CSQ field."""
        field_lower = field_name.lower()
        
        # Essential fields for ML training
        essential_fields = [
            'consequence', 'gene', 'clin_sig', 'sift', 'polyphen',
            'cds_position', 'protein_position', 'amino_acids'
        ]
        
        # High importance fields
        high_fields = [
            'impact', 'symbol', 'exon', 'intron', 'canonical',
            'domains', 'transcription_factors'
        ]
        
        # Artifact fields
        artifact_fields = [
            'variant_class', 'symbol_source', 'hgnc_id', 'flags'
        ]
        
        if any(essential in field_lower for essential in essential_fields):
            return 'ESSENTIAL'
        elif any(high in field_lower for high in high_fields):
            return 'HIGH'
        elif any(artifact in field_lower for artifact in artifact_fields):
            return 'ARTIFACT'
        else:
            return 'LOW'
    
    def detect_concatenation(self, values):
        """Detect concatenation patterns in field values."""
        concatenation_stats = {}
        
        for pattern_name, pattern_regex in self.concat_patterns.items():
            matches = []
            for value in values:
                if value and isinstance(value, str) and re.search(pattern_regex, value):
                    matches.append(value)
            
            if matches:
                concatenation_stats[pattern_name] = {
                    'count': len(matches),
                    'percentage': (len(matches) / len(values)) * 100,
                    'examples': matches[:3]  # Store first 3 examples
                }
        
        return concatenation_stats
    
    def analyze_csq_field(self, field_name, field_index):
        """Analyze a specific CSQ field for data quality issues."""
        print(f"Analyzing CSQ field {field_index+1}/{len(self.csq_fields)}: {field_name}")
        
        # Extract all values for this field
        field_values = []
        non_null_values = []
        
        for variant in self.vcf_data:
            if field_index < len(variant['csq_values']):
                value = variant['csq_values'][field_index]
                field_values.append(value)
                
                if value and value != '.' and value != '':
                    non_null_values.append(value)
            else:
                field_values.append(None)
        
        # Basic statistics
        total_values = len(field_values)
        non_null_count = len(non_null_values)
        null_count = total_values - non_null_count
        unique_values = len(set(non_null_values)) if non_null_values else 0
        
        # Data type detection
        data_type = self._detect_data_type(non_null_values)
        
        # Concatenation analysis
        concatenation_patterns = self.detect_concatenation(non_null_values)
        has_concatenation = bool(concatenation_patterns)
        
        # Value length statistics
        value_lengths = [len(str(v)) for v in non_null_values if v]
        length_stats = {}
        if value_lengths:
            length_stats = {
                'min': min(value_lengths),
                'max': max(value_lengths),
                'mean': np.mean(value_lengths),
                'median': np.median(value_lengths)
            }
        
        # Clinical relevance
        clinical_relevance = self.get_clinical_relevance(field_name)
        clinical_score = self.clinical_scores[clinical_relevance]
        
        # Pipeline source attribution
        pipeline_source = self._determine_pipeline_source(field_name, has_concatenation)
        
        # ML readiness
        ml_ready = self._is_ml_ready(field_name, has_concatenation, clinical_relevance)
        
        # Sample values
        sample_values = non_null_values[:5] if non_null_values else []
        
        return {
            'analysis': {
                'total_values': total_values,
                'non_null_values': str(non_null_count),
                'null_values': str(null_count),
                'unique_values': unique_values,
                'data_type': data_type,
                'concatenation': concatenation_patterns,
                'value_lengths': length_stats,
                'sample_values': [str(v)[:100] for v in sample_values]  # Truncate long values
            },
            'clinical_relevance': clinical_relevance,
            'clinical_score': clinical_score,
            'pipeline_source': pipeline_source,
            'has_concatenation': has_concatenation,
            'ml_ready': ml_ready
        }
    
    def _detect_data_type(self, values):
        """Detect the data type of field values."""
        if not values:
            return 'unknown'
        
        # Check if all values are numeric
        numeric_count = 0
        for value in values[:100]:  # Sample first 100 values
            try:
                float(value)
                numeric_count += 1
            except (ValueError, TypeError):
                pass
        
        if numeric_count > len(values[:100]) * 0.8:
            return 'numeric'
        else:
            return 'object'
    
    def _determine_pipeline_source(self, field_name, has_concatenation):
        """Determine which pipeline stage introduced the field."""
        field_lower = field_name.lower()
        
        # VEP-specific fields
        vep_fields = [
            'consequence', 'impact', 'symbol', 'gene', 'feature_type',
            'feature', 'biotype', 'exon', 'intron', 'hgvsc', 'hgvsp',
            'cdna_position', 'cds_position', 'protein_position', 'amino_acids',
            'codons', 'existing_variation', 'distance', 'strand', 'flags',
            'symbol_source', 'hgnc_id', 'canonical', 'ccds', 'ensp',
            'swissprot', 'trembl', 'uniparc', 'uniprot_isoform',
            'sift', 'polyphen', 'domains', 'clin_sig', 'somatic',
            'pubmed', 'var_synonyms', 'af', 'afr_af', 'amr_af', 'eas_af',
            'eur_af', 'sas_af', 'gnomade_af', 'gnomade_afr_af',
            'gnomade_amr_af', 'gnomade_asj_af', 'gnomade_eas_af',
            'gnomade_fin_af', 'gnomade_mid_af', 'gnomade_nfe_af',
            'gnomade_remaining_af', 'gnomade_sas_af'
        ]
        
        if any(vep_field in field_lower for vep_field in vep_fields):
            if has_concatenation:
                return 'VEP_ANNOTATION'
            else:
                return 'VEP_CLEAN'
        else:
            return 'VEP_ANNOTATION'
    
    def _is_ml_ready(self, field_name, has_concatenation, clinical_relevance):
        """Determine if field is ready for ML training."""
        # Field is ML-ready if it's clinically relevant and has no concatenation
        return (clinical_relevance in ['ESSENTIAL', 'HIGH'] and 
                not has_concatenation and 
                clinical_relevance != 'ARTIFACT')
    
    def run_comprehensive_audit(self):
        """Execute the complete VCF data quality audit."""
        print("\n🔍 STARTING COMPREHENSIVE VCF AUDIT")
        print("=" * 50)
        
        # Parse VCF header
        self.parse_vcf_header()
        
        # Load VCF data
        total_variants, annotated_variants = self.load_vcf_data()
        
        # Analyze each CSQ field
        for i, field_name in enumerate(self.csq_fields):
            field_analysis = self.analyze_csq_field(field_name, i)
            self.audit_results['column_analysis'][field_name] = field_analysis
        
        # Generate summary statistics
        self._generate_summary_stats(total_variants, annotated_variants)
        
        # Generate concatenation analysis
        self._generate_concatenation_analysis()
        
        # Generate clinical scoring
        self._generate_clinical_scoring()
        
        # Generate pipeline attribution
        self._generate_pipeline_attribution()
        
        # Generate recommendations
        self._generate_recommendations()
        
        print(f"\n✅ Audit completed: {self._count_concatenated_fields()}/{len(self.csq_fields)} fields have concatenation issues")
    
    def _generate_summary_stats(self, total_variants, annotated_variants):
        """Generate overall summary statistics."""
        concatenated_fields = self._count_concatenated_fields()
        
        self.audit_results['summary'] = {
            'total_columns': len(self.csq_fields),
            'concatenated_columns': concatenated_fields,
            'concatenated_percentage': (concatenated_fields / len(self.csq_fields)) * 100,
            'total_variants': total_variants,
            'annotated_variants': annotated_variants,
            'annotation_rate': (annotated_variants / total_variants) * 100 if total_variants > 0 else 0,
            'audit_timestamp': datetime.now().isoformat()
        }
    
    def _count_concatenated_fields(self):
        """Count fields with concatenation issues."""
        count = 0
        for field_data in self.audit_results['column_analysis'].values():
            if field_data['has_concatenation']:
                count += 1
        return count
    
    def _generate_concatenation_analysis(self):
        """Generate concatenation analysis."""
        patterns_by_frequency = Counter()
        sources_with_concatenation = Counter()
        
        for field_name, field_data in self.audit_results['column_analysis'].items():
            if field_data['has_concatenation']:
                # Count concatenation patterns
                for pattern in field_data['analysis']['concatenation'].keys():
                    patterns_by_frequency[pattern] += 1
                
                # Count sources with concatenation
                sources_with_concatenation[field_data['pipeline_source']] += 1
        
        self.audit_results['concatenation_analysis'] = {
            'patterns_by_frequency': dict(patterns_by_frequency),
            'sources_with_concatenation': dict(sources_with_concatenation),
            'most_problematic_patterns': patterns_by_frequency.most_common()
        }
    
    def _generate_clinical_scoring(self):
        """Generate clinical relevance scoring analysis."""
        ranking = defaultdict(list)
        essential_issues = 0
        ml_ready_count = 0
        
        for field_name, field_data in self.audit_results['column_analysis'].items():
            relevance = field_data['clinical_relevance']
            ranking[relevance].append({
                'column': field_name,
                'has_concatenation': field_data['has_concatenation'],
                'ml_ready': field_data['ml_ready']
            })
            
            if relevance == 'ESSENTIAL' and field_data['has_concatenation']:
                essential_issues += 1
            
            if field_data['ml_ready']:
                ml_ready_count += 1
        
        self.audit_results['clinical_scoring'] = {
            'ranking_by_importance': dict(ranking),
            'essential_columns_with_issues': essential_issues,
            'total_ml_ready_columns': ml_ready_count
        }
    
    def _generate_pipeline_attribution(self):
        """Generate pipeline source attribution analysis."""
        issues_by_source = defaultdict(list)
        source_counts = Counter()
        
        for field_name, field_data in self.audit_results['column_analysis'].items():
            source = field_data['pipeline_source']
            if field_data['has_concatenation']:
                issues_by_source[source].append(field_name)
            source_counts[source] += 1
        
        self.audit_results['pipeline_attribution'] = {
            'issues_by_source': dict(issues_by_source),
            'vep_introduced_issues': len(issues_by_source.get('VEP_ANNOTATION', [])),
            'source_dataset_issues': len(issues_by_source.get('SOURCE_DATASET', [])),
            'merge_stage_issues': len(issues_by_source.get('MERGE_STAGE', []))
        }
    
    def _generate_recommendations(self):
        """Generate actionable recommendations."""
        recommendations = []
        
        # Critical recommendation for essential fields
        essential_issues = []
        for field_name, field_data in self.audit_results['column_analysis'].items():
            if (field_data['clinical_relevance'] == 'ESSENTIAL' and 
                field_data['has_concatenation']):
                essential_issues.append(field_name)
        
        if essential_issues:
            recommendations.append({
                'priority': 'CRITICAL',
                'action': 'Fix concatenation in essential fields',
                'columns': essential_issues,
                'rationale': 'These fields are critical for ML training but have concatenation issues'
            })
        
        # High recommendation for VEP issues
        vep_issues = self.audit_results['pipeline_attribution']['issues_by_source'].get('VEP_ANNOTATION', [])
        if len(vep_issues) >= 5:
            recommendations.append({
                'priority': 'HIGH',
                'action': 'Optimize VEP parameters to reduce concatenation',
                'columns': vep_issues[:10],  # Show first 10
                'rationale': f'VEP introduced concatenation in {len(vep_issues)} fields'
            })
        
        self.audit_results['recommendations'] = recommendations
    
    def save_results(self, output_dir):
        """Save audit results to multiple formats."""
        os.makedirs(output_dir, exist_ok=True)
        
        # Save JSON results
        json_file = os.path.join(output_dir, 'comprehensive_audit_results.json')
        with open(json_file, 'w') as f:
            json.dump(self.audit_results, f, indent=2, default=str)
        
        # Save text report
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
            f.write("PROSTATE CANCER VCF DATA QUALITY AUDIT REPORT\n")
            f.write("=" * 50 + "\n\n")
            
            # Summary
            summary = self.audit_results['summary']
            f.write("AUDIT SUMMARY\n")
            f.write(f"Timestamp: {summary['audit_timestamp']}\n")
            f.write(f"Total CSQ fields analyzed: {summary['total_columns']}\n")
            f.write(f"Fields with concatenation: {summary['concatenated_columns']} ({summary['concatenated_percentage']:.1f}%)\n")
            f.write(f"Total variants: {summary['total_variants']:,}\n")
            f.write(f"Annotated variants: {summary['annotated_variants']:,}\n")
            f.write(f"Annotation rate: {summary['annotation_rate']:.1f}%\n\n")
            
            # Recommendations
            f.write("PRIORITY RECOMMENDATIONS\n")
            f.write("-" * 25 + "\n")
            for i, rec in enumerate(self.audit_results['recommendations'], 1):
                f.write(f"{i}. [{rec['priority']}] {rec['action']}\n")
                f.write(f"   Rationale: {rec['rationale']}\n")
                f.write(f"   Affected fields: {len(rec['columns'])}\n\n")
            
            # Pipeline attribution
            attr = self.audit_results['pipeline_attribution']
            f.write("PIPELINE SOURCE ATTRIBUTION\n")
            f.write("-" * 30 + "\n")
            for source, fields in attr['issues_by_source'].items():
                f.write(f"{source}: {len(fields)} fields\n")
            f.write("\n")
            
            # Top problematic fields
            f.write("TOP 10 MOST PROBLEMATIC FIELDS\n")
            f.write("-" * 35 + "\n")
            
            # Sort fields by concatenation severity
            problematic_fields = []
            for field_name, data in self.audit_results['column_analysis'].items():
                if data['has_concatenation']:
                    total_concat = sum(info['count'] for info in data['analysis']['concatenation'].values())
                    problematic_fields.append((field_name, data, total_concat))
            
            problematic_fields.sort(key=lambda x: x[2], reverse=True)
            
            for i, (field_name, data, total_concat) in enumerate(problematic_fields[:10], 1):
                f.write(f" {i:2d}. {field_name} ({data['clinical_relevance']})\n")
                for pattern, info in data['analysis']['concatenation'].items():
                    f.write(f"     {pattern}: {info['count']} values ({info['percentage']:.1f}%)\n")
                f.write("\n")
    
    def _generate_csv_summary(self, filename):
        """Generate CSV summary for easy analysis."""
        rows = []
        for field_name, data in self.audit_results['column_analysis'].items():
            row = {
                'column_name': field_name,
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
    # Get file paths from environment variables (set by bash script)
    vcf_file = os.environ.get('VCF_INPUT_FILE', 
                             "/u/aa107/uiuc-cancer-research/data/processed/vep/vep_annotated.vcf")
    enhanced_csv = os.environ.get('ENHANCED_CSV_FILE', None)
    analysis_mode = os.environ.get('ANALYSIS_MODE', 'VCF')
    
    output_dir = "/u/aa107/uiuc-cancer-research/results/validation/phase1_audit"
    
    print("🧬 PROSTATE CANCER VCF DATA QUALITY AUDITOR")
    print("Phase 1: Comprehensive CSQ Field Analysis")
    print("=" * 50)
    print(f"📁 VCF Input: {vcf_file}")
    print(f"📊 Analysis Mode: {analysis_mode}")
    print("")
    
    # Initialize auditor
    auditor = ProstateCancerVCFAuditor(vcf_file, enhanced_csv)
    
    # Run comprehensive audit
    auditor.run_comprehensive_audit()
    
    # Save results
    auditor.save_results(output_dir)
    
    print("\n🎯 VCF AUDIT COMPLETED SUCCESSFULLY!")
    print(f"Results saved to: {output_dir}")
    print("\nNext steps:")
    print("1. Review generated reports")
    print("2. Compare VCF-level vs CSV-level concatenation")
    print("3. Identify VCF-to-CSV conversion issues")
    print("4. Implement targeted fixes")


if __name__ == "__main__":
    main()