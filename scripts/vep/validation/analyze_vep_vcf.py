#!/usr/bin/env python3
"""
Comprehensive VEP VCF Analysis Script
Analyzes VEP annotated VCF file to understand available features for TabNet training

Location: /u/aa107/uiuc-cancer-research/scripts/vep/validation/analyze_vep_vcf.py
Usage: python analyze_vep_vcf.py
"""

import re
import gzip
from collections import Counter, defaultdict
from pathlib import Path

def open_vcf_file(vcf_path):
    """Open VCF file (handles both .vcf and .vcf.gz)"""
    if str(vcf_path).endswith('.gz'):
        return gzip.open(vcf_path, 'rt')
    else:
        return open(vcf_path, 'r')

def parse_csq_header(header_line):
    """Extract CSQ field names from VEP header"""
    # Example: ##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations...Format: Allele|Consequence|IMPACT|...">
    match = re.search(r'Format: ([^"]+)', header_line)
    if match:
        return match.group(1).split('|')
    return []

def analyze_csq_annotation(csq_string, csq_fields):
    """Parse a single CSQ annotation string"""
    if not csq_string or not csq_fields:
        return {}
    
    # Split CSQ annotation by commas (multiple transcripts)
    transcripts = csq_string.split(',')
    
    # Analyze first transcript (usually canonical)
    if transcripts:
        values = transcripts[0].split('|')
        # Create field->value mapping
        annotation = {}
        for i, field in enumerate(csq_fields):
            if i < len(values):
                annotation[field] = values[i] if values[i] else None
            else:
                annotation[field] = None
        return annotation
    return {}

def analyze_vcf_file(vcf_path, report_path=None):
    """Comprehensive VCF analysis"""
    
    # Setup output - both console and file
    def output_line(text=""):
        print(text)
        if report_path:
            with open(report_path, 'a', encoding='utf-8') as f:
                f.write(text + '\n')
    
    # Clear report file if it exists
    if report_path:
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write("")  # Clear file
    
    output_line(f"🔍 Analyzing VEP VCF: {vcf_path}")
    output_line("=" * 70)
    
    # Initialize counters
    total_variants = 0
    annotated_variants = 0
    csq_fields = []
    consequences = Counter()
    impacts = Counter()
    genes = Counter()
    functional_scores = defaultdict(int)
    
    # Sample annotations for inspection
    sample_annotations = []
    
    try:
        with open_vcf_file(vcf_path) as vcf:
            for line_num, line in enumerate(vcf, 1):
                line = line.strip()
                
                # Parse header
                if line.startswith('##'):
                    if 'ID=CSQ' in line:
                        csq_fields = parse_csq_header(line)
                        output_line(f"📋 CSQ Fields Found: {len(csq_fields)} fields")
                        output_line(f"🔗 CSQ Format: {' | '.join(csq_fields[:10])}{'...' if len(csq_fields) > 10 else ''}")
                        output_line()
                
                elif line.startswith('#CHROM'):
                    output_line(f"📊 Starting variant analysis...")
                    output_line()
                
                # Parse variants
                elif line and not line.startswith('#'):
                    total_variants += 1
                    
                    # Split VCF line
                    fields = line.split('\t')
                    if len(fields) >= 8:
                        info_field = fields[7]
                        
                        # Extract CSQ annotation
                        csq_match = re.search(r'CSQ=([^;]+)', info_field)
                        if csq_match:
                            annotated_variants += 1
                            csq_string = csq_match.group(1)
                            
                            # Parse annotation
                            annotation = analyze_csq_annotation(csq_string, csq_fields)
                            
                            if annotation:
                                # Collect statistics
                                if annotation.get('Consequence'):
                                    consequences[annotation['Consequence']] += 1
                                
                                if annotation.get('IMPACT'):
                                    impacts[annotation['IMPACT']] += 1
                                
                                if annotation.get('SYMBOL'):
                                    genes[annotation['SYMBOL']] += 1
                                
                                # Check for functional scores
                                for field, value in annotation.items():
                                    if value and value != '.':
                                        if any(score_type in field.upper() for score_type in ['SIFT', 'POLYPHEN', 'CADD', 'REVEL', 'METALR']):
                                            functional_scores[field] += 1
                                
                                # Collect samples
                                if len(sample_annotations) < 5:
                                    sample_annotations.append({
                                        'variant': f"{fields[0]}:{fields[1]} {fields[3]}>{fields[4]}",
                                        'annotation': annotation
                                    })
                    
                    # Progress indicator
                    if total_variants % 50000 == 0:
                        output_line(f"📈 Processed {total_variants:,} variants...")
                
                # Safety limit for very large files
                if total_variants >= 200000:
                    output_line(f"⚠️  Analysis limited to first {total_variants:,} variants")
                    break
    
    except Exception as e:
        output_line(f"❌ Error reading VCF file: {e}")
        return
    
    # Generate comprehensive report
    output_line("\n" + "=" * 70)
    output_line("📊 VEP ANNOTATION ANALYSIS REPORT")
    output_line("=" * 70)
    
    # Basic statistics
    output_line(f"🔢 BASIC STATISTICS:")
    output_line(f"   Total variants processed: {total_variants:,}")
    output_line(f"   Annotated variants: {annotated_variants:,}")
    output_line(f"   Annotation rate: {(annotated_variants/total_variants*100) if total_variants > 0 else 0:.1f}%")
    output_line(f"   CSQ fields available: {len(csq_fields)}")
    output_line()
    
    # Functional scores availability
    output_line(f"🧬 FUNCTIONAL PREDICTION SCORES:")
    if functional_scores:
        for score_type, count in sorted(functional_scores.items()):
            percentage = (count / annotated_variants * 100) if annotated_variants > 0 else 0
            output_line(f"   ✅ {score_type}: {count:,} variants ({percentage:.1f}%)")
    else:
        output_line("   ❌ No functional prediction scores detected")
    output_line()
    
    # Impact distribution
    output_line(f"🎯 VARIANT IMPACT DISTRIBUTION:")
    total_impacts = sum(impacts.values())
    for impact, count in impacts.most_common():
        percentage = (count / total_impacts * 100) if total_impacts > 0 else 0
        output_line(f"   {impact}: {count:,} ({percentage:.1f}%)")
    output_line()
    
    # Top consequences
    output_line(f"⚡ TOP VARIANT CONSEQUENCES:")
    for consequence, count in consequences.most_common(10):
        percentage = (count / annotated_variants * 100) if annotated_variants > 0 else 0
        output_line(f"   {consequence}: {count:,} ({percentage:.1f}%)")
    output_line()
    
    # Top genes
    output_line(f"🧬 TOP AFFECTED GENES:")
    for gene, count in genes.most_common(10):
        output_line(f"   {gene}: {count:,} variants")
    output_line()
    
    # Sample annotations
    output_line(f"🔬 SAMPLE ANNOTATIONS:")
    for i, sample in enumerate(sample_annotations, 1):
        output_line(f"\n   Sample {i}: {sample['variant']}")
        # Show key fields
        key_fields = ['Consequence', 'IMPACT', 'SYMBOL', 'Feature', 'HGVSp']
        for field in key_fields:
            value = sample['annotation'].get(field, 'N/A')
            if value:
                output_line(f"      {field}: {value}")
    output_line()
    
    # TabNet readiness assessment
    output_line("🤖 TABNET READINESS ASSESSMENT:")
    
    # Feature count estimation
    non_empty_fields = len([f for f in csq_fields if any(sample['annotation'].get(f) for sample in sample_annotations)])
    output_line(f"   Estimated usable features: ~{non_empty_fields}")
    
    # Data volume
    if annotated_variants >= 100000:
        output_line(f"   ✅ Data volume: Excellent ({annotated_variants:,} variants)")
    elif annotated_variants >= 50000:
        output_line(f"   ✅ Data volume: Good ({annotated_variants:,} variants)")
    else:
        output_line(f"   ⚠️  Data volume: Limited ({annotated_variants:,} variants)")
    
    # Functional scores
    if len(functional_scores) >= 5:
        output_line(f"   ✅ Functional scores: Rich ({len(functional_scores)} types)")
    elif len(functional_scores) >= 2:
        output_line(f"   ⚠️  Functional scores: Basic ({len(functional_scores)} types)")
    else:
        output_line(f"   ❌ Functional scores: Limited ({len(functional_scores)} types)")
    
    # Impact diversity
    if len(impacts) >= 3:
        output_line(f"   ✅ Impact diversity: Good ({len(impacts)} levels)")
    else:
        output_line(f"   ⚠️  Impact diversity: Limited ({len(impacts)} levels)")
    
    # Overall assessment
    output_line(f"\n🎯 OVERALL ASSESSMENT:")
    score = 0
    if annotated_variants >= 100000: score += 2
    elif annotated_variants >= 50000: score += 1
    
    if len(functional_scores) >= 5: score += 2
    elif len(functional_scores) >= 2: score += 1
    
    if non_empty_fields >= 60: score += 2
    elif non_empty_fields >= 40: score += 1
    
    if score >= 5:
        output_line("   🚀 EXCELLENT - Ready for TabNet training")
    elif score >= 3:
        output_line("   ✅ GOOD - Ready with minor enhancements")
    elif score >= 1:
        output_line("   ⚠️  FAIR - Usable but needs improvement")
    else:
        output_line("   ❌ POOR - Significant issues need resolution")
    
    output_line("\n" + "=" * 70)
    output_line("Analysis complete! 🎉")
    output_line("=" * 70)

def main():
    """Main analysis function"""
    print("🧬 VEP VCF Comprehensive Analysis Tool")
    print("📍 For TabNet Prostate Cancer Classification")
    print()
    
    # File paths
    vcf_path = Path("/u/aa107/uiuc-cancer-research/data/processed/vep/vep_annotated.vcf")
    report_path = Path("/u/aa107/uiuc-cancer-research/data/processed/vep/vep_analysis_report.txt")
    
    # Check if file exists
    if not vcf_path.exists():
        print(f"❌ VCF file not found: {vcf_path}")
        print("Please check the file path and try again.")
        return
    
    # Create output directory if needed
    report_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Get file size
    file_size_mb = vcf_path.stat().st_size / (1024 * 1024)
    print(f"📁 VCF Input: {vcf_path}")
    print(f"📏 Size: {file_size_mb:.1f} MB")
    print(f"📄 Report Output: {report_path}")
    print()
    
    # Run analysis
    analyze_vcf_file(vcf_path, report_path)
    
    print(f"\n✅ Analysis complete!")
    print(f"📄 Full report saved to: {report_path}")
    print(f"🔍 Review the report for detailed TabNet readiness assessment")

if __name__ == "__main__":
    main()