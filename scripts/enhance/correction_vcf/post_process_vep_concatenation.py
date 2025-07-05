#!/usr/bin/env python3
"""
VEP Post-Processing Concatenation Correction Script
Fixes concatenated fields in VEP annotated VCF files for TabNet ML training

Location: /u/aa107/uiuc-cancer-research/scripts/enhance/correction_vcf/post_process_vep_concatenation.py
Usage: python post_process_vep_concatenation.py

Author: PhD Research Student, University of Illinois
Contact: aa107@illinois.edu
"""

import os
import re
import sys
from datetime import datetime
from collections import defaultdict, Counter

# Field-specific cleaning rules
FIELD_CLEANING_RULES = {
    1: {  # Consequence - Use severity ranking
        'name': 'Consequence',
        'separators': ['&', '/'],
        'method': 'severity_ranking',
        'severity_table': 'CONSEQUENCE_SEVERITY'
    },
    33: {  # DOMAINS - Keep first domain only (most significant)
        'name': 'DOMAINS', 
        'separators': ['&'],
        'method': 'keep_first'
    },
    50: {  # CLIN_SIG - Use severity ranking
        'name': 'CLIN_SIG',
        'separators': ['&', '/'],
        'method': 'severity_ranking', 
        'severity_table': 'CLIN_SIG_SEVERITY'
    },
    53: {  # PUBMED - Count total publications
        'name': 'PUBMED',
        'separators': ['&'],
        'method': 'count_values'
    },
    54: {  # VAR_SYNONYMS - Keep shortest synonym
        'name': 'VAR_SYNONYMS',
        'separators': ['&'],
        'method': 'keep_shortest'
    }
}

# Severity ranking tables
CONSEQUENCE_SEVERITY = {
    'transcript_ablation': 10,
    'splice_acceptor_variant': 9,
    'splice_donor_variant': 9, 
    'stop_gained': 8,
    'frameshift_variant': 8,
    'stop_lost': 7,
    'start_lost': 7,
    'missense_variant': 4,
    'splice_region_variant': 3,
    'synonymous_variant': 2,
    'intron_variant': 1,
    'upstream_gene_variant': 0,
    'downstream_gene_variant': 0
}

CLIN_SIG_SEVERITY = {
    'pathogenic': 5,
    'likely_pathogenic': 4,
    'uncertain_significance': 3,
    'likely_benign': 2,
    'benign': 1,
    'not_provided': 0
}

class VEPConcatenationCleaner:
    """Cleans concatenated fields in VEP VCF files"""
    
    def __init__(self):
        self.stats = {
            'total_lines': 0,
            'processed_lines': 0,
            'cleaned_lines': 0,
            'field_changes': defaultdict(int),
            'errors': 0
        }
        
    def log_debug(self, message, level="INFO"):
        """Debug logging with timestamp"""
        timestamp = datetime.now().strftime("%H:%M:%S")
        print(f"[{timestamp}] {level}: {message}")
        
    def get_highest_severity_value(self, values, severity_table_name):
        """Select value with highest severity score"""
        severity_table = globals()[severity_table_name]
        
        max_severity = -1
        best_value = values[0] if values else ''
        
        self.log_debug(f"Evaluating severity for values: {values[:3]}{'...' if len(values) > 3 else ''}", "DEBUG")
        
        for value in values:
            severity = severity_table.get(value.lower(), -1)
            self.log_debug(f"  {value} -> severity: {severity}", "DEBUG")
            if severity > max_severity:
                max_severity = severity
                best_value = value
                
        self.log_debug(f"Selected highest severity: {best_value} (score: {max_severity})", "DEBUG")
        return best_value
    
    def apply_cleaning_method(self, value, rules, field_name):
        """Apply specific cleaning method based on field rules"""
        # Split by all possible separators
        parts = [value]
        for sep in rules['separators']:
            new_parts = []
            for part in parts:
                new_parts.extend(part.split(sep))
            parts = new_parts
        
        # Remove empty parts
        parts = [p.strip() for p in parts if p.strip()]
        
        if not parts:
            self.log_debug(f"No valid parts found for {field_name}: {value}", "WARN")
            return ''
        
        method = rules['method']
        original_count = len(parts)
        
        self.log_debug(f"Cleaning {field_name}: {original_count} concatenated values using {method}", "DEBUG")
        
        if method == 'severity_ranking':
            result = self.get_highest_severity_value(parts, rules['severity_table'])
        elif method == 'keep_first':
            result = parts[0]
            self.log_debug(f"Keeping first value: {result}", "DEBUG")
        elif method == 'count_values':
            result = str(len(parts))
            self.log_debug(f"Counting values: {original_count} -> {result}", "DEBUG")
        elif method == 'keep_shortest':
            result = min(parts, key=len)
            self.log_debug(f"Keeping shortest: {result} (length: {len(result)})", "DEBUG")
        else:
            result = value  # Fallback
            self.log_debug(f"Unknown method {method}, keeping original", "WARN")
        
        if result != value:
            self.stats['field_changes'][field_name] += 1
        
        return result
    
    def clean_concatenated_fields(self, fields):
        """Clean concatenated values in CSQ fields using universal strategy"""
        cleaned_fields = fields.copy()
        changes_made = False
        
        for position, rules in FIELD_CLEANING_RULES.items():
            if position < len(fields):
                original_value = fields[position]
                field_name = rules['name']
                
                # Check if field contains concatenation
                has_concatenation = any(sep in original_value for sep in rules['separators'])
                
                if has_concatenation and original_value.strip():
                    self.log_debug(f"Found concatenation in {field_name} (pos {position}): {original_value[:50]}{'...' if len(original_value) > 50 else ''}")
                    
                    cleaned_value = self.apply_cleaning_method(original_value, rules, field_name)
                    cleaned_fields[position] = cleaned_value
                    changes_made = True
                    
                    self.log_debug(f"  BEFORE: {original_value}")
                    self.log_debug(f"  AFTER:  {cleaned_value}")
                    
        return cleaned_fields, changes_made
    
    def process_vcf_line(self, line):
        """Process single VCF variant line"""
        if line.startswith('#'):
            return line, False  # Keep header lines unchanged
        
        cols = line.strip().split('\t')
        if len(cols) < 8:
            return line, False
        
        info_field = cols[7]  # Column 8 contains INFO
        
        # Extract CSQ field
        csq_match = re.search(r'CSQ=([^;]+)', info_field)
        if not csq_match:
            return line, False
        
        csq_data = csq_match.group(1)
        
        # Process multiple transcripts (comma-separated)
        transcripts = csq_data.split(',')
        cleaned_transcripts = []
        line_changed = False
        
        for i, transcript in enumerate(transcripts):
            fields = transcript.split('|')
            if len(fields) >= 60:  # Ensure complete CSQ annotation
                cleaned_fields, transcript_changed = self.clean_concatenated_fields(fields)
                cleaned_transcripts.append('|'.join(cleaned_fields))
                if transcript_changed:
                    line_changed = True
                    self.log_debug(f"Modified transcript {i+1}/{len(transcripts)}")
            else:
                cleaned_transcripts.append(transcript)
                self.log_debug(f"Skipping incomplete transcript {i+1} (only {len(fields)} fields)", "WARN")
        
        if line_changed:
            # Reconstruct INFO field
            new_csq = ','.join(cleaned_transcripts)
            new_info = re.sub(r'CSQ=[^;]+', f'CSQ={new_csq}', info_field)
            cols[7] = new_info
            return '\t'.join(cols), True
        
        return line, False
    
    def process_file(self, input_file, output_file):
        """Main file processing function"""
        self.log_debug("🧬 VEP Concatenation Post-Processing Started")
        self.log_debug(f"📁 Input: {input_file}")
        self.log_debug(f"📁 Output: {output_file}")
        
        # Validate input file
        if not os.path.exists(input_file):
            self.log_debug(f"❌ Input file not found: {input_file}", "ERROR")
            return False
        
        # Get file size for progress tracking
        file_size = os.path.getsize(input_file) / (1024 * 1024)  # MB
        self.log_debug(f"📏 File size: {file_size:.1f} MB")
        
        try:
            with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
                for line_num, line in enumerate(infile, 1):
                    self.stats['total_lines'] += 1
                    
                    # Progress reporting
                    if line_num % 10000 == 0:
                        self.log_debug(f"📊 Processed {line_num:,} lines... ({self.stats['cleaned_lines']:,} modified)")
                    
                    try:
                        cleaned_line, was_changed = self.process_vcf_line(line)
                        outfile.write(cleaned_line + '\n')
                        
                        self.stats['processed_lines'] += 1
                        if was_changed:
                            self.stats['cleaned_lines'] += 1
                            
                    except Exception as e:
                        self.log_debug(f"Error processing line {line_num}: {str(e)}", "ERROR")
                        outfile.write(line)  # Write original line on error
                        self.stats['errors'] += 1
            
            self.log_debug("✅ Processing complete!")
            self.print_summary()
            return True
            
        except Exception as e:
            self.log_debug(f"❌ Critical error: {str(e)}", "ERROR")
            return False
    
    def print_summary(self):
        """Print processing summary statistics"""
        self.log_debug("\n" + "=" * 60)
        self.log_debug("📈 PROCESSING SUMMARY")
        self.log_debug("=" * 60)
        self.log_debug(f"   📊 Total lines: {self.stats['total_lines']:,}")
        self.log_debug(f"   🔧 Lines processed: {self.stats['processed_lines']:,}")
        self.log_debug(f"   ✨ Lines modified: {self.stats['cleaned_lines']:,}")
        self.log_debug(f"   ⚠️  Errors: {self.stats['errors']:,}")
        
        if self.stats['processed_lines'] > 0:
            mod_rate = (self.stats['cleaned_lines'] / self.stats['processed_lines']) * 100
            self.log_debug(f"   📈 Modification rate: {mod_rate:.2f}%")
        
        self.log_debug("\n🔧 FIELD-SPECIFIC CHANGES:")
        for field_name, count in self.stats['field_changes'].items():
            self.log_debug(f"   {field_name}: {count:,} changes")
        
        self.log_debug("\n🎯 NEXT STEPS:")
        self.log_debug("   1. Run validation script to verify improvements")
        self.log_debug("   2. Check that concatenation rates dropped below 2%")
        self.log_debug("   3. Proceed with TabNet training if validation passes")
        self.log_debug("=" * 60)


def main():
    """Main execution function"""
    input_file = "/u/aa107/uiuc-cancer-research/data/processed/vep/vep_annotated.vcf"
    output_file = "/u/aa107/uiuc-cancer-research/data/processed/vep/vep_annotated_clean.vcf"
    
    # Create output directory if needed
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # Initialize cleaner and process file
    cleaner = VEPConcatenationCleaner()
    success = cleaner.process_file(input_file, output_file)
    
    if success:
        print("\n🎉 VEP POST-PROCESSING COMPLETED SUCCESSFULLY!")
        print(f"✅ Clean VCF saved to: {output_file}")
        return 0
    else:
        print("\n❌ VEP POST-PROCESSING FAILED!")
        return 1


if __name__ == "__main__":
    sys.exit(main())