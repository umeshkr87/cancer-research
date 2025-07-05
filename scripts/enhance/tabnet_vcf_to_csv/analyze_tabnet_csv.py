#!/usr/bin/env python3
"""
TabNet CSV Feature Analyzer
Analyzes the converted CSV file to understand feature distribution and TabNet readiness

Location: /u/aa107/uiuc-cancer-research/scripts/enhance/analyze_tabnet_csv.py
Usage: python analyze_tabnet_csv.py
"""

import pandas as pd
import numpy as np
from pathlib import Path
import argparse

def analyze_csv_features(csv_path, output_report=None):
    """Analyze CSV features for TabNet training readiness"""
    
    def output_line(text=""):
        print(text)
        if output_report:
            with open(output_report, 'a', encoding='utf-8') as f:
                f.write(text + '\n')
    
    # Clear report file
    if output_report:
        with open(output_report, 'w', encoding='utf-8') as f:
            f.write("")
    
    output_line("🔍 TabNet CSV Feature Analysis")
    output_line("=" * 50)
    
    try:
        # Load CSV
        output_line(f"📁 Loading CSV: {csv_path}")
        df = pd.read_csv(csv_path)
        output_line(f"✅ Loaded successfully!")
        output_line()
        
        # Basic statistics
        output_line("📊 BASIC DATASET STATISTICS")
        output_line("-" * 30)
        output_line(f"Total variants: {len(df):,}")
        output_line(f"Total features: {len(df.columns)}")
        output_line(f"Memory usage: {df.memory_usage(deep=True).sum() / (1024**2):.1f} MB")
        output_line()
        
        # Feature type analysis
        numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
        text_cols = df.select_dtypes(include=[object]).columns.tolist()
        
        output_line("🔢 FEATURE TYPE BREAKDOWN")
        output_line("-" * 30)
        output_line(f"Numeric features: {len(numeric_cols)}")
        output_line(f"Text features: {len(text_cols)}")
        output_line()
        
        # Missing data analysis
        output_line("❓ MISSING DATA ANALYSIS")
        output_line("-" * 30)
        missing_stats = df.isnull().sum()
        missing_pct = (missing_stats / len(df) * 100).round(1)
        
        # Show features with significant missing data
        high_missing = missing_pct[missing_pct > 20].sort_values(ascending=False)
        if len(high_missing) > 0:
            output_line("Features with >20% missing data:")
            for col, pct in high_missing.head(10).items():
                output_line(f"   {col}: {pct}%")
        else:
            output_line("✅ No features with >20% missing data")
        
        output_line(f"Overall missing data: {df.isnull().sum().sum():,} cells ({(df.isnull().sum().sum() / df.size * 100):.1f}%)")
        output_line()
        
        # Key feature analysis
        output_line("🧬 KEY FEATURE ANALYSIS")
        output_line("-" * 30)
        
        # Basic variant info
        if 'variant_type' in df.columns:
            output_line("Variant type distribution:")
            for vtype, count in df['variant_type'].value_counts().items():
                pct = count / len(df) * 100
                output_line(f"   {vtype}: {count:,} ({pct:.1f}%)")
            output_line()
        
        # Impact distribution
        if 'IMPACT' in df.columns:
            output_line("Variant impact distribution:")
            for impact, count in df['IMPACT'].value_counts().items():
                pct = count / len(df) * 100
                output_line(f"   {impact}: {count:,} ({pct:.1f}%)")
            output_line()
        
        # Functional scores availability
        output_line("🔬 FUNCTIONAL SCORES AVAILABILITY")
        output_line("-" * 30)
        
        functional_features = ['sift_score', 'sift_prediction', 'polyphen_score', 'polyphen_prediction']
        for feature in functional_features:
            if feature in df.columns:
                non_null = df[feature].notna().sum()
                pct = non_null / len(df) * 100
                output_line(f"   {feature}: {non_null:,} ({pct:.1f}%)")
        output_line()
        
        # Gene analysis
        if 'SYMBOL' in df.columns:
            output_line("🧬 TOP AFFECTED GENES")
            output_line("-" * 30)
            top_genes = df['SYMBOL'].value_counts().head(10)
            for gene, count in top_genes.items():
                if pd.notna(gene):
                    output_line(f"   {gene}: {count:,} variants")
            output_line()
        
        # Therapeutic pathway analysis
        output_line("💊 THERAPEUTIC PATHWAY INDICATORS")
        output_line("-" * 30)
        
        pathway_features = ['dna_repair_pathway', 'mismatch_repair_pathway', 'hormone_pathway', 'is_important_gene']
        for feature in pathway_features:
            if feature in df.columns:
                positive_count = (df[feature] == 1).sum()
                pct = positive_count / len(df) * 100
                output_line(f"   {feature}: {positive_count:,} variants ({pct:.1f}%)")
        output_line()
        
        # TabNet feature readiness
        output_line("🤖 TABNET TRAINING READINESS")
        output_line("-" * 30)
        
        # Count usable features (numeric or categorical with reasonable categories)
        usable_numeric = len([col for col in numeric_cols if df[col].notna().sum() > len(df) * 0.1])
        usable_categorical = len([col for col in text_cols if df[col].nunique() <= 50 and df[col].notna().sum() > len(df) * 0.1])
        total_usable = usable_numeric + usable_categorical
        
        output_line(f"Usable numeric features: {usable_numeric}")
        output_line(f"Usable categorical features: {usable_categorical}")
        output_line(f"Total usable features: {total_usable}")
        output_line()
        
        # Data quality score
        quality_score = 0
        
        # Volume score
        if len(df) >= 100000:
            quality_score += 2
            volume_status = "✅ Excellent"
        elif len(df) >= 50000:
            quality_score += 1
            volume_status = "✅ Good"
        else:
            volume_status = "⚠️  Limited"
        
        # Feature score
        if total_usable >= 60:
            quality_score += 2
            feature_status = "✅ Rich"
        elif total_usable >= 40:
            quality_score += 1
            feature_status = "✅ Adequate"
        else:
            feature_status = "⚠️  Limited"
        
        # Missing data score
        overall_missing_pct = df.isnull().sum().sum() / df.size * 100
        if overall_missing_pct < 10:
            quality_score += 2
            missing_status = "✅ Low"
        elif overall_missing_pct < 30:
            quality_score += 1
            missing_status = "⚠️  Moderate"
        else:
            missing_status = "❌ High"
        
        output_line("📋 QUALITY ASSESSMENT")
        output_line("-" * 30)
        output_line(f"Data volume: {volume_status} ({len(df):,} variants)")
        output_line(f"Feature richness: {feature_status} ({total_usable} usable features)")
        output_line(f"Missing data: {missing_status} ({overall_missing_pct:.1f}%)")
        output_line()
        
        # Overall readiness
        if quality_score >= 5:
            overall_status = "🚀 EXCELLENT - Ready for TabNet training"
        elif quality_score >= 3:
            overall_status = "✅ GOOD - Ready with minor preprocessing"
        elif quality_score >= 1:
            overall_status = "⚠️  FAIR - Needs improvement"
        else:
            overall_status = "❌ POOR - Major issues"
        
        output_line(f"🎯 OVERALL READINESS: {overall_status}")
        output_line()
        
        # Recommendations
        output_line("💡 RECOMMENDATIONS")
        output_line("-" * 30)
        
        if quality_score >= 5:
            output_line("✅ Proceed with TabNet training")
            output_line("✅ Configure 6-8 decision steps for interpretability")
            output_line("✅ Use attention mechanisms for clinical explainability")
        elif quality_score >= 3:
            output_line("🔧 Minor preprocessing recommended:")
            if overall_missing_pct > 10:
                output_line("   - Handle missing data (imputation or removal)")
            if total_usable < 50:
                output_line("   - Consider feature engineering for more predictors")
            output_line("✅ Then proceed with TabNet training")
        else:
            output_line("⚠️  Address data quality issues before training:")
            if len(df) < 50000:
                output_line("   - Increase dataset size if possible")
            if total_usable < 30:
                output_line("   - Engineer more features from available data")
            if overall_missing_pct > 50:
                output_line("   - Significant missing data cleanup needed")
        
        output_line()
        output_line("=" * 50)
        output_line("Analysis complete! 🎉")
        output_line("=" * 50)
        
        return True
        
    except Exception as e:
        output_line(f"❌ Error analyzing CSV: {e}")
        return False

def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='Analyze TabNet CSV features')
    parser.add_argument('--input', '-i',
                       default='/u/aa107/uiuc-cancer-research/data/processed/tabnet_csv/prostate_variants_tabnet.csv',
                       help='Input CSV file path')
    parser.add_argument('--report', '-r',
                       default='/u/aa107/uiuc-cancer-research/data/processed/tabnet_csv/feature_analysis_report.txt',
                       help='Output report file path')
    
    args = parser.parse_args()
    
    csv_path = Path(args.input)
    report_path = Path(args.report) if args.report else None
    
    if not csv_path.exists():
        print(f"❌ CSV file not found: {csv_path}")
        print("Please run VCF to CSV conversion first")
        return False
    
    # Create report directory
    if report_path:
        report_path.parent.mkdir(parents=True, exist_ok=True)
    
    print("🔍 Analyzing TabNet CSV features...")
    success = analyze_csv_features(csv_path, report_path)
    
    if success and report_path:
        print(f"📄 Detailed report saved to: {report_path}")
    
    return success

if __name__ == "__main__":
    main()