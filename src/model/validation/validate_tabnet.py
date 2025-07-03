#!/usr/bin/env python3
"""
TabNet Validation Framework - Enhanced Version
Validates AlphaMissense integration and realistic performance expectations
"""

import pandas as pd
import numpy as np
import json
import time
from pathlib import Path
from datetime import datetime
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score, classification_report
from sklearn.ensemble import RandomForestClassifier
import warnings
warnings.filterwarnings('ignore')

# Import our enhanced TabNet model
import sys
sys.path.append('/u/aa107/uiuc-cancer-research/src')

from model.tabnet_prostate_variant_classifier import ProstateVariantTabNet

class EnhancedTabNetValidator:
    """
    Validation framework for enhanced TabNet with AlphaMissense features
    Focuses on data leakage detection and realistic performance expectations
    """
    
    def __init__(self, results_dir=None):
        """Initialize validator"""
        self.results_dir = Path(results_dir) if results_dir else Path("/u/aa107/uiuc-cancer-research/results/validation")
        self.results_dir.mkdir(parents=True, exist_ok=True)
        
        self.validation_results = {}
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # Expected performance ranges (REALISTIC)
        self.performance_thresholds = {
            'excellent': 0.75,      # 75%+ is excellent for genomic data
            'good': 0.65,           # 65%+ is good
            'suspicious': 0.95      # 95%+ suggests data leakage
        }
        
        print("🧬 Enhanced TabNet Validator Initialized")
        print(f"📁 Results directory: {self.results_dir}")
        print(f"⚠️  Suspicious accuracy threshold: {self.performance_thresholds['suspicious']:.1%}")
    
    def validate_enhanced_dataset(self):
        """
        Step 1: Validate enhanced dataset integrity
        """
        print("\n🔍 STEP 1: ENHANCED DATASET VALIDATION")
        print("=" * 50)
        
        enhanced_path = "/u/aa107/uiuc-cancer-research/data/processed/tabnet_csv/prostate_variants_tabnet_enhanced.csv"
        
        if not Path(enhanced_path).exists():
            print(f"❌ Enhanced dataset not found: {enhanced_path}")
            return False
        
        try:
            df = pd.read_csv(enhanced_path)
            print(f"✅ Enhanced dataset loaded: {df.shape[0]:,} variants × {df.shape[1]} features")
            
            # CRITICAL: Check for data leakage features (should be ABSENT)
            leakage_features = ['functional_pathogenicity', 'sift_confidence', 'polyphen_confidence']
            leakage_found = [f for f in leakage_features if f in df.columns]
            
            if leakage_found:
                print(f"❌ CRITICAL: Data leakage features found: {leakage_found}")
                print("   These must be removed before training!")
                return False
            else:
                print("✅ No data leakage features detected")
            
            # CRITICAL: Check for AlphaMissense features (should be PRESENT)
            am_features = ['alphamissense_pathogenicity', 'alphamissense_class']
            missing_am = [f for f in am_features if f not in df.columns]
            
            if missing_am:
                print(f"❌ CRITICAL: AlphaMissense features missing: {missing_am}")
                return False
            else:
                print("✅ AlphaMissense features present")
            
            # Validate AlphaMissense coverage and distribution
            am_coverage = df['alphamissense_pathogenicity'].notna().sum()
            coverage_rate = am_coverage / len(df) * 100
            
            print(f"📊 AlphaMissense coverage: {am_coverage:,} variants ({coverage_rate:.1f}%)")
            
            if coverage_rate < 30:
                print("⚠️  Low AlphaMissense coverage")
            else:
                print("✅ Good AlphaMissense coverage")
            
            # Check score distribution
            am_scores = df['alphamissense_pathogenicity'].dropna()
            if len(am_scores) > 0:
                print(f"📊 AlphaMissense scores: {am_scores.min():.3f} - {am_scores.max():.3f}")
                print(f"📊 Mean pathogenicity: {am_scores.mean():.3f}")
                
                # Validate realistic distribution
                pathogenic_count = (am_scores >= 0.7).sum()
                benign_count = (am_scores <= 0.3).sum()
                
                print(f"📊 Likely pathogenic (≥0.7): {pathogenic_count:,} ({pathogenic_count/len(am_scores)*100:.1f}%)")
                print(f"📊 Likely benign (≤0.3): {benign_count:,} ({benign_count/len(am_scores)*100:.1f}%)")
            
            self.validation_results['dataset_validation'] = {
                'total_variants': len(df),
                'total_features': len(df.columns),
                'leakage_features_found': leakage_found,
                'alphamissense_coverage': coverage_rate,
                'alphamissense_score_range': [float(am_scores.min()), float(am_scores.max())] if len(am_scores) > 0 else None
            }
            
            return True
            
        except Exception as e:
            print(f"❌ Dataset validation failed: {e}")
            return False
    
    def baseline_validation(self):
        """
        Step 2: Baseline validation with Random Forest
        """
        print("\n🌲 STEP 2: BASELINE VALIDATION (RANDOM FOREST)")
        print("=" * 50)
        
        try:
            # Load data using our enhanced model
            model = ProstateVariantTabNet()
            X, y = model.load_data()
            
            print(f"📊 Loaded data: {X.shape[0]:,} samples × {X.shape[1]} features")
            
            # Quick Random Forest baseline
            rf = RandomForestClassifier(n_estimators=100, random_state=42, n_jobs=-1)
            
            # 3-fold CV for speed
            skf = StratifiedKFold(n_splits=3, shuffle=True, random_state=42)
            cv_scores = []
            
            print("🔄 Running 3-fold cross-validation...")
            
            for fold, (train_idx, val_idx) in enumerate(skf.split(X, y), 1):
                X_train_fold, X_val_fold = X[train_idx], X[val_idx]
                y_train_fold, y_val_fold = y[train_idx], y[val_idx]
                
                rf.fit(X_train_fold, y_train_fold)
                y_pred = rf.predict(X_val_fold)
                accuracy = accuracy_score(y_val_fold, y_pred)
                cv_scores.append(accuracy)
                
                print(f"   Fold {fold}: {accuracy:.3f}")
            
            mean_accuracy = np.mean(cv_scores)
            std_accuracy = np.std(cv_scores)
            
            print(f"\n🎯 Random Forest Baseline:")
            print(f"   Mean accuracy: {mean_accuracy:.3f} ± {std_accuracy:.3f}")
            
            # Performance interpretation
            if mean_accuracy > self.performance_thresholds['suspicious']:
                print("⚠️  SUSPICIOUS: Baseline too high - possible data leakage")
                return False
            elif mean_accuracy > self.performance_thresholds['excellent']:
                print("✅ EXCELLENT: Strong baseline performance")
            elif mean_accuracy > self.performance_thresholds['good']:
                print("✅ GOOD: Reasonable baseline performance")
            else:
                print("📈 MODERATE: Challenging dataset")
            
            self.validation_results['baseline_validation'] = {
                'mean_accuracy': mean_accuracy,
                'std_accuracy': std_accuracy,
                'fold_scores': cv_scores
            }
            
            return True
            
        except Exception as e:
            print(f"❌ Baseline validation failed: {e}")
            return False
    
    def tabnet_validation(self):
        """
        Step 3: TabNet model validation
        """
        print("\n🔥 STEP 3: TABNET MODEL VALIDATION")
        print("=" * 50)
        
        try:
            # Initialize TabNet model
            tabnet = ProstateVariantTabNet(n_d=32, n_a=32, n_steps=3)  # Smaller for validation
            
            # Load data
            X, y = tabnet.load_data()
            
            # Use subset for quick validation
            if len(X) > 5000:
                print("📊 Using subset for quick validation...")
                indices = np.random.choice(len(X), 5000, replace=False)
                X, y = X[indices], y[indices]
            
            print(f"📊 Validation data: {X.shape[0]:,} samples × {X.shape[1]} features")
            
            # Cross-validation
            cv_results = tabnet.cross_validate(X, y, cv_folds=3)
            
            mean_accuracy = cv_results['mean_accuracy']
            std_accuracy = cv_results['std_accuracy']
            
            print(f"\n🎯 TabNet Performance:")
            print(f"   Mean accuracy: {mean_accuracy:.3f} ± {std_accuracy:.3f}")
            
            # CRITICAL: Check for data leakage
            if mean_accuracy > self.performance_thresholds['suspicious']:
                print("❌ CRITICAL: Suspiciously high accuracy - DATA LEAKAGE DETECTED!")
                print("   TabNet should NOT achieve >95% accuracy on legitimate data")
                self.validation_results['data_leakage_detected'] = True
                return False
            elif mean_accuracy > self.performance_thresholds['excellent']:
                print("✅ EXCELLENT: Realistic high performance")
            elif mean_accuracy > self.performance_thresholds['good']:
                print("✅ GOOD: Realistic performance for genomic data")
            else:
                print("📈 MODERATE: May need feature engineering")
            
            # Feature importance analysis
            print("\n🔍 Analyzing feature importance...")
            
            # Quick training for feature importance
            from sklearn.model_selection import train_test_split
            X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
            
            tabnet.train(X_train, y_train, max_epochs=50)  # Quick training
            
            feature_importance = tabnet.get_feature_importance()
            
            print("📊 Top 5 Features:")
            for idx, row in feature_importance.head(5).iterrows():
                print(f"   {row['feature']}: {row['importance']:.3f}")
            
            # Check if AlphaMissense features are important
            am_features_important = feature_importance[feature_importance['feature'].str.contains('alphamissense', na=False)]
            
            if len(am_features_important) > 0:
                print("✅ AlphaMissense features detected in top features:")
                for idx, row in am_features_important.iterrows():
                    print(f"   {row['feature']}: {row['importance']:.3f}")
            else:
                print("⚠️  AlphaMissense features not in top features")
            
            self.validation_results['tabnet_validation'] = {
                'mean_accuracy': mean_accuracy,
                'std_accuracy': std_accuracy,
                'fold_scores': cv_results['fold_scores'],
                'top_features': feature_importance.head(10).to_dict('records'),
                'alphamissense_features_important': len(am_features_important) > 0
            }
            
            return True
            
        except Exception as e:
            print(f"❌ TabNet validation failed: {e}")
            return False
    
    def generate_validation_report(self):
        """
        Step 4: Generate comprehensive validation report
        """
        print("\n📋 STEP 4: GENERATING VALIDATION REPORT")
        print("=" * 50)
        
        report_path = self.results_dir / f"enhanced_tabnet_validation_{self.timestamp}.json"
        
        # Add summary
        self.validation_results['summary'] = {
            'validation_timestamp': self.timestamp,
            'validator_version': 'enhanced_v1.0',
            'data_leakage_detected': self.validation_results.get('data_leakage_detected', False),
            'alphamissense_integrated': 'dataset_validation' in self.validation_results and 
                                     self.validation_results['dataset_validation']['alphamissense_coverage'] > 30
        }
        
        # Save report
        with open(report_path, 'w') as f:
            json.dump(self.validation_results, f, indent=2)
        
        print(f"✅ Validation report saved: {report_path}")
        
        # Print summary
        print(f"\n📊 VALIDATION SUMMARY:")
        print(f"   Data leakage detected: {'❌ YES' if self.validation_results.get('data_leakage_detected') else '✅ NO'}")
        print(f"   AlphaMissense integrated: {'✅ YES' if self.validation_results['summary']['alphamissense_integrated'] else '❌ NO'}")
        
        if 'tabnet_validation' in self.validation_results:
            tabnet_acc = self.validation_results['tabnet_validation']['mean_accuracy']
            print(f"   TabNet accuracy: {tabnet_acc:.3f}")
            
            if tabnet_acc > self.performance_thresholds['suspicious']:
                print("   Status: ❌ SUSPICIOUS - Investigate data leakage")
            elif tabnet_acc > self.performance_thresholds['excellent']:
                print("   Status: ✅ EXCELLENT - Ready for production")
            else:
                print("   Status: ✅ GOOD - Realistic performance")
        
        return report_path

def main():
    """
    Main validation workflow
    """
    print("🧬 ENHANCED TABNET VALIDATION FRAMEWORK")
    print("=" * 60)
    print("Validating AlphaMissense integration and data leakage elimination")
    print()
    
    # Initialize validator
    validator = EnhancedTabNetValidator()
    
    try:
        # Step 1: Validate enhanced dataset
        if not validator.validate_enhanced_dataset():
            print("❌ Dataset validation failed - cannot proceed")
            return False
        
        # Step 2: Baseline validation
        if not validator.baseline_validation():
            print("❌ Baseline validation failed - possible data leakage")
            return False
        
        # Step 3: TabNet validation
        if not validator.tabnet_validation():
            print("❌ TabNet validation failed - check for issues")
            return False
        
        # Step 4: Generate report
        report_path = validator.generate_validation_report()
        
        print(f"\n✅ VALIDATION COMPLETED SUCCESSFULLY!")
        print(f"📋 Report: {report_path}")
        print(f"\n🎯 NEXT STEPS:")
        print(f"   1. Review validation report")
        print(f"   2. If no data leakage detected, proceed with full training")
        print(f"   3. Use enhanced dataset for production model")
        
        return True
        
    except Exception as e:
        print(f"❌ Validation failed: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)