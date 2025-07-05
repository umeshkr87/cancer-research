#!/usr/bin/env python3
"""
Enhanced TabNet Validation Framework
Fixed data leakage issue - uses same logic as working test_environment.py
"""

import pandas as pd
import numpy as np
import sys
import json
import traceback
from pathlib import Path
from datetime import datetime
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.preprocessing import StandardScaler

# Add project root to path
sys.path.append('/u/aa107/uiuc-cancer-research/src')

class EnhancedTabNetValidator:
    def __init__(self):
        self.project_dir = Path("/u/aa107/uiuc-cancer-research")
        self.results_dir = self.project_dir / "results/validation"
        self.results_dir.mkdir(parents=True, exist_ok=True)
        
        # Performance thresholds
        self.performance_thresholds = {
            'excellent': 0.85,
            'good': 0.75,
            'acceptable': 0.65,
            'suspicious': 0.95  # Above this indicates data leakage
        }
        
        self.validation_results = {
            'timestamp': datetime.now().isoformat(),
            'data_leakage_detected': False,
            'alphamissense_integrated': False,
            'summary': {}
        }
    
    def enhanced_dataset_validation(self):
        """
        Step 1: Enhanced dataset validation using SAME logic as test_environment.py
        """
        print("\n🔍 STEP 1: ENHANCED DATASET VALIDATION")
        print("=" * 50)
        
        try:
            # Use clean dataset path
            clean_path = "/u/aa107/uiuc-cancer-research/data/processed/tabnet_csv/prostate_variants_tabnet_clean.csv"
            
            if not Path(clean_path).exists():
                print(f"❌ Clean dataset not found: {clean_path}")
                return False
            
            # Load dataset
            df = pd.read_csv(clean_path, low_memory=False)
            print(f"✅ Enhanced dataset loaded: {df.shape[0]:,} variants × {df.shape[1]} features")
            
            # CRITICAL: Check for data leakage features
            leakage_features = ['sift_prediction', 'polyphen_prediction', 'functional_pathogenicity']
            leakage_found = [f for f in leakage_features if f in df.columns]
            
            if leakage_found:
                print(f"❌ DATA LEAKAGE DETECTED: {leakage_found}")
                self.validation_results['data_leakage_detected'] = True
                return False
            else:
                print("✅ No data leakage features detected")
            
            # Check AlphaMissense integration
            am_features = ['alphamissense_pathogenicity', 'alphamissense_class']
            missing_am = [f for f in am_features if f not in df.columns]
            
            if missing_am:
                print(f"❌ AlphaMissense features missing: {missing_am}")
                return False
            else:
                print("✅ AlphaMissense features present")
                am_coverage = df['alphamissense_pathogenicity'].notna().sum()
                coverage_pct = am_coverage / len(df) * 100
                print(f"📊 AlphaMissense coverage: {am_coverage:,} variants ({coverage_pct:.1f}%)")
                
                if coverage_pct >= 30:
                    print("✅ Good AlphaMissense coverage")
                    self.validation_results['alphamissense_integrated'] = True
                else:
                    print("⚠️  Low AlphaMissense coverage")
                
                # Show score distribution
                am_scores = df['alphamissense_pathogenicity'].dropna()
                print(f"📊 AlphaMissense scores: {am_scores.min():.3f} - {am_scores.max():.3f}")
                print(f"📊 Mean pathogenicity: {am_scores.mean():.3f}")
                
                # Show pathogenicity distribution
                high_path = (am_scores >= 0.7).sum()
                low_path = (am_scores <= 0.3).sum()
                print(f"📊 Likely pathogenic (≥0.7): {high_path:,} ({high_path/len(am_scores)*100:.1f}%)")
                print(f"📊 Likely benign (≤0.3): {low_path:,} ({low_path/len(am_scores)*100:.1f}%)")
            
            self.validation_results['dataset_validation'] = {
                'shape': df.shape,
                'data_leakage_detected': len(leakage_found) > 0,
                'alphamissense_coverage': coverage_pct
            }
            
            return True
            
        except Exception as e:
            print(f"❌ Dataset validation failed: {e}")
            traceback.print_exc()
            return False
    
    def load_data_with_proper_features(self):
        """
        Load data using SAME logic as working test_environment.py
        """
        # Import custom model to use its working logic
        try:
            from model.tabnet_prostate_variant_classifier import ProstateVariantTabNet
            
            # Initialize model
            model = ProstateVariantTabNet()
            
            # Use model's working data loading logic
            X, y = model.load_data()
            
            print(f"📊 Loaded data: {X.shape[0]:,} samples × {X.shape[1]} features")
            print(f"🎯 Target distribution:")
            target_counts = pd.Series(y).value_counts()
            for target, count in target_counts.items():
                pct = count / len(y) * 100
                print(f"   {target}: {count:,} ({pct:.1f}%)")
            
            return X, y, model.feature_names
            
        except Exception as e:
            print(f"❌ Data loading failed: {e}")
            traceback.print_exc()
            return None, None, None
    
    def baseline_validation(self):
        """
        Step 2: Baseline validation using SAME features as test_environment.py
        """
        print("\n🌲 STEP 2: BASELINE VALIDATION (RANDOM FOREST)")
        print("=" * 50)
        
        try:
            # Load data using working logic
            X, y, feature_names = self.load_data_with_proper_features()
            
            if X is None:
                print("❌ Could not load data")
                return False
            
            print(f"📊 Features used: {feature_names}")
            
            # Use subset for quick validation
            if len(X) > 10000:
                print("📊 Using subset for quick validation...")
                indices = np.random.choice(len(X), 10000, replace=False)
                X_subset = X[indices]
                y_subset = y[indices]
            else:
                X_subset = X
                y_subset = y
            
            # Cross-validation
            print(f"🔄 Running 3-fold cross-validation...")
            skf = StratifiedKFold(n_splits=3, shuffle=True, random_state=42)
            
            cv_scores = []
            for fold, (train_idx, val_idx) in enumerate(skf.split(X_subset, y_subset), 1):
                X_train, X_val = X_subset[train_idx], X_subset[val_idx]
                y_train, y_val = y_subset[train_idx], y_subset[val_idx]
                
                # Scale features
                scaler = StandardScaler()
                X_train_scaled = scaler.fit_transform(X_train)
                X_val_scaled = scaler.transform(X_val)
                
                # Train Random Forest
                rf = RandomForestClassifier(n_estimators=50, max_depth=10, random_state=42)
                rf.fit(X_train_scaled, y_train)
                
                # Predict
                y_pred = rf.predict(X_val_scaled)
                accuracy = accuracy_score(y_val, y_pred)
                cv_scores.append(accuracy)
                
                print(f"   Fold {fold}: {accuracy:.3f}")
            
            mean_accuracy = np.mean(cv_scores)
            std_accuracy = np.std(cv_scores)
            
            print(f"\n🎯 Random Forest Baseline:")
            print(f"   Mean accuracy: {mean_accuracy:.3f} ± {std_accuracy:.3f}")
            
            # Check for data leakage
            if mean_accuracy > self.performance_thresholds['suspicious']:
                print("❌ SUSPICIOUS: Baseline too high - possible data leakage")
                self.validation_results['data_leakage_detected'] = True
                return False
            elif mean_accuracy > self.performance_thresholds['good']:
                print("✅ EXCELLENT: Good baseline performance")
            elif mean_accuracy > self.performance_thresholds['acceptable']:
                print("✅ GOOD: Acceptable baseline performance")
            else:
                print("✅ MODERATE: Expected for complex genomic data")
            
            self.validation_results['baseline_validation'] = {
                'mean_accuracy': mean_accuracy,
                'std_accuracy': std_accuracy,
                'fold_scores': cv_scores
            }
            
            return True
            
        except Exception as e:
            print(f"❌ Baseline validation failed: {e}")
            traceback.print_exc()
            return False
    
    def tabnet_validation(self):
        """
        Step 3: TabNet model validation using SAME logic as test_environment.py
        """
        print("\n🔥 STEP 3: TABNET MODEL VALIDATION")
        print("=" * 50)
        
        try:
            from model.tabnet_prostate_variant_classifier import ProstateVariantTabNet
            
            # Initialize TabNet model (smaller for validation)
            tabnet = ProstateVariantTabNet(n_d=32, n_a=32, n_steps=3)
            
            # Load data using working logic
            X, y, feature_names = self.load_data_with_proper_features()
            
            if X is None:
                print("❌ Could not load data")
                return False
            
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
                self.validation_results['data_leakage_detected'] = True
                return False
            elif mean_accuracy > self.performance_thresholds['excellent']:
                print("⚠️  HIGH: Check features - might indicate leakage")
            elif mean_accuracy > self.performance_thresholds['good']:
                print("✅ EXCELLENT: Target performance achieved")
            elif mean_accuracy > self.performance_thresholds['acceptable']:
                print("✅ GOOD: Realistic performance")
            else:
                print("✅ MODERATE: Expected for validation subset")
            
            self.validation_results['tabnet_validation'] = {
                'mean_accuracy': mean_accuracy,
                'std_accuracy': std_accuracy,
                'fold_scores': cv_results['fold_scores']
            }
            
            return True
            
        except Exception as e:
            print(f"❌ TabNet validation failed: {e}")
            traceback.print_exc()
            return False
    
    def generate_report(self):
        """
        Generate validation report
        """
        print("\n📋 GENERATING VALIDATION REPORT")
        print("=" * 40)
        
        # Update summary
        self.validation_results['summary'] = {
            'data_leakage_detected': self.validation_results['data_leakage_detected'],
            'alphamissense_integrated': self.validation_results['alphamissense_integrated'],
            'validation_passed': not self.validation_results['data_leakage_detected']
        }
        
        # Save report
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        report_file = self.results_dir / f"enhanced_tabnet_validation_{timestamp}.json"
        
        with open(report_file, 'w') as f:
            json.dump(self.validation_results, f, indent=2)
        
        print(f"✅ Report saved: {report_file}")
        
        # Print summary
        print(f"\n📊 VALIDATION SUMMARY:")
        print(f"   Data leakage detected: {'❌ YES' if self.validation_results['data_leakage_detected'] else '✅ NO'}")
        print(f"   AlphaMissense integrated: {'✅ YES' if self.validation_results['alphamissense_integrated'] else '❌ NO'}")
        
        if 'baseline_validation' in self.validation_results:
            baseline_acc = self.validation_results['baseline_validation']['mean_accuracy']
            print(f"   Baseline accuracy: {baseline_acc:.3f}")
        
        if 'tabnet_validation' in self.validation_results:
            tabnet_acc = self.validation_results['tabnet_validation']['mean_accuracy']
            print(f"   TabNet accuracy: {tabnet_acc:.3f}")
        
        return report_file

def main():
    """
    Main validation pipeline
    """
    print("🧬 ENHANCED TABNET VALIDATION FRAMEWORK")
    print("=" * 60)
    print("Fixed data leakage issue - uses same logic as test_environment.py")
    print()
    
    # Initialize validator
    validator = EnhancedTabNetValidator()
    
    print("🧬 Enhanced TabNet Validator Initialized")
    print(f"📁 Results directory: {validator.results_dir}")
    print(f"⚠️  Suspicious accuracy threshold: {validator.performance_thresholds['suspicious']*100:.1f}%")
    
    # Step 1: Dataset validation
    if not validator.enhanced_dataset_validation():
        print("❌ Dataset validation failed - aborting")
        return False
    
    # Step 2: Baseline validation
    if not validator.baseline_validation():
        print("❌ Baseline validation failed - possible data leakage")
        return False
    
    # Step 3: TabNet validation
    if not validator.tabnet_validation():
        print("❌ TabNet validation failed")
        return False
    
    # Generate report
    report_file = validator.generate_report()
    
    print(f"\n🎉 VALIDATION COMPLETED SUCCESSFULLY!")
    print(f"📋 Report: {report_file}")
    
    if not validator.validation_results['data_leakage_detected']:
        print("✅ NO DATA LEAKAGE DETECTED - Ready for production training!")
        print("\n🎯 Next steps:")
        print("   1. Run full training: python src/model/tabnet_prostate_variant_classifier.py")
        print("   2. Or submit cluster job for full training")
    else:
        print("❌ DATA LEAKAGE DETECTED - Fix required before training")
    
    return True

if __name__ == "__main__":
    try:
        success = main()
        sys.exit(0 if success else 1)
    except Exception as e:
        print(f"❌ Validation failed: {e}")
        traceback.print_exc()
        sys.exit(1)