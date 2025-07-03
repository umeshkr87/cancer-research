#!/usr/bin/env python3
"""
Test Environment for Enhanced TabNet Prostate Cancer Classification
Validates AlphaMissense integration and absence of data leakage
"""

import pandas as pd
import numpy as np
import sys
import traceback
from pathlib import Path

# Add project root to path
sys.path.append('/u/aa107/uiuc-cancer-research/src')

def test_enhanced_dataset():
    """Test the enhanced dataset with AlphaMissense features"""
    print("\n🧬 Testing Enhanced Dataset...")
    
    try:
        enhanced_path = "/u/aa107/uiuc-cancer-research/data/processed/tabnet_csv/prostate_variants_tabnet_enhanced.csv"
        
        if not Path(enhanced_path).exists():
            print(f"  ❌ Enhanced dataset not found: {enhanced_path}")
            print(f"  💡 Run: bash scripts/enhance/functional_enhancement/run_functional_imputation.sh")
            return False
        
        df = pd.read_csv(enhanced_path)
        print(f"  ✅ Enhanced dataset loaded: {df.shape[0]:,} variants, {df.shape[1]} features")
        
        # CRITICAL: Check no data leakage features
        leakage_features = ['functional_pathogenicity', 'sift_confidence', 'polyphen_confidence']
        leakage_found = [f for f in leakage_features if f in df.columns]
        
        if leakage_found:
            print(f"  ❌ CRITICAL: Data leakage features found: {leakage_found}")
            return False
        else:
            print(f"  ✅ No data leakage features detected")
        
        # CRITICAL: Check AlphaMissense features present
        alphamissense_features = ['alphamissense_pathogenicity', 'alphamissense_class']
        missing_am = [f for f in alphamissense_features if f not in df.columns]
        
        if missing_am:
            print(f"  ❌ CRITICAL: AlphaMissense features missing: {missing_am}")
            return False
        else:
            print(f"  ✅ AlphaMissense features present")
        
        # Check AlphaMissense coverage
        am_coverage = df['alphamissense_pathogenicity'].notna().sum()
        coverage_rate = am_coverage / len(df) * 100
        print(f"  📊 AlphaMissense coverage: {am_coverage:,} variants ({coverage_rate:.1f}%)")
        
        if coverage_rate < 30:
            print(f"  ⚠️  Low coverage - expected ~43%")
        else:
            print(f"  ✅ Good coverage rate")
        
        # Check AlphaMissense score distribution
        am_scores = df['alphamissense_pathogenicity'].dropna()
        if len(am_scores) > 0:
            print(f"  📊 AlphaMissense score range: {am_scores.min():.3f} - {am_scores.max():.3f}")
            print(f"  📊 Mean pathogenicity: {am_scores.mean():.3f}")
            
            # Check for realistic distribution
            if am_scores.min() < 0 or am_scores.max() > 1:
                print(f"  ❌ Invalid score range - should be 0-1")
                return False
            else:
                print(f"  ✅ Valid score range")
        
        return True
        
    except Exception as e:
        print(f"  ❌ Enhanced dataset test failed: {e}")
        traceback.print_exc()
        return False

def test_pytorch_tabnet():
    """Test PyTorch TabNet installation"""
    print("\n🔥 Testing PyTorch TabNet...")
    
    try:
        import torch
        print(f"  ✅ PyTorch version: {torch.__version__}")
        
        # Check CUDA availability
        if torch.cuda.is_available():
            print(f"  ✅ CUDA available: {torch.cuda.get_device_name(0)}")
        else:
            print(f"  ⚠️  CUDA not available - will use CPU")
        
        from pytorch_tabnet.tab_model import TabNetClassifier
        print(f"  ✅ TabNet imported successfully")
        
        # Test TabNet initialization
        model = TabNetClassifier(n_d=8, n_a=8, n_steps=3)
        print(f"  ✅ TabNet model initialized")
        
        return True
        
    except ImportError as e:
        print(f"  ❌ TabNet import failed: {e}")
        print(f"  💡 Install with: pip install pytorch-tabnet")
        return False
    except Exception as e:
        print(f"  ❌ TabNet test failed: {e}")
        return False

def test_sklearn_dependencies():
    """Test scikit-learn dependencies"""
    print("\n🔬 Testing Scikit-learn Dependencies...")
    
    try:
        from sklearn.model_selection import train_test_split, StratifiedKFold
        from sklearn.preprocessing import LabelEncoder
        from sklearn.metrics import accuracy_score, classification_report
        print(f"  ✅ All sklearn imports successful")
        
        # Test basic functionality
        X = np.random.rand(100, 5)
        y = np.random.choice(['A', 'B', 'C'], 100)
        
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
        
        le = LabelEncoder()
        y_encoded = le.fit_transform(y)
        
        print(f"  ✅ Basic sklearn functionality works")
        return True
        
    except Exception as e:
        print(f"  ❌ Sklearn test failed: {e}")
        return False

def test_custom_tabnet_model():
    """Test our custom TabNet model"""
    print("\n🧬 Testing Custom TabNet Model...")
    
    try:
        # Import custom model
        from model.tabnet_prostate_variant_classifier import ProstateVariantTabNet
        print(f"  ✅ Custom TabNet model imported")
        
        # Test initialization
        model = ProstateVariantTabNet(n_d=32, n_a=32, n_steps=3)
        print(f"  ✅ Model initialized")
        
        # Test feature groups
        print(f"  📊 Feature groups: {len(model.feature_groups)}")
        for group_name, features in model.feature_groups.items():
            print(f"    {group_name}: {len(features)} features")
        
        # Check AlphaMissense features in groups
        all_features = []
        for features in model.feature_groups.values():
            all_features.extend(features)
        
        am_features_in_groups = [f for f in all_features if 'alphamissense' in f]
        if am_features_in_groups:
            print(f"  ✅ AlphaMissense features in groups: {am_features_in_groups}")
        else:
            print(f"  ❌ No AlphaMissense features in feature groups")
            return False
        
        # Check banned features
        print(f"  🚫 Banned features: {model.banned_features}")
        
        return True
        
    except ImportError as e:
        print(f"  ❌ Custom model import failed: {e}")
        print(f"  💡 Check file path and syntax")
        return False
    except Exception as e:
        print(f"  ❌ Custom model test failed: {e}")
        traceback.print_exc()
        return False

def test_data_loading():
    """Test data loading with custom model"""
    print("\n📁 Testing Data Loading...")
    
    try:
        from model.tabnet_prostate_variant_classifier import ProstateVariantTabNet
        
        model = ProstateVariantTabNet()
        
        # Test data loading (should work with enhanced dataset)
        X, y = model.load_data()
        
        print(f"  ✅ Data loaded successfully")
        print(f"  📊 Shape: {X.shape}")
        print(f"  🎯 Classes: {np.unique(y)}")
        print(f"  📝 Features: {len(model.feature_names)}")
        
        # Check for AlphaMissense features in selected features
        am_features = [f for f in model.feature_names if 'alphamissense' in f]
        if am_features:
            print(f"  ✅ AlphaMissense features selected: {am_features}")
        else:
            print(f"  ⚠️  No AlphaMissense features in selected features")
        
        return True
        
    except Exception as e:
        print(f"  ❌ Data loading failed: {e}")
        traceback.print_exc()
        return False

def test_quick_training():
    """Test quick training on small subset"""
    print("\n🚀 Testing Quick Training...")
    
    try:
        from model.tabnet_prostate_variant_classifier import ProstateVariantTabNet
        from sklearn.model_selection import train_test_split
        
        model = ProstateVariantTabNet(n_d=16, n_a=16, n_steps=2)  # Small for testing
        
        # Load data
        X, y = model.load_data()
        
        # Use small subset for quick test
        if len(X) > 1000:
            indices = np.random.choice(len(X), 1000, replace=False)
            X, y = X[indices], y[indices]
        
        # Split data
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
        
        print(f"  📊 Training on {len(X_train)} samples...")
        
        # Quick training
        model.train(X_train, y_train, max_epochs=10)  # Very few epochs for testing
        
        # Quick evaluation
        accuracy = model.evaluate(X_test, y_test)
        
        print(f"  🎯 Test accuracy: {accuracy:.3f}")
        
        # Check for suspicious accuracy
        if accuracy > 0.95:
            print(f"  ⚠️  WARNING: Suspiciously high accuracy - check for data leakage")
        else:
            print(f"  ✅ Realistic accuracy for quick test")
        
        return True
        
    except Exception as e:
        print(f"  ❌ Quick training failed: {e}")
        traceback.print_exc()
        return False

def run_all_tests():
    """Run all environment tests"""
    print("🧪 ENHANCED TABNET ENVIRONMENT TESTING")
    print("=" * 50)
    print("Testing AlphaMissense integration and data leakage elimination")
    print()
    
    tests = [
        ("Enhanced Dataset", test_enhanced_dataset),
        ("PyTorch TabNet", test_pytorch_tabnet),
        ("Sklearn Dependencies", test_sklearn_dependencies),
        ("Custom TabNet Model", test_custom_tabnet_model),
        ("Data Loading", test_data_loading),
        ("Quick Training", test_quick_training),
    ]
    
    results = []
    for test_name, test_func in tests:
        try:
            result = test_func()
            results.append((test_name, result))
        except Exception as e:
            print(f"❌ {test_name} test crashed: {e}")
            results.append((test_name, False))
    
    # Summary
    print(f"\n📊 TEST SUMMARY")
    print("=" * 30)
    
    passed = 0
    total = len(results)
    
    for test_name, result in results:
        status = "✅ PASS" if result else "❌ FAIL"
        print(f"{status} {test_name}")
        if result:
            passed += 1
    
    print(f"\n🎯 OVERALL: {passed}/{total} tests passed")
    
    if passed == total:
        print("✅ ALL TESTS PASSED - Environment ready for enhanced TabNet training!")
        return True
    else:
        print("❌ Some tests failed - fix issues before proceeding")
        return False

if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)