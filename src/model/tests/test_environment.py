#!/usr/bin/env python3
"""
Complete TabNet Environment Validation
Tests: imports, GPU access, TabNet functionality, project structure
"""

import sys
import os
import traceback

def test_python_version():
    """Test Python version"""
    print("Testing Python version...")
    print(f"Python version: {sys.version}")
    if sys.version_info >= (3, 11):
        print("✅ Python version is 3.11+")
        return True
    else:
        print("⚠️  Python version is older than 3.11")
        return False

def test_basic_imports():
    """Test all required imports"""
    print("\nTesting basic imports...")
    required_packages = [
        ('numpy', 'np'),
        ('pandas', 'pd'),
        ('sklearn', None),
        ('matplotlib.pyplot', 'plt'),
        ('seaborn', 'sns')
    ]
    
    success = True
    for package, alias in required_packages:
        try:
            if alias:
                exec(f"import {package} as {alias}")
            else:
                exec(f"import {package}")
            print(f"  ✅ {package}")
        except ImportError as e:
            print(f"  ❌ {package}: {e}")
            success = False
    
    return success

def test_pytorch():
    """Test PyTorch installation"""
    print("\nTesting PyTorch...")
    try:
        import torch
        import torchvision
        print(f"  ✅ PyTorch version: {torch.__version__}")
        print(f"  ✅ Torchvision version: {torchvision.__version__}")
        print(f"  ✅ CUDA available: {torch.cuda.is_available()}")
        if torch.cuda.is_available():
            print(f"  ✅ CUDA version: {torch.version.cuda}")
        return True
    except ImportError as e:
        print(f"  ❌ PyTorch import failed: {e}")
        return False

def test_pytorch_tabnet():
    """Test PyTorch TabNet import and basic functionality"""
    print("\nTesting PyTorch TabNet...")
    try:
        from pytorch_tabnet.tab_model import TabNetClassifier
        from pytorch_tabnet.tab_model import TabNetRegressor
        print("  ✅ TabNet import successful")
        
        # Test basic initialization
        model = TabNetClassifier(device_name='cpu')  # Use CPU for testing
        print("  ✅ TabNet classifier initialization successful")
        return True
    except ImportError as e:
        print(f"  ❌ TabNet import failed: {e}")
        return False
    except Exception as e:
        print(f"  ❌ TabNet initialization failed: {e}")
        return False

def test_gpu_access():
    """Test GPU availability and details"""
    print("\nTesting GPU access...")
    try:
        import torch
        
        if not torch.cuda.is_available():
            print("  ⚠️  CUDA not available (this is normal on login nodes)")
            print("  ℹ️  GPU access will be tested on compute nodes")
            return True
        
        print(f"  ✅ CUDA available: {torch.cuda.is_available()}")
        print(f"  ✅ CUDA version: {torch.version.cuda}")
        print(f"  ✅ GPU count: {torch.cuda.device_count()}")
        
        for i in range(torch.cuda.device_count()):
            gpu_name = torch.cuda.get_device_name(i)
            print(f"  ✅ GPU {i}: {gpu_name}")
            
            # Test basic GPU operations
            device = torch.device(f'cuda:{i}')
            test_tensor = torch.randn(10, 10).to(device)
            result = torch.mm(test_tensor, test_tensor)
            print(f"  ✅ GPU {i} computation test passed")
        
        return True
        
    except Exception as e:
        print(f"  ❌ GPU test failed: {e}")
        traceback.print_exc()
        return False

def test_tabnet_functionality():
    """Test TabNet with synthetic data"""
    print("\nTesting TabNet functionality...")
    try:
        from pytorch_tabnet.tab_model import TabNetClassifier
        import numpy as np
        import torch
        
        # Create synthetic data matching our expected structure
        n_samples = 200
        n_features = 80  # Our expected feature count
        n_classes = 5    # Our classification targets
        
        # Generate synthetic data
        X_synthetic = np.random.randn(n_samples, n_features).astype(np.float32)
        y_synthetic = np.random.randint(0, n_classes, n_samples)
        
        print(f"  ✅ Synthetic data created: {X_synthetic.shape}, {n_classes} classes")
        
        # Split data
        split_idx = int(0.8 * n_samples)
        X_train, X_val = X_synthetic[:split_idx], X_synthetic[split_idx:]
        y_train, y_val = y_synthetic[:split_idx], y_synthetic[split_idx:]
        
        # Initialize TabNet with our expected parameters
        device = 'cuda' if torch.cuda.is_available() else 'cpu'
        model = TabNetClassifier(
            n_d=32, n_a=32, n_steps=3,  # Smaller for testing
            gamma=1.3, lambda_sparse=1e-3,
            device_name=device,
            verbose=0  # Suppress training output for testing
        )
        
        print(f"  ✅ TabNet model initialized on {device}")
        
        # Test training (minimal epochs)
        model.fit(
            X_train, y_train,
            eval_set=[(X_val, y_val)],
            max_epochs=3,  # Very short for testing
            patience=3,
            batch_size=32,
            virtual_batch_size=16
        )
        
        print("  ✅ TabNet training completed")
        
        # Test prediction
        predictions = model.predict(X_val)
        probabilities = model.predict_proba(X_val)
        
        print(f"  ✅ Predictions shape: {predictions.shape}")
        print(f"  ✅ Probabilities shape: {probabilities.shape}")
        
        # Test feature importance
        feature_importance = model.feature_importances_
        print(f"  ✅ Feature importance shape: {feature_importance.shape}")
        
        # Test attention/explanation
        explain_matrix, masks = model.explain(X_val)
        print(f"  ✅ Explanation matrix shape: {explain_matrix.shape}")
        print(f"  ✅ Number of attention masks: {len(masks)}")
        
        return True
        
    except Exception as e:
        print(f"  ❌ TabNet functionality test failed: {e}")
        traceback.print_exc()
        return False

def test_project_structure():
    """Test project file structure"""
    print("\nTesting project structure...")
    
    expected_files = [
        'src/model/tabnet_prostate_variant_classifier.py',
        'requirements.txt',
        'config/pipeline_config.py',
        'src/data/preprocessing/tcga_prad_loader.py'
    ]
    
    success = True
    for file_path in expected_files:
        if os.path.exists(file_path):
            print(f"  ✅ {file_path}")
        else:
            print(f"  ❌ {file_path} (missing)")
            success = False
    
    return success

def test_tabnet_model_import():
    """Test importing our custom TabNet model"""
    print("\nTesting custom TabNet model import...")
    try:
        sys.path.append('src/model')
        from tabnet_prostate_variant_classifier import ProstateVariantTabNet
        
        # Test initialization
        model = ProstateVariantTabNet()
        print("  ✅ Custom TabNet model imported successfully")
        print(f"  ✅ Feature groups defined: {list(model.feature_groups.keys())}")
        
        return True
        
    except ImportError as e:
        print(f"  ❌ Custom model import failed: {e}")
        return False
    except Exception as e:
        print(f"  ❌ Custom model test failed: {e}")
        return False

def main():
    """Run comprehensive environment validation"""
    print("=" * 70)
    print("COMPLETE TABNET ENVIRONMENT VALIDATION")
    print("=" * 70)
    
    tests = [
        ("Python Version", test_python_version),
        ("Basic Imports", test_basic_imports),
        ("PyTorch", test_pytorch),
        ("PyTorch TabNet", test_pytorch_tabnet),
        ("GPU Access", test_gpu_access),
        ("TabNet Functionality", test_tabnet_functionality),
        ("Project Structure", test_project_structure),
        ("Custom TabNet Model", test_tabnet_model_import)
    ]
    
    results = []
    for test_name, test_func in tests:
        print(f"\n{'-' * 50}")
        print(f"RUNNING: {test_name}")
        print(f"{'-' * 50}")
        try:
            result = test_func()
            results.append((test_name, result))
        except Exception as e:
            print(f"❌ {test_name} failed with exception: {e}")
            results.append((test_name, False))
    
    # Summary
    print(f"\n{'=' * 70}")
    print("ENVIRONMENT TEST SUMMARY")
    print(f"{'=' * 70}")
    
    passed = sum(1 for _, result in results if result)
    total = len(results)
    
    for test_name, result in results:
        status = "✅ PASS" if result else "❌ FAIL"
        print(f"{test_name:.<50} {status}")
    
    print(f"\nOVERALL: {passed}/{total} tests passed")
    
    if passed == total:
        print("\n🎉 ALL TESTS PASSED - Environment ready for TabNet development!")
        print("✅ Ready to proceed with synthetic data generation")
        print("✅ Ready for H100 GPU training")
    else:
        print(f"\n⚠️  {total - passed} tests failed - check issues above")
        
    return passed == total

if __name__ == "__main__":
    main()