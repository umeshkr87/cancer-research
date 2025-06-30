#!/usr/bin/env python3
"""
Updated Environment Tests for TabNet Prostate Cancer Variant Classification
Tests real dataset loading, confidence features, and H100 GPU functionality
"""

import sys
import os
import traceback
import subprocess
from pathlib import Path

# Setup environment before importing anything else
def setup_environment():
    """Setup conda environment and install packages if needed"""
    print("🔧 Setting up environment...")
    
    # Check if we're already in the right environment
    conda_env = os.environ.get('CONDA_DEFAULT_ENV', '')
    if conda_env == 'tabnet-prostate':
        print(f"✅ Already in conda environment: {conda_env}")
        return True
    
    # Try to activate conda environment
    try:
        # First try to initialize conda if not done
        conda_init_cmd = [
            'bash', '-c', 
            'module load anaconda3 2>/dev/null; '
            'eval "$(conda shell.bash hook)"; '
            'conda activate tabnet-prostate; '
            'python -c "import sys; print(sys.executable)"'
        ]
        
        result = subprocess.run(conda_init_cmd, capture_output=True, text=True, timeout=30)
        if result.returncode == 0:
            python_path = result.stdout.strip()
            print(f"✅ Conda environment available: {python_path}")
            
            # Re-execute this script with the correct python
            if python_path != sys.executable:
                print("🔄 Re-executing with conda python...")
                os.execv(python_path, [python_path] + sys.argv)
            return True
            
    except Exception as e:
        print(f"⚠️  Conda environment setup failed: {e}")
    
    # Fallback: install packages in current environment
    print("🔧 Installing required packages in current environment...")
    packages = ['torch', 'pandas', 'numpy', 'scikit-learn', 'pytorch-tabnet']
    
    for package in packages:
        try:
            __import__(package.replace('-', '_'))
            print(f"✅ {package} already available")
        except ImportError:
            print(f"📦 Installing {package}...")
            try:
                subprocess.check_call([sys.executable, '-m', 'pip', 'install', package, '--quiet'])
                print(f"✅ {package} installed")
            except subprocess.CalledProcessError:
                print(f"❌ Failed to install {package}")
                return False
    
    return True

# Setup environment first
if not setup_environment():
    print("❌ Environment setup failed")
    sys.exit(1)

# Now import packages
try:
    import pandas as pd
    import numpy as np
    import torch
except ImportError as e:
    print(f"❌ Import failed even after setup: {e}")
    print("💡 Try running: module load anaconda3 && conda activate tabnet-prostate")
    sys.exit(1)

# Add project root to path
sys.path.append('/u/aa107/uiuc-cancer-research/src/model')
sys.path.append('/u/aa107/uiuc-cancer-research')

def test_python_version():
    """Test Python version compatibility"""
    print("Testing Python version...")
    version = sys.version_info
    
    if version.major == 3 and version.minor >= 8:
        print(f"  ✅ Python {version.major}.{version.minor}.{version.micro}")
        return True
    else:
        print(f"  ❌ Python {version.major}.{version.minor}.{version.micro} (requires 3.8+)")
        return False

def test_core_imports():
    """Test essential package imports"""
    print("\nTesting core package imports...")
    
    packages = {
        'pandas': 'pd',
        'numpy': 'np', 
        'torch': 'torch',
        'sklearn': 'sklearn',
        'matplotlib': 'plt'
    }
    
    success = True
    for package, alias in packages.items():
        try:
            if package == 'pandas':
                import pandas as pd
                print(f"  ✅ pandas {pd.__version__}")
            elif package == 'numpy':
                import numpy as np
                print(f"  ✅ numpy {np.__version__}")
            elif package == 'torch':
                import torch
                print(f"  ✅ pytorch {torch.__version__}")
            elif package == 'sklearn':
                import sklearn
                print(f"  ✅ scikit-learn {sklearn.__version__}")
            elif package == 'matplotlib':
                import matplotlib.pyplot as plt
                import matplotlib
                print(f"  ✅ matplotlib {matplotlib.__version__}")
        except ImportError as e:
            print(f"  ❌ {package} import failed: {e}")
            success = False
    
    return success

def test_pytorch_tabnet():
    """Test PyTorch TabNet import and basic functionality"""
    print("\nTesting PyTorch TabNet...")
    
    try:
        from pytorch_tabnet.tab_model import TabNetClassifier
        print("  ✅ pytorch-tabnet imported successfully")
        
        # Test basic instantiation
        model = TabNetClassifier(n_d=8, n_a=8, n_steps=3, verbose=0)
        print("  ✅ TabNetClassifier instantiated")
        
        return True
        
    except ImportError as e:
        print(f"  ❌ pytorch-tabnet import failed: {e}")
        print("  🔧 Installing pytorch-tabnet...")
        try:
            subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'pytorch-tabnet', '--quiet'])
            from pytorch_tabnet.tab_model import TabNetClassifier
            print("  ✅ pytorch-tabnet installed and imported successfully")
            return True
        except Exception as install_error:
            print(f"  ❌ Installation failed: {install_error}")
            return False
    except Exception as e:
        print(f"  ❌ TabNet test failed: {e}")
        return False

def test_gpu_access():
    """Test GPU availability and H100 detection"""
    print("\nTesting GPU access...")
    
    try:
        import torch
        
        if torch.cuda.is_available():
            device_count = torch.cuda.device_count()
            device_name = torch.cuda.get_device_name(0)
            device_props = torch.cuda.get_device_properties(0)
            memory_gb = device_props.total_memory / 1e9
            
            print(f"  ✅ CUDA available: {torch.version.cuda}")
            print(f"  ✅ GPU devices: {device_count}")
            print(f"  ✅ Primary GPU: {device_name}")
            print(f"  ✅ GPU memory: {memory_gb:.1f}GB")
            
            # Check if it's H100 or equivalent
            if "H100" in device_name or "A100" in device_name or memory_gb > 40:
                print("  🚀 High-performance GPU detected (H100/A100-class)")
                return True
            else:
                print("  ⚠️  Standard GPU detected (will work but slower)")
                return True
        else:
            print("  ❌ CUDA not available")
            print("  💡 Will fall back to CPU training")
            return False
            
    except Exception as e:
        print(f"  ❌ GPU test failed: {e}")
        return False

def test_dataset_loading():
    """Test loading the actual imputed prostate cancer dataset"""
    print("\nTesting dataset loading...")
    
    dataset_path = "/u/aa107/uiuc-cancer-research/data/processed/tabnet_csv/prostate_variants_tabnet_imputed.csv"
    
    try:
        # Check if file exists
        if not os.path.exists(dataset_path):
            print(f"  ❌ Dataset not found: {dataset_path}")
            print("  💡 Ensure the imputed dataset has been generated")
            return False
        
        # Load dataset
        print(f"  📁 Loading: {os.path.basename(dataset_path)}")
        df = pd.read_csv(dataset_path)
        
        print(f"  ✅ Dataset loaded: {len(df):,} variants")
        print(f"  ✅ Features available: {df.shape[1]}")
        
        # Check for critical confidence features from imputation breakthrough
        confidence_features = ['sift_confidence', 'polyphen_confidence', 'functional_pathogenicity']
        missing_confidence = [f for f in confidence_features if f not in df.columns]
        
        if missing_confidence:
            print(f"  ⚠️  Missing confidence features: {missing_confidence}")
        else:
            print("  ✅ All confidence features present")
        
        # Check functional scores coverage
        sift_coverage = df['sift_score'].notna().sum() / len(df) * 100
        polyphen_coverage = df['polyphen_score'].notna().sum() / len(df) * 100
        
        print(f"  📊 SIFT coverage: {sift_coverage:.1f}%")
        print(f"  📊 PolyPhen coverage: {polyphen_coverage:.1f}%")
        
        # Expect ~100% coverage after imputation
        if sift_coverage > 95 and polyphen_coverage > 95:
            print("  ✅ Excellent functional score coverage (imputation successful)")
        else:
            print("  ⚠️  Functional score coverage below expected (check imputation)")
        
        return True
        
    except Exception as e:
        print(f"  ❌ Dataset loading failed: {e}")
        traceback.print_exc()
        return False

def test_custom_tabnet_model():
    """Test our custom TabNet model implementation"""
    print("\nTesting custom TabNet model...")
    
    try:
        # First ensure config can be imported
        try:
            from config.pipeline_config import PipelineConfig
            print("  ✅ Pipeline config imported")
        except ImportError as e:
            print(f"  ⚠️  Pipeline config import failed: {e}")
            print("  💡 This is expected if running outside project structure")
        
        try:
            from tabnet_prostate_variant_classifier import ProstateVariantTabNet
            print("  ✅ Custom TabNet model imported")
        except ImportError as e:
            print(f"  ❌ Custom model import failed: {e}")
            print("  💡 Ensure you're running from the project root directory")
            return False
        
        # Test model initialization
        model = ProstateVariantTabNet(n_d=32, n_a=32, n_steps=3)  # Small for testing
        print("  ✅ Model initialized successfully")
        
        # Test feature groups
        feature_groups = model.feature_groups
        expected_groups = ['functional_scores', 'genomic_context', 'pathway_indicators', 'clinical_impact']
        
        for group in expected_groups:
            if group in feature_groups:
                print(f"  ✅ Feature group '{group}': {len(feature_groups[group])} features")
            else:
                print(f"  ❌ Missing feature group: {group}")
        
        return True
        
    except Exception as e:
        print(f"  ❌ Custom model test failed: {e}")
        traceback.print_exc()
        return False

def test_tabnet_with_real_data():
    """Test TabNet training with real dataset (quick test)"""
    print("\nTesting TabNet with real data (quick test)...")
    
    try:
        from tabnet_prostate_variant_classifier import ProstateVariantTabNet
        
        # Initialize model with small parameters for quick test
        model = ProstateVariantTabNet(n_d=16, n_a=16, n_steps=2)
        
        # Try to load real data
        X, y = model.load_data()
        print(f"  ✅ Data loaded: {X.shape}")
        
        # Quick train on subset
        subset_size = min(1000, len(X))  # Use small subset for testing
        X_subset = X[:subset_size]
        y_subset = y[:subset_size]
        
        print(f"  🧪 Testing on {subset_size} variants...")
        
        # Split data
        from sklearn.model_selection import train_test_split
        X_train, X_val, y_train, y_val = train_test_split(
            X_subset, y_subset, test_size=0.3, stratify=y_subset, random_state=42
        )
        
        # Quick training (1 epoch)
        print("  🚀 Quick training test...")
        device = 'cuda' if torch.cuda.is_available() else 'cpu'
        
        from pytorch_tabnet.tab_model import TabNetClassifier
        quick_model = TabNetClassifier(
            n_d=16, n_a=16, n_steps=2,
            device_name=device, verbose=0
        )
        
        quick_model.fit(
            X_train, y_train,
            eval_set=[(X_val, y_val)],
            max_epochs=2,  # Very quick test
            batch_size=64,
            virtual_batch_size=32
        )
        
        # Test prediction
        predictions = quick_model.predict(X_val)
        print(f"  ✅ Predictions generated: {predictions.shape}")
        
        # Test attention
        explain_matrix, masks = quick_model.explain(X_val[:5])
        print(f"  ✅ Attention analysis: {explain_matrix.shape}")
        
        return True
        
    except Exception as e:
        print(f"  ❌ Real data test failed: {e}")
        traceback.print_exc()
        return False

def test_project_structure():
    """Test project directory structure"""
    print("\nTesting project structure...")
    
    base_path = "/u/aa107/uiuc-cancer-research"
    
    expected_structure = {
        'data/processed/tabnet_csv/': ['prostate_variants_tabnet_imputed.csv'],
        'src/model/': ['tabnet_prostate_variant_classifier.py'],
        'config/': ['pipeline_config.py'],
        'models/': [],  # Will be created
        'results/': [],  # Will be created
        'logs/': []  # Will be created
    }
    
    success = True
    for directory, files in expected_structure.items():
        dir_path = os.path.join(base_path, directory)
        
        if os.path.exists(dir_path):
            print(f"  ✅ Directory: {directory}")
            
            for file in files:
                file_path = os.path.join(dir_path, file)
                if os.path.exists(file_path):
                    size_mb = os.path.getsize(file_path) / 1e6
                    print(f"    ✅ {file} ({size_mb:.1f}MB)")
                else:
                    print(f"    ❌ {file} (missing)")
                    success = False
        else:
            print(f"  ⚠️  Directory missing: {directory}")
            # Create directory if it's for outputs
            if directory in ['models/', 'results/', 'logs/']:
                os.makedirs(dir_path, exist_ok=True)
                print(f"    ✅ Created: {directory}")
    
    return success

def test_configuration():
    """Test configuration loading and validation"""
    print("\nTesting configuration...")
    
    try:
        from config.pipeline_config import config
        
        # Test basic config access
        print(f"  ✅ Project root: {config.PROJECT_ROOT}")
        print(f"  ✅ Target accuracy: {config.TARGET_ACCURACY}")
        print(f"  ✅ TabNet steps: {config.DEFAULT_TABNET_PARAMS['n_steps']}")
        
        # Test GPU config
        gpu_config = config.get_gpu_config()
        print(f"  ✅ GPU config: {gpu_config['device']}")
        
        # Test directory creation
        config.create_directories()
        print("  ✅ Directories created")
        
        # Test file validation
        files_exist = config.validate_data_files()
        if files_exist:
            print("  ✅ All required files present")
        else:
            print("  ⚠️  Some required files missing")
        
        return True
        
    except Exception as e:
        print(f"  ❌ Configuration test failed: {e}")
        traceback.print_exc()
        return False

def main():
    """Run complete environment validation"""
    print("=" * 70)
    print("TABNET PROSTATE CANCER - COMPLETE ENVIRONMENT VALIDATION")
    print("Updated for Week 5 with real dataset and confidence features")
    print("=" * 70)
    
    tests = [
        ("Python Version", test_python_version),
        ("Core Imports", test_core_imports), 
        ("PyTorch TabNet", test_pytorch_tabnet),
        ("GPU Access", test_gpu_access),
        ("Dataset Loading", test_dataset_loading),
        ("Project Structure", test_project_structure),
        ("Configuration", test_configuration),
        ("Custom TabNet Model", test_custom_tabnet_model),
        ("Real Data Training", test_tabnet_with_real_data)
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
    
    # Final summary
    print(f"\n{'=' * 70}")
    print("ENVIRONMENT VALIDATION SUMMARY")
    print(f"{'=' * 70}")
    
    passed = sum(1 for _, result in results if result)
    total = len(results)
    
    for test_name, result in results:
        status = "✅ PASS" if result else "❌ FAIL"
        print(f"{test_name:.<40} {status}")
    
    print(f"\nOVERALL RESULT: {passed}/{total} tests passed")
    
    if passed == total:
        print("\n🎉 ALL TESTS PASSED!")
        print("✅ Environment ready for Week 5 TabNet optimization")
        print("✅ H100 GPU ready for hyperparameter search")
        print("✅ Real dataset with confidence features loaded")
        print("✅ Ready to achieve 80-85% accuracy target")
        print("\n🚀 Next steps:")
        print("   1. Run: python src/model/validate_tabnet.py")
        print("   2. Run: python scripts/optimization/hyperparameter_optimization.py")
    elif passed >= total * 0.8:
        print(f"\n⚠️  {total - passed} tests failed, but core functionality works")
        print("✅ Can proceed with optimization (fix issues in parallel)")
        print("\n🔧 Common fixes:")
        if not any(result for name, result in results if "GPU" in name):
            print("   - GPU issues: Submit jobs to GPU partition instead")
        if not any(result for name, result in results if "Dataset" in name):
            print("   - Dataset issues: Check if imputed dataset exists")
        print("   - Missing packages: Run 'pip install pytorch-tabnet pandas numpy'")
    else:
        print(f"\n❌ {total - passed} critical tests failed")
        print("🔧 Fix issues before proceeding with optimization")
        print("\n💡 Quick fixes:")
        print("   1. Ensure you're in project root: cd /u/aa107/uiuc-cancer-research")
        print("   2. Load environment: module load anaconda3")
        print("   3. Create environment: conda create -n tabnet-prostate python=3.11 -y")
        print("   4. Activate environment: conda activate tabnet-prostate")
        print("   5. Install packages: pip install pytorch-tabnet pandas numpy scikit-learn")
    
    print(f"\n{'=' * 70}")
    
    return passed >= total * 0.8

if __name__ == "__main__":
    print("🧬 Starting TabNet Environment Validation...")
    print("🔧 Note: This script will auto-install missing packages if needed")
    print()
    
    success = main()
    exit(0 if success else 1)