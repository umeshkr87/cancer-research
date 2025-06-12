# TabNet Prostate Cancer Variant Classification

Interpretable deep learning framework for prostate cancer genomic variant classification using TabNet architecture with attention mechanisms.

## Dependencies

### System Requirements
- **Platform**: University of Illinois Campus Cluster
- **GPU**: H100 nodes (IllinoisComputes-GPU partition)
- **Account**: aa107-ic allocation
- **Python**: 3.11+

### Environment Setup
```bash
# Load anaconda and create environment
module load anaconda3
conda create -n tabnet-prostate python=3.11 -y
conda activate tabnet-prostate

# Install PyTorch with CUDA support
conda install pytorch torchvision torchaudio pytorch-cuda=11.8 -c pytorch -c nvidia -y

# Install remaining dependencies
pip install -r requirements.txt
```

## Validation

### Test Environment
```bash
# Quick local test
python src/model/tests/test_environment.py

# Complete GPU test
cd src/model/tests
sbatch test_gpu_environment.sbatch

# Automated test with cleanup
./run_tabnet_tests.sh
```

### Success Criteria
- 8/8 tests pass
- H100 GPU detected (NVIDIA A100-80GB)
- TabNet trains on synthetic data
- Attention mechanisms functional

## Model Usage

```python
from src.model.tabnet_prostate_variant_classifier import ProstateVariantTabNet

# Initialize model
model = ProstateVariantTabNet(n_d=64, n_a=64, n_steps=6)

# Train with your data
model.train(X_train, y_train, X_val, y_val)

# Get predictions with interpretability
results = model.predict_with_explanation(X_test, feature_names)
```

## Project Structure

```
├── requirements.txt              # Dependencies
├── src/model/
│   ├── tabnet_prostate_variant_classifier.py    # Main TabNet model
│   └── tests/                    # Environment validation
├── src/data/preprocessing/       # Data processing modules
├── config/pipeline_config.py    # Configuration
└── notebooks/                   # Analysis notebooks
```

## Next Steps

1. Generate synthetic data (80 features, 5 classes)
2. Test with real TCGA-PRAD + COSMIC + ClinVar data
3. Optimize attention patterns for clinical interpretability
4. Benchmark H100 performance

## Contact

**Team**: Abraham Arellano Tavara (aa107@illinois.edu)  
**Mentor**: Jathurshan Pradeepkumar  
**Support**: help@campuscluster.illinois.edu