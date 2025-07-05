#!/usr/bin/env python3
"""
Pipeline Configuration for TabNet Prostate Cancer Variant Classification
Centralized configuration management for reproducible experiments
"""

import os
from pathlib import Path

class PipelineConfig:
    """
    Configuration class for TabNet prostate cancer classification pipeline
    """
    
    # =====================
    # PROJECT PATHS
    # =====================
    PROJECT_ROOT = Path("/u/aa107/uiuc-cancer-research")
    DATA_DIR = PROJECT_ROOT / "data"
    PROCESSED_DATA_DIR = DATA_DIR / "processed" / "tabnet_csv"
    MODELS_DIR = PROJECT_ROOT / "models"
    RESULTS_DIR = PROJECT_ROOT / "results"
    LOGS_DIR = PROJECT_ROOT / "logs"
    
    # Input data files
    IMPUTED_DATASET = PROCESSED_DATA_DIR / "prostate_variants_tabnet_imputed.csv"
    FEATURE_ANALYSIS_REPORT = PROCESSED_DATA_DIR / "feature_analysis_report.txt"
    IMPUTATION_REPORT = PROCESSED_DATA_DIR / "imputation_report.txt"
    
    # =====================
    # TABNET HYPERPARAMETERS
    # =====================
    
    # Performance target based on project goals
    TARGET_ACCURACY = 0.82  # 82% target (80-85% range)
    COMPETITIVE_BASELINE = 0.93  # XGBoost ROC baseline to beat
    
    # TabNet Architecture - Default Configuration
    DEFAULT_TABNET_PARAMS = {
        'n_d': 64,                    # Decision prediction layer width
        'n_a': 64,                    # Attention embedding dimension
        'n_steps': 6,                 # Number of decision steps
        'gamma': 1.3,                 # Feature reusage parameter
        'lambda_sparse': 1e-3,        # Sparsity regularization
        'mask_type': 'entmax',        # Attention mask type
        'optimizer_fn': 'Adam',       # Optimizer
        'optimizer_params': {'lr': 2e-2},  # Learning rate
        'scheduler_fn': None,         # Learning rate scheduler
        'scheduler_params': {},
        'seed': 42                    # Random seed for reproducibility
    }
    
    # Hyperparameter Optimization Ranges for H100 GPU Search
    HYPERPARAMETER_GRID = {
        'n_d': [32, 64, 128],         # Decision layer width
        'n_a': [32, 64, 128],         # Attention dimension
        'n_steps': [6, 7, 8],         # Decision steps (interpretability vs performance)
        'gamma': [1.0, 1.3, 1.5],     # Feature reuse strength
        'lambda_sparse': [1e-4, 1e-3, 1e-2],  # Sparsity regularization
        'learning_rate': [1e-3, 2e-3, 1e-2, 2e-2]  # Learning rates
    }
    
    # =====================
    # TRAINING CONFIGURATION
    # =====================
    TRAINING_PARAMS = {
        'max_epochs': 200,            # Maximum training epochs
        'patience': 15,               # Early stopping patience
        'batch_size': 1024,           # Batch size (optimized for H100)
        'virtual_batch_size': 128,    # Virtual batch size for TabNet
        'num_workers': 4,             # Data loading workers
        'drop_last': False,           # Keep last incomplete batch
        'eval_metric': 'accuracy',    # Primary evaluation metric
        'test_size': 0.2,             # Train/test split ratio
        'validation_size': 0.2,       # Validation set size
        'stratify': True,             # Stratified sampling
        'random_state': 42            # Random seed
    }
    
    # =====================
    # CLINICAL FEATURES
    # =====================
    
    # Critical features from imputation breakthrough
    CONFIDENCE_FEATURES = [
        'sift_confidence',            # Quality of SIFT score (1.0=observed, 0.5-0.7=imputed)
        'polyphen_confidence',        # Quality of PolyPhen score
        'functional_pathogenicity'    # Composite pathogenicity score
    ]
    
    # Core functional prediction scores
    FUNCTIONAL_FEATURES = [
        'sift_score',                 # SIFT deleteriousness score
        'polyphen_score',             # PolyPhen pathogenicity score
        'cadd_phred',                 # CADD scaled score
        'conservation_score'          # Conservation across species
    ]
    
    # Prostate-specific pathway indicators
    PATHWAY_FEATURES = [
        'dna_repair_pathway',         # DNA repair genes (BRCA1/2, ATM, etc.)
        'mismatch_repair_pathway',    # MMR genes (MLH1, MSH2, etc.)
        'hormone_pathway',            # Androgen receptor pathway
        'is_important_gene'           # Important prostate cancer genes
    ]
    
    # Clinical interpretability groups for attention analysis
    CLINICAL_FEATURE_GROUPS = {
        'PARP_inhibitor_relevant': [
            'dna_repair_pathway', 'brca1_interaction', 'brca2_interaction',
            'atm_interaction', 'palb2_interaction', 'chek2_interaction'
        ],
        'hormone_therapy_relevant': [
            'hormone_pathway', 'ar_pathway_score', 'hormone_sensitivity',
            'adt_response', 'enzalutamide_response'
        ],
        'functional_predictions': [
            'sift_score', 'polyphen_score', 'cadd_phred', 'conservation_score'
        ],
        'confidence_indicators': [
            'sift_confidence', 'polyphen_confidence', 'functional_pathogenicity'
        ]
    }
    
    # =====================
    # TARGET CLASSIFICATION
    # =====================
    
    # 5-class variant classification system
    TARGET_CLASSES = [
        'Actionable Pathogenic',      # Direct therapeutic implications
        'Likely Actionable',          # Probable therapeutic relevance
        'VUS',                        # Uncertain significance
        'Likely Benign',              # Probably not therapeutically relevant
        'Benign'                      # No therapeutic implications
    ]
    
    # Clinical significance mapping
    CLINICAL_SIGNIFICANCE_MAP = {
        'pathogenic': 'Actionable Pathogenic',
        'likely_pathogenic': 'Likely Actionable',
        'uncertain_significance': 'VUS',
        'likely_benign': 'Likely Benign',
        'benign': 'Benign'
    }
    
    # =====================
    # VALIDATION CONFIGURATION
    # =====================
    VALIDATION_PARAMS = {
        'cv_folds': 5,                # 5-fold cross-validation
        'cv_strategy': 'stratified',  # Stratified K-fold
        'metrics': [                  # Evaluation metrics
            'accuracy',
            'precision_macro',
            'recall_macro',
            'f1_macro',
            'roc_auc_ovr'             # One-vs-rest ROC AUC
        ],
        'confidence_interval': 0.95,  # Statistical confidence level
        'bootstrap_samples': 1000     # Bootstrap samples for CI
    }
    
    # =====================
    # H100 GPU OPTIMIZATION
    # =====================
    GPU_CONFIG = {
        'device': 'cuda',             # GPU device
        'gpu_partition': 'IllinoisComputes-GPU',  # Campus cluster partition
        'max_gpu_hours': 40,          # Week 5 allocation
        'optimization_budget': 36,    # Hours for hyperparameter search
        'validation_budget': 4,       # Hours for final validation
        'parallel_jobs': 2,           # Parallel optimization jobs
        'memory_limit': '80GB',       # H100 memory limit
        'batch_optimization': True    # Enable batch optimization
    }
    
    # Grid search configuration
    OPTIMIZATION_CONFIG = {
        'search_strategy': 'grid',    # Grid search for systematic exploration
        'total_configurations': 108, # 3*3*3*3*3*4 combinations
        'early_stopping': True,      # Stop poor configurations early
        'save_best_n': 5,           # Save top 5 models
        'evaluation_time': 20,      # Minutes per configuration
        'target_metric': 'accuracy'  # Primary optimization metric
    }
    
    # =====================
    # LOGGING AND OUTPUT
    # =====================
    LOGGING_CONFIG = {
        'log_level': 'INFO',
        'log_file': 'tabnet_training.log',
        'experiment_tracking': True,
        'save_attention_patterns': True,
        'save_feature_importance': True,
        'save_predictions': True
    }
    
    OUTPUT_CONFIG = {
        'model_save_format': 'pytorch',
        'results_format': 'json',
        'plots_format': 'png',
        'save_intermediate': True,
        'compress_outputs': True
    }
    
    # =====================
    # INTERPRETABILITY
    # =====================
    INTERPRETABILITY_CONFIG = {
        'attention_analysis': True,
        'feature_importance_analysis': True,
        'clinical_case_studies': 15,     # Number of detailed case studies
        'attention_visualization': True,
        'pathway_analysis': True,
        'expert_validation_ready': True
    }
    
    @classmethod
    def create_directories(cls):
        """Create necessary project directories"""
        directories = [
            cls.MODELS_DIR,
            cls.RESULTS_DIR,
            cls.LOGS_DIR,
            cls.PROCESSED_DATA_DIR
        ]
        
        for directory in directories:
            directory.mkdir(parents=True, exist_ok=True)
            print(f"✅ Directory ready: {directory}")
    
    @classmethod
    def validate_data_files(cls):
        """Validate that required data files exist"""
        required_files = [
            cls.IMPUTED_DATASET,
            cls.FEATURE_ANALYSIS_REPORT,
            cls.IMPUTATION_REPORT
        ]
        
        missing_files = []
        for file_path in required_files:
            if not file_path.exists():
                missing_files.append(file_path)
            else:
                print(f"✅ Data file found: {file_path.name}")
        
        if missing_files:
            print("❌ Missing required files:")
            for file_path in missing_files:
                print(f"   - {file_path}")
            return False
        
        return True
    
    @classmethod
    def get_gpu_config(cls):
        """Get GPU configuration for H100 optimization"""
        import torch
        
        if torch.cuda.is_available():
            gpu_name = torch.cuda.get_device_name(0)
            gpu_memory = torch.cuda.get_device_properties(0).total_memory / 1e9
            
            config = cls.GPU_CONFIG.copy()
            config['detected_gpu'] = gpu_name
            config['available_memory'] = f"{gpu_memory:.1f}GB"
            config['cuda_version'] = torch.version.cuda
            
            print(f"✅ GPU detected: {gpu_name}")
            print(f"✅ Memory available: {gpu_memory:.1f}GB")
            
            return config
        else:
            print("⚠️  No GPU detected - falling back to CPU")
            return {'device': 'cpu', 'available_memory': 'N/A'}
    
    @classmethod
    def summary(cls):
        """Print configuration summary"""
        print("=" * 60)
        print("TABNET PROSTATE CANCER CLASSIFICATION CONFIG")
        print("=" * 60)
        print(f"Project Root: {cls.PROJECT_ROOT}")
        print(f"Dataset: {cls.IMPUTED_DATASET.name}")
        print(f"Target Accuracy: {cls.TARGET_ACCURACY:.1%}")
        print(f"TabNet Steps: {cls.DEFAULT_TABNET_PARAMS['n_steps']}")
        print(f"CV Folds: {cls.VALIDATION_PARAMS['cv_folds']}")
        print(f"H100 GPU Hours: {cls.GPU_CONFIG['max_gpu_hours']}")
        print(f"Grid Search Size: {cls.OPTIMIZATION_CONFIG['total_configurations']}")
        print("=" * 60)

# Create global config instance
config = PipelineConfig()

if __name__ == "__main__":
    # Test configuration
    config.summary()
    config.create_directories()
    config.validate_data_files()
    config.get_gpu_config()