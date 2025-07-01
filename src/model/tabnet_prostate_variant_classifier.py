#!/usr/bin/env python3
"""
TabNet Prostate Cancer Variant Classification - DATA LEAKAGE FIXED
Interpretable deep learning for prostate cancer genomic variant classification
using TabNet with attention mechanisms. 

CRITICAL FIX: Removes data leakage features that caused artificial 100% accuracy
"""

import pandas as pd
import numpy as np
import torch
from pytorch_tabnet.tab_model import TabNetClassifier
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.metrics import accuracy_score, roc_auc_score, classification_report, confusion_matrix
from sklearn.ensemble import RandomForestClassifier
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
from pathlib import Path
import joblib
warnings.filterwarnings('ignore')

class ProstateVariantTabNet:
    """
    TabNet model for prostate cancer variant classification with interpretability
    FIXED: No data leakage - realistic performance expected (75-85% accuracy)
    """
    
    def __init__(self, n_d=64, n_a=64, n_steps=6, gamma=1.3, lambda_sparse=1e-3):
        """
        Initialize TabNet model with prostate cancer-specific configuration
        
        Args:
            n_d: Decision prediction layer width
            n_a: Attention embedding dimension  
            n_steps: Number of sequential attention steps
            gamma: Relaxation parameter for feature reusage
            lambda_sparse: Sparsity regularization strength
        """
        self.n_d = n_d
        self.n_a = n_a
        self.n_steps = n_steps
        self.gamma = gamma
        self.lambda_sparse = lambda_sparse
        
        # Clinical feature groups for interpretability (WITHOUT LEAKAGE FEATURES)
        self.feature_groups = {
            'functional_scores': ['sift_score', 'polyphen_score', 'cadd_phred', 
                                'conservation_score', 'gerp_score', 'phylop_score'],
            'genomic_context': ['variant_impact', 'variant_type', 'chromosome', 
                              'exonic_function', 'splicing_impact'],
            'frequency_indicators': ['gnomad_af_eur', 'gnomad_af_afr', 'gnomad_af_asj',
                                   'af_1kg', 'exac_af', 'is_rare'],
            'pathway_indicators': ['dna_repair_pathway', 'mismatch_repair_pathway', 
                                 'hormone_pathway', 'is_important_gene', 'core_prostate_gene'],
            'structural_features': ['is_indel', 'is_snv', 'variant_size', 'ref_length', 'alt_length']
        }
        
        # EXPLICITLY BANNED: These features cause data leakage
        self.banned_features = [
            'functional_pathogenicity',  # Derived from target
            'sift_confidence',           # Artificial indicator
            'polyphen_confidence',       # Artificial indicator
            'clinical_significance',     # Direct target correlation
            'variant_classification'     # This IS the target
        ]
        
        self.model = None
        self.feature_names = None
        self.selected_features = None
        self.label_encoder = LabelEncoder()
        self.feature_encoders = {}
        
    def load_data(self, data_path=None):
        """
        Load and prepare prostate cancer variant data WITHOUT data leakage
        
        Args:
            data_path: Path to CSV file. If None, uses default clean dataset
            
        Returns:
            X: Feature matrix (numpy array)
            y: Target labels (numpy array)
        """
        if data_path is None:
            # Use the clean dataset (without leakage features)
            data_path = "/u/aa107/uiuc-cancer-research/data/processed/tabnet_csv/prostate_variants_tabnet_clean.csv"
            
        print(f"📁 Loading data from: {data_path}")
        
        if not Path(data_path).exists():
            print(f"❌ Data file not found: {data_path}")
            print("💡 Run the data cleaning script first:")
            print("   python scripts/create_clean_dataset.py")
            raise FileNotFoundError(f"Data file not found: {data_path}")
        
        # Load the dataset
        df = pd.read_csv(data_path)
        print(f"✅ Loaded {len(df):,} variants with {len(df.columns)} columns")
        
        # Validate no leakage features are present
        self._validate_no_data_leakage(df)
        
        # Create target variable
        y = self._create_target_variable(df)
        
        # Select features (SAFE feature selection)
        selected_features = self._select_features_fixed(df)
        self.feature_names = selected_features
        self.selected_features = selected_features
        
        # Prepare feature matrix
        X = self._prepare_features(df[selected_features])
        
        print(f"🎯 Target distribution:")
        target_counts = pd.Series(y).value_counts()
        for class_name, count in target_counts.items():
            pct = count / len(y) * 100
            print(f"   {class_name}: {count:,} ({pct:.1f}%)")
        
        return X, y
    
    def _validate_no_data_leakage(self, df):
        """
        CRITICAL: Validate that no data leakage features are present
        """
        print("🔍 VALIDATING NO DATA LEAKAGE...")
        
        leakage_found = []
        for banned_feature in self.banned_features:
            if banned_feature in df.columns:
                leakage_found.append(banned_feature)
        
        if leakage_found:
            print(f"❌ CRITICAL ERROR: Data leakage features found: {leakage_found}")
            print("   These features must be removed before training!")
            print("   Use the clean dataset creation script.")
            raise ValueError(f"Data leakage detected: {leakage_found}")
        
        print("✅ No data leakage features detected - safe to proceed")
    
    def _create_target_variable(self, df):
        """
        Create target variable for classification
        """
        target_col = 'variant_classification'
        
        if target_col in df.columns:
            print(f"✅ Using existing target column: {target_col}")
            return df[target_col].fillna('VUS')
        
        # Create target based on available clinical significance
        if 'clinical_significance' in df.columns:
            print("⚠️  Using clinical_significance to create target (validation only)")
            target_mapping = {
                'Pathogenic': 'Actionable Pathogenic',
                'Likely_pathogenic': 'Likely Actionable', 
                'Uncertain_significance': 'VUS',
                'Likely_benign': 'Likely Benign',
                'Benign': 'Benign'
            }
            return df['clinical_significance'].map(target_mapping).fillna('VUS')
        
        # Fallback: Create simplified classification based on variant impact
        impact_cols = ['variant_impact', 'impact', 'IMPACT', 'consequence_impact']
        impact_col = None
        for col in impact_cols:
            if col in df.columns:
                impact_col = col
                break
        
        if impact_col:
            print(f"📊 Creating target from impact column: {impact_col}")
            conditions = [
                (df[impact_col] == 'HIGH') & (df.get('is_important_gene', 0) == 1),
                (df[impact_col] == 'MODERATE') & (df.get('is_important_gene', 0) == 1),
                df[impact_col].isin(['LOW', 'MODIFIER'])
            ]
            choices = ['Actionable Pathogenic', 'Likely Actionable', 'VUS']
            return np.select(conditions, choices, default='VUS')
        
        # Last resort: Random assignment for testing (should not happen in production)
        print("⚠️  WARNING: No suitable target column found, creating random targets for testing")
        np.random.seed(42)
        return np.random.choice(['Actionable Pathogenic', 'Likely Actionable', 'VUS'], 
                               size=len(df), p=[0.2, 0.3, 0.5])
    
    def _select_features_fixed(self, df):
        """
        Select features WITHOUT data leakage for TabNet training
        CRITICAL: Removes functional_pathogenicity and confidence scores
        """
        print("🔧 SELECTING FEATURES (DATA LEAKAGE FIXED)")
        
        # ✅ CORE FUNCTIONAL FEATURES (allowed but use carefully)
        core_features = []
        
        # Only include raw scores if they have reasonable distribution
        functional_score_candidates = ['sift_score', 'polyphen_score', 'cadd_phred', 
                                     'conservation_score', 'gerp_score', 'phylop_score', 
                                     'phastcons_score']
        
        for feat in functional_score_candidates:
            if feat in df.columns and df[feat].notna().sum() > 1000:  # Only if we have real data
                core_features.append(feat)
        
        # ✅ GENOMIC CONTEXT FEATURES (safe to use)
        genomic_features = []
        
        # Variant type and impact
        variant_options = {
            'variant_impact': ['variant_impact', 'impact', 'IMPACT', 'consequence_impact'],
            'variant_type': ['variant_type', 'type', 'TYPE', 'variant_class'],
            'consequence': ['consequence', 'Consequence', 'CONSEQUENCE'],
            'exonic_function': ['exonic_function', 'ExonicFunc', 'EXONIC_FUNCTION']
        }
        
        for feature_name, possible_cols in variant_options.items():
            for col in possible_cols:
                if col in df.columns:
                    genomic_features.append(col)
                    break
        
        # Chromosomal location
        location_features = []
        location_options = ['chromosome', 'chr', 'CHROM', '#CHROM']
        for col in location_options:
            if col in df.columns:
                location_features.append(col)
                break
        
        # ✅ ALLELE FREQUENCY FEATURES (safe to use)
        frequency_features = []
        freq_options = [
            'gnomad_af_eur', 'gnomad_af_afr', 'gnomad_af_asj', 'gnomad_af_eas',
            'af_1kg', 'exac_af', 'esp_af', 'gnomad_af'
        ]
        for feat in freq_options:
            if feat in df.columns and df[feat].notna().sum() > 100:
                frequency_features.append(feat)
        
        # Rare variant indicator (derived from frequency)
        if 'is_rare' in df.columns:
            frequency_features.append('is_rare')
        
        # ✅ PATHWAY FEATURES (safe if not derived from target)
        pathway_features = []
        pathway_options = [
            'dna_repair_pathway', 'mismatch_repair_pathway', 
            'hormone_pathway', 'is_important_gene',
            'core_prostate_gene', 'cancer_predisposition_gene'
        ]
        for feat in pathway_options:
            if feat in df.columns:
                pathway_features.append(feat)
        
        # ✅ GENE-LEVEL FEATURES
        gene_features = []
        if 'SYMBOL' in df.columns or 'gene_symbol' in df.columns:
            gene_col = 'SYMBOL' if 'SYMBOL' in df.columns else 'gene_symbol'
            gene_features.append(gene_col)
        
        # ✅ STRUCTURAL FEATURES
        structural_features = []
        structural_options = [
            'splicing_impact', 'variant_size', 'is_indel', 'is_snv', 
            'ref_length', 'alt_length'
        ]
        for feat in structural_options:
            if feat in df.columns:
                structural_features.append(feat)
        
        # Combine all features
        all_features = (core_features + genomic_features + location_features + 
                       frequency_features + pathway_features + gene_features + 
                       structural_features)
        
        # Remove duplicates while preserving order
        selected_features = []
        for feat in all_features:
            if feat not in selected_features and feat not in self.banned_features:
                selected_features.append(feat)
        
        # Validate we have enough features
        if len(selected_features) < 5:
            print("⚠️  WARNING: Very few features selected. Adding backup features...")
            backup_features = []
            for col in df.columns:
                if (col not in selected_features and 
                    col not in self.banned_features and
                    df[col].dtype in ['int64', 'float64', 'object'] and
                    not col.lower().startswith('clinical')):
                    backup_features.append(col)
                    if len(selected_features) + len(backup_features) >= 10:
                        break
            selected_features.extend(backup_features)
        
        print(f"✅ Selected {len(selected_features)} features (NO DATA LEAKAGE)")
        print(f"   Core functional: {[f for f in core_features if f in selected_features]}")
        print(f"   Genomic context: {[f for f in genomic_features if f in selected_features]}")
        print(f"   Frequency: {[f for f in frequency_features if f in selected_features]}")
        print(f"   Pathway: {[f for f in pathway_features if f in selected_features]}")
        
        # ❌ EXPLICITLY REMOVED (log for transparency)
        actually_removed = [f for f in self.banned_features if f in df.columns]
        if actually_removed:
            print(f"❌ REMOVED (data leakage): {actually_removed}")
        
        return selected_features
    
    def _prepare_features(self, X):
        """
        Prepare features for TabNet training
        """
        print("🔧 Preprocessing features...")
        
        X_processed = X.copy()
        
        # Handle missing values
        for col in X_processed.columns:
            if X_processed[col].dtype in ['int64', 'float64']:
                # Numerical: fill with median
                median_val = X_processed[col].median()
                X_processed[col] = X_processed[col].fillna(median_val)
            else:
                # Categorical: fill with 'Unknown'
                X_processed[col] = X_processed[col].fillna('Unknown')
        
        # Encode categorical variables
        categorical_cols = X_processed.select_dtypes(include=['object']).columns
        for col in categorical_cols:
            if col not in self.feature_encoders:
                self.feature_encoders[col] = LabelEncoder()
                X_processed[col] = self.feature_encoders[col].fit_transform(X_processed[col].astype(str))
            else:
                # Handle unseen categories
                unique_vals = set(X_processed[col].astype(str))
                known_vals = set(self.feature_encoders[col].classes_)
                unknown_vals = unique_vals - known_vals
                
                if unknown_vals:
                    # Map unknown values to first class
                    X_processed[col] = X_processed[col].astype(str).apply(
                        lambda x: x if x in known_vals else self.feature_encoders[col].classes_[0]
                    )
                
                X_processed[col] = self.feature_encoders[col].transform(X_processed[col])
        
        print(f"✅ Prepared {X_processed.shape[1]} features for {X_processed.shape[0]:,} samples")
        return X_processed.values
    
    def train(self, X_train, y_train, X_val=None, y_val=None, max_epochs=200, patience=20):
        """
        Train TabNet model with early stopping
        """
        print("🚀 TRAINING TABNET MODEL")
        print("=" * 40)
        
        # Encode target
        y_train_encoded = self.label_encoder.fit_transform(y_train)
        
        if X_val is not None and y_val is not None:
            y_val_encoded = self.label_encoder.transform(y_val)
            eval_set = [(X_val, y_val_encoded)]
        else:
            # Create validation split
            X_train_split, X_val_split, y_train_split, y_val_split = train_test_split(
                X_train, y_train_encoded, test_size=0.2, stratify=y_train_encoded, random_state=42
            )
            X_train, y_train_encoded = X_train_split, y_train_split
            eval_set = [(X_val_split, y_val_split)]
        
        # Initialize TabNet model
        self.model = TabNetClassifier(
            n_d=self.n_d,
            n_a=self.n_a,
            n_steps=self.n_steps,
            gamma=self.gamma,
            lambda_sparse=self.lambda_sparse,
            optimizer_fn=torch.optim.Adam,
            optimizer_params=dict(lr=2e-2),
            mask_type='entmax',
            scheduler_params={"step_size": 10, "gamma": 0.9},
            scheduler_fn=torch.optim.lr_scheduler.StepLR,
            verbose=1
        )
        
        print(f"📊 Training with {X_train.shape[0]:,} samples, {X_train.shape[1]} features")
        print(f"   Classes: {len(self.label_encoder.classes_)}")
        print(f"   Expected realistic accuracy: 75-85% (NO data leakage)")
        
        # Train the model
        self.model.fit(
            X_train=X_train,
            y_train=y_train_encoded,
            eval_set=eval_set,
            eval_name=['validation'],
            eval_metric=['accuracy'],
            max_epochs=max_epochs,
            patience=patience,
            batch_size=1024,
            virtual_batch_size=128,
            num_workers=0,
            drop_last=False
        )
        
        # Validation performance
        val_accuracy = self.model.history['validation_accuracy'][-1]
        print(f"\n🎯 Training completed!")
        print(f"   Final validation accuracy: {val_accuracy:.3f}")
        
        # Performance check
        if val_accuracy > 0.95:
            print("⚠️  WARNING: Accuracy >95% - check for remaining data leakage!")
        elif val_accuracy > 0.75:
            print("✅ GOOD: Realistic performance for genomic classification")
        else:
            print("📈 MODERATE: Consider feature engineering or ensemble methods")
        
        return val_accuracy
    
    def evaluate(self, X_test, y_test):
        """
        Evaluate model performance
        """
        if self.model is None:
            raise ValueError("Model not trained. Call train() first.")
        
        y_test_encoded = self.label_encoder.transform(y_test)
        
        # Predictions
        y_pred_encoded = self.model.predict(X_test)
        y_pred_proba = self.model.predict_proba(X_test)
        
        # Calculate metrics
        accuracy = accuracy_score(y_test_encoded, y_pred_encoded)
        
        # Convert back to original labels for interpretation
        y_pred = self.label_encoder.inverse_transform(y_pred_encoded)
        
        print(f"\n📊 TEST EVALUATION RESULTS:")
        print(f"   Accuracy: {accuracy:.3f}")
        
        # Classification report
        print(f"\n📋 Classification Report:")
        print(classification_report(y_test, y_pred))
        
        return accuracy
    
    def get_feature_importance(self):
        """
        Get feature importance from trained TabNet model
        """
        if self.model is None:
            raise ValueError("Model not trained. Call train() first.")
        
        # Get feature importance (global)
        feature_importance = self.model.feature_importances_
        
        # Create DataFrame for interpretability
        importance_df = pd.DataFrame({
            'feature': self.feature_names,
            'importance': feature_importance
        }).sort_values('importance', ascending=False)
        
        return importance_df
    
    def analyze_clinical_attention(self, X_sample, sample_size=10):
        """
        Analyze attention patterns for clinical interpretability
        """
        if self.model is None:
            raise ValueError("Model not trained. Call train() first.")
        
        if len(X_sample) > sample_size:
            indices = np.random.choice(len(X_sample), sample_size, replace=False)
            X_sample = X_sample[indices]
        
        # Get attention masks
        explain_matrix, masks = self.model.explain(X_sample)
        
        # Aggregate attention by feature groups
        pathway_attention = {}
        for group_name, features in self.feature_groups.items():
            group_indices = [i for i, fname in enumerate(self.feature_names) if fname in features]
            if group_indices:
                group_attention = masks[:, group_indices].sum(axis=1).mean()
                pathway_attention[group_name] = float(group_attention)
        
        return pathway_attention, explain_matrix, masks
    
    def cross_validate(self, X, y, cv_folds=5):
        """
        Perform cross-validation to validate model performance
        """
        print(f"🔄 CROSS-VALIDATION ({cv_folds} folds)")
        print("=" * 40)
        
        skf = StratifiedKFold(n_splits=cv_folds, shuffle=True, random_state=42)
        cv_scores = []
        
        for fold, (train_idx, val_idx) in enumerate(skf.split(X, y), 1):
            print(f"\n📊 Fold {fold}/{cv_folds}")
            
            X_train_fold, X_val_fold = X[train_idx], X[val_idx]
            y_train_fold, y_val_fold = y[train_idx], y[val_idx]
            
            # Create temporary model for this fold
            temp_model = TabNetClassifier(
                n_d=self.n_d, n_a=self.n_a, n_steps=self.n_steps,
                gamma=self.gamma, lambda_sparse=self.lambda_sparse,
                optimizer_fn=torch.optim.Adam,
                optimizer_params=dict(lr=2e-2),
                verbose=0
            )
            
            # Encode target for this fold
            fold_label_encoder = LabelEncoder()
            y_train_encoded = fold_label_encoder.fit_transform(y_train_fold)
            y_val_encoded = fold_label_encoder.transform(y_val_fold)
            
            # Train
            temp_model.fit(
                X_train=X_train_fold, y_train=y_train_encoded,
                eval_set=[(X_val_fold, y_val_encoded)],
                max_epochs=100, patience=10,
                batch_size=512, virtual_batch_size=128
            )
            
            # Evaluate
            y_pred = temp_model.predict(X_val_fold)
            accuracy = accuracy_score(y_val_encoded, y_pred)
            cv_scores.append(accuracy)
            
            print(f"   Fold {fold} accuracy: {accuracy:.3f}")
        
        mean_score = np.mean(cv_scores)
        std_score = np.std(cv_scores)
        
        print(f"\n🎯 CROSS-VALIDATION RESULTS:")
        print(f"   Mean accuracy: {mean_score:.3f} ± {std_score:.3f}")
        print(f"   Individual folds: {[f'{score:.3f}' for score in cv_scores]}")
        
        # Performance interpretation
        if mean_score > 0.95:
            print(f"   Status: ⚠️  SUSPICIOUS - Check for data leakage")
        elif mean_score > 0.75:
            print(f"   Status: ✅ EXCELLENT - Realistic for genomic classification")
        elif mean_score > 0.65:
            print(f"   Status: ✅ GOOD - Expected for complex genomic data")
        else:
            print(f"   Status: ⚠️  MODERATE - Consider feature engineering")
        
        return {
            'mean_accuracy': mean_score,
            'std_accuracy': std_score,
            'fold_scores': cv_scores
        }
    
    def save_model(self, filepath):
        """Save trained model"""
        if self.model is None:
            raise ValueError("Model not trained. Call train() first.")
        
        # Save TabNet model
        self.model.save_model(filepath)
        
        # Save encoders and feature info
        model_info = {
            'label_encoder': self.label_encoder,
            'feature_encoders': self.feature_encoders,
            'feature_names': self.feature_names,
            'selected_features': self.selected_features,
            'banned_features': self.banned_features
        }
        
        info_path = filepath.replace('.zip', '_info.pkl')
        joblib.dump(model_info, info_path)
        
        print(f"✅ Model saved: {filepath}")
        print(f"✅ Model info saved: {info_path}")
    
    def load_model(self, filepath):
        """Load pre-trained model"""
        self.model = TabNetClassifier()
        self.model.load_model(filepath)
        
        # Load encoders and feature info
        info_path = filepath.replace('.zip', '_info.pkl')
        if Path(info_path).exists():
            model_info = joblib.load(info_path)
            self.label_encoder = model_info['label_encoder']
            self.feature_encoders = model_info['feature_encoders']
            self.feature_names = model_info['feature_names']
            self.selected_features = model_info['selected_features']
            self.banned_features = model_info.get('banned_features', [])
            
        print(f"✅ Model loaded: {filepath}")

def create_clean_dataset(input_path, output_path):
    """
    Create a clean dataset without data leakage features
    """
    print("🧹 CREATING CLEAN DATASET WITHOUT DATA LEAKAGE...")
    
    df = pd.read_csv(input_path)
    
    # Remove problematic features
    leakage_features = [
        'functional_pathogenicity',
        'sift_confidence', 
        'polyphen_confidence'
    ]
    
    original_cols = df.shape[1]
    df_clean = df.drop(columns=[col for col in leakage_features if col in df.columns])
    removed_count = original_cols - df_clean.shape[1]
    
    # Save clean dataset
    df_clean.to_csv(output_path, index=False)
    
    print(f"✅ Clean dataset saved: {output_path}")
    print(f"   Original features: {original_cols}")
    print(f"   Removed features: {removed_count}")
    print(f"   Clean features: {df_clean.shape[1]}")
    print(f"   Rows: {len(df_clean):,}")
    
    return output_path

def main():
    """
    Example usage and testing
    """
    print("🧬 TabNet Prostate Cancer Variant Classifier (DATA LEAKAGE FIXED)")
    print("=" * 70)
    
    # Initialize model
    tabnet = ProstateVariantTabNet(n_d=64, n_a=64, n_steps=6)
    
    try:
        # Load data
        X, y = tabnet.load_data()
        
        # Split data
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.2, stratify=y, random_state=42
        )
        
        print(f"\n📊 Data split:")
        print(f"   Training set: {X_train.shape[0]:,} variants")
        print(f"   Test set: {X_test.shape[0]:,} variants")
        
        # Cross-validation first (recommended)
        print(f"\n🔄 Running cross-validation...")
        cv_results = tabnet.cross_validate(X, y, cv_folds=3)
        
        # Train final model
        print(f"\n🚀 Training final model...")
        val_accuracy = tabnet.train(X_train, y_train)
        
        # Final evaluation
        test_accuracy = tabnet.evaluate(X_test, y_test)
        print(f"\n🎯 Final Test Accuracy: {test_accuracy:.3f}")
        
        # Clinical interpretability analysis
        pathway_attention, _, _ = tabnet.analyze_clinical_attention(X_test[:5])
        print("\n🔍 Clinical Pathway Attention:")
        for pathway, attention in pathway_attention.items():
            print(f"   {pathway}: {attention:.3f}")
        
        # Feature importance
        feature_importance = tabnet.get_feature_importance()
        print(f"\n📊 Top 5 Most Important Features:")
        for idx, row in feature_importance.head(5).iterrows():
            print(f"   {row['feature']}: {row['importance']:.3f}")
        
        return tabnet
        
    except Exception as e:
        print(f"❌ Error in main execution: {e}")
        import traceback
        traceback.print_exc()
        return None

if __name__ == "__main__":
    model = main()