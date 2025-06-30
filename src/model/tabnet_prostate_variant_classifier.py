#!/usr/bin/env python3
"""
TabNet Prostate Cancer Variant Classification
Interpretable deep learning for prostate cancer genomic variant classification
using TabNet with attention mechanisms and functional score confidence.
"""

import pandas as pd
import numpy as np
import torch
from pytorch_tabnet.tab_model import TabNetClassifier
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.metrics import accuracy_score, roc_auc_score, classification_report
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

class ProstateVariantTabNet:
    """
    TabNet model for prostate cancer variant classification with interpretability
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
        
        # Clinical feature groups for interpretability
        self.feature_groups = {
            'functional_scores': ['sift_score', 'polyphen_score', 'cadd_phred', 
                                'sift_confidence', 'polyphen_confidence', 'functional_pathogenicity'],
            'genomic_context': ['conservation_score', 'gnomad_af_eur', 'gnomad_af_afr', 
                              'regulatory_region', 'variant_type', 'chromosome'],
            'pathway_indicators': ['dna_repair_pathway', 'mismatch_repair_pathway', 
                                 'hormone_pathway', 'is_important_gene'],
            'clinical_impact': ['variant_impact', 'exonic_function', 'splicing_impact']
        }
        
        self.model = None
        self.feature_names = None
        self.label_encoder = LabelEncoder()
        self.scaler = StandardScaler()
        
    def load_data(self, csv_path="/u/aa107/uiuc-cancer-research/data/processed/tabnet_csv/prostate_variants_tabnet_imputed.csv"):
        """
        Load and preprocess the enhanced prostate cancer variant dataset
        
        Args:
            csv_path: Path to the imputed dataset CSV
            
        Returns:
            X: Feature matrix
            y: Target labels
        """
        print(f"Loading dataset from: {csv_path}")
        
        try:
            df = pd.read_csv(csv_path)
            print(f"✅ Loaded {len(df):,} variants with {df.shape[1]} features")
            
            # Define target variable (adjust based on actual column name)
            target_col = 'variant_classification'
            if target_col not in df.columns:
                # Fallback: create target from clinical significance
                df[target_col] = self._create_target_from_clinical_significance(df)
            
            # Select features for TabNet
            feature_cols = self._select_features(df)
            X = df[feature_cols].copy()
            y = df[target_col].copy()
            
            # Preprocess features
            X = self._preprocess_features(X)
            y = self.label_encoder.fit_transform(y)
            
            self.feature_names = feature_cols
            print(f"✅ Preprocessed to {X.shape[1]} features, {len(np.unique(y))} classes")
            
            return X, y
            
        except Exception as e:
            print(f"❌ Error loading data: {e}")
            raise
    
    def _create_target_from_clinical_significance(self, df):
        """Create 5-class target from clinical significance if not available"""
        target_mapping = {
            'pathogenic': 'Actionable Pathogenic',
            'likely_pathogenic': 'Likely Actionable', 
            'uncertain_significance': 'VUS',
            'likely_benign': 'Likely Benign',
            'benign': 'Benign'
        }
        
        # Use clinical_significance if available
        if 'clinical_significance' in df.columns:
            return df['clinical_significance'].map(target_mapping).fillna('VUS')
        
        # Check for various impact column names in your dataset
        impact_cols = ['variant_impact', 'impact', 'IMPACT', 'consequence_impact', 'vep_impact']
        impact_col = None
        for col in impact_cols:
            if col in df.columns:
                impact_col = col
                break
        
        if impact_col:
            # Create target based on available impact column
            conditions = [
                (df[impact_col] == 'HIGH') & (df.get('is_important_gene', 0) == 1),
                (df[impact_col] == 'MODERATE') & (df.get('is_important_gene', 0) == 1),
                df[impact_col].isin(['LOW', 'MODIFIER'])
            ]
            choices = ['Actionable Pathogenic', 'Likely Actionable', 'VUS']
            return np.select(conditions, choices, default='VUS')
        else:
            # Fallback: create simplified classification based on functional scores
            print("  ⚠️  No impact column found, using functional scores for classification")
            conditions = [
                (df['sift_score'] <= 0.05) & (df['polyphen_score'] >= 0.85),
                (df['sift_score'] <= 0.05) | (df['polyphen_score'] >= 0.85),
                (df['sift_score'] > 0.05) & (df['polyphen_score'] < 0.85)
            ]
            choices = ['Actionable Pathogenic', 'Likely Actionable', 'VUS']
            return np.select(conditions, choices, default='VUS')
    
    def _select_features(self, df):
        """Select relevant features for TabNet training"""
        # Core functional features (must have)
        core_features = ['sift_score', 'polyphen_score', 'cadd_phred']
        
        # Confidence features (from imputation breakthrough)
        confidence_features = ['sift_confidence', 'polyphen_confidence', 'functional_pathogenicity']
        
        # Pathway indicators
        pathway_features = ['dna_repair_pathway', 'mismatch_repair_pathway', 'hormone_pathway', 'is_important_gene']
        
        # Additional genomic context - check multiple possible column names
        context_feature_options = {
            'variant_impact': ['variant_impact', 'impact', 'IMPACT', 'consequence_impact', 'vep_impact'],
            'variant_type': ['variant_type', 'type', 'TYPE', 'variant_class', 'allele'],
            'chromosome': ['chromosome', 'chr', 'CHROM', '#CHROM', 'seqname'],
            'gnomad_af': ['gnomad_af_eur', 'gnomad_af', 'AF', 'af', 'frequency']
        }
        
        context_features = []
        for feature_name, possible_cols in context_feature_options.items():
            for col in possible_cols:
                if col in df.columns:
                    context_features.append(col)
                    break  # Take first match
        
        # Combine all features that exist in the dataset
        all_candidate_features = core_features + confidence_features + pathway_features + context_features
        available_features = [f for f in all_candidate_features if f in df.columns]
        
        # If we don't have enough features, add some backup features
        if len(available_features) < 10:
            backup_features = []
            for col in df.columns:
                if col not in available_features and df[col].dtype in ['int64', 'float64', 'object']:
                    backup_features.append(col)
                    if len(available_features) + len(backup_features) >= 15:
                        break
            available_features.extend(backup_features)
        
        print(f"✅ Selected {len(available_features)} features from {len(all_candidate_features)} candidates")
        print(f"   Core features: {[f for f in core_features if f in available_features]}")
        print(f"   Confidence features: {[f for f in confidence_features if f in available_features]}")
        
        return available_features
    
    def _preprocess_features(self, X):
        """Preprocess features for TabNet"""
        # Handle categorical variables
        categorical_cols = X.select_dtypes(include=['object']).columns
        for col in categorical_cols:
            le = LabelEncoder()
            X[col] = le.fit_transform(X[col].astype(str))
        
        # Fill missing values with median
        X = X.fillna(X.median())
        
        # Scale numerical features
        numerical_cols = X.select_dtypes(include=[np.number]).columns
        X[numerical_cols] = self.scaler.fit_transform(X[numerical_cols])
        
        return X.values
    
    def train(self, X_train, y_train, X_val=None, y_val=None, max_epochs=200, patience=15):
        """
        Train TabNet model with early stopping
        
        Args:
            X_train, y_train: Training data
            X_val, y_val: Validation data (optional)
            max_epochs: Maximum training epochs
            patience: Early stopping patience
        """
        print("🚀 Starting TabNet training...")
        
        # Create validation set if not provided
        if X_val is None or y_val is None:
            X_train, X_val, y_train, y_val = train_test_split(
                X_train, y_train, test_size=0.2, stratify=y_train, random_state=42
            )
        
        # Initialize TabNet
        device = 'cuda' if torch.cuda.is_available() else 'cpu'
        print(f"✅ Using device: {device}")
        
        self.model = TabNetClassifier(
            n_d=self.n_d,
            n_a=self.n_a,
            n_steps=self.n_steps,
            gamma=self.gamma,
            lambda_sparse=self.lambda_sparse,
            optimizer_fn=torch.optim.Adam,
            optimizer_params=dict(lr=2e-2),
            mask_type='entmax',
            device_name=device,
            verbose=1
        )
        
        # Train model
        self.model.fit(
            X_train, y_train,
            eval_set=[(X_val, y_val)],
            max_epochs=max_epochs,
            patience=patience,
            batch_size=1024,
            virtual_batch_size=128,
            num_workers=0,
            drop_last=False
        )
        
        print("✅ Training completed!")
        
        # Evaluate on validation set
        val_accuracy = self.evaluate(X_val, y_val)
        print(f"✅ Validation Accuracy: {val_accuracy:.3f}")
        
        return val_accuracy
    
    def predict(self, X):
        """Make predictions"""
        if self.model is None:
            raise ValueError("Model not trained. Call train() first.")
        return self.model.predict(X)
    
    def predict_proba(self, X):
        """Get prediction probabilities"""
        if self.model is None:
            raise ValueError("Model not trained. Call train() first.")
        return self.model.predict_proba(X)
    
    def evaluate(self, X_test, y_test):
        """Evaluate model performance"""
        predictions = self.predict(X_test)
        accuracy = accuracy_score(y_test, predictions)
        return accuracy
    
    def get_feature_importance(self):
        """Get feature importance scores"""
        if self.model is None:
            raise ValueError("Model not trained. Call train() first.")
        return self.model.feature_importances_
    
    def get_attention_masks(self, X):
        """Get attention masks for interpretability"""
        if self.model is None:
            raise ValueError("Model not trained. Call train() first.")
        
        explain_matrix, masks = self.model.explain(X)
        return explain_matrix, masks
    
    def analyze_clinical_attention(self, X, sample_indices=None):
        """
        Analyze attention patterns for clinical interpretability
        Focus on PARP inhibitor and hormone therapy relevant features
        """
        if sample_indices is None:
            sample_indices = [0]  # Analyze first sample by default
        
        explain_matrix, masks = self.get_attention_masks(X[sample_indices])
        
        # Group attention by clinical pathways
        pathway_attention = {}
        for group_name, feature_list in self.feature_groups.items():
            if self.feature_names:
                # Find indices of features in this group
                group_indices = [i for i, fname in enumerate(self.feature_names) if fname in feature_list]
                if group_indices:
                    pathway_attention[group_name] = explain_matrix[0, group_indices].sum()
        
        return pathway_attention, explain_matrix, masks
    
    def save_model(self, filepath):
        """Save trained model"""
        if self.model is None:
            raise ValueError("Model not trained. Call train() first.")
        self.model.save_model(filepath)
    
    def load_model(self, filepath):
        """Load pre-trained model"""
        self.model = TabNetClassifier()
        self.model.load_model(filepath)

def main():
    """
    Example usage and testing
    """
    print("🧬 TabNet Prostate Cancer Variant Classifier")
    print("=" * 60)
    
    # Initialize model
    tabnet = ProstateVariantTabNet(n_d=64, n_a=64, n_steps=6)
    
    try:
        # Load data
        X, y = tabnet.load_data()
        
        # Split data
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.2, stratify=y, random_state=42
        )
        
        print(f"Training set: {X_train.shape[0]:,} variants")
        print(f"Test set: {X_test.shape[0]:,} variants")
        
        # Train model
        val_accuracy = tabnet.train(X_train, y_train)
        
        # Final evaluation
        test_accuracy = tabnet.evaluate(X_test, y_test)
        print(f"🎯 Final Test Accuracy: {test_accuracy:.3f}")
        
        # Clinical interpretability analysis
        pathway_attention, _, _ = tabnet.analyze_clinical_attention(X_test[:5])
        print("\n🔍 Clinical Pathway Attention:")
        for pathway, attention in pathway_attention.items():
            print(f"  {pathway}: {attention:.3f}")
        
        # Feature importance
        feature_importance = tabnet.get_feature_importance()
        print(f"\n📊 Top feature importance: {feature_importance.max():.3f}")
        
        return tabnet
        
    except Exception as e:
        print(f"❌ Error in main execution: {e}")
        return None

if __name__ == "__main__":
    model = main()