import torch
import torch.nn as nn
import pandas as pd
import numpy as np
from pytorch_tabnet.tab_model import TabNetClassifier
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.metrics import accuracy_score, classification_report
import matplotlib.pyplot as plt
import seaborn as sns

class ProstateVariantTabNet:
    def __init__(self, n_d=64, n_a=64, n_steps=6, gamma=1.3, lambda_sparse=1e-3):
        """
        TabNet for Prostate Cancer Variant Classification
        
        Args:
            n_d: Width of the decision prediction layer
            n_a: Width of the attention embedding for each mask
            n_steps: Number of decision steps
            gamma: Coefficient for feature reusage in the masks
            lambda_sparse: Sparsity regularization coefficient
        """
        self.n_d = n_d
        self.n_a = n_a  
        self.n_steps = n_steps
        self.gamma = gamma
        self.lambda_sparse = lambda_sparse
        
        self.model = None
        self.feature_scaler = StandardScaler()
        self.label_encoder = LabelEncoder()
        
        # Feature groups for interpretability
        self.feature_groups = {
            'genomic_context': list(range(0, 20)),      # Conservation, functional predictions
            'prostate_biology': list(range(20, 45)),    # AR pathway, DNA repair
            'therapeutic_context': list(range(45, 65)), # PARP inhibitor, hormone therapy
            'clinical_annotations': list(range(65, 80)) # Tumor stage, PSA, survival
        }
        
    def load_merged_datasets(self, tcga_path, cosmic_path, clinvar_path):
        """
        Load and merge the three primary datasets
        Your colleague's merged dataset should have these columns
        """
        # Placeholder for merged dataset loading
        # This will be replaced with your colleague's merger output
        print("Loading merged datasets...")
        
        # Expected merged dataset structure:
        expected_columns = [
            # Genomic Context (20 features)
            'phylop_score', 'gerp_score', 'cadd_score', 'sift_score', 'polyphen_score',
            'gnomad_af_total', 'gnomad_af_eur', 'gnomad_af_afr', 'conservation_score',
            'regulatory_region', 'exonic_function', 'splicing_impact', 'lof_prediction',
            'missense_prediction', 'variant_type', 'chromosome', 'position', 'ref_alt',
            'gene_symbol', 'transcript_id',
            
            # Prostate Biology (25 features) 
            'ar_pathway_score', 'dna_repair_score', 'pi3k_pathway_score', 'pten_impact',
            'brca1_interaction', 'brca2_interaction', 'atm_interaction', 'palb2_interaction',
            'chek2_interaction', 'rad51_interaction', 'prostate_expression_level',
            'tissue_specificity', 'hormone_sensitivity', 'age_correlation', 'gleason_correlation',
            'psa_correlation', 'metastasis_correlation', 'recurrence_score', 'pathway_disruption',
            'protein_domain_impact', 'structural_impact', 'binding_affinity_change',
            'enzymatic_activity', 'cellular_localization', 'protein_stability',
            
            # Therapeutic Context (20 features)
            'parp_inhibitor_score', 'olaparib_response', 'rucaparib_response', 'niraparib_response',
            'hormone_therapy_score', 'adt_response', 'enzalutamide_response', 'abiraterone_response',
            'immunotherapy_score', 'msi_status', 'tmb_score', 'neoantigen_load',
            'pd1_expression', 'pdl1_expression', 'clinical_trial_eligibility', 'drug_target_score',
            'resistance_prediction', 'combination_therapy_score', 'biomarker_status', 'actionability_score',
            
            # Clinical Annotations (15 features)
            'tumor_stage', 'gleason_grade', 'psa_level', 'age_diagnosis', 'family_history',
            'treatment_response', 'survival_months', 'progression_status', 'metastasis_sites',
            'prior_treatments', 'comorbidities', 'performance_status', 'risk_stratification',
            'followup_time', 'outcome_status',
            
            # Target variable
            'variant_classification'  # 5 classes: Actionable Pathogenic, Likely Actionable, VUS, Likely Benign, Benign
        ]
        
        return expected_columns
    
    def preprocess_features(self, df):
        """
        Preprocess the merged dataset for TabNet training
        """
        # Separate features and target
        feature_cols = [col for col in df.columns if col != 'variant_classification']
        X = df[feature_cols].copy()
        y = df['variant_classification'].copy()
        
        # Handle missing values
        X = X.fillna(X.median())
        
        # Encode categorical features
        categorical_cols = X.select_dtypes(include=['object']).columns
        for col in categorical_cols:
            le = LabelEncoder()
            X[col] = le.fit_transform(X[col].astype(str))
        
        # Scale features
        X_scaled = self.feature_scaler.fit_transform(X)
        
        # Encode target labels
        y_encoded = self.label_encoder.fit_transform(y)
        
        return X_scaled, y_encoded, feature_cols
    
    def build_model(self):
        """
        Initialize TabNet model with prostate cancer-specific configuration
        """
        self.model = TabNetClassifier(
            n_d=self.n_d,
            n_a=self.n_a,
            n_steps=self.n_steps,
            gamma=self.gamma,
            lambda_sparse=self.lambda_sparse,
            optimizer_fn=torch.optim.Adam,
            optimizer_params=dict(lr=2e-2, weight_decay=1e-5),
            mask_type='entmax',  # Sparse attention masks
            scheduler_params=dict(step_size=50, gamma=0.9),
            scheduler_fn=torch.optim.lr_scheduler.StepLR,
            verbose=1,
            device_name='cuda' if torch.cuda.is_available() else 'cpu'
        )
        
        return self.model
    
    def train(self, X_train, y_train, X_val, y_val, max_epochs=200, patience=20):
        """
        Train TabNet model with early stopping
        """
        if self.model is None:
            self.build_model()
            
        self.model.fit(
            X_train=X_train, y_train=y_train,
            eval_set=[(X_val, y_val)],
            eval_name=['val'],
            max_epochs=max_epochs,
            patience=patience,
            batch_size=512,  # Optimized for H100 GPU
            virtual_batch_size=128,
            num_workers=4,
            drop_last=False
        )
        
        return self.model
    
    def get_feature_importance(self):
        """
        Extract feature importance from trained TabNet
        """
        if self.model is None:
            raise ValueError("Model must be trained first")
            
        importance = self.model.feature_importances_
        return importance
    
    def get_attention_masks(self, X):
        """
        Extract attention masks for interpretability analysis
        """
        if self.model is None:
            raise ValueError("Model must be trained first")
            
        # Get attention masks from each decision step
        explain_matrix, masks = self.model.explain(X)
        
        return explain_matrix, masks
    
    def analyze_clinical_patterns(self, X, feature_names, patient_ids=None):
        """
        Analyze attention patterns for clinical interpretability
        """
        explain_matrix, masks = self.get_attention_masks(X)
        
        # Group attention by feature categories
        group_attention = {}
        for group_name, indices in self.feature_groups.items():
            group_attention[group_name] = explain_matrix[:, indices].sum(axis=1)
        
        # Create interpretability report
        interpretation = {
            'feature_importance': self.get_feature_importance(),
            'group_attention': group_attention,
            'decision_masks': masks,
            'feature_groups': self.feature_groups
        }
        
        return interpretation
    
    def predict_with_explanation(self, X, feature_names):
        """
        Generate predictions with clinical explanations
        """
        if self.model is None:
            raise ValueError("Model must be trained first")
            
        # Get predictions
        predictions = self.model.predict(X)
        probabilities = self.model.predict_proba(X)
        
        # Get explanations
        interpretation = self.analyze_clinical_patterns(X, feature_names)
        
        # Convert to clinical labels
        clinical_labels = self.label_encoder.inverse_transform(predictions)
        
        results = {
            'predictions': clinical_labels,
            'probabilities': probabilities,
            'interpretation': interpretation
        }
        
        return results
    
    def visualize_attention(self, X, feature_names, sample_idx=0):
        """
        Visualize attention patterns for a specific sample
        """
        explain_matrix, masks = self.get_attention_masks(X[sample_idx:sample_idx+1])
        
        # Create attention heatmap
        plt.figure(figsize=(15, 8))
        
        # Plot feature importance
        plt.subplot(2, 2, 1)
        feature_imp = self.get_feature_importance()
        top_features = np.argsort(feature_imp)[-20:]
        plt.barh(range(len(top_features)), feature_imp[top_features])
        plt.yticks(range(len(top_features)), [feature_names[i] for i in top_features])
        plt.title('Top 20 Feature Importance')
        plt.xlabel('Importance Score')
        
        # Plot attention by group
        plt.subplot(2, 2, 2)
        group_attention = {}
        for group_name, indices in self.feature_groups.items():
            group_attention[group_name] = explain_matrix[0, indices].sum()
        
        plt.bar(group_attention.keys(), group_attention.values())
        plt.title('Attention by Feature Group')
        plt.ylabel('Attention Weight')
        plt.xticks(rotation=45)
        
        # Plot decision step masks
        plt.subplot(2, 1, 2)
        plt.imshow(masks[0].T, aspect='auto', cmap='Blues')
        plt.title('Decision Step Attention Masks')
        plt.xlabel('Decision Steps')
        plt.ylabel('Features')
        plt.colorbar()
        
        plt.tight_layout()
        plt.show()
        
        return plt.gcf()

# Example usage workflow
def main():
    """
    Example workflow for training TabNet on prostate cancer variants
    """
    # Initialize model
    tabnet_model = ProstateVariantTabNet(
        n_d=64,
        n_a=64, 
        n_steps=6,
        gamma=1.3,
        lambda_sparse=1e-3
    )
    
    # Load your colleague's merged dataset
    # df = pd.read_csv('merged_prostate_variants.csv')
    
    # For now, create placeholder for testing structure
    print("TabNet Prostate Cancer Variant Classifier initialized")
    print("Ready to train once merged datasets are available")
    print(f"Expected feature groups: {tabnet_model.feature_groups}")
    
    return tabnet_model

if __name__ == "__main__":
    model = main()