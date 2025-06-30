#!/usr/bin/env python3
"""
TabNet Validation Framework for Prostate Cancer Variant Classification
Comprehensive validation with clinical interpretability analysis
"""

import pandas as pd
import numpy as np
import json
import time
from pathlib import Path
from datetime import datetime
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score, f1_score,
    roc_auc_score, classification_report, confusion_matrix
)
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

# Import our custom classes
import sys
sys.path.append('/u/aa107/uiuc-cancer-research/src/model')
sys.path.append('/u/aa107/uiuc-cancer-research')

from tabnet_prostate_variant_classifier import ProstateVariantTabNet
from config.pipeline_config import config

class TabNetValidator:
    """
    Comprehensive validation framework for TabNet prostate cancer classification
    Includes clinical interpretability analysis and statistical significance testing
    """
    
    def __init__(self, results_dir=None):
        """
        Initialize validator
        
        Args:
            results_dir: Directory to save validation results
        """
        self.results_dir = Path(results_dir) if results_dir else config.RESULTS_DIR
        self.results_dir.mkdir(parents=True, exist_ok=True)
        
        self.validation_results = {}
        self.clinical_analysis = {}
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
    def run_cross_validation(self, model, X, y, cv_folds=5):
        """
        Run stratified k-fold cross-validation
        
        Args:
            model: ProstateVariantTabNet instance
            X: Feature matrix
            y: Target labels
            cv_folds: Number of CV folds
            
        Returns:
            cv_results: Dictionary with CV metrics
        """
        print(f"🔄 Running {cv_folds}-fold cross-validation...")
        
        cv = StratifiedKFold(n_splits=cv_folds, shuffle=True, random_state=42)
        
        cv_results = {
            'fold_accuracies': [],
            'fold_precisions': [],
            'fold_recalls': [],
            'fold_f1s': [],
            'fold_aucs': [],
            'attention_patterns': [],
            'feature_importances': []
        }
        
        for fold, (train_idx, val_idx) in enumerate(cv.split(X, y)):
            print(f"  📊 Processing fold {fold + 1}/{cv_folds}...")
            
            # Split data
            X_train_fold, X_val_fold = X[train_idx], X[val_idx]
            y_train_fold, y_val_fold = y[train_idx], y[val_idx]
            
            # Create fresh model for each fold
            fold_model = ProstateVariantTabNet(
                n_d=model.n_d,
                n_a=model.n_a,
                n_steps=model.n_steps,
                gamma=model.gamma,
                lambda_sparse=model.lambda_sparse
            )
            
            # Train model
            fold_model.train(X_train_fold, y_train_fold, X_val_fold, y_val_fold)
            
            # Evaluate
            y_pred = fold_model.predict(X_val_fold)
            y_proba = fold_model.predict_proba(X_val_fold)
            
            # Calculate metrics
            accuracy = accuracy_score(y_val_fold, y_pred)
            precision = precision_score(y_val_fold, y_pred, average='macro', zero_division=0)
            recall = recall_score(y_val_fold, y_pred, average='macro', zero_division=0)
            f1 = f1_score(y_val_fold, y_pred, average='macro', zero_division=0)
            
            # ROC AUC for multiclass
            try:
                auc = roc_auc_score(y_val_fold, y_proba, multi_class='ovr', average='macro')
            except:
                auc = 0.0
            
            cv_results['fold_accuracies'].append(accuracy)
            cv_results['fold_precisions'].append(precision)
            cv_results['fold_recalls'].append(recall)
            cv_results['fold_f1s'].append(f1)
            cv_results['fold_aucs'].append(auc)
            
            # Collect interpretability data
            feature_importance = fold_model.get_feature_importance()
            cv_results['feature_importances'].append(feature_importance)
            
            # Analyze attention patterns for clinical interpretability
            attention_analysis, _, _ = fold_model.analyze_clinical_attention(X_val_fold[:100])
            cv_results['attention_patterns'].append(attention_analysis)
            
            print(f"    ✅ Fold {fold + 1} accuracy: {accuracy:.3f}")
        
        # Calculate summary statistics
        cv_results['mean_accuracy'] = np.mean(cv_results['fold_accuracies'])
        cv_results['std_accuracy'] = np.std(cv_results['fold_accuracies'])
        cv_results['mean_precision'] = np.mean(cv_results['fold_precisions'])
        cv_results['mean_recall'] = np.mean(cv_results['fold_recalls'])
        cv_results['mean_f1'] = np.mean(cv_results['fold_f1s'])
        cv_results['mean_auc'] = np.mean(cv_results['fold_aucs'])
        
        # Calculate confidence intervals
        cv_results['accuracy_ci'] = self._calculate_confidence_interval(cv_results['fold_accuracies'])
        
        print(f"✅ Cross-validation completed!")
        print(f"   Mean accuracy: {cv_results['mean_accuracy']:.3f} ± {cv_results['std_accuracy']:.3f}")
        print(f"   95% CI: [{cv_results['accuracy_ci'][0]:.3f}, {cv_results['accuracy_ci'][1]:.3f}]")
        
        return cv_results
    
    def _calculate_confidence_interval(self, scores, confidence=0.95):
        """Calculate confidence interval for metric scores"""
        n = len(scores)
        mean = np.mean(scores)
        std_err = stats.sem(scores)
        margin = std_err * stats.t.ppf((1 + confidence) / 2, n - 1)
        return [mean - margin, mean + margin]
    
    def analyze_clinical_interpretability(self, model, X, y, n_samples=100):
        """
        Analyze clinical interpretability of TabNet attention patterns
        Focus on PARP inhibitor and hormone therapy decision support
        """
        print("🔍 Analyzing clinical interpretability...")
        
        # Select diverse samples for analysis
        sample_indices = np.random.choice(len(X), min(n_samples, len(X)), replace=False)
        X_samples = X[sample_indices]
        y_samples = y[sample_indices]
        
        # Get attention patterns
        pathway_attention, explain_matrix, masks = model.analyze_clinical_attention(X_samples)
        
        clinical_analysis = {
            'pathway_attention_summary': pathway_attention,
            'attention_stability': {},
            'feature_importance_ranking': {},
            'clinical_relevance_scores': {}
        }
        
        # Analyze attention stability across samples
        all_attentions = []
        for i in range(min(50, len(X_samples))):  # Analyze subset for stability
            sample_attention, _, _ = model.analyze_clinical_attention(X_samples[i:i+1])
            all_attentions.append(sample_attention)
        
        # Calculate attention stability (coefficient of variation)
        for pathway in pathway_attention.keys():
            pathway_scores = [att.get(pathway, 0) for att in all_attentions]
            mean_score = np.mean(pathway_scores)
            std_score = np.std(pathway_scores)
            cv = std_score / mean_score if mean_score > 0 else 1.0
            clinical_analysis['attention_stability'][pathway] = {
                'mean': mean_score,
                'std': std_score,
                'coefficient_of_variation': cv,
                'stability_score': 1 - min(cv, 1.0)  # Higher is more stable
            }
        
        # Feature importance analysis
        feature_importance = model.get_feature_importance()
        if model.feature_names:
            importance_dict = dict(zip(model.feature_names, feature_importance))
            # Sort by importance
            sorted_features = sorted(importance_dict.items(), key=lambda x: x[1], reverse=True)
            clinical_analysis['feature_importance_ranking'] = dict(sorted_features[:20])
        
        # Clinical relevance scoring
        clinical_analysis['clinical_relevance_scores'] = self._score_clinical_relevance(
            pathway_attention, clinical_analysis['attention_stability']
        )
        
        print("✅ Clinical interpretability analysis completed")
        
        self.clinical_analysis = clinical_analysis
        return clinical_analysis
    
    def _score_clinical_relevance(self, pathway_attention, stability_analysis):
        """
        Score clinical relevance of attention patterns
        Based on alignment with known prostate cancer biology
        """
        relevance_weights = {
            'functional_scores': 0.30,      # Critical for pathogenicity assessment
            'pathway_indicators': 0.25,     # Important for therapeutic decisions
            'genomic_context': 0.20,        # Important for variant interpretation
            'clinical_impact': 0.25         # Important for clinical translation
        }
        
        clinical_scores = {}
        for pathway, weight in relevance_weights.items():
            attention_score = pathway_attention.get(pathway, 0)
            stability_score = stability_analysis.get(pathway, {}).get('stability_score', 0)
            
            # Combined score: attention strength × stability × clinical weight
            clinical_score = attention_score * stability_score * weight
            clinical_scores[pathway] = {
                'attention': attention_score,
                'stability': stability_score,
                'weight': weight,
                'clinical_score': clinical_score
            }
        
        # Overall clinical relevance score
        total_clinical_score = sum(scores['clinical_score'] for scores in clinical_scores.values())
        clinical_scores['overall_clinical_relevance'] = total_clinical_score
        
        return clinical_scores
    
    def evaluate_performance_targets(self, cv_results):
        """
        Evaluate if model meets performance targets
        """
        print("🎯 Evaluating performance targets...")
        
        target_accuracy = config.TARGET_ACCURACY  # 82%
        mean_accuracy = cv_results['mean_accuracy']
        accuracy_ci = cv_results['accuracy_ci']
        
        performance_assessment = {
            'target_accuracy': target_accuracy,
            'achieved_accuracy': mean_accuracy,
            'accuracy_ci_lower': accuracy_ci[0],
            'accuracy_ci_upper': accuracy_ci[1],
            'meets_target': mean_accuracy >= target_accuracy,
            'target_in_ci': accuracy_ci[0] <= target_accuracy <= accuracy_ci[1],
            'performance_gap': mean_accuracy - target_accuracy,
            'statistical_significance': self._test_significance_vs_target(
                cv_results['fold_accuracies'], target_accuracy
            )
        }
        
        if performance_assessment['meets_target']:
            print(f"✅ Target achieved: {mean_accuracy:.3f} >= {target_accuracy:.3f}")
        else:
            gap = target_accuracy - mean_accuracy
            print(f"⚠️  Target missed by {gap:.3f}: {mean_accuracy:.3f} < {target_accuracy:.3f}")
        
        return performance_assessment
    
    def _test_significance_vs_target(self, scores, target, alpha=0.05):
        """Test if achieved performance is significantly different from target"""
        t_stat, p_value = stats.ttest_1samp(scores, target)
        return {
            't_statistic': t_stat,
            'p_value': p_value,
            'is_significant': p_value < alpha,
            'better_than_target': t_stat > 0 and p_value < alpha
        }
    
    def generate_comprehensive_report(self, model, X, y):
        """
        Generate comprehensive validation report
        """
        print("📊 Generating comprehensive validation report...")
        
        report_data = {
            'validation_info': {
                'timestamp': self.timestamp,
                'dataset_size': len(X),
                'n_features': X.shape[1],
                'n_classes': len(np.unique(y)),
                'target_accuracy': config.TARGET_ACCURACY,
                'model_config': {
                    'n_d': model.n_d,
                    'n_a': model.n_a,
                    'n_steps': model.n_steps,
                    'gamma': model.gamma,
                    'lambda_sparse': model.lambda_sparse
                }
            }
        }
        
        # Run cross-validation
        cv_results = self.run_cross_validation(model, X, y)
        report_data['cross_validation'] = cv_results
        
        # Performance target evaluation
        performance_assessment = self.evaluate_performance_targets(cv_results)
        report_data['performance_assessment'] = performance_assessment
        
        # Clinical interpretability analysis
        clinical_analysis = self.analyze_clinical_interpretability(model, X, y)
        report_data['clinical_interpretability'] = clinical_analysis
        
        # Save detailed report
        report_file = self.results_dir / f"tabnet_validation_report_{self.timestamp}.json"
        with open(report_file, 'w') as f:
            json.dump(report_data, f, indent=2, default=str)
        
        # Generate summary
        self._generate_summary_report(report_data)
        
        print(f"✅ Comprehensive report saved: {report_file}")
        
        return report_data
    
    def _generate_summary_report(self, report_data):
        """Generate human-readable summary report"""
        summary_file = self.results_dir / f"tabnet_summary_{self.timestamp}.txt"
        
        with open(summary_file, 'w') as f:
            f.write("=" * 70 + "\n")
            f.write("TABNET PROSTATE CANCER VALIDATION SUMMARY\n")
            f.write("=" * 70 + "\n\n")
            
            # Basic info
            info = report_data['validation_info']
            f.write(f"Validation Date: {info['timestamp']}\n")
            f.write(f"Dataset: {info['dataset_size']:,} variants, {info['n_features']} features\n")
            f.write(f"Classes: {info['n_classes']}\n")
            f.write(f"Target Accuracy: {info['target_accuracy']:.1%}\n\n")
            
            # Model configuration
            config_info = info['model_config']
            f.write("Model Configuration:\n")
            f.write(f"  Decision Steps: {config_info['n_steps']}\n")
            f.write(f"  Attention Dimension: {config_info['n_a']}\n")
            f.write(f"  Network Dimension: {config_info['n_d']}\n\n")
            
            # Performance results
            cv = report_data['cross_validation']
            f.write("Cross-Validation Results (5-fold):\n")
            f.write(f"  Accuracy: {cv['mean_accuracy']:.3f} ± {cv['std_accuracy']:.3f}\n")
            f.write(f"  95% CI: [{cv['accuracy_ci'][0]:.3f}, {cv['accuracy_ci'][1]:.3f}]\n")
            f.write(f"  Precision: {cv['mean_precision']:.3f}\n")
            f.write(f"  Recall: {cv['mean_recall']:.3f}\n")
            f.write(f"  F1-Score: {cv['mean_f1']:.3f}\n")
            f.write(f"  ROC AUC: {cv['mean_auc']:.3f}\n\n")
            
            # Performance assessment
            perf = report_data['performance_assessment']
            f.write("Performance Assessment:\n")
            status = "✅ ACHIEVED" if perf['meets_target'] else "❌ MISSED"
            f.write(f"  Target Status: {status}\n")
            f.write(f"  Gap to Target: {perf['performance_gap']:+.3f}\n\n")
            
            # Clinical interpretability
            clinical = report_data['clinical_interpretability']
            f.write("Clinical Interpretability:\n")
            relevance = clinical['clinical_relevance_scores']['overall_clinical_relevance']
            f.write(f"  Overall Clinical Relevance: {relevance:.3f}\n")
            
            f.write("\n  Pathway Attention Analysis:\n")
            for pathway, scores in clinical['clinical_relevance_scores'].items():
                if pathway != 'overall_clinical_relevance':
                    f.write(f"    {pathway}: {scores['clinical_score']:.3f}\n")
            
            f.write("\n" + "=" * 70 + "\n")
            
        print(f"✅ Summary report saved: {summary_file}")
    
    def quick_validation(self, model_params=None, max_samples=5000):
        """
        Quick validation for development/testing
        Uses subset of data for faster iteration
        """
        print("⚡ Running quick validation...")
        
        # Use default params if not provided
        if model_params is None:
            model_params = {
                'n_d': 32, 'n_a': 32, 'n_steps': 3,
                'gamma': 1.3, 'lambda_sparse': 1e-3
            }
        
        # Load data
        model = ProstateVariantTabNet(**model_params)
        X, y = model.load_data()
        
        # Use subset for quick validation
        if len(X) > max_samples:
            indices = np.random.choice(len(X), max_samples, replace=False)
            X, y = X[indices], y[indices]
            print(f"  Using subset: {len(X):,} variants")
        
        # Quick 3-fold CV
        cv_results = self.run_cross_validation(model, X, y, cv_folds=3)
        
        # Quick performance check
        mean_accuracy = cv_results['mean_accuracy']
        target_met = mean_accuracy >= config.TARGET_ACCURACY
        
        print(f"\n⚡ Quick Validation Results:")
        print(f"   Accuracy: {mean_accuracy:.3f}")
        print(f"   Target: {'✅ MET' if target_met else '❌ MISSED'}")
        
        return cv_results

def main():
    """
    Main validation workflow
    """
    print("🧬 TabNet Prostate Cancer Validation Framework")
    print("=" * 60)
    
    # Initialize validator
    validator = TabNetValidator()
    
    try:
        # Load model and data
        print("📁 Loading model and data...")
        model = ProstateVariantTabNet(
            n_d=config.DEFAULT_TABNET_PARAMS['n_d'],
            n_a=config.DEFAULT_TABNET_PARAMS['n_a'],
            n_steps=config.DEFAULT_TABNET_PARAMS['n_steps']
        )
        
        X, y = model.load_data()
        print(f"✅ Loaded {len(X):,} variants with {X.shape[1]} features")
        
        # Generate comprehensive validation report
        validation_report = validator.generate_comprehensive_report(model, X, y)
        
        # Print key results
        cv_results = validation_report['cross_validation']
        performance = validation_report['performance_assessment']
        
        print(f"\n🎯 FINAL RESULTS:")
        print(f"   Accuracy: {cv_results['mean_accuracy']:.3f} ± {cv_results['std_accuracy']:.3f}")
        print(f"   Target: {performance['target_accuracy']:.3f}")
        print(f"   Status: {'✅ ACHIEVED' if performance['meets_target'] else '❌ MISSED'}")
        
        return validation_report
        
    except Exception as e:
        print(f"❌ Validation failed: {e}")
        import traceback
        traceback.print_exc()
        return None

if __name__ == "__main__":
    results = main()