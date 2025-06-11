"""
TCGA-PRAD specific data preprocessor.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Union, Dict, Any

from .base_preprocessor import BasePreprocessor
from .tcga_prad_loader import (
    load_clinical_data,
    load_biospecimen_data,
    process_image_data,
    merge_datasets
)
from src.utils.logger import setup_logger

logger = setup_logger("tcga_prad_preprocessor")

class TCGAPRADPreprocessor(BasePreprocessor):
    def __init__(self, config: Dict[str, Any]):
        super().__init__(config)
        self.logger = logger
    
    def load_data(self, path: Union[str, Path]) -> pd.DataFrame:
        """
        Load TCGA-PRAD data from the specified paths.
        
        Args:
            path: Dictionary containing paths to different data sources
        
        Returns:
            pd.DataFrame: Merged dataset containing clinical, biospecimen, and DICOM data
        """
        self.logger.info("Starting TCGA-PRAD data loading")
        
        # Load clinical data
        clinical_df = load_clinical_data(str(self.config["RAW_DATA"]["TCGA_PRAD"]["clinical"]))
        
        # Load biospecimen data
        biospecimen_df = load_biospecimen_data(str(self.config["RAW_DATA"]["TCGA_PRAD"]["biospecimen"]))
        
        # Process DICOM images
        dicom_df = process_image_data(
            str(self.config["RAW_DATA"]["TCGA_PRAD"]["dicom_root"]),
            self.config["PREPROCESSING_CONFIG"]
        )
        
        # Merge all datasets
        merged_df = merge_datasets(
            clinical_df,
            biospecimen_df,
            dicom_df,
            self.config
        )
        
        return merged_df
    
    def clean_data(self, data: pd.DataFrame) -> pd.DataFrame:
        """
        Clean TCGA-PRAD specific data.
        
        Args:
            data: Input DataFrame
            
        Returns:
            pd.DataFrame: Cleaned data
        """
        self.logger.info("Starting data cleaning")
        
        # First, handle any columns that might contain lists
        for col in data.columns:
            if data[col].apply(lambda x: isinstance(x, list)).any():
                self.logger.info(f"Converting list column {col} to string representation")
                data[col] = data[col].apply(lambda x: str(x) if isinstance(x, list) else x)
        
        # Remove duplicates
        data = data.drop_duplicates()
        self.logger.info(f"Removed duplicates. Remaining rows: {len(data)}")
        
        # Identify column types after list conversion
        numeric_cols = []
        categorical_cols = []
        
        for col in data.columns:
            # Skip columns that should not be processed
            if col in ['submitter_id', 'dicom_file', 'dicom_path']:
                continue
                
            # Try to convert to numeric
            try:
                pd.to_numeric(data[col], errors='raise')
                numeric_cols.append(col)
            except (ValueError, TypeError):
                categorical_cols.append(col)
        
        self.logger.info(f"Identified {len(numeric_cols)} numeric columns and {len(categorical_cols)} categorical columns")
        
        # Handle missing values in numeric columns
        for col in numeric_cols:
            median_val = pd.to_numeric(data[col], errors='coerce').median()
            if pd.isna(median_val):
                self.logger.warning(f"Column {col} has all missing/invalid numeric values")
                data[col] = 0  # or any other appropriate default value
            else:
                data[col] = pd.to_numeric(data[col], errors='coerce').fillna(median_val)
        
        # Handle missing values in categorical columns
        for col in categorical_cols:
            mode_val = data[col].mode().iloc[0] if not data[col].mode().empty else "Unknown"
            data[col] = data[col].fillna(mode_val)
        
        self.logger.info("Completed missing value imputation")
        
        return data
    
    def transform_data(self, data: pd.DataFrame) -> pd.DataFrame:
        """
        Transform TCGA-PRAD data.
        
        Args:
            data: Input DataFrame
            
        Returns:
            pd.DataFrame: Transformed data
        """
        self.logger.info("Starting data transformation")
        
        # Only transform numeric columns that we identified as truly numeric
        numeric_cols = []
        for col in data.columns:
            try:
                if col not in ['submitter_id', 'dicom_file', 'dicom_path']:
                    pd.to_numeric(data[col], errors='raise')
                    numeric_cols.append(col)
            except (ValueError, TypeError):
                continue
        
        # Log transform numeric columns if needed
        for col in numeric_cols:
            try:
                min_val = pd.to_numeric(data[col], errors='coerce').min()
                if pd.notna(min_val) and min_val > 0:
                    data[col] = np.log1p(pd.to_numeric(data[col], errors='coerce'))
                    self.logger.debug(f"Applied log transformation to column: {col}")
            except Exception as e:
                self.logger.warning(f"Could not transform column {col}: {str(e)}")
        
        self.logger.info("Completed data transformation")
        return data
    
    def validate_data(self, data: pd.DataFrame) -> bool:
        """
        Validate TCGA-PRAD data.
        
        Args:
            data: Input DataFrame
            
        Returns:
            bool: True if validation passes, False otherwise
        """
        self.logger.info("Starting data validation")
        
        # Check for missing values
        missing_ratio = data.isnull().sum().max() / len(data)
        if missing_ratio > self.config["VALIDATION_CONFIG"]["missing_threshold"]:
            self.logger.error(f"Missing value ratio {missing_ratio:.2f} exceeds threshold")
            return False
        
        # Check for high correlation in numeric features
        numeric_cols = []
        for col in data.columns:
            try:
                if col not in ['submitter_id', 'dicom_file', 'dicom_path']:
                    pd.to_numeric(data[col], errors='raise')
                    numeric_cols.append(col)
            except (ValueError, TypeError):
                continue
        
        if len(numeric_cols) > 1:
            try:
                numeric_data = data[numeric_cols].apply(pd.to_numeric, errors='coerce')
                corr_matrix = numeric_data.corr().abs()
                
                # Count highly correlated pairs (excluding self-correlations)
                high_corr_pairs = 0
                for i in range(len(numeric_cols)):
                    for j in range(i + 1, len(numeric_cols)):
                        if corr_matrix.iloc[i, j] > self.config["VALIDATION_CONFIG"]["correlation_threshold"]:
                            high_corr_pairs += 1
                            self.logger.warning(f"High correlation ({corr_matrix.iloc[i, j]:.3f}) between {numeric_cols[i]} and {numeric_cols[j]}")
                
                if high_corr_pairs > self.config["VALIDATION_CONFIG"]["max_correlated_pairs"]:
                    self.logger.error(f"Found {high_corr_pairs} highly correlated feature pairs (threshold: {self.config['VALIDATION_CONFIG']['max_correlated_pairs']})")
                    return False
                else:
                    self.logger.info(f"Found {high_corr_pairs} highly correlated feature pairs (within acceptable limit)")
                
            except Exception as e:
                self.logger.warning(f"Could not compute correlations: {str(e)}")
        
        self.logger.info("Data validation passed")
        return True 