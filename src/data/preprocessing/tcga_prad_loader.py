"""
Data loader module for TCGA-PRAD data.
"""

import json
import os
from pathlib import Path
import pandas as pd
import pydicom
from typing import Dict, List, Optional

from src.utils.logger import setup_logger

logger = setup_logger("tcga_prad_loader")

def load_clinical_data(clinical_file_path: str) -> pd.DataFrame:
    """
    Load and process clinical data from JSON file
    
    Parameters:
        clinical_file_path: Path to clinical JSON file
    Returns:
        DataFrame containing processed clinical data
    """
    try:
        with open(clinical_file_path) as f:
            clinical_data = json.load(f)
        
        clinical_df = pd.json_normalize(clinical_data)
        clinical_df["submitter_id"] = clinical_df["submitter_id"].str.upper()
        logger.info(f"Loaded {len(clinical_df)} clinical records")
        return clinical_df
    except Exception as e:
        logger.error(f"Error loading clinical data: {e}")
        raise

def load_biospecimen_data(biospecimen_file_path: str) -> pd.DataFrame:
    """
    Load and process biospecimen data from JSON file
    
    Parameters:
        biospecimen_file_path: Path to biospecimen JSON file
    Returns:
        DataFrame containing processed biospecimen data
    """
    try:
        with open(biospecimen_file_path) as f:
            biospecimen_data = json.load(f)
        
        biospecimen_df = pd.json_normalize(biospecimen_data)
        biospecimen_df["submitter_id"] = biospecimen_df["submitter_id"].str.upper()
        logger.info(f"Loaded {len(biospecimen_df)} biospecimen records")
        return biospecimen_df
    except Exception as e:
        logger.error(f"Error loading biospecimen data: {e}")
        raise

def process_dicom_file(dcm_path: str, config: Dict) -> Optional[Dict]:
    """
    Process a single DICOM file and extract relevant metadata.
    
    The function extracts key metadata from DICOM files including:
    - PatientID (maps to submitter_id in clinical data)
    - StudyDate (date of the imaging study)
    - Modality (type of imaging, e.g., CT, MRI)
    - StudyInstanceUID (unique identifier for the study)
    - SeriesInstanceUID (unique identifier for the series within study)
    
    Parameters:
        dcm_path: Path to DICOM file
        config: Configuration dictionary containing DICOM processing parameters
    Returns:
        Dictionary containing DICOM metadata or None if processing fails
    """
    try:
        # Read DICOM file without loading pixel data for efficiency
        dcm = pydicom.dcmread(dcm_path, stop_before_pixels=config["dicom"]["stop_before_pixels"])
        
        def safe_get_value(dcm, field):
            """Safely get DICOM field value, handling MultiValue fields"""
            try:
                value = getattr(dcm, field, None)
                if value is None:
                    return None
                # Convert MultiValue to string or list
                if hasattr(value, 'VM') and value.VM > 1:
                    return list(value)
                return str(value)
            except Exception:
                return None
        
        # Extract required fields
        metadata = {
            "submitter_id": safe_get_value(dcm, "PatientID"),  # Links to clinical data
            "dicom_file": os.path.basename(dcm_path),
            "dicom_path": dcm_path,
            "StudyDate": safe_get_value(dcm, "StudyDate"),
            "Modality": safe_get_value(dcm, "Modality"),
            "StudyInstanceUID": safe_get_value(dcm, "StudyInstanceUID"),
            "SeriesInstanceUID": safe_get_value(dcm, "SeriesInstanceUID"),
            # Additional useful DICOM tags
            "SliceThickness": safe_get_value(dcm, "SliceThickness"),
            "ImageType": safe_get_value(dcm, "ImageType"),
            "PixelSpacing": safe_get_value(dcm, "PixelSpacing"),
            "StudyDescription": safe_get_value(dcm, "StudyDescription")
        }
            
        return metadata
    except Exception as e:
        logger.warning(f"Error reading {dcm_path}: {e}")
        return None

def process_image_data(dicom_root: str, config: Dict) -> pd.DataFrame:
    """
    Process all DICOM files in a directory and create a structured dataset.
    
    This function:
    1. Recursively finds all DICOM files in the directory
    2. Extracts metadata from each file
    3. Creates a DataFrame with one row per DICOM file
    4. Standardizes patient IDs to match clinical data format
    
    The resulting DataFrame can be merged with clinical data using submitter_id
    
    Parameters:
        dicom_root: Root directory containing DICOM files
        config: Configuration dictionary
    Returns:
        DataFrame containing processed DICOM metadata
    """
    dicom_records = []
    total_files = 0
    processed_files = 0
    
    logger.info(f"Starting DICOM processing from directory: {dicom_root}")
    
    # Recursively process all DICOM files
    for dcm_path in Path(dicom_root).rglob("*.dcm"):
        total_files += 1
        record = process_dicom_file(str(dcm_path), config)
        if record:
            dicom_records.append(record)
            processed_files += 1
            
        if total_files % 1000 == 0:
            logger.info(f"Processed {total_files} DICOM files...")
    
    dicom_df = pd.DataFrame(dicom_records)
    
    if not dicom_df.empty:
        # Standardize patient IDs to match clinical data format
        dicom_df["submitter_id"] = dicom_df["submitter_id"].str.upper()
        
        # Add summary statistics
        logger.info(f"Total DICOM files found: {total_files}")
        logger.info(f"Successfully processed files: {processed_files}")
        logger.info(f"Number of unique patients: {dicom_df['submitter_id'].nunique()}")
        logger.info(f"Number of unique studies: {dicom_df['StudyInstanceUID'].nunique()}")
        logger.info(f"Modalities found: {dicom_df['Modality'].unique().tolist()}")
    else:
        logger.warning("No DICOM files were successfully processed")
        
    return dicom_df

def merge_datasets(clinical_df: pd.DataFrame, 
                  biospecimen_df: pd.DataFrame, 
                  dicom_df: pd.DataFrame,
                  config: Dict) -> pd.DataFrame:
    """
    Merge clinical, biospecimen, and DICOM data.
    
    The merge process:
    1. First merges clinical and biospecimen data using submitter_id
    2. Then merges with DICOM metadata using submitter_id
    
    This creates a comprehensive dataset where:
    - Each row represents a unique DICOM image
    - Clinical and biospecimen data is repeated for all images of the same patient
    - Relationships can be analyzed between clinical features and imaging data
    
    Parameters:
        clinical_df: Clinical data DataFrame
        biospecimen_df: Biospecimen data DataFrame
        dicom_df: DICOM metadata DataFrame
        config: Configuration dictionary
    Returns:
        Merged DataFrame containing all data
    """
    # Validate input data
    if len(clinical_df) < config["VALIDATION_CONFIG"]["min_clinical_records"]:
        raise ValueError(f"Insufficient clinical records: {len(clinical_df)}")
    
    if len(dicom_df) < config["VALIDATION_CONFIG"]["min_dicom_files"]:
        raise ValueError(f"Insufficient DICOM files: {len(dicom_df)}")
    
    # Merge clinical and biospecimen data
    full_clinical_df = pd.merge(clinical_df, biospecimen_df, 
                               on="submitter_id", how="outer")
    logger.info(f"Combined clinical and biospecimen: {len(full_clinical_df)} records")
    logger.info(f"Number of unique patients in clinical data: {full_clinical_df['submitter_id'].nunique()}")
    
    # Merge with DICOM data
    merged_df = pd.merge(full_clinical_df, dicom_df, 
                        on="submitter_id", how="inner")
    
    # Log merge statistics
    logger.info(f"Final merged dataset: {len(merged_df)} rows")
    logger.info(f"Number of patients with both clinical and imaging data: {merged_df['submitter_id'].nunique()}")
    logger.info(f"Average number of images per patient: {len(merged_df) / merged_df['submitter_id'].nunique():.2f}")
    
    return merged_df 