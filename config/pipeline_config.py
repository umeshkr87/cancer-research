"""
Configuration settings for the data processing pipeline.
"""

from pathlib import Path

# Base paths
BASE_DIR = Path(__file__).parent.parent
DATA_DIR = BASE_DIR / "data"

# Raw data paths
RAW_DATA = {
    "TCGA_PRAD": {
        "clinical": DATA_DIR / "raw" / "TCGA-PRAD" / "clinical" / "clinical.project-tcga-prad.json",
        "biospecimen": DATA_DIR / "raw" / "TCGA-PRAD" / "clinical" / "biospecimen.project-tcga-prad.json",
        "dicom_root": DATA_DIR / "raw" / "TCGA-PRAD" / "images" / "manifest-LKlXeErz8879128910887416605"
    },
    "COSMIC": DATA_DIR / "raw" / "COSMIC",
    "CLINVAR": DATA_DIR / "raw" / "clinvar"
}

# Processed data paths
PROCESSED_DATA = DATA_DIR / "processed"
INTERIM_DATA = DATA_DIR / "interim"

# Data processing parameters
PREPROCESSING_CONFIG = {
    "random_seed": 42,
    "test_size": 0.2,
    "validation_size": 0.2,
    "dicom": {
        "stop_before_pixels": True,  # Only read DICOM metadata, not pixel data
        "required_fields": ["PatientID", "StudyDate", "Modality", "StudyInstanceUID", "SeriesInstanceUID"]
    }
}

# Feature engineering parameters
FEATURE_CONFIG = {
    "numeric_features": [],  # Will be populated during processing
    "categorical_features": [],  # Will be populated during processing
    "text_features": []  # Will be populated during processing
}

# Validation thresholds
VALIDATION_CONFIG = {
    "missing_threshold": 0.8,  # Increased threshold as clinical data might have more missing values
    "correlation_threshold": 0.98,  # Increased threshold to allow for naturally correlated medical features
    "min_clinical_records": 1,  # Minimum number of clinical records required
    "min_dicom_files": 1,  # Minimum number of DICOM files required
    "max_correlated_pairs": 10  # Maximum number of allowed highly correlated pairs
} 