import json
import os
from pathlib import Path
import pandas as pd
import pydicom
from typing import Dict, List, Optional

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
        print(f"✅ Loaded {len(clinical_df)} clinical records.")
        return clinical_df
    except Exception as e:
        print(f"❌ Error loading clinical data: {e}")
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
        print(f"✅ Loaded {len(biospecimen_df)} biospecimen records.")
        return biospecimen_df
    except Exception as e:
        print(f"❌ Error loading biospecimen data: {e}")
        raise

def process_dicom_file(dcm_path: str) -> Optional[Dict]:
    """
    Process a single DICOM file and extract relevant metadata
    
    Parameters:
        dcm_path: Path to DICOM file
    Returns:
        Dictionary containing DICOM metadata or None if processing fails
    """
    try:
        dcm = pydicom.dcmread(dcm_path, stop_before_pixels=True)
        return {
            "submitter_id": getattr(dcm, "PatientID", None),
            "dicom_file": os.path.basename(dcm_path),
            "dicom_path": dcm_path,
            "StudyDate": getattr(dcm, "StudyDate", None),
            "Modality": getattr(dcm, "Modality", None),
            "StudyInstanceUID": getattr(dcm, "StudyInstanceUID", None),
            "SeriesInstanceUID": getattr(dcm, "SeriesInstanceUID", None),
        }
    except Exception as e:
        print(f"⚠️ Error reading {dcm_path}: {e}")
        return None

def process_image_data(dicom_root: str) -> pd.DataFrame:
    """
    Process all DICOM files in a directory
    
    Parameters:
        dicom_root: Root directory containing DICOM files
    Returns:
        DataFrame containing processed DICOM metadata
    """
    dicom_records = []
    
    for dcm_path in Path(dicom_root).rglob("*.dcm"):
        record = process_dicom_file(str(dcm_path))
        if record:
            dicom_records.append(record)
    
    dicom_df = pd.DataFrame(dicom_records)
    if not dicom_df.empty:
        dicom_df["submitter_id"] = dicom_df["submitter_id"].str.upper()
    print(f"✅ Extracted metadata from {len(dicom_df)} DICOM files.")
    return dicom_df

def merge_datasets(clinical_df: pd.DataFrame, 
                  biospecimen_df: pd.DataFrame, 
                  dicom_df: pd.DataFrame) -> pd.DataFrame:
    """
    Merge clinical, biospecimen, and DICOM data
    
    Parameters:
        clinical_df: Clinical data DataFrame
        biospecimen_df: Biospecimen data DataFrame
        dicom_df: DICOM metadata DataFrame
    Returns:
        Merged DataFrame containing all data
    """
    full_clinical_df = pd.merge(clinical_df, biospecimen_df, 
                               on="submitter_id", how="outer")
    print(f"✅ Combined clinical and biospecimen: {len(full_clinical_df)} records.")
    
    merged_df = pd.merge(full_clinical_df, dicom_df, 
                        on="submitter_id", how="inner")
    print(f"✅ Final merged dataset: {len(merged_df)} rows.")
    return merged_df

def save_merged_data(merged_df: pd.DataFrame, output_path: str) -> None:
    """
    Save the merged dataset to a CSV file
    """
    merged_df.to_csv(output_path, index=False)
    print(f"📁 Saved merged dataset to: {output_path}")

def main():
    # Define project paths
    project_root = Path(__file__).parent.parent
    data_dir = project_root / "data"
    
    # Input paths
    clinical_path = data_dir / "raw/clinical/clinical.project-tcga-prad.json"
    biospecimen_path = data_dir / "raw/clinical/biospecimen.project-tcga-prad.json"
    dicom_root = data_dir / "raw/images/manifest-LKlXeErz8879128910887416605/TCGA-PRAD"
    
    # Output paths
    output_dir = data_dir / "processed/merged_datasets"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "merged_clinical_dicom_dataset.csv"

    try:
        # Load data
        clinical_df = load_clinical_data(str(clinical_path))
        biospecimen_df = load_biospecimen_data(str(biospecimen_path))
        dicom_df = process_image_data(str(dicom_root))
        
        # Merge all data
        merged_df = merge_datasets(clinical_df, biospecimen_df, dicom_df)
        
        # Save results
        merged_df.to_csv(output_path, index=False)
        print(f"📁 Saved merged dataset to: {output_path}")
        
        # Preview results
        print("\nPreview of merged dataset:")
        print(merged_df.head())
        
    except Exception as e:
        print(f"❌ Error in main processing: {e}")
        raise

if __name__ == "__main__":
    main()
