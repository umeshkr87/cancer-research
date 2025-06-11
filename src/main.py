"""
Main script to run the data processing pipeline.
"""

import sys
from pathlib import Path
from typing import Dict, Any

# Add src to Python path
sys.path.append(str(Path(__file__).parent.parent))

from src.data.preprocessing.tcga_prad_preprocessor import TCGAPRADPreprocessor
from src.utils.logger import setup_logger
from config.pipeline_config import (
    RAW_DATA,
    PROCESSED_DATA,
    PREPROCESSING_CONFIG,
    VALIDATION_CONFIG
)

def process_tcga_prad(config: Dict[str, Any]) -> None:
    """
    Process TCGA-PRAD data.
    
    Args:
        config: Configuration dictionary
    """
    logger = setup_logger("tcga_prad_pipeline")
    logger.info("Starting TCGA-PRAD data processing")
    
    try:
        # Initialize preprocessor
        preprocessor = TCGAPRADPreprocessor(config)
        
        # Create output directory if it doesn't exist
        output_dir = PROCESSED_DATA / "tcga_prad"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Process data
        output_path = output_dir / "merged_clinical_dicom_dataset.csv"
        
        # Note: We pass None as input_path because the preprocessor will use
        # the paths from the config dictionary
        preprocessor.process(None, output_path)
        
        logger.info(f"Successfully processed TCGA-PRAD data. Output saved to {output_path}")
        
    except Exception as e:
        logger.error(f"Error processing TCGA-PRAD data: {str(e)}")
        raise

def main():
    """Main pipeline execution function."""
    logger = setup_logger("main_pipeline")
    logger.info("Starting data processing pipeline")
    
    # Create necessary directories
    PROCESSED_DATA.mkdir(parents=True, exist_ok=True)
    
    # Combine all configuration
    config = {
        "RAW_DATA": RAW_DATA,
        "PREPROCESSING_CONFIG": PREPROCESSING_CONFIG,
        "VALIDATION_CONFIG": VALIDATION_CONFIG
    }
    
    try:
        # Process TCGA-PRAD data
        process_tcga_prad(config)
        
        # Add processing for COSMIC and ClinVar data here
        # TODO: Implement COSMIC and ClinVar processing
        
        logger.info("Pipeline completed successfully")
        
    except Exception as e:
        logger.error(f"Pipeline failed: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main() 