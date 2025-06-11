"""
Base preprocessor class for data preprocessing.
"""

from abc import ABC, abstractmethod
import pandas as pd
from pathlib import Path
from typing import Union, Dict, Any

class BasePreprocessor(ABC):
    def __init__(self, config: Dict[str, Any]):
        self.config = config
        
    @abstractmethod
    def load_data(self, path: Union[str, Path]) -> pd.DataFrame:
        """Load data from the specified path."""
        pass
    
    @abstractmethod
    def clean_data(self, data: pd.DataFrame) -> pd.DataFrame:
        """Clean the loaded data."""
        pass
    
    @abstractmethod
    def transform_data(self, data: pd.DataFrame) -> pd.DataFrame:
        """Transform the cleaned data."""
        pass
    
    @abstractmethod
    def validate_data(self, data: pd.DataFrame) -> bool:
        """Validate the transformed data."""
        pass
    
    def process(self, input_path: Union[str, Path], output_path: Union[str, Path]) -> None:
        """
        Execute the complete preprocessing pipeline.
        
        Args:
            input_path: Path to input data
            output_path: Path to save processed data
        """
        # Load data
        data = self.load_data(input_path)
        
        # Clean data
        data = self.clean_data(data)
        
        # Transform data
        data = self.transform_data(data)
        
        # Validate data
        if not self.validate_data(data):
            raise ValueError("Data validation failed")
        
        # Save processed data
        data.to_csv(output_path, index=False) 