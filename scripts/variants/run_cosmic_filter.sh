#!/bin/bash
set -e

echo "🚀 Automated COSMIC Prostate Cancer Filter"
echo "========================================="

# Get script directory (works from anywhere)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

echo "Project root: $PROJECT_ROOT"

# Load environment
module load anaconda3
source activate tabnet-prostate

# Install packages if needed
conda install pandas numpy -y -q

# Navigate to script directory
cd "$SCRIPT_DIR"

# Run the filtering script
python filter_cosmic_prostate.py

echo "✅ Done! Check: $PROJECT_ROOT/data/raw/variants/cosmic_prostate.csv"