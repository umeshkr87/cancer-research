#!/bin/bash
# Simple VEP submission script
# Usage: bash submit_vep.sh

SCRIPT_DIR="/u/aa107/uiuc-cancer-research/scripts/vep"

# Make script executable
chmod +x ${SCRIPT_DIR}/run_vep_annotation.sh

# Submit to cluster
echo "Submitting VEP annotation job..."
sbatch ${SCRIPT_DIR}/run_vep_annotation.sh

echo "Job submitted. Monitor with: squeue -u $USER"
echo "Check logs in: ${SCRIPT_DIR}/vep_*.out"