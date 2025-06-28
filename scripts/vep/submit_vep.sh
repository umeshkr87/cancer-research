#!/bin/bash
# VEP submission script - Phase 3 optimized
# Usage: bash submit_vep.sh

SCRIPT_DIR="/u/aa107/uiuc-cancer-research/scripts/vep"

# Make script executable
chmod +x ${SCRIPT_DIR}/run_vep_annotation.sh

echo "🎯 Submitting Phase 3 optimized VEP annotation job..."
echo "Target: Improve success rate from 91.1% to 95%+"

# Submit to cluster
sbatch ${SCRIPT_DIR}/run_vep_annotation.sh

echo "✅ Job submitted to IllinoisComputes partition"
echo "📊 Monitor with: squeue -u $USER"
echo "📋 Check logs in: ${SCRIPT_DIR}/vep_*.out"
echo "📁 Results will be in: /u/aa107/uiuc-cancer-research/data/processed/vep/"