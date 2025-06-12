#!/bin/bash
#SBATCH --job-name=tabnet_env_test
#SBATCH --partition=IllinoisComputes-GPU
#SBATCH --account=aa107-ic
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --time=00:20:00
#SBATCH --output=env_test_%j.out
#SBATCH --error=env_test_%j.err

echo "Starting TabNet Environment Test on GPU Node"
echo "============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Date: $(date)"
echo "GPU: $CUDA_VISIBLE_DEVICES"

# Navigate to project root
cd $SLURM_SUBMIT_DIR/../../../
echo "Working directory: $(pwd)"

# Load environment with explicit conda initialization
echo "Loading anaconda module..."
module load anaconda3

# Initialize conda for bash (required for compute nodes)
echo "Initializing conda..."
eval "$(conda shell.bash hook)"

echo "Activating tabnet-prostate environment..."
conda activate tabnet-prostate

# Verify environment activation
echo "Python location: $(which python)"
echo "Python version: $(python --version)"
echo "Conda environment: $CONDA_DEFAULT_ENV"

# Additional verification
python -c "import sys; print(f'Python executable: {sys.executable}')"

# Run comprehensive test
echo ""
echo "Running comprehensive environment test..."
echo "========================================="
python src/model/tests/test_environment.py

echo ""
echo "Environment test completed at $(date)"
echo "Check output above for results"