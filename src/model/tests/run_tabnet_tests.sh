#!/bin/bash
echo "Starting TabNet environment validation..."

# Activate environment
module load anaconda3
conda activate tabnet-prostate

# Run local tests
echo "Running local tests..."
python src/model/tests/test_environment.py

# Submit GPU tests
echo "Submitting GPU tests..."
cd src/model/tests
JOB_ID=$(sbatch test_gpu_environment.sbatch | cut -d ' ' -f 4)
echo "GPU test job submitted: $JOB_ID"

# Wait for completion and show results
echo "Waiting for GPU test completion..."
while [ $(squeue -j $JOB_ID -h | wc -l) -gt 0 ]; do
    sleep 30
done

echo "GPU test completed. Results:"
cat env_test_${JOB_ID}.out

# Clean up output files
rm -f env_test_${JOB_ID}.out env_test_${JOB_ID}.err
echo "Output files cleaned up."