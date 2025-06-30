#!/bin/bash
#SBATCH --job-name=tabnet_validate
#SBATCH --partition=IllinoisComputes-GPU
#SBATCH --account=aa107-ic
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --time=02:00:00
#SBATCH --output=validation_%j.out

cd /u/aa107/uiuc-cancer-research
module load anaconda3
eval "$(conda shell.bash hook)"
conda activate tabnet-prostate

python src/model/validation/validate_tabnet.py
