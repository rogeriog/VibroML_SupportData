#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --account=htforft
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gpus=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=8:00:00
#SBATCH --job-name=uq_large
#SBATCH --output=uq_large_%j.out
#SBATCH --error=uq_large_%j.err

# ==============================================================================
# UQ run for large structures (>150 atoms) on gpu partition (8hr wall time)
# Only processes structures 37, 38, 39, 44
# ==============================================================================

set -eo pipefail

module purge
source ~/.bashrc
module load EasyBuild/2022a CUDA/11.7.0 cuDNN/8.4.1.50-CUDA-11.7.0

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MP_START_METHOD=spawn
export PYTHONUNBUFFERED=1
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

# Only run large structures (skip small ones)
export VIBROML_SKIP_SMALL=1
export VIBROML_SKIP_LARGE=0

export PYTHONPATH=/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/01_codebase/VibroML:$PYTHONPATH

conda activate /gpfs/scratch/acad/htforft/rgouvea/conda_envs/orb_env

cd /gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/02_runs_and_data/pet_phonons

echo "=== UQ Comprehensive Run: LARGE structures only (gpu) ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
date

python run_uq_validation_comprehensive.py

echo "done"
date
