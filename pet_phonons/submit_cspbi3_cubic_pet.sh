#!/bin/bash
#SBATCH --partition=debug-gpu
#SBATCH --account=htforft
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gpus=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=1:59:00
#SBATCH --job-name=csk_cubic_pet
#SBATCH --output=csk_cubic_pet_%j.out
#SBATCH --error=csk_cubic_pet_%j.err

# ==============================================================================
# PET-MAD UQ + I_curv for the new CsPbI3 cubic Pm-3m aristotype.
# MACE I_curv was already computed by the previous job (15390049).
# This script runs the PET-MAD ensemble UQ workflow and appends the result
# to the existing uq_validation_summary.csv.
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

WORKDIR="/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/02_runs_and_data/pet_phonons"
cd "$WORKDIR"

ORB_ENV="/gpfs/scratch/acad/htforft/rgouvea/conda_envs/orb_env"

export PYTHONPATH=/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/01_codebase/VibroML:$PYTHONPATH
export VIBROML_PYTHON="$ORB_ENV/bin/python"

echo "============================================"
echo "CsPbI3 cubic Pm-3m — PET-MAD UQ + I_curv"
echo "Job ID: $SLURM_JOB_ID"
echo "Date: $(date)"
echo "============================================"

"$ORB_ENV/bin/python" run_cspbi3_cubic_uq.py

echo ""
echo "All done: $(date)"
