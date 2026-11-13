#!/bin/bash
#SBATCH --partition=debug-gpu
#SBATCH --account=htforft
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gpus=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=2:00:00
#SBATCH --job-name=uq_validation

# ============================================================
# PET-MAD UQ Validation Batch
# Runs 105 structures with resume logic.
# Resubmit to continue from last unfinished structure.
#
# Uses orb_env (has upet/metatrain) + PYTHONPATH for VibroML.
# ============================================================

module purge
source ~/.bashrc
module load EasyBuild/2022a CUDA/11.7.0 cuDNN/8.4.1.50-CUDA-11.7.0

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MP_START_METHOD=spawn
export PYTHONUNBUFFERED=1
export PYTHONPATH=/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/01_codebase/VibroML:$PYTHONPATH

export VIBROML_PYTHON=/gpfs/scratch/acad/htforft/rgouvea/conda_envs/orb_env/bin/python

echo "=== UQ Validation Batch ==="
echo "Date: $(date)"
echo "Node: $SLURM_NODELIST"
echo "Job ID: $SLURM_JOB_ID"
echo "CPUs: $SLURM_CPUS_PER_TASK, GPUs: $SLURM_GPUS_ON_NODE"
echo "Python: $VIBROML_PYTHON"
echo ""

cd /gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/02_runs_and_data/pet_phonons

$VIBROML_PYTHON run_uq_validation_batch.py

echo ""
echo "=== Batch complete ==="
echo "Date: $(date)"
echo "Results: $(wc -l < uq_validation_summary.csv 2>/dev/null || echo 0) lines in CSV"
