#!/bin/bash
#SBATCH --job-name=pet_phonons_uq
#SBATCH --time=1:30:00
#SBATCH --output=reproduce_slurm_%j.log
#SBATCH --partition=debug-gpu
#SBATCH --gpus=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --account=htforft

module purge
source ~/.bashrc
module load EasyBuild/2022a CUDA/11.7.0 cuDNN/8.4.1.50-CUDA-11.7.0

conda activate /gpfs/scratch/acad/htforft/rgouvea/conda_envs/orb_env

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK
export PYTHONUNBUFFERED=1
export MP_START_METHOD=spawn

echo "=== PET-MAD phonon UQ reproduction ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node:   $SLURMD_NODENAME"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null || echo "(no GPU info)"
date

cd /gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/pet_phonons

python reproduce_pet_phonons.py

echo "=== Done ==="
date
