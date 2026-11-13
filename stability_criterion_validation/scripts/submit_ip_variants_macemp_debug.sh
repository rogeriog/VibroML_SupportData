#!/bin/bash
# PORTABILITY: Update WORKDIR, CONDA_ENV, and VIBROML_CODEBASE paths
# below to match your cluster layout before submitting.
# Export VIBROML_CODEBASE or edit the script to point to your VibroML checkout.
#SBATCH --partition=debug-gpu
#SBATCH --account=htforft
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gpus=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=0:30:00
#SBATCH --job-name=ip_test_macemp
#SBATCH --output=ip_test_macemp_%j.out
#SBATCH --error=ip_test_macemp_%j.err

module purge
source ~/.bashrc
module load EasyBuild/2022a CUDA/11.7.0 cuDNN/8.4.1.50-CUDA-11.7.0

conda activate /gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/04_envs_and_archives/vibroml_env

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MP_START_METHOD=spawn
export PYTHONUNBUFFERED=1

cd /gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/02_runs_and_data/pet_phonons

# Test with 2 structures only
python compute_ip_variants.py \
    --engine macemp \
    --model-name medium \
    --structures uq_validation_structures.json \
    --output ip_test_macemp.json \
    --q_mesh 6 6 6 \
    --device cuda \
    --max_structures 2
