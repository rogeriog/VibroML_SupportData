#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --account=htforft
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gpus=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=8:00:00
#SBATCH --job-name=sc222_uma
#SBATCH --output=sc222_uma_%j.out
#SBATCH --error=sc222_uma_%j.err

module purge
source ~/.bashrc
module load EasyBuild/2022a CUDA/11.7.0 cuDNN/8.4.1.50-CUDA-11.7.0
conda activate /gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/04_envs_and_archives/uma_env

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MP_START_METHOD=spawn
export PYTHONUNBUFFERED=1
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
export VIBROML_FAIRCHEM_MODELS=/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/01_codebase/VibroML/fairchem_models

cd /gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/02_runs_and_data/pet_phonons

# Recompute the two zero-dispersion ground-truth structures with a REAL supercell.
# idx=34 (Bi2Sn2O7 Fd-3m aristotype, GT positive) and idx=39 (Cc, GT negative)
# previously used supercell 1x1x1 -> q-independent dynamical matrix -> no dispersion.
python compute_ip_variants.py \
    --engine uma  \
    --structures uq_validation_structures.json \
    --only-idx 34 39 \
    --supercell-override 2 \
    --output ip_variants_uma_sc222.json \
    --freqs-dir ip_variants_freqs_uma_sc222 \
    --q_mesh 6 6 6 \
    --device cuda
