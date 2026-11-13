#!/bin/bash
#SBATCH --partition=batch
#SBATCH --account=htforft
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=24:00:00
#SBATCH --job-name=sc222_esen39_cpu
#SBATCH --output=sc222_esen39_cpu_%j.out
#SBATCH --error=sc222_esen39_cpu_%j.err

module purge
source ~/.bashrc
module load EasyBuild/2022a CUDA/11.7.0 cuDNN/8.4.1.50-CUDA-11.7.0
conda activate /gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/04_envs_and_archives/esen_env

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MP_START_METHOD=spawn
export PYTHONUNBUFFERED=1
export VIBROML_FAIRCHEM_MODELS=/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/01_codebase/VibroML/fairchem_models

# Hide GPU to force eSEN to use CPU (1408 atoms OOMs even on 80GB A100)
export CUDA_VISIBLE_DEVICES=""

cd /gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/02_runs_and_data/pet_phonons

# eSEN idx=39 only, on CPU (1408 atoms OOMs even on 80GB GPU)
python compute_ip_variants.py \
    --engine esen \
    --structures uq_validation_structures.json \
    --only-idx 39 \
    --supercell-override 2 \
    --output ip_variants_esen_sc222_idx39.json \
    --freqs-dir ip_variants_freqs_esen_sc222_idx39 \
    --q_mesh 6 6 6 \
    --device cpu
