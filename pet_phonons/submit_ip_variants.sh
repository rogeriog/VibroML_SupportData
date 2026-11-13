#!/bin/bash
#SBATCH --partition=debug-gpu
#SBATCH --account=htforft
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gpus=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=1:59:00
#SBATCH --job-name=ip_variants
#SBATCH --output=ip_variants_%j.out
#SBATCH --error=ip_variants_%j.err

# ==============================================================================
# Compute I_p variants for all benchmark structures using BOTH MACE and PET-MAD
# on a uniform (6,6,6) q-mesh. Saves raw frequencies and I_p for p = 0.5..inf.
#
# Step 1: MACE-OMAT (vibroml_env)
# Step 2: PET-MAD (orb_env)
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

VIBROML_ENV="/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/04_envs_and_archives/vibroml_env"
ORB_ENV="/gpfs/scratch/acad/htforft/rgouvea/conda_envs/orb_env"

export PYTHONPATH=/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/01_codebase/VibroML:$PYTHONPATH

echo "============================================"
echo "I_p variants computation (MACE + PET)"
echo "Job ID: $SLURM_JOB_ID"
echo "Date: $(date)"
echo "============================================"

# --- Step 1: MACE ---
echo ""
echo "=== Step 1: MACE-OMAT ==="
conda activate "$VIBROML_ENV"

"$VIBROML_ENV/bin/python" compute_ip_variants.py \
    --engine mace \
    --structures uq_validation_structures.json \
    --output ip_variants_mace.json \
    --q_mesh 6 6 6 \
    --delta 0.05

conda deactivate

# --- Step 2: PET-MAD ---
echo ""
echo "=== Step 2: PET-MAD ==="
export VIBROML_PYTHON="$ORB_ENV/bin/python"

"$ORB_ENV/bin/python" compute_ip_variants.py \
    --engine pet \
    --structures uq_validation_structures.json \
    --output ip_variants_pet.json \
    --q_mesh 6 6 6 \
    --delta 0.05

echo ""
echo "All done: $(date)"
