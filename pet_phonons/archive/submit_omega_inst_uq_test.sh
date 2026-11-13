#!/bin/bash
#SBATCH --partition=debug-gpu
#SBATCH --account=htforft
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gpus=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=1:59:00
#SBATCH --job-name=omega_uq_test
#SBATCH --output=omega_uq_test_%j.out
#SBATCH --error=omega_uq_test_%j.err

set -eo pipefail

module purge
source ~/.bashrc
module load EasyBuild/2022a CUDA/11.7.0 cuDNN/8.4.1.50-CUDA-11.7.0

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MP_START_METHOD=spawn
export PYTHONUNBUFFERED=1

WORKDIR="/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/02_runs_and_data/pet_phonons"
cd "$WORKDIR"

ORB_ENV="/gpfs/scratch/acad/htforft/rgouvea/conda_envs/orb_env"

echo "============================================"
echo "PET-MAD UQ test (1 structure: LiF simple cubic)"
echo "Job ID: $SLURM_JOB_ID"
echo "Date: $(date)"
echo "============================================"

"$ORB_ENV/bin/python" compute_omega_inst.py \
    --engine pet \
    --model_name pet-mad-xs \
    --structures test_uq_single.json \
    --output test_uq_single.json.json \
    --csv_output test_uq_single.csv \
    --q_mesh 6 6 6 \
    --device cuda \
    --fmax 1e-3 \
    --symmetry_constrained \
    --phonon_uq

echo ""
echo "Done: $(date)"
echo ""
echo "=== CSV contents ==="
cat test_uq_single.csv
