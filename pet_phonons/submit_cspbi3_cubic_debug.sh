#!/bin/bash
#SBATCH --partition=debug-gpu
#SBATCH --account=htforft
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gpus=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=1:59:00
#SBATCH --job-name=csk_cubic_uq
#SBATCH --output=csk_cubic_uq_%j.out
#SBATCH --error=csk_cubic_uq_%j.err

# ==============================================================================
# Compute all metrics for the missing CsPbI3 cubic Pm-3m aristotype:
#   1. MACE-OMAT I_curv on uniform q-mesh (vibroml_env)
#   2. PET-MAD UQ + I_curv via comprehensive script (orb_env)
#
# The comprehensive script has resume logic and will skip all completed
# structures, processing only the new CsPbI3 cubic (idx=106).
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
echo "CsPbI3 cubic Pm-3m — full metric computation"
echo "Job ID: $SLURM_JOB_ID"
echo "Date: $(date)"
echo "============================================"

# --- Step 1: MACE-OMAT I_curv ---
echo ""
echo "=== Step 1: MACE-OMAT I_curv ==="
conda activate "$VIBROML_ENV"

MACE_OUTPUT="omega_inst_cspbi3_cubic_mace.json"
"$VIBROML_ENV/bin/python" compute_omega_inst.py \
    --engine mace \
    --model_name medium-omat-0 \
    --structures cspbi3_cubic_only.json \
    --output "$MACE_OUTPUT" \
    --q_mesh 8 8 8 \
    --delta 0.05 \
    --fmax 1e-3 \
    --symmetry_constrained

echo "MACE I_curv done. Results:"
cat "$MACE_OUTPUT"

conda deactivate

# --- Step 2: PET-MAD UQ + I_curv (comprehensive script) ---
echo ""
echo "=== Step 2: PET-MAD UQ + I_curv (comprehensive) ==="
export VIBROML_PYTHON="$ORB_ENV/bin/python"
export VIBROML_SKIP_LARGE=0

# The comprehensive script will resume and only process the new structure
"$ORB_ENV/bin/python" run_uq_validation_comprehensive.py

echo ""
echo "All done: $(date)"
