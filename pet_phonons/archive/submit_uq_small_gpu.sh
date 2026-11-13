#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --account=htforft
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gpus=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=4:00:00
#SBATCH --job-name=uq_comprehensive
#SBATCH --output=uq_comprehensive_%j.out
#SBATCH --error=uq_comprehensive_%j.err

# ==============================================================================
# Comprehensive PET-MAD UQ validation run:
#   - Symmetry-constrained relaxation (same-structure comparison with MACE)
#   - vibroml.main --phonon-uq (band structure plots, npz, UQ summary)
#   - I_curv computation on uniform q-mesh
#   - All results saved to uq_validation_summary.csv
#
# Resume logic: skips structures already in the CSV.
# Resubmit this script until all 105 structures are completed.
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

# Skip large structures (>150 atoms) on gpu to avoid timeouts
# They are handled separately by submit_uq_large.sh on gpu partition
export VIBROML_SKIP_LARGE=1

WORKDIR="/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/02_runs_and_data/pet_phonons"
cd "$WORKDIR"

ORB_ENV="/gpfs/scratch/acad/htforft/rgouvea/conda_envs/orb_env"

# VibroML needs to be on PYTHONPATH for the subprocess call to vibroml.main
export PYTHONPATH=/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/01_codebase/VibroML:$PYTHONPATH
export VIBROML_PYTHON="$ORB_ENV/bin/python"

echo "============================================"
echo "Comprehensive UQ validation run"
echo "Job ID: $SLURM_JOB_ID"
echo "Date: $(date)"
echo "Python: $VIBROML_PYTHON"
echo "PYTHONPATH: $PYTHONPATH"
echo "============================================"

# Run the comprehensive batch script
"$ORB_ENV/bin/python" run_uq_validation_comprehensive.py

echo ""
echo "All done: $(date)"
