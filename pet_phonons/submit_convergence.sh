#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --account=htforft
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gpus=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=1:00:00
#SBATCH --job-name=conv_sweep
#SBATCH --output=convergence_sweep_%j.out
#SBATCH --error=convergence_sweep_%j.err

# ==============================================================================
# Convergence sweep: phonon stability metrics (soft_frac, I_curv, omega_min)
# vs q-mesh size (4x4x4, 6x6x6, 8x8x8) for 4 representative high-symmetry
# aristotype structures (LiF, CsPbI3, HfO2, Bi2Sn2O7) using MACE-OMAT.
#
# Outputs:
#   convergence_sweep.json  - full results
#   convergence_sweep.csv   - tabular results
#   figures_5eng/convergence_sweep.{pdf,png} - convergence plot
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

export PYTHONPATH=/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/01_codebase/VibroML:$PYTHONPATH

echo "============================================"
echo "Convergence sweep (MACE-OMAT)"
echo "Job ID: $SLURM_JOB_ID"
echo "Date:   $(date)"
echo "Node:   $SLURM_JOB_NODELIST"
echo "GPUs:   $SLURM_GPUS_ON_NODE  CPUs: $SLURM_CPUS_PER_TASK"
echo "============================================"

conda activate "$VIBROML_ENV"

"$VIBROML_ENV/bin/python" compute_convergence_sweep.py \
    --engine mace \
    --model-name medium-omat-0 \
    --device cuda \
    --output convergence_sweep.json

echo ""
echo "All done: $(date)"
