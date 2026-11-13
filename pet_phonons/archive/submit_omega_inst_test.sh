#!/bin/bash
#SBATCH --partition=debug-gpu
#SBATCH --account=htforft
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gpus=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=1:30:00
#SBATCH --job-name=omega_inst_test
#SBATCH --output=omega_inst_test_%j.out
#SBATCH --error=omega_inst_test_%j.err

# ==============================================================================
# Test run: compute I_curv (mean negative curvature) for 3 test structures
# using both MACE-OMAT and PET-MAD on a uniform q-point mesh.
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

WORKDIR="/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/02_runs_and_data/pet_phonons"
cd "$WORKDIR"

STRUCTURES="test_omega_inst_structures.json"
Q_MESH="6 6 6"   # uniform q-mesh; 6x6x6 = 216 q-points (no Gamma, even mesh)

# Environment paths
VIBROML_ENV="/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/04_envs_and_archives/vibroml_env"
ORB_ENV="/gpfs/scratch/acad/htforft/rgouvea/conda_envs/orb_env"

echo "============================================"
echo "Omega_inst test run"
echo "Job ID: $SLURM_JOB_ID"
echo "Date: $(date)"
echo "Structures: $STRUCTURES"
echo "Q-mesh: $Q_MESH"
echo "============================================"

# --- MACE-OMAT ---
echo ""
echo "========== MACE-OMAT =========="
echo "Started: $(date)"
"$VIBROML_ENV/bin/python" compute_omega_inst.py \
    --engine mace \
    --model_name medium-omat-0 \
    --structures "$STRUCTURES" \
    --output test_omega_inst_mace.json \
    --q_mesh $Q_MESH \
    --device cuda \
    --fmax 1e-3

echo "MACE finished: $(date)"

# --- PET-MAD ---
echo ""
echo "========== PET-MAD =========="
echo "Started: $(date)"
"$ORB_ENV/bin/python" compute_omega_inst.py \
    --engine pet \
    --model_name pet-mad-xs \
    --structures "$STRUCTURES" \
    --output test_omega_inst_pet.json \
    --q_mesh $Q_MESH \
    --device cuda \
    --fmax 1e-3

echo "PET-MAD finished: $(date)"

# --- Summary ---
echo ""
echo "============================================"
echo "SUMMARY"
echo "============================================"
echo ""
echo "--- MACE-OMAT ---"
cat test_omega_inst_mace.json | "$ORB_ENV/bin/python" -c "
import json, sys
data = json.load(sys.stdin)
print(f'{\"compound\":25s} {\"I_curv(THz^2)\":>14s} {\"omega_min(THz)\":>15s} {\"soft_frac\":>10s} {\"n_imag\":>8s}')
for r in data:
    print(f'{r[\"compound\"]:25s} {r[\"i_curv\"]:14.4f} {r[\"omega_min\"]:15.4f} {r[\"soft_frac\"]:10.4f} {r[\"n_imag\"]:8d}')
"
echo ""
echo "--- PET-MAD ---"
cat test_omega_inst_pet.json | "$ORB_ENV/bin/python" -c "
import json, sys
data = json.load(sys.stdin)
print(f'{\"compound\":25s} {\"I_curv(THz^2)\":>14s} {\"omega_min(THz)\":>15s} {\"soft_frac\":>10s} {\"n_imag\":>8s}')
for r in data:
    print(f'{r[\"compound\"]:25s} {r[\"i_curv\"]:14.4f} {r[\"omega_min\"]:15.4f} {r[\"soft_frac\"]:10.4f} {r[\"n_imag\"]:8d}')
"

echo ""
echo "All done: $(date)"
