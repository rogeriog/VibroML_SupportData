#!/bin/bash
#SBATCH --partition=debug-gpu
#SBATCH --account=htforft
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gpus=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=1:59:00
#SBATCH --job-name=omega_inst_full
#SBATCH --output=omega_inst_full_%j.out
#SBATCH --error=omega_inst_full_%j.err

# ==============================================================================
# Full run: compute I_curv + PET-MAD ensemble UQ for all 105 structures
#
# MACE-OMAT: already completed (omega_inst_full_mace.json has 105 results)
#   - unconstrained relaxation, I_curv on uniform q-mesh
#
# PET-MAD: NEW run with:
#   - symmetry-constrained relaxation (same-structure comparison)
#   - ensemble UQ (128 members): max_std, std_at_min, min_minus_3std
#   - I_curv from ensemble mean frequencies
#
# Resume logic: skips structures already in the output JSON/CSV.
# Resubmit this script until all structures are completed.
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

STRUCTURES="uq_validation_structures.json"
Q_MESH="6 6 6"

ORB_ENV="/gpfs/scratch/acad/htforft/rgouvea/conda_envs/orb_env"

PET_OUT="omega_inst_full_pet_symconstrained.json"
PET_CSV="omega_inst_full_pet_symconstrained.csv"

echo "============================================"
echo "Omega_inst full run: PET-MAD UQ + sym-constrained"
echo "Job ID: $SLURM_JOB_ID"
echo "Date: $(date)"
echo "Q-mesh: $Q_MESH"
echo "PET-MAD: symmetry-constrained relaxation + ensemble UQ"
echo "============================================"

# --- PET-MAD (symmetry-constrained + UQ) ---
echo ""
echo "========== PET-MAD (sym-constrained + UQ) =========="
echo "Started: $(date)"
"$ORB_ENV/bin/python" compute_omega_inst.py \
    --engine pet \
    --model_name pet-mad-xs \
    --structures "$STRUCTURES" \
    --output "$PET_OUT" \
    --csv_output "$PET_CSV" \
    --q_mesh $Q_MESH \
    --device cuda \
    --fmax 1e-3 \
    --symmetry_constrained \
    --phonon_uq

echo "PET-MAD finished: $(date)"

# --- Summary ---
echo ""
echo "============================================"
echo "SUMMARY"
echo "============================================"
"$ORB_ENV/bin/python" -c "
import json, os

def load(path):
    if os.path.exists(path):
        with open(path) as f:
            return {r['compound'] + '_' + str(r.get('idx',-1)): r for r in json.load(f)}
    return {}

mace = load('omega_inst_full_mace.json')
pet = load('$PET_OUT')

print(f'MACE results: {len(mace)}')
print(f'PET results: {len(pet)}')
print()

if pet:
    print(f'{\"compound\":25s} {\"cat\":15s} {\"MACE I_curv\":>12s} {\"PET I_curv\":>12s} {\"PET max_std\":>12s} {\"PET min-3s\":>12s}')
    print('-'*100)
    for key in sorted(pet.keys()):
        p = pet[key]
        m = mace.get(key, {})
        print(f'{p[\"compound\"]:25s} {p.get(\"category\",\"?\"):15s} '
              f'{m.get(\"i_curv\",0):12.4f} {p[\"i_curv\"]:12.4f} '
              f'{p.get(\"max_std\",0):12.4f} {p.get(\"min_minus_3std\",0):12.4f}')
"

echo ""
echo "All done: $(date)"
