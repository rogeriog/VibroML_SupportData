#!/bin/bash
#SBATCH --partition=debug-gpu
#SBATCH --account=htforft
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gpus=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=1:30:00
#SBATCH --job-name=omega_inst_symtest
#SBATCH --output=omega_inst_symtest_%j.out
#SBATCH --error=omega_inst_symtest_%j.err

# ==============================================================================
# Test: symmetry-constrained PET-MAD relaxation + I_curv computation
# on 3 test structures. Compares with unconstrained PET-MAD results.
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
Q_MESH="6 6 6"

ORB_ENV="/gpfs/scratch/acad/htforft/rgouvea/conda_envs/orb_env"

echo "============================================"
echo "Omega_inst symmetry-constrained PET-MAD test"
echo "Job ID: $SLURM_JOB_ID"
echo "Date: $(date)"
echo "============================================"

# --- PET-MAD with symmetry-constrained relaxation ---
echo ""
echo "========== PET-MAD (symmetry-constrained) =========="
echo "Started: $(date)"
"$ORB_ENV/bin/python" compute_omega_inst.py \
    --engine pet \
    --model_name pet-mad-xs \
    --structures "$STRUCTURES" \
    --output test_omega_inst_pet_symconstrained.json \
    --q_mesh $Q_MESH \
    --device cuda \
    --fmax 1e-3 \
    --symmetry_constrained

echo "PET-MAD (symconstrained) finished: $(date)"

# --- Summary comparison ---
echo ""
echo "============================================"
echo "COMPARISON: unconstrained vs symmetry-constrained PET-MAD"
echo "============================================"
echo ""
"$ORB_ENV/bin/python" -c "
import json

with open('test_omega_inst_pet.json') as f:
    unc = {r['compound']: r for r in json.load(f)}
with open('test_omega_inst_pet_symconstrained.json') as f:
    sym = {r['compound']: r for r in json.load(f)}
with open('test_omega_inst_mace.json') as f:
    mace = {r['compound']: r for r in json.load(f)}

print(f'{\"compound\":25s} | {\"MACE I_curv\":>12s} | {\"PET unc\":>12s} | {\"PET sym\":>12s} | {\"MACE wmin\":>10s} | {\"PET unc wmin\":>12s} | {\"PET sym wmin\":>12s}')
print('-'*120)
for c in mace:
    m = mace[c]
    u = unc.get(c, {})
    s = sym.get(c, {})
    print(f'{c:25s} | {m[\"i_curv\"]:12.4f} | {u.get(\"i_curv\",0):12.4f} | {s.get(\"i_curv\",0):12.4f} | {m[\"omega_min\"]:10.4f} | {u.get(\"omega_min\",0):12.4f} | {s.get(\"omega_min\",0):12.4f}')
"

echo ""
echo "All done: $(date)"
