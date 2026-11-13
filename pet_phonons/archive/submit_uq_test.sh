#!/bin/bash
#SBATCH --partition=debug-gpu
#SBATCH --account=htforft
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gpus=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=0:30:00
#SBATCH --job-name=vibroml_uq_test

module purge
source ~/.bashrc
module load EasyBuild/2022a CUDA/11.7.0 cuDNN/8.4.1.50-CUDA-11.7.0

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MP_START_METHOD=spawn
export PYTHONUNBUFFERED=1
export PYTHONPATH=/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/01_codebase/VibroML:$PYTHONPATH

PYTHON=/gpfs/scratch/acad/htforft/rgouvea/conda_envs/orb_env/bin/python
SCRIPT_DIR=/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/02_runs_and_data/pet_phonons

echo "=== VibroML PET-MAD UQ Test (single structure) ==="
echo "Job ID: $SLURM_JOB_ID  |  Node: $(hostname)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null || echo "No GPU visible"
date

# Run just the first structure by creating a temp manifest with only entry 1
$PYTHON -c "
import json
with open('${SCRIPT_DIR}/uq_validation_structures.json') as f:
    s = json.load(f)
with open('${SCRIPT_DIR}/uq_test_single.json', 'w') as f:
    json.dump([s[0]], f)
print(f'Test structure: {s[0][\"compound\"]} {s[0][\"classification\"]} {s[0][\"mace_freq_thz\"]:.3f} THz')
"

# Temporarily point the batch script at the single-entry manifest
$PYTHON -c "
import sys
sys.path.insert(0, '${SCRIPT_DIR}')
import run_uq_validation_batch as batch
batch.STRUCTURES_JSON = '${SCRIPT_DIR}/uq_test_single.json'
batch.OUTPUT_DIR = '${SCRIPT_DIR}/uq_validation_runs'
batch.SUMMARY_CSV = '${SCRIPT_DIR}/uq_test_summary.csv'
batch.main()
"

echo "=== Done ==="
date
