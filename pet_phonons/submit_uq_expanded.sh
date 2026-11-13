#!/bin/bash
#SBATCH --job-name=uq_expanded
#SBATCH --time=8:00:00
#SBATCH --output=/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/02_runs_and_data/pet_phonons/logs/uq_expanded_%j.out
#SBATCH --partition=gpu
#SBATCH --gpus=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --account=htforft

module purge
source ~/.bashrc
module load EasyBuild/2022a CUDA/11.7.0 cuDNN/8.4.1.50-CUDA-11.7.0

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MP_START_METHOD=spawn
export PYTHONUNBUFFERED=1
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

export UQ_STRUCTURES_JSON=/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/02_runs_and_data/pet_phonons/uq_validation_structures_expanded.json
export VIBROML_SKIP_LARGE=1

ORB_ENV="/gpfs/scratch/acad/htforft/rgouvea/conda_envs/orb_env"
export PYTHONPATH=/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/01_codebase/VibroML:$PYTHONPATH
export VIBROML_PYTHON="$ORB_ENV/bin/python"

cd /gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/02_runs_and_data/pet_phonons

echo "=== UQ Expanded Run (40 new quaternary) ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Python: $VIBROML_PYTHON"
echo "PYTHONPATH: $PYTHONPATH"
date

conda activate $ORB_ENV

python run_uq_validation_comprehensive.py

echo "done"
date
