#!/bin/bash
#SBATCH --job-name=vibroml_neb
#SBATCH --time=00:30:00
#SBATCH --output=log_phononsneb.txt
#SBATCH --partition=debug-gpu 
#SBATCH --nodes=1                                
#SBATCH --ntasks-per-node=1                      
#SBATCH --cpus-per-task=8                        # Allocate 8 cores to match the 1/4 node ratio
#SBATCH --mem=60G                                # Reduced from 120G to meet the 60GB per GPU limit
#SBATCH --gpus=1                            
#SBATCH --account=htforft                        
#SBATCH --array=0           

module purge

source ~/.bashrc

module load EasyBuild/2022a CUDA/11.7.0 cuDNN/8.4.1.50-CUDA-11.7.0

conda activate /gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/vibroml_env
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Keep this to prevent deadlock
export MP_START_METHOD=spawn
export PYTHONUNBUFFERED=1
echo "start"
date
VIBRO_DIR="/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/VibroML"
CIF_DIR="./interp_images"
export PYTHONPATH=$PYTHONPATH:$VIBRO_DIR

# Initialize conda functions in bash
source $(conda info --base)/etc/profile.d/conda.sh

# 1. Run MACE-OMAT (Default)
echo "------------------------------------------"
echo "Running MACE-OMAT (medium-omat-0)..."
conda activate /gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/vibroml_env
python eval_energy_forces.py mace_omat "$CIF_DIR/*.cif"

# 2. Run MACE-MP (Foundational medium)
echo "------------------------------------------"
echo "Running MACE-MP (foundation medium)..."
python eval_energy_forces.py mace_mp "$CIF_DIR/*.cif"

# 3. Run eSEN
echo "------------------------------------------"
echo "Running eSEN..."
conda activate /gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/esen_env
python eval_energy_forces.py esen "$CIF_DIR/*.cif"

# 4. Run UMA
echo "------------------------------------------"
echo "Running UMA..."
conda activate /gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/uma_env
python eval_energy_forces.py uma "$CIF_DIR/*.cif"

# 5. Final Merge
echo "------------------------------------------"
echo "Merging all data..."
python - <<EOF
import pandas as pd
engines = ['mace_omat', 'mace_mp', 'esen', 'uma']
dfs = [pd.read_csv(f'results_{e}.csv') for e in engines]

final = dfs[0]
for df in dfs[1:]:
    final = final.merge(df, on='filename')

final.to_csv('final_benchmarking_results.csv', index=False)
print("COMPLETED: final_benchmarking_results.csv")
EOF