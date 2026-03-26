#!/bin/bash
#SBATCH --job-name=vibro_mp_evals
#SBATCH --output=mp_eval_batch_%j.log
#SBATCH --partition=debug-gpu 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --gpus=1
#SBATCH --time=2:00:00
#SBATCH --account=htforft

# 1. Clean Environment to prevent NumPy path pollution
module purge

source ~/.bashrc

module load EasyBuild/2022a CUDA/11.7.0 cuDNN/8.4.1.50-CUDA-11.7.0

conda activate /gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/vibroml_env
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Keep this to prevent deadlock
export MP_START_METHOD=spawn
# 3. Paths
ROOT_CIF_DIR="/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/cifs_to_eval"

FMAX=0.0005
FMAX=0.00001
MODELS=("mace")
MAX_ATOMS=1000
MAX_ATOMS=2000


echo "Starting dynamic evaluation at $(date)"

# 4. Loop through CIFs
find "$ROOT_CIF_DIR" -name "*.cif" | while read -r cif_path; do
    
    # Filter out internal result files
    # if [[ "$cif_path" == *"_relaxed"* ]] || [[ "$cif_path" == *"_unconverged"* ]] || [[ "$cif_path" == *"_primitive"* ]] || [[ "$cif_path" == *"_conventional"* ]]; then
    #     continue
    # fi

    # --- Dynamic Supercell Logic ---
    # Use Python to count atoms in primitive and find max N x N x N
    SC_DIM=$(python -c "
import ase.io
try:
    atoms = ase.io.read('$cif_path')
    n_prim = len(atoms)
    # Find the largest N such that n_prim * N^3 <= $MAX_ATOMS
    # We use cubic expansion N x N x N for simplicity in phonon sampling
    n_factor = int(($MAX_ATOMS / n_prim)**(1/3))
    n_factor = max(1, n_factor) # Ensure at least 1x1x1
    print(f'{n_factor},{n_factor},{n_factor}')
except:
    print('2,2,2') # Fallback
")
    # ------------------------------

    folder_name=$(basename "$(dirname "$cif_path")")
    file_base=$(basename "$cif_path" .cif)

    for model in "${MODELS[@]}"; do
        PREFIX="EVAL_${folder_name}_${model}"
        
        # --- NEW SKIP LOGIC START ---
        # 1. Identify all folders that match this calculation's naming pattern
        # The vibroml output format is: PREFIX_Filename_ESEN_phonon_output_Timestamp
        # So we search for: ${PREFIX}_${file_base}_*
        
        candidates=( "${PREFIX}_${file_base}_"* )
        
        already_done=false
        
        # Check if the glob found actual directories
        if [ -e "${candidates[0]}" ]; then
            for dir in "${candidates[@]}"; do
                # Check if this specific directory contains the success marker
                if [ -d "$dir" ] && ls "$dir"/*phonon_run_summary.txt 1> /dev/null 2>&1; then
                    already_done=true
                    # echo "--> SKIP: Found completed run for $file_base in $dir"
                    break
                fi
            done
        fi

        if [ "$already_done" = true ]; then
            continue
        fi
        # --- NEW SKIP LOGIC END ---

        echo "----------------------------------------------------------------"
        echo "Processing: $file_base"
        echo "Primitive Atoms: $(python -c "import ase.io; print(len(ase.io.read('$cif_path')))")"
        echo "Target Supercell: $SC_DIM"
        echo "Model: $model"
        echo "----------------------------------------------------------------"


        vibroml --cif "$cif_path" \
            --engine "$model" \
            --model_name "medium" \
            --fmax "$FMAX" \
            --supercell "$SC_DIM" \
            --output-prefix "$PREFIX" 

    done
done

echo "Finished at $(date)"