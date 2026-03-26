#!/bin/bash
#SBATCH --job-name=batch_md_mace
#SBATCH --output=md_mace_batch_%j.log
#SBATCH --partition=gpu 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --gpus=1
#SBATCH --time=20:00:00
#SBATCH --account=htforft

# 1. Clean Environment & Load Modules
module purge
source ~/.bashrc
module load EasyBuild/2022a CUDA/11.7.0 cuDNN/8.4.1.50-CUDA-11.7.0

# 2. Activate the UMA environment (from the md_stability script)
conda activate /gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/vibroml_env

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK
# Keep this to prevent deadlock (from md script)
export MP_START_METHOD=spawn
export PYTHONUNBUFFERED=1

# 3. Configuration
# !!! CHANGE THIS PATH TO YOUR TARGET CIF FOLDER !!!
ROOT_CIF_DIR="/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/cifs_to_eval"

# MD Stability Settings
ENGINE="mace"
TEMP=300
TIME_PS=10
SC_SIZE="2x2x2"

echo "Starting MD Stability batch evaluation at $(date)"

# 4. Loop through CIFs
find "$ROOT_CIF_DIR" -name "*.cif" | while read -r cif_path; do

    folder_name=$(basename "$(dirname "$cif_path")")
    file_base=$(basename "$cif_path" .cif)

    # Define a prefix for the output folder
    # Format: MD_STABILITY_{ParentFolder}_{Engine}
    PREFIX="MD_${folder_name}_${ENGINE}"
    
    # --- SKIP LOGIC ---
    # We look for folders matching: PREFIX_Filename_Engine_*
    # Inside, we check for: md_stability_analysis/*_md_stability_report.txt
    
    candidates=( "${PREFIX}_${file_base}_"* )
    already_done=false
    
    # Check if the glob found actual directories
    if [ -e "${candidates[0]}" ]; then
        for dir in "${candidates[@]}"; do
            # Check if directory exists and contains the specific report file
            if [ -d "$dir" ] && ls "$dir"/md_stability_analysis/*_md_stability_report.txt 1> /dev/null 2>&1; then
                already_done=true
                echo "--> SKIP: Found completed MD report for $file_base in $dir"
                break
            fi
        done
    fi

    if [ "$already_done" = true ]; then
        continue
    fi
    # ------------------

    echo "----------------------------------------------------------------"
    echo "Processing: $file_base"
    echo "Folder: $folder_name"
    echo "Engine: $ENGINE | Supercell: $SC_SIZE"
    echo "----------------------------------------------------------------"

    # Execute vibroml with MD Stability parameters
    vibroml --cif "$cif_path" \
        --method md_stability \
        --engine "$ENGINE" \
        --model_name "medium" \
        --temp "$TEMP" \
        --time "$TIME_PS" \
        --supercell-size "$SC_SIZE" \
        --output-prefix "$PREFIX"

done

echo "Finished at $(date)"