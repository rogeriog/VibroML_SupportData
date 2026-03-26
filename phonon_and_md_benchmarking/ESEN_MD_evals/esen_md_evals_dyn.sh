#!/bin/bash

# ================= CONFIGURATION =================
# Path to your CIFs
ROOT_CIF_DIR="/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/cifs_to_eval"
# Path to your conda environment (ESEN environment)
ENV_PATH="/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/esen_env"

ENGINE="esen"
TEMP=300
TIME_PS=10
# The limit to prevent OOM/timeouts (same as UMA)
MAX_ATOMS=400
# =================================================

# 1. Activate environment to access Python for atom counting
source ~/.bashrc
module load EasyBuild/2022a CUDA/11.7.0 cuDNN/8.4.1.50-CUDA-11.7.0
conda activate "$ENV_PATH"

echo "Starting MD Job Dispatcher for ESEN..."
echo "Max atoms allowed per simulation: $MAX_ATOMS"

# Create a logs directory to keep things clean
mkdir -p batch_logs_esen

find "$ROOT_CIF_DIR" -name "*.cif" | sort | while read -r cif_path; do

    folder_name=$(basename "$(dirname "$cif_path")")
    file_base=$(basename "$cif_path" .cif)
    PREFIX="MD_${folder_name}_${ENGINE}"

    # --- SKIP LOGIC: Check if already done ---
    if ls "${PREFIX}_${file_base}_"*/md_stability_analysis/*_md_stability_report.txt 1> /dev/null 2>&1; then
        # echo "--> SKIP: Found completed report for $file_base"
        continue
    fi

    # --- DYNAMIC SUPERCELL CALCULATION ---
    # Use python to get the number of atoms in the primitive cell
    PRIM_ATOMS=$(python -c "from ase.io import read; print(len(read('$cif_path')))" 2>/dev/null)
    
    # Priority list of supercells to try
    SC_SIZE="1x1x1" # Default fallback
    
    ATOMS_222=$(( PRIM_ATOMS * 8 ))
    ATOMS_221=$(( PRIM_ATOMS * 4 ))
    ATOMS_211=$(( PRIM_ATOMS * 2 ))
    ATOMS_111=$(( PRIM_ATOMS * 1 ))

    if [ "$ATOMS_222" -le "$MAX_ATOMS" ]; then
        SC_SIZE="2x2x2"
        echo "  [NORMAL] $file_base: 2x2x2 ($ATOMS_222 atoms) fits limit."
    elif [ "$ATOMS_221" -le "$MAX_ATOMS" ]; then
        SC_SIZE="2x2x1"
        echo "  [ADJUST] $file_base: 2x2x2 too large ($ATOMS_222). Using 2x2x1 ($ATOMS_221 atoms)."
    elif [ "$ATOMS_211" -le "$MAX_ATOMS" ]; then
        SC_SIZE="2x1x1"
        echo "  [ADJUST] $file_base: 2x2x1 too large ($ATOMS_221). Using 2x1x1 ($ATOMS_211 atoms)."
    else
        SC_SIZE="1x1x1"
        echo "  [ADJUST] $file_base: 2x1x1 too large ($ATOMS_211). Using 1x1x1 ($ATOMS_111 atoms)."
    fi

    if [ "$ATOMS_111" -gt "$MAX_ATOMS" ]; then
        echo "  [WARNING] $file_base: Even 1x1x1 ($ATOMS_111 atoms) exceeds limit. Submitting anyway."
    fi

    # --- JOB SUBMISSION ---
    JOB_SCRIPT="temp_submit_esen_${file_base}.slurm"
    
    cat <<EOF > "$JOB_SCRIPT"
#!/bin/bash
#SBATCH --job-name=MD_ESEN_${file_base}
#SBATCH --output=batch_logs_esen/md_esen_${file_base}_%j.log
#SBATCH --partition=gpu 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --gpus=1
#SBATCH --time=10:00:00  
#SBATCH --account=htforft

# Load Modules & Env
module purge
source ~/.bashrc
module load EasyBuild/2022a CUDA/11.7.0 cuDNN/8.4.1.50-CUDA-11.7.0
conda activate $ENV_PATH

# Performance Settings
export OMP_NUM_THREADS=\$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=\$SLURM_CPUS_PER_TASK
export MP_START_METHOD=spawn
export PYTHONUNBUFFERED=1

# --- CHANGED: FIX FOR CUDA OOM ---
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

echo "Processing: $file_base"
echo "Atoms in primitive: $PRIM_ATOMS"
echo "Selected Supercell: $SC_SIZE"

# Run VibroML
vibroml --cif "$cif_path" \\
    --method md_stability \\
    --engine "$ENGINE" \\
    --temp "$TEMP" \\
    --time "$TIME_PS" \\
    --supercell-size "$SC_SIZE" \\
    --output-prefix "$PREFIX"

EOF

    # Submit the job and delete the temp file
    sbatch "$JOB_SCRIPT"
    rm "$JOB_SCRIPT"
    
    sleep 1

done

echo "All ESEN jobs submitted."