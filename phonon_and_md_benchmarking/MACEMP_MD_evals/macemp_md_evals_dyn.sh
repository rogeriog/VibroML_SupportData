#!/bin/bash

# ================= CONFIGURATION =================
# Path to your CIFs
ROOT_CIF_DIR="/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/cifs_to_eval"
# Path to your conda environment
ENV_PATH="/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/vibroml_env"

ENGINE="mace"
TEMP=300
TIME_PS=10
# The limit you requested to prevent CUDA OOM
MAX_ATOMS=800
# =================================================

# 1. Activate environment to access Python for atom counting
source ~/.bashrc
module load EasyBuild/2022a CUDA/11.7.0 cuDNN/8.4.1.50-CUDA-11.7.0
conda activate "$ENV_PATH"

echo "Starting MD Job Dispatcher..."
echo "Max atoms allowed per simulation: $MAX_ATOMS"

# Create a logs directory to keep things clean
mkdir -p batch_logs

find "$ROOT_CIF_DIR" -name "*.cif" | sort | while read -r cif_path; do

    folder_name=$(basename "$(dirname "$cif_path")")
    file_base=$(basename "$cif_path" .cif)
    PREFIX="MD_${folder_name}_${ENGINE}"

    # --- SKIP LOGIC: Check if already done ---
    # Looks for the specific report file in any matching output folder
   
    if ls "${PREFIX}_${file_base}_"*/md_stability_analysis/*_md_stability_report.txt 1> /dev/null 2>&1; then
        echo "--> SKIP: Found completed report for $file_base"
        continue
    fi
    
    # Also skip if graph_break marker exists (indicates unstable but detected)
    if [ -f "${PREFIX}_${file_base}_graph_break_marker.txt" ]; then
        echo "--> SKIP: Found graph_break marker for $file_base (unstable detected)"
        continue
    fi

    # --- DYNAMIC SUPERCELL CALCULATION ---
    # Use python to get the number of atoms in the primitive cell
    PRIM_ATOMS=$(python -c "from ase.io import read; print(len(read('$cif_path')))" 2>/dev/null)
    
    # Priority list of supercells to try
    # 2x2x2 (Preferred) -> 2x2x1 -> 2x1x1 -> 1x1x1 (Fallback)
    SC_SIZE="1x1x1" # Default fallback
    
    # Calculate atom counts for different configurations
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

    # Check if even 1x1x1 exceeds the limit (unlikely but possible for huge unit cells)
    if [ "$ATOMS_111" -gt "$MAX_ATOMS" ]; then
        echo "  [WARNING] $file_base: Even 1x1x1 ($ATOMS_111 atoms) exceeds limit. Submitting anyway but expect OOM."
    fi

    # --- JOB SUBMISSION ---
    # Write a temporary SLURM script for THIS specific structure
    JOB_SCRIPT="temp_submit_${file_base}.slurm"
    LOG_FILE="batch_logs/md_${file_base}_%j.log"
    
    cat <<EOF > "$JOB_SCRIPT"
#!/bin/bash
#SBATCH --job-name=MD_${file_base}
#SBATCH --output=${LOG_FILE}
#SBATCH --partition=gpu 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --gpus=1
#SBATCH --time=6:00:00
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
# Fix for CUDA OOM fragmentation
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

echo "Processing: $file_base"
echo "Atoms in primitive: $PRIM_ATOMS"
echo "Selected Supercell: $SC_SIZE"

# Run VibroML with MACEMP (medium model)
vibroml --cif "$cif_path" \
    --method md_stability \
    --engine "$ENGINE" \
    --model_name "medium" \
    --temp "$TEMP" \
    --time "$TIME_PS" \
    --supercell-size "$SC_SIZE" \
    --output-prefix "$PREFIX"

# --- GRAPH BREAK DETECTION ---
# Check if MD failed due to graph break (atoms exceeding model cutoff)
# This indicates the structure is unstable, not an error to rerun
# The log file might not be immediately available with %j replaced, so check generically
sleep 2
for logfile in batch_logs/md_${file_base}_*.log; do
    if [ -f "$logfile" ]; then
        if grep -q "No edges found in input system" "$logfile" 2>/dev/null || \
           grep -q "atoms are farther apart than the radius cutoff" "$logfile" 2>/dev/null; then
            echo "[GRAPH_BREAK] Detected unstable structure (atoms exceeded model cutoff)"
            echo "Graph break detected at $(date)" > "${PREFIX}_${file_base}_graph_break_marker.txt"
            echo "Composition: $file_base" >> "${PREFIX}_${file_base}_graph_break_marker.txt"
            echo "This run should NOT be rerun - the material is dynamically unstable"
        fi
        break
    fi
done

EOF

    # Submit the job and delete the temp file
    sbatch "$JOB_SCRIPT"
    rm "$JOB_SCRIPT"
    
    # Optional: Sleep briefly to be kind to the scheduler
    sleep 1

done

echo "All jobs submitted."
