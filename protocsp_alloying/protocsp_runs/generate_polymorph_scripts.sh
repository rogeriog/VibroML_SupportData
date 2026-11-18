#!/bin/bash

# List of polymorph CIF files
polymorphs=("Cs2KInI6_top1_I42m.cif" "Cs2KInI6_top2_P-31m.cif" "Cs2KInI6_top3_C2c.cif" "Cs2KInI6_top4_P-3.cif" "Cs2KInI6_top5_P-1.cif" "Cs2KInI6_top6_Fmm2.cif" "Cs2KInI6_top7_P2_1c.cif")

# Base script template
base_script="#!/bin/bash
#SBATCH --job-name=protocsp_Cs2KInI6_POLY
#SBATCH --time=2:00:00
#SBATCH --output=log_cs2kini6_protocsp_POLY.txt
#SBATCH --partition=gpu
#SBATCH --gpus=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --account=htforft

# 1. Clean Environment to prevent NumPy path pollution
module purge

source ~/.bashrc

module load EasyBuild/2022a CUDA/11.7.0 cuDNN/8.4.1.50-CUDA-11.7.0

conda activate /gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/vibroml_env
export OMP_NUM_THREADS=\$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=\$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=\$SLURM_CPUS_PER_TASK
export PYTHONUNBUFFERED=1

unset TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD
# Keep this to prevent deadlock
export MP_START_METHOD=spawn

echo \"start\"
date

# x = 0.0 (Pure K on B'-site)
mkdir -p cs2kini6_phase_diagram_POLY/x_0p0/y_{0p00,0p17,0p33,0p50,0p67}
python protocsp/main.py \"Cs2KInI6\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir cs2kini6_phase_diagram_POLY/x_0p0/y_0p00 --symmetrize --mlip > cs2kini6_phase_diagram_POLY/x_0p0/y_0p00/log.log
python protocsp/main.py \"Cs2KInI5Br\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir cs2kini6_phase_diagram_POLY/x_0p0/y_0p17 --symmetrize --mlip > cs2kini6_phase_diagram_POLY/x_0p0/y_0p17/log.log
python protocsp/main.py \"Cs2KInI4Br2\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir cs2kini6_phase_diagram_POLY/x_0p0/y_0p33 --symmetrize --mlip > cs2kini6_phase_diagram_POLY/x_0p0/y_0p33/log.log
python protocsp/main.py \"Cs2KInI3Br3\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir cs2kini6_phase_diagram_POLY/x_0p0/y_0p50 --symmetrize --mlip > cs2kini6_phase_diagram_POLY/x_0p0/y_0p50/log.log
python protocsp/main.py \"Cs2KInI2Br4\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir cs2kini6_phase_diagram_POLY/x_0p0/y_0p67 --symmetrize --mlip > cs2kini6_phase_diagram_POLY/x_0p0/y_0p67/log.log

# x = 0.1 (10% Na substitution)
mkdir -p cs2kini6_phase_diagram_POLY/x_0p1/y_{0p00,0p17,0p33,0p50,0p67}
python protocsp/main.py \"Cs2K0.9Na0.1InI6\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir cs2kini6_phase_diagram_POLY/x_0p1/y_0p00 --symmetrize --mlip > cs2kini6_phase_diagram_POLY/x_0p1/y_0p00/log.log
python protocsp/main.py \"Cs2K0.9Na0.1InI5Br\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir cs2kini6_phase_diagram_POLY/x_0p1/y_0p17 --symmetrize --mlip > cs2kini6_phase_diagram_POLY/x_0p1/y_0p17/log.log
python protocsp/main.py \"Cs2K0.9Na0.1InI4Br2\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir cs2kini6_phase_diagram_POLY/x_0p1/y_0p33 --symmetrize --mlip > cs2kini6_phase_diagram_POLY/x_0p1/y_0p33/log.log
python protocsp/main.py \"Cs2K0.9Na0.1InI3Br3\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir cs2kini6_phase_diagram_POLY/x_0p1/y_0p50 --symmetrize --mlip > cs2kini6_phase_diagram_POLY/x_0p1/y_0p50/log.log
python protocsp/main.py \"Cs2K0.9Na0.1InI2Br4\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir cs2kini6_phase_diagram_POLY/x_0p1/y_0p67 --symmetrize --mlip > cs2kini6_phase_diagram_POLY/x_0p1/y_0p67/log.log

# x = 0.2 (20% Na substitution)
mkdir -p cs2kini6_phase_diagram_POLY/x_0p2/y_{0p00,0p17,0p33,0p50,0p67}
python protocsp/main.py \"Cs2K0.8Na0.2InI6\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir cs2kini6_phase_diagram_POLY/x_0p2/y_0p00 --symmetrize --mlip > cs2kini6_phase_diagram_POLY/x_0p2/y_0p00/log.log
python protocsp/main.py \"Cs2K0.8Na0.2InI5Br\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir cs2kini6_phase_diagram_POLY/x_0p2/y_0p17 --symmetrize --mlip > cs2kini6_phase_diagram_POLY/x_0p2/y_0p17/log.log
python protocsp/main.py \"Cs2K0.8Na0.2InI4Br2\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir cs2kini6_phase_diagram_POLY/x_0p2/y_0p33 --symmetrize --mlip > cs2kini6_phase_diagram_POLY/x_0p2/y_0p33/log.log
python protocsp/main.py \"Cs2K0.8Na0.2InI3Br3\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir cs2kini6_phase_diagram_POLY/x_0p2/y_0p50 --symmetrize --mlip > cs2kini6_phase_diagram_POLY/x_0p2/y_0p50/log.log
python protocsp/main.py \"Cs2K0.8Na0.2InI2Br4\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir cs2kini6_phase_diagram_POLY/x_0p2/y_0p67 --symmetrize --mlip > cs2kini6_phase_diagram_POLY/x_0p2/y_0p67/log.log

# x = 0.3 (30% Na substitution)
mkdir -p cs2kini6_phase_diagram_POLY/x_0p3/y_{0p00,0p17,0p33,0p50,0p67}
python protocsp/main.py \"Cs2K0.7Na0.3InI6\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir cs2kini6_phase_diagram_POLY/x_0p3/y_0p00 --symmetrize --mlip > cs2kini6_phase_diagram_POLY/x_0p3/y_0p00/log.log
python protocsp/main.py \"Cs2K0.7Na0.3InI5Br\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir cs2kini6_phase_diagram_POLY/x_0p3/y_0p17 --symmetrize --mlip > cs2kini6_phase_diagram_POLY/x_0p3/y_0p17/log.log
python protocsp/main.py \"Cs2K0.7Na0.3InI4Br2\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir cs2kini6_phase_diagram_POLY/x_0p3/y_0p33 --symmetrize --mlip > cs2kini6_phase_diagram_POLY/x_0p3/y_0p33/log.log
python protocsp/main.py \"Cs2K0.7Na0.3InI3Br3\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir cs2kini6_phase_diagram_POLY/x_0p3/y_0p50 --symmetrize --mlip > cs2kini6_phase_diagram_POLY/x_0p3/y_0p50/log.log
python protocsp/main.py \"Cs2K0.7Na0.3InI2Br4\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir cs2kini6_phase_diagram_POLY/x_0p3/y_0p67 --symmetrize --mlip > cs2kini6_phase_diagram_POLY/x_0p3/y_0p67/log.log

# x = 0.4 (40% Na substitution)
mkdir -p cs2kini6_phase_diagram_POLY/x_0p4/y_{0p00,0p17,0p33,0p50,0p67}
python protocsp/main.py \"Cs2K0.6Na0.4InI6\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir cs2kini6_phase_diagram_POLY/x_0p4/y_0p00 --symmetrize --mlip > cs2kini6_phase_diagram_POLY/x_0p4/y_0p00/log.log
python protocsp/main.py \"Cs2K0.6Na0.4InI5Br\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir cs2kini6_phase_diagram_POLY/x_0p4/y_0p17 --symmetrize --mlip > cs2kini6_phase_diagram_POLY/x_0p4/y_0p17/log.log
python protocsp/main.py \"Cs2K0.6Na0.4InI4Br2\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir cs2kini6_phase_diagram_POLY/x_0p4/y_0p33 --symmetrize --mlip > cs2kini6_phase_diagram_POLY/x_0p4/y_0p33/log.log
python protocsp/main.py \"Cs2K0.6Na0.4InI3Br3\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir cs2kini6_phase_diagram_POLY/x_0p4/y_0p50 --symmetrize --mlip > cs2kini6_phase_diagram_POLY/x_0p4/y_0p50/log.log
python protocsp/main.py \"Cs2K0.6Na0.4InI2Br4\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir cs2kini6_phase_diagram_POLY/x_0p4/y_0p67 --symmetrize --mlip > cs2kini6_phase_diagram_POLY/x_0p4/y_0p67/log.log

echo \"done\"
date"

for i in "${!polymorphs[@]}"; do
    poly_num=$((i+1))
    poly_name="${polymorphs[$i]}"
    script_name="cg_submit_gpu6_Cs2KInI6_top${poly_num}.sh"
    echo "Generating $script_name for $poly_name"
    echo "$base_script" | sed "s/POLY/top${poly_num}/g" | sed "s/BASE_CIF/$poly_name/g" > "$script_name"
    chmod +x "$script_name"
done

echo "All scripts generated."