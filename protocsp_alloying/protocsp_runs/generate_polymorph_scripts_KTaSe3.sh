#!/bin/bash

# List of 'top' polymorph CIF files
polymorphs_top=(
  "KTaSe3_top1_P2_1_c_0p00meV.cif"
  "KTaSe3_top2_P1_2p77meV.cif"
  "KTaSe3_top3_P-1_4p40meV.cif"
  "KTaSe3_top4_P1_5p48meV.cif"
  "KTaSe3_top5_P1_6p16meV.cif"
  "KTaSe3_top6_P2_1_m_7p42meV.cif"
  "KTaSe3_top7_P1_8p23meV.cif"
)

# List of 'top_2d' (unique) polymorph CIF files
polymorphs_2d=(
  "KTaSe3_unique_23.cif"
  "KTaSe3_unique_27.cif"
  "KTaSe3_unique_29.cif"
)

# Base script template
base_script="#!/bin/bash
#SBATCH --job-name=protocsp_KTaSe3_POLY
#SBATCH --time=2:00:00
#SBATCH --output=log_ktase3_protocsp_POLY.txt
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

# ==========================================
# COMMENTED OUT: OLD RUNS THAT ALREADY RAN
# ==========================================

# x = 0.0 (Pure K and Ta)
# mkdir -p ktase3_phase_diagram_POLY/x_0p0/y_{0p0,0p1,0p2,0p3}
# python protocsp/main.py \"KTaSe3\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir ktase3_phase_diagram_POLY/x_0p0/y_0p0 --symmetrize --mlip > ktase3_phase_diagram_POLY/x_0p0/y_0p0/log.log
# python protocsp/main.py \"KTaSe2.7O0.3\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir ktase3_phase_diagram_POLY/x_0p0/y_0p1 --symmetrize --mlip > ktase3_phase_diagram_POLY/x_0p0/y_0p1/log.log
# python protocsp/main.py \"KTaSe2.4O0.6\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir ktase3_phase_diagram_POLY/x_0p0/y_0p2 --symmetrize --mlip > ktase3_phase_diagram_POLY/x_0p0/y_0p2/log.log
# python protocsp/main.py \"KTaSe2.1O0.9\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir ktase3_phase_diagram_POLY/x_0p0/y_0p3 --symmetrize --mlip > ktase3_phase_diagram_POLY/x_0p0/y_0p3/log.log

# x = 0.1 (10% Ba and Zr substitution)
# mkdir -p ktase3_phase_diagram_POLY/x_0p1/y_{0p0,0p1,0p2,0p3}
# python protocsp/main.py \"K0.9Ba0.1Ta0.9Zr0.1Se3\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir ktase3_phase_diagram_POLY/x_0p1/y_0p0 --symmetrize --mlip > ktase3_phase_diagram_POLY/x_0p1/y_0p0/log.log
# python protocsp/main.py \"K0.9Ba0.1Ta0.9Zr0.1Se2.7O0.3\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir ktase3_phase_diagram_POLY/x_0p1/y_0p1 --symmetrize --mlip > ktase3_phase_diagram_POLY/x_0p1/y_0p1/log.log
# python protocsp/main.py \"K0.9Ba0.1Ta0.9Zr0.1Se2.4O0.6\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir ktase3_phase_diagram_POLY/x_0p1/y_0p2 --symmetrize --mlip > ktase3_phase_diagram_POLY/x_0p1/y_0p2/log.log
# python protocsp/main.py \"K0.9Ba0.1Ta0.9Zr0.1Se2.1O0.9\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir ktase3_phase_diagram_POLY/x_0p1/y_0p3 --symmetrize --mlip > ktase3_phase_diagram_POLY/x_0p1/y_0p3/log.log

# x = 0.2 (20% Ba and Zr substitution)
# mkdir -p ktase3_phase_diagram_POLY/x_0p2/y_{0p0,0p1,0p2,0p3}
# python protocsp/main.py \"K0.8Ba0.2Ta0.8Zr0.2Se3\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir ktase3_phase_diagram_POLY/x_0p2/y_0p0 --symmetrize --mlip > ktase3_phase_diagram_POLY/x_0p2/y_0p0/log.log
# python protocsp/main.py \"K0.8Ba0.2Ta0.8Zr0.2Se2.7O0.3\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir ktase3_phase_diagram_POLY/x_0p2/y_0p1 --symmetrize --mlip > ktase3_phase_diagram_POLY/x_0p2/y_0p1/log.log
# python protocsp/main.py \"K0.8Ba0.2Ta0.8Zr0.2Se2.4O0.6\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir ktase3_phase_diagram_POLY/x_0p2/y_0p2 --symmetrize --mlip > ktase3_phase_diagram_POLY/x_0p2/y_0p2/log.log
# python protocsp/main.py \"K0.8Ba0.2Ta0.8Zr0.2Se2.1O0.9\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir ktase3_phase_diagram_POLY/x_0p2/y_0p3 --symmetrize --mlip > ktase3_phase_diagram_POLY/x_0p2/y_0p3/log.log

# x = 0.3 (30% Ba and Zr substitution)
# mkdir -p ktase3_phase_diagram_POLY/x_0p3/y_{0p0,0p1,0p2,0p3}
# python protocsp/main.py \"K0.7Ba0.3Ta0.7Zr0.3Se3\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir ktase3_phase_diagram_POLY/x_0p3/y_0p0 --symmetrize --mlip > ktase3_phase_diagram_POLY/x_0p3/y_0p0/log.log
# python protocsp/main.py \"K0.7Ba0.3Ta0.7Zr0.3Se2.7O0.3\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir ktase3_phase_diagram_POLY/x_0p3/y_0p1 --symmetrize --mlip > ktase3_phase_diagram_POLY/x_0p3/y_0p1/log.log
# python protocsp/main.py \"K0.7Ba0.3Ta0.7Zr0.3Se2.4O0.6\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir ktase3_phase_diagram_POLY/x_0p3/y_0p2 --symmetrize --mlip > ktase3_phase_diagram_POLY/x_0p3/y_0p2/log.log
# python protocsp/main.py \"K0.7Ba0.3Ta0.7Zr0.3Se2.1O0.9\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir ktase3_phase_diagram_POLY/x_0p3/y_0p3 --symmetrize --mlip > ktase3_phase_diagram_POLY/x_0p3/y_0p3/log.log


# ==========================================
# NEW RUNS: O = 1.2 and Ba/Zr = 0.4
# ==========================================

# Adding O = 1.2 (y = 0.4) for previous x values
mkdir -p ktase3_phase_diagram_POLY/x_0p0/y_0p4
python protocsp/main.py \"KTaSe1.8O1.2\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir ktase3_phase_diagram_POLY/x_0p0/y_0p4 --symmetrize --mlip > ktase3_phase_diagram_POLY/x_0p0/y_0p4/log.log

mkdir -p ktase3_phase_diagram_POLY/x_0p1/y_0p4
python protocsp/main.py \"K0.9Ba0.1Ta0.9Zr0.1Se1.8O1.2\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir ktase3_phase_diagram_POLY/x_0p1/y_0p4 --symmetrize --mlip > ktase3_phase_diagram_POLY/x_0p1/y_0p4/log.log

mkdir -p ktase3_phase_diagram_POLY/x_0p2/y_0p4
python protocsp/main.py \"K0.8Ba0.2Ta0.8Zr0.2Se1.8O1.2\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir ktase3_phase_diagram_POLY/x_0p2/y_0p4 --symmetrize --mlip > ktase3_phase_diagram_POLY/x_0p2/y_0p4/log.log

mkdir -p ktase3_phase_diagram_POLY/x_0p3/y_0p4
python protocsp/main.py \"K0.7Ba0.3Ta0.7Zr0.3Se1.8O1.2\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir ktase3_phase_diagram_POLY/x_0p3/y_0p4 --symmetrize --mlip > ktase3_phase_diagram_POLY/x_0p3/y_0p4/log.log

# x = 0.4 (40% Ba and Zr substitution) - all O values
mkdir -p ktase3_phase_diagram_POLY/x_0p4/y_{0p0,0p1,0p2,0p3,0p4}
python protocsp/main.py \"K0.6Ba0.4Ta0.6Zr0.4Se3\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir ktase3_phase_diagram_POLY/x_0p4/y_0p0 --symmetrize --mlip > ktase3_phase_diagram_POLY/x_0p4/y_0p0/log.log
python protocsp/main.py \"K0.6Ba0.4Ta0.6Zr0.4Se2.7O0.3\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir ktase3_phase_diagram_POLY/x_0p4/y_0p1 --symmetrize --mlip > ktase3_phase_diagram_POLY/x_0p4/y_0p1/log.log
python protocsp/main.py \"K0.6Ba0.4Ta0.6Zr0.4Se2.4O0.6\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir ktase3_phase_diagram_POLY/x_0p4/y_0p2 --symmetrize --mlip > ktase3_phase_diagram_POLY/x_0p4/y_0p2/log.log
python protocsp/main.py \"K0.6Ba0.4Ta0.6Zr0.4Se2.1O0.9\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir ktase3_phase_diagram_POLY/x_0p4/y_0p3 --symmetrize --mlip > ktase3_phase_diagram_POLY/x_0p4/y_0p3/log.log
python protocsp/main.py \"K0.6Ba0.4Ta0.6Zr0.4Se1.8O1.2\" --base-cif BASE_CIF --top-k 20 --save-cif --output-dir ktase3_phase_diagram_POLY/x_0p4/y_0p4 --symmetrize --mlip > ktase3_phase_diagram_POLY/x_0p4/y_0p4/log.log

echo \"done\"
date"

# 1. Generate scripts for 'top' polymorphs
echo "Generating scripts for 'top' polymorphs..."
for i in "${!polymorphs_top[@]}"; do
    poly_num=$((i+1))
    poly_name="${polymorphs_top[$i]}"
    script_name="cg_submit_gpu5_KTaSe3_top${poly_num}.sh"
    echo "  -> $script_name for $poly_name"
    echo "$base_script" | sed "s/POLY/top${poly_num}/g" | sed "s/BASE_CIF/$poly_name/g" > "$script_name"
    chmod +x "$script_name"
done

# 2. Generate scripts for 'top_2d' (unique) polymorphs
echo "Generating scripts for 'top_2d' polymorphs..."
for i in "${!polymorphs_2d[@]}"; do
    poly_num=$((i+1))
    poly_name="${polymorphs_2d[$i]}"
    script_name="cg_submit_gpu5_KTaSe3_top_2d_${poly_num}.sh"
    echo "  -> $script_name for $poly_name"
    echo "$base_script" | sed "s/POLY/2d_top${poly_num}/g" | sed "s/BASE_CIF/$poly_name/g" > "$script_name"
    chmod +x "$script_name"
done

echo "All 10 scripts generated successfully!"