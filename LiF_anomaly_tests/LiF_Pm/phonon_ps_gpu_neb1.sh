#!/bin/bash
#SBATCH --job-name=vibroml_neb
#SBATCH --time=2:00:00
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

vibroml --cif LiFsimplecubic_sg_P63mc.cif \
    --final_cif OPTRANDOM_sg_Pm_top_4_iterfinal_sample4_LiFsimplecubic_energy_m11p8712_primitive_freqp511p2674THz.cif \
    --method neb \
    --engine mace \
    --neb_num_images 7 \
    --neb_spring_constant 3.0 \
    --neb_max_iterations 2000 \
    --neb_force_tolerance 0.25

echo "done"
date

