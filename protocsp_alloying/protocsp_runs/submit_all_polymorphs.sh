#!/bin/bash

# Submit all polymorph jobs
# for i in {1..7}; do
#     script="cg_submit_gpu6_Cs2KInI6_top${i}.sh"
#     if [ -f "$script" ]; then
#         echo "Submitting $script"
#         sbatch "$script"
#     else
#         echo "Warning: $script not found"
#     fi
# done




for i in {1..3}; do
    script="cg_submit_gpu5_KTaSe3_top_2d_${i}.sh"
    if [ -f "$script" ]; then
        echo "Submitting $script"
        sbatch "$script"
    else
        echo "Warning: $script not found"
    fi
done


# Submit all polymorph jobs
for i in {1..7}; do
    script="cg_submit_gpu5_KTaSe3_top${i}.sh"
    if [ -f "$script" ]; then
        echo "Submitting $script"
        sbatch "$script"
    else
        echo "Warning: $script not found"
    fi
done

echo "All polymorph jobs submitted."