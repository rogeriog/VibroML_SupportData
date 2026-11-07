#!/bin/bash
# Run all 63 VASP single points for the B2 NEB image verification.
# Usage (on a machine with vasp_std): bash run_all_vasp.sh [ncores]
# Produces OUTCAR in each image directory; afterwards run decompose_forces.py.

NCORES=${1:-24}
cd "$(dirname "$0")/dft_inputs/neb_images_to_gs" || exit 1

for path_dir in */; do
    for img_dir in "${path_dir}"*/; do
        if [ -f "${img_dir}OUTCAR" ] && grep -q "reached required accuracy" "${img_dir}OUTCAR"; then
            echo "skip (done): $img_dir"
            continue
        fi
        echo "running: $img_dir"
        (cd "$img_dir" && mpirun -np "$NCORES" vasp_std > vasp_stdout.txt 2>&1)
    done
done
echo "All images complete."
