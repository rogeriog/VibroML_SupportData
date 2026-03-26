import os
import sys
import glob
import numpy as np
import pandas as pd
from ase.io import read

# Ensure we can import vibroml from the root directory
sys.path.append('/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/VibroML')
from vibroml.utils.structure_utils import initialize_calculator

def run_worker():
    if len(sys.argv) < 3:
        print("Usage: python worker_eval.py <engine_key> <cif_pattern>")
        sys.exit(1)
    
    engine_key = sys.argv[1] # e.g., 'mace_omat', 'mace_mp', 'esen', 'uma'
    cif_pattern = sys.argv[2]
    
    # Model configuration mapping
    # MACE engines use 'mace' engine type but different model_names
    config = {
        "mace_omat": {"engine": "mace", "model": "medium-omat-0", "path": None},
        "mace_mp":   {"engine": "mace", "model": "medium",        "path": None},
        "esen":      {"engine": "esen", "model": None,            "path": "/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/VibroML/fairchem_models/esen_30m_omat.pt"},
        "uma":       {"engine": "uma",  "model": None,            "path": "/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/VibroML/fairchem_models/uma-m-1p1.pt"}
    }
    
    if engine_key not in config:
        print(f"Error: Engine key {engine_key} not recognized.")
        sys.exit(1)
        
    setup = config[engine_key]

    # Initialize calculator
    calc = initialize_calculator(
        engine=setup["engine"], 
        model_name=setup["model"],
        checkpoint_model_path=setup["path"]
    )
    
    if calc is None:
        print(f"Failed to initialize {engine_key}")
        sys.exit(1)

    cif_files = sorted(glob.glob(cif_pattern))
    results = []

    print(f"Engine: {engine_key} | Evaluating {len(cif_files)} structures...")

    for cif in cif_files:
        try:
            atoms = read(cif)
            atoms.calc = calc
            
            # Divide total energy by number of atoms for eV/atom
            total_energy = atoms.get_potential_energy()
            energy_per_atom = total_energy / len(atoms)
            
            forces = atoms.get_forces()
            max_force = np.max(np.linalg.norm(forces, axis=1))
            
            results.append({
                "filename": os.path.basename(cif),
                f"{engine_key}_energy_ev_per_atom": energy_per_atom,
                f"{engine_key}_max_force": max_force
            })
        except Exception as e:
            print(f"Error processing {cif}: {e}")

    # Save results to a temporary CSV
    df = pd.DataFrame(results)
    output_name = f"results_{engine_key}.csv"
    df.to_csv(output_name, index=False)
    print(f"Saved results to {output_name}")

if __name__ == "__main__":
    run_worker()