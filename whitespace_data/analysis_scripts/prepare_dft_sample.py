"""Prepare VASP single-point inputs for the white-space DFT verification
(GMR item B): a random subsample of GA winners vs their starting structures.

Samples 200 systems (seed 42) uniformly from the 2,649 white-space systems
that have a deposited top-1 GA structure. For each system writes two run
dirs: <name>_initial (Alexandria input) and <name>_winner (top-1 GA
structure). INCAR/KPOINTS conventions match the earlier B2 NEB statics
package (ENCUT 520, ALGO Fast, pymatgen grid density 3000).

POTCARs are NOT written (no PP library on this machine) — generate them on
the VASP machine with make_potcars.py.
"""
import os
import random
import sqlite3
import json
from pymatgen.core import Structure
from pymatgen.io.vasp import Poscar, Incar, Kpoints

DB = "/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/02_runs_and_data/" \
     "VibroML_SupportData/whitespace_data/vibroml_database.db"
CIFROOT = "/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/" \
          "02_runs_and_data/candidates_lowcoverage_cif"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "dft_inputs", "whitespace_winners")
N_SAMPLE = 200
SEED = 42

INCAR_DICT = {
    "ALGO": "Fast", "EDIFF": 1e-6, "ENAUG": 1360, "ENCUT": 520,
    "GGA": "Pe", "IBRION": -1, "ICHARG": 2, "ISIF": 0, "ISMEAR": 0,
    "SIGMA": 0.05, "ISPIN": 2, "LASPH": True, "LREAL": False, "NSW": 0,
    "PREC": "Accurate", "LWAVE": False, "LCHARG": False, "NCORE": 4,
}

db = sqlite3.connect(DB)
pool = db.execute(
    """SELECT m.id, m.material_name, m.dataset, m.cif_path, g.cif_path,
              m.initial_relaxed_energy_per_atom_eV, g.energy_per_atom_eV
       FROM ga_structures g JOIN materials m ON g.material_id = m.id
       WHERE g.is_final_top = 1 AND g.rank = 1""").fetchall()

# keep only systems where both CIFs are actually deposited
pool = [r for r in pool
        if os.path.exists(os.path.join(CIFROOT, r[3]))
        and os.path.exists(os.path.join(CIFROOT, r[4]))]

random.seed(SEED)
sample = random.sample(pool, N_SAMPLE)
random.shuffle(sample)

manifest = []
missing = []
elements = set()
for mid, name, dataset, cif_init, cif_win, e_init, e_win in sample:
    p_init = os.path.join(CIFROOT, cif_init)
    p_win = os.path.join(CIFROOT, cif_win)
    if not (os.path.exists(p_init) and os.path.exists(p_win)):
        missing.append(name)
        continue
    for tag, path in (("initial", p_init), ("winner", p_win)):
        s = Structure.from_file(path)
        run = f"{dataset}_{name}_{tag}"
        rd = os.path.join(OUT, run)
        os.makedirs(rd, exist_ok=True)
        Poscar(s).write_file(os.path.join(rd, "POSCAR"))
        Incar(INCAR_DICT).write_file(os.path.join(rd, "INCAR"))
        Kpoints.automatic_density(s, 3000).write_file(
            os.path.join(rd, "KPOINTS"))
        elements.update(s.composition.elements)
    manifest.append({
        "run_dir": f"{dataset}_{name}",
        "material": name, "dataset": dataset,
        "e_mace_init": e_init, "e_mace_winner": e_win,
        "dE_mace": (e_win - e_init) if e_init is not None else None,
        "elements": sorted(str(el) for el in
                           Structure.from_file(p_init).composition.elements),
    })

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "sample_manifest.json"), "w") as f:
    json.dump(manifest, f, indent=1)

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "elements_needed.txt"), "w") as f:
    f.write("\n".join(sorted(e.symbol for e in elements)) + "\n")

print(f"pool={len(pool)}  sampled={len(sample)}  written={len(manifest)}  "
      f"missing_cifs={len(missing)}")
print(f"elements needed: {len(elements)}")
if missing:
    print("missing:", missing[:10])
