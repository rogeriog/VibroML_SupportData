#!/usr/bin/env python3
"""
Run the full UQ workflow for a single structure (the new CsPbI3 cubic).
Bypasses the resume logic of the comprehensive script by using a separate
structures JSON and a separate output CSV.
"""
import os, sys, json, time

# Reuse functions from the comprehensive script
sys.path.insert(0, os.path.dirname(__file__))
from run_uq_validation_comprehensive import (
    run_single_structure, save_summary_csv, load_mace_icurv
)

WORKDIR = os.path.dirname(os.path.abspath(__file__))
STRUCTURES_JSON = os.path.join(WORKDIR, "cspbi3_cubic_only.json")
SUMMARY_CSV = os.path.join(WORKDIR, "uq_validation_summary.csv")
MACE_ICURV_JSON = os.path.join(WORKDIR, "omega_inst_cspbi3_cubic_mace.json")

# Load the new structure
with open(STRUCTURES_JSON) as f:
    structures = json.load(f)

entry = structures[0]
# Use idx=106 (position in the full uq_validation_structures.json)
# so it integrates with the existing CSV
IDX = 106

print(f"Processing CsPbI3 cubic Pm-3m (idx={IDX})")
print(f"  filepath: {entry['filepath']}")
print(f"  natoms: {entry['natoms']}")
print(f"  supercell: {entry['supercell']}")

# Load MACE I_curv from the just-computed computation
mace_icurv = {}
if os.path.exists(MACE_ICURV_JSON):
    with open(MACE_ICURV_JSON) as f:
        mace_data = json.load(f)
    if mace_data:
        r = mace_data[0]
        mace_icurv = {
            "mace_i_curv": r.get("i_curv"),
            "mace_omega_inst": r.get("omega_inst"),
            "mace_i_curv_omega_min": r.get("omega_min"),
            "mace_i_curv_soft_frac": r.get("soft_frac"),
            "mace_i_curv_n_imag": r.get("n_imag"),
            "mace_i_curv_n_modes": r.get("n_modes"),
        }
        print(f"  MACE I_curv: {mace_icurv['mace_i_curv']:.4f} THz²")
else:
    print("  WARNING: MACE I_curv JSON not found")

# Load PET-MAD calculator
print("\nLoading PET-MAD calculator...")
from upet.calculator import UPETCalculator
pet_calc = UPETCalculator(model="pet-mad-xs", device="cuda")

# Run the full workflow
result = run_single_structure(entry, IDX, 1, pet_calc, mace_icurv)

if result:
    # Append to the existing summary CSV
    # Load existing results
    import csv
    existing = []
    if os.path.exists(SUMMARY_CSV):
        with open(SUMMARY_CSV) as f:
            reader = csv.DictReader(f)
            existing = list(reader)

    # Remove any stale entry with the same idx
    existing = [r for r in existing if int(r.get('idx', 0)) != IDX]

    # Add new result
    existing.append(result)
    save_summary_csv(existing)
    print(f"\nResult appended to {SUMMARY_CSV}")
    print(f"\nKey metrics:")
    print(f"  MACE I_curv:     {result.get('mace_i_curv', '?')} THz²")
    print(f"  PET I_curv:      {result.get('pet_i_curv', '?')} THz²")
    print(f"  PET omega_min:   {result.get('pet_omega_min', '?')} THz")
    print(f"  PET std_at_min:  {result.get('pet_std_at_min_thz', '?')} THz")
    print(f"  PET max_std:     {result.get('pet_max_std_thz', '?')} THz")
else:
    print("FAILED to compute result")
    sys.exit(1)
