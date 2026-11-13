#!/usr/bin/env python3
"""Build the UQ validation manifest: 43 benchmark + 2 ProtoCSP + 50 quaternary + 10 quinary."""
import os, re, json, random
import numpy as np

random.seed(42)

base_remediation = "/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/02_runs_and_data/examples_VibroML/VibroML_SupportData/remediation_results"

# Reference CIFs
ref_cifs = {
    "LiF": os.path.join(base_remediation, "LiF_simple_cubic", "LiFsimplecubic.cif"),
    "CsPbI3": os.path.join(base_remediation, "CsPbI3", "CsPbI3.cif"),
    "HfO2": os.path.join(base_remediation, "HfO2", "HfO2_Fm3m_GA_phonon_output_20250826-161333", "initial_relaxation_for_single_run", "HfO2_Fm3m_relaxed.cif"),
    "Bi2Sn2O7": os.path.join(base_remediation, "Bi2Sn2O7", "Sn2Bi2O7_MACE_TRADITIONAL_ALL_phonon_output_20260123-012716", "initial_relaxation_for_single_run", "Sn2Bi2O7_relaxed.cif"),
}

# Load pre-computed data
with open("/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/02_runs_and_data/pet_phonons/benchmark_all_phonon_analyses.json") as f:
    all_entries = json.load(f)
with open("/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/02_runs_and_data/pet_phonons/tableS1_entries.json") as f:
    table_s1 = json.load(f)
with open("/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/02_runs_and_data/pet_phonons/tableS1_matched_by_phonon.json") as f:
    matched = json.load(f)

MANIFEST = []
idx = 0
used_cif_paths = set()

def count_atoms(cif_path):
    with open(cif_path) as f:
        content = f.read()
    lines = content.split('\n')
    in_loop = False
    atom_count = 0
    for line in lines:
        if '_atom_site_fract_x' in line:
            in_loop = True
            continue
        if in_loop and line.strip() and not line.startswith('_') and not line.startswith('loop_'):
            atom_count += 1
        elif in_loop and (line.startswith('loop_') or line.startswith('_')):
            break
    return atom_count if atom_count > 0 else 2

def get_mace_stats_from_dir(pd_path):
    """Get ω_min and soft_frac from a phonon output directory."""
    omega_min = None
    soft_frac = None
    soft_frac_15 = None
    natoms = 0
    
    for f in os.listdir(pd_path):
        if f.endswith("_phonon_run_summary.txt"):
            try:
                with open(os.path.join(pd_path, f)) as fh:
                    for line in fh:
                        if "Most Negative Frequency (overall):" in line:
                            omega_min = float(line.split(":")[-1].replace("THz", "").strip())
                        if "Total Number of Atoms:" in line:
                            natoms = int(line.split(":")[-1].strip())
            except:
                pass
            break
    
    bs_files = [f for f in os.listdir(pd_path) if f.startswith("band_structure_energies")]
    if bs_files:
        try:
            data = np.loadtxt(os.path.join(pd_path, bs_files[0]))
            soft_frac = float((data < -0.1).sum() / data.size)
            soft_frac_15 = float((data < -1.5).sum() / data.size)
        except:
            pass
    
    return omega_min, soft_frac, soft_frac_15, natoms

# 1. Add matched benchmark polymorphs
for m in matched:
    if m["cif_path"] and m["cif_path"] not in used_cif_paths:
        idx += 1
        used_cif_paths.add(m["cif_path"])
        # Determine supercell based on natoms
        sc = 2 if m["natoms"] and m["natoms"] <= 20 else 1
        MANIFEST.append({
            "order": idx,
            "category": "benchmark",
            "compound": m["tableS1"]["compound"],
            "space_group": m["tableS1"]["space_group"],
            "filename": os.path.basename(m["cif_path"]),
            "filepath": m["cif_path"],
            "natoms": m["natoms"],
            "supercell": sc,
            "mace_omega_min": m["mace_omega_min"],
            "mace_soft_frac": m["mace_soft_frac"],
            "mace_soft_frac_below_1p5": m["mace_soft_frac_below_1p5"],
        })

# 2. Add reference structures
ref_to_add = {
    "LiF": {"sg": "Pm\\bar{3}m", "omega_min": -6.07, "soft_frac": 0.2883, "sf_15": 0.2883},
    "CsPbI3": {"sg": "Pm\\bar{3}m", "omega_min": -0.74, "soft_frac": 0.2387, "sf_15": 0.2387},
    "Bi2Sn2O7": {"sg": "Fd\\bar{3}m", "omega_min": -4.15, "soft_frac": 0.0669, "sf_15": 0.0669},
}

for comp, info in ref_to_add.items():
    cif = ref_cifs[comp]
    if cif and os.path.exists(cif) and cif not in used_cif_paths:
        idx += 1
        used_cif_paths.add(cif)
        natoms = count_atoms(cif)
        MANIFEST.append({
            "order": idx,
            "category": "benchmark",
            "compound": comp,
            "space_group": info["sg"],
            "filename": os.path.basename(cif),
            "filepath": cif,
            "natoms": natoms,
            "supercell": 2 if natoms <= 20 else 1,
            "mace_omega_min": info["omega_min"],
            "mace_soft_frac": info["soft_frac"],
            "mace_soft_frac_below_1p5": info["sf_15"],
        })

# 3. Add remaining unmatched Bi2Sn2O7 entries by energy
ref_e_bi2sn2o7 = -6.0317
bi2sn2o7_unmatched = [
    {"sg": "P2_12_12_1", "rel_e": 20.6},
    {"sg": "P1", "rel_e": 0.2},
    {"sg": "P1", "rel_e": -12.2},
    {"sg": "Cc", "rel_e": -12.6},
    {"sg": "Pna2_1", "rel_e": -12.7},
    {"sg": "P2_1", "rel_e": -13.1},
    {"sg": "P2_12_12_1", "rel_e": -13.8},
    {"sg": "Cc", "rel_e": -14.1},
]

bi2sn2o7_entries = [e for e in all_entries if e["compound"] == "Bi2Sn2O7" and e["cif_path"] and e["cif_path"] not in used_cif_paths]

for target in bi2sn2o7_unmatched:
    target_e = ref_e_bi2sn2o7 + target["rel_e"] / 1000.0
    candidates = [e for e in bi2sn2o7_entries if e["cif_path"] not in used_cif_paths and e["energy_per_atom"] is not None]
    if candidates:
        best = min(candidates, key=lambda x: abs(x["energy_per_atom"] - target_e))
        if abs(best["energy_per_atom"] - target_e) * 1000 < 10:
            idx += 1
            used_cif_paths.add(best["cif_path"])
            MANIFEST.append({
                "order": idx,
                "category": "benchmark",
                "compound": "Bi2Sn2O7",
                "space_group": target["sg"],
                "filename": os.path.basename(best["cif_path"]),
                "filepath": best["cif_path"],
                "natoms": best["natoms"],
                "supercell": 1,  # Bi2Sn2O7 is large
                "mace_omega_min": best["mace_omega_min"],
                "mace_soft_frac": best["mace_soft_frac"],
                "mace_soft_frac_below_1p5": best["mace_soft_frac_below_1p5"],
            })

# 4. Add LiF Pm anomalous structure
ref_e_lif = -4.8212
target_e = ref_e_lif + (-7347.9) / 1000.0
lif_entries = [e for e in all_entries if e["compound"] == "LiF" and e["cif_path"] and e["cif_path"] not in used_cif_paths]
if lif_entries:
    candidates = [e for e in lif_entries if e["energy_per_atom"] is not None]
    if candidates:
        best = min(candidates, key=lambda x: abs(x["energy_per_atom"] - target_e))
        if abs(best["energy_per_atom"] - target_e) * 1000 < 200:
            idx += 1
            used_cif_paths.add(best["cif_path"])
            MANIFEST.append({
                "order": idx,
                "category": "benchmark",
                "compound": "LiF",
                "space_group": "Pm",
                "filename": os.path.basename(best["cif_path"]),
                "filepath": best["cif_path"],
                "natoms": best["natoms"],
                "supercell": 2,
                "mace_omega_min": best["mace_omega_min"],
                "mace_soft_frac": best["mace_soft_frac"],
                "mace_soft_frac_below_1p5": best["mace_soft_frac_below_1p5"],
            })

print(f"Benchmark structures: {len(MANIFEST)}")

# 5. Add ProtoCSP final validated structures
protocsp_base = {
    "Cs2K0.9Na0.1InI3Br3": "/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/02_runs_and_data/ProtoCSP_runs/cs2kini6",
    "K0.6Ba0.4Ta0.6Zr0.4Se2.1O0.9": "/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/02_runs_and_data/ProtoCSP_runs/ktase3",
}

for compound, pbase in protocsp_base.items():
    for f in os.listdir(pbase):
        if f.endswith(".cif") and compound in f:
            cif_path = os.path.join(pbase, f)
            if cif_path not in used_cif_paths:
                idx += 1
                used_cif_paths.add(cif_path)
                natoms = count_atoms(cif_path)
                MANIFEST.append({
                    "order": idx,
                    "category": "protocsp_final",
                    "compound": compound,
                    "space_group": "unknown",
                    "filename": f,
                    "filepath": cif_path,
                    "natoms": natoms,
                    "supercell": 1,
                    "mace_omega_min": 0.0,
                    "mace_soft_frac": 0.0,
                    "mace_soft_frac_below_1p5": 0.0,
                })
                print(f"Added ProtoCSP: {compound} ({natoms} atoms)")
            break

print(f"After ProtoCSP: {len(MANIFEST)}")

# 6. Add 50 random quaternary abandoned
quat_base = "/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/02_runs_and_data/candidates_lowcoverage_cif/quaternary_abandoned_all"
quat_cifs = sorted([f for f in os.listdir(quat_base) if f.endswith(".cif")])
quat_sample = random.sample(quat_cifs, min(50, len(quat_cifs)))

for cif_fname in quat_sample:
    cif_path = os.path.join(quat_base, cif_fname)
    base_name = cif_fname.replace(".cif", "")
    phonon_dirs = [d for d in os.listdir(quat_base) if d.startswith(base_name) and "_MACE_GA_phonon_output_" in d]
    
    omega_min, soft_frac, soft_frac_15, natoms = None, None, None, 0
    if phonon_dirs:
        omega_min, soft_frac, soft_frac_15, natoms = get_mace_stats_from_dir(
            os.path.join(quat_base, phonon_dirs[0]))
    
    idx += 1
    MANIFEST.append({
        "order": idx,
        "category": "quaternary_abandoned",
        "compound": cif_fname.split("_")[0],
        "space_group": "unknown",
        "filename": cif_fname,
        "filepath": cif_path,
        "natoms": natoms,
        "supercell": 2 if natoms <= 20 else 1,
        "mace_omega_min": omega_min,
        "mace_soft_frac": soft_frac,
        "mace_soft_frac_below_1p5": soft_frac_15,
    })

print(f"After quaternary: {len(MANIFEST)}")

# 7. Add 10 random quinary abandoned
quin_base = "/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/02_runs_and_data/candidates_lowcoverage_cif/quintenary_abandoned_all"
quin_cifs = sorted([f for f in os.listdir(quin_base) if f.endswith(".cif")])
quin_sample = random.sample(quin_cifs, min(10, len(quin_cifs)))

for cif_fname in quin_sample:
    cif_path = os.path.join(quin_base, cif_fname)
    base_name = cif_fname.replace(".cif", "")
    phonon_dirs = [d for d in os.listdir(quin_base) if d.startswith(base_name) and "_MACE_GA_phonon_output_" in d]
    
    omega_min, soft_frac, soft_frac_15, natoms = None, None, None, 0
    if phonon_dirs:
        omega_min, soft_frac, soft_frac_15, natoms = get_mace_stats_from_dir(
            os.path.join(quin_base, phonon_dirs[0]))
    
    idx += 1
    MANIFEST.append({
        "order": idx,
        "category": "quinary_abandoned",
        "compound": cif_fname.split("_")[0],
        "space_group": "unknown",
        "filename": cif_fname,
        "filepath": cif_path,
        "natoms": natoms,
        "supercell": 2 if natoms <= 20 else 1,
        "mace_omega_min": omega_min,
        "mace_soft_frac": soft_frac,
        "mace_soft_frac_below_1p5": soft_frac_15,
    })

print(f"After quinary: {len(MANIFEST)}")

# Summary
print(f"\n=== FINAL MANIFEST ===")
print(f"Total: {len(MANIFEST)} structures")
cats = {}
for m in MANIFEST:
    cats[m["category"]] = cats.get(m["category"], 0) + 1
for cat, count in sorted(cats.items()):
    print(f"  {cat}: {count}")

print(f"\nWith MACE stats: {sum(1 for m in MANIFEST if m['mace_omega_min'] is not None)}")
print(f"Without MACE stats: {sum(1 for m in MANIFEST if m['mace_omega_min'] is None)}")

omegas = [m["mace_omega_min"] for m in MANIFEST if m["mace_omega_min"] is not None]
print(f"\nω_min range: {min(omegas):.3f} to {max(omegas):.3f} THz")
print(f"  Stable (ω > -1.5): {sum(1 for o in omegas if o > -1.5)}")
print(f"  Quasi (-3 < ω <= -1.5): {sum(1 for o in omegas if -3 < o <= -1.5)}")
print(f"  Unstable (ω <= -3): {sum(1 for o in omegas if o <= -3)}")

# Save
out_path = "/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/02_runs_and_data/pet_phonons/uq_validation_structures.json"
with open(out_path, "w") as f:
    json.dump(MANIFEST, f, indent=2)
print(f"\nSaved to {out_path}")
