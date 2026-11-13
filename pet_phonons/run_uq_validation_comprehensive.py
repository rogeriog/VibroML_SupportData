#!/usr/bin/env python3
"""
Comprehensive PET-MAD UQ validation for VibroML stability metrics.

For each structure, this script:
  1. Performs symmetry-constrained relaxation with PET-MAD (preserving
     the initial space group, so MACE and PET-MAD evaluate the SAME structure)
  2. Saves the relaxed CIF
  3. Calls vibroml.main --no-relax --phonon-uq on the relaxed CIF, which
     generates:
       - Phonon band structure plots (ensemble + mean±std) as PNG
       - UQ summary text file (min_freq, max_std, std_at_min, min_minus_3std)
       - npz data file (mean_freqs, std_freqs per q-point per band)
       - initial_settings.json
  4. Computes I_curv (mean negative curvature, THz^2) on a uniform q-mesh
     using PET-MAD mean forces
  5. Copies original MACE phonon data (band structure, plots) alongside
  6. Saves all results incrementally to a summary CSV

Resume logic: detects completed structures from existing output dirs and CSV,
skipping them on re-submission.

Usage:
  python run_uq_validation_comprehensive.py [--dry-run]
"""

import os, sys, json, subprocess, time, csv, shutil, re
import numpy as np
from pathlib import Path

STRUCTURES_JSON = os.environ.get(
    "UQ_STRUCTURES_JSON",
    os.path.join(os.path.dirname(__file__), "uq_validation_structures.json"))
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "uq_validation_runs")
SUMMARY_CSV = os.path.join(os.path.dirname(__file__), "uq_validation_summary.csv")
MACE_ICURV_JSON = os.path.join(os.path.dirname(__file__), "omega_inst_full_mace.json")

# Use orb_env Python (has upet/metatrain/uqphonon); VibroML is on PYTHONPATH
VIBROML_PYTHON = os.environ.get("VIBROML_PYTHON", sys.executable)
VIBROML_MODULE = "vibroml.main"

# Q-mesh for I_curv computation (uniform)
I_CURV_Q_MESH = (6, 6, 6)
ACOUSTIC_THRESHOLD = 0.05  # THz

EV_TO_THZ = 4.135667696e-3


def get_space_group(atoms, symprec=0.01):
    """Get the space group number and symbol."""
    try:
        import spglib
        cell = (atoms.cell.array, atoms.get_scaled_positions(), atoms.numbers)
        sg = spglib.get_spacegroup(cell, symprec=symprec)
        ds = spglib.get_symmetry_dataset(cell, symprec=symprec)
        return ds["number"], str(sg)
    except Exception:
        return -1, "unknown"


def symmetry_constrained_relax(cif_path, calc, fmax=1e-3, max_steps=200, symprec=0.01):
    """Load CIF, relax with symmetry constraint, save relaxed CIF, return path."""
    from ase.io import read, write
    from ase.optimize import FIRE
    from ase.constraints import FixSymmetry

    atoms = read(cif_path)
    sg_no, sg_sym = get_space_group(atoms, symprec)

    atoms.calc = calc
    constraint = FixSymmetry(atoms, symprec=symprec,
                             adjust_positions=True, adjust_cell=True)
    atoms.set_constraint(constraint)

    try:
        opt = FIRE(atoms, logfile=None)
        opt.run(fmax=fmax, steps=max_steps)
        print(f"    Sym-constrained relaxed: fmax={atoms.get_forces().max():.6f} eV/A")
    except Exception as e:
        print(f"    Relaxation failed: {e}, using original structure")

    # Remove constraint before saving
    atoms.set_constraint()

    sg_no_r, sg_sym_r = get_space_group(atoms, symprec)

    # Save relaxed CIF
    relaxed_path = cif_path.replace(".cif", "_symrelaxed.cif")
    write(relaxed_path, atoms)

    return relaxed_path, sg_no, sg_sym, sg_no_r, sg_sym_r, len(atoms)


def compute_i_curv(cif_path, calc, supercell, q_mesh=I_CURV_Q_MESH, delta=0.05):
    """Compute I_curv on a uniform q-mesh using ASE Phonons with mean forces."""
    from ase.io import read
    from ase.phonons import Phonons
    from ase.dft.kpoints import monkhorst_pack

    atoms = read(cif_path)
    atoms.calc = calc

    ph_dir = "phonon_temp_i_curv"
    if os.path.exists(ph_dir):
        shutil.rmtree(ph_dir)

    ph = Phonons(atoms, calc, supercell=supercell, delta=delta, name=ph_dir)
    ph.run()
    ph.read()

    qpoints = monkhorst_pack(q_mesh)
    energies_eV = ph.band_structure(qpoints, modes=False, verbose=False)
    freqs_THz = energies_eV / EV_TO_THZ

    gamma_mask_q = np.all(np.abs(qpoints) < 1e-10, axis=1)
    gamma_indices = np.where(gamma_mask_q)[0]

    ph.clean()

    # Compute I_curv
    freqs = freqs_THz.ravel()
    omega_min = float(min(freqs.min(), 0.0))

    exclude_mask = np.zeros_like(freqs_THz, dtype=bool)
    if len(gamma_indices) > 0:
        for gi in gamma_indices:
            exclude_mask[gi, :3] = True
    near_zero_mask = np.abs(freqs_THz) < ACOUSTIC_THRESHOLD
    exclude_mask = exclude_mask | near_zero_mask

    kept = freqs_THz[~exclude_mask]
    n_modes = len(kept)
    n_excluded = int(exclude_mask.sum())

    if n_modes == 0:
        return {"i_curv": 0.0, "omega_inst": 0.0, "omega_min": omega_min,
                "soft_frac": 0.0, "n_modes": 0, "n_imag": 0, "n_excluded": n_excluded}

    n_imag = int((kept < -ACOUSTIC_THRESHOLD).sum())
    soft_frac = float(n_imag) / n_modes

    neg_curvature = np.where(kept < 0, kept ** 2, 0.0)
    i_curv = float(neg_curvature.mean())
    omega_inst = float(np.sqrt(i_curv))

    return {"i_curv": i_curv, "omega_inst": omega_inst, "omega_min": omega_min,
            "soft_frac": soft_frac, "n_modes": n_modes, "n_imag": n_imag,
            "n_excluded": n_excluded}


def copy_original_phonon_data(cif_path, out_dir, entry):
    """Copy original MACE phonon band structure data and plots alongside the PET UQ results."""
    fname = os.path.basename(cif_path)
    m = re.search(r'energy_(m\d+p\d+)', fname)
    energy_tag = m.group(1) if m else ""

    category = entry.get("category", "")
    if category in ("quaternary_abandoned", "quinary_abandoned"):
        base_dir = os.path.dirname(cif_path)
        base_name = fname.replace(".cif", "")
        phonon_dirs = [d for d in os.listdir(base_dir)
                       if d.startswith(base_name) and "_MACE_GA_phonon_output_" in d]
        if not phonon_dirs:
            return
        phonon_dir = os.path.join(base_dir, phonon_dirs[0])
        src_dir = phonon_dir
    else:
        cif_dir = os.path.dirname(cif_path)
        if "final_structures" in cif_dir:
            phonon_dir = os.path.dirname(cif_dir)
        else:
            phonon_dir = cif_dir

        if not os.path.isdir(phonon_dir):
            return
        analysis_dirs = [d for d in os.listdir(phonon_dir)
                         if d.startswith("final_phonon_analysis_") and energy_tag in d]
        if analysis_dirs:
            src_dir = os.path.join(phonon_dir, analysis_dirs[0])
        else:
            src_dir = phonon_dir

    if not os.path.isdir(src_dir):
        return

    dest_dir = os.path.join(out_dir, "original_mace_phonon_data")
    os.makedirs(dest_dir, exist_ok=True)

    for f in os.listdir(src_dir):
        if any(f.endswith(ext) for ext in (".txt", ".png", ".svg", ".json")):
            if "band_structure" in f or "dos" in f or "phonon_run_summary" in f \
               or "special_point" in f or "k_point" in f or f.endswith(".png") \
               or f.endswith(".svg"):
                try:
                    shutil.copy2(os.path.join(src_dir, f), os.path.join(dest_dir, f))
                except:
                    pass


def parse_summary(path):
    """Parse a PET-MAD UQ summary file."""
    data = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if "Most negative mean frequency:" in line:
                data["min_freq"] = float(line.split(":")[-1].replace("THz", "").strip())
            elif "Max std across all q-points:" in line:
                data["max_std"] = float(line.split(":")[-1].replace("THz", "").strip())
            elif "Std at min-frequency q-point:" in line:
                data["std_at_min"] = float(line.split(":")[-1].replace("THz", "").strip())
            elif "Min freq minus 3*std:" in line:
                data["min_minus_3std"] = float(line.split(":")[-1].replace("THz", "").strip())
            elif "Dynamically stable:" in line:
                data["stable"] = "Yes" in line
            elif "Ensemble members:" in line:
                data["n_ensemble"] = int(line.split(":")[-1].strip())
    return data


def run_single_structure(entry, idx, total, pet_calc, mace_icurv):
    """Run comprehensive PET-MAD UQ + I_curv on a single structure."""
    cif_path = entry["filepath"]
    compound = entry["compound"]
    category = entry.get("category", "unknown")
    sc = entry["supercell"]
    name = entry.get("filename", os.path.basename(cif_path)).replace("_conventional", "").replace(".cif", "")
    prefix = f"{idx:03d}_{compound}_{category[:4]}_{name[:40]}"

    out_dir = os.path.join(OUTPUT_DIR, prefix)
    os.makedirs(out_dir, exist_ok=True)

    # Copy original MACE phonon data
    copy_original_phonon_data(cif_path, out_dir, entry)

    mace_omega = entry.get("mace_omega_min")
    mace_str = f"ω={mace_omega:.3f}" if mace_omega is not None else "ω=?"
    print(f"\n[{idx}/{total}] {compound} ({category}) {mace_str} THz, "
          f"{entry.get('natoms', '?')} atoms, sc={sc}x{sc}x{sc}")
    print(f"  Output: {out_dir}")

    t_start = time.time()

    # --- Step 1: Symmetry-constrained relaxation ---
    print(f"  Step 1: Symmetry-constrained relaxation with PET-MAD...")
    try:
        relaxed_path, sg_no, sg_sym, sg_no_r, sg_sym_r, natoms = symmetry_constrained_relax(
            cif_path, pet_calc, fmax=1e-3
        )
    except Exception as e:
        print(f"  RELAXATION FAILED: {e}")
        import traceback
        traceback.print_exc()
        return None

    # --- Step 2: Call vibroml.main for UQ phonon plots/npz/summary ---
    print(f"  Step 2: Running vibroml.main --phonon-uq (plots + UQ stats)...")
    cmd = [
        VIBROML_PYTHON, "-m", VIBROML_MODULE,
        "--cif", relaxed_path,
        "--engine", "pet",
        "--model_name", "pet-mad-xs",
        "--supercell", str(sc),
        "--delta", "0.05",
        "--fmax", "1e-4",
        "--units", "THz",
        "--phonon-uq",
        "--no-relax",  # we already relaxed with symmetry constraint
        "--relaxation-patience", "15",
    ]

    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=3600,
            cwd=out_dir,
            env={**os.environ, "PYTHONUNBUFFERED": "1"}
        )

        if result.returncode != 0:
            print(f"  VIBROML FAILED (exit {result.returncode})")
            print(f"  stderr: {result.stderr[:500]}")
            return None

    except subprocess.TimeoutExpired:
        print(f"  VIBROML TIMEOUT")
        return None
    except Exception as e:
        print(f"  VIBROML ERROR: {e}")
        return None

    # Find vibroml output directory
    vibroml_output_dirs = [d for d in os.listdir(out_dir)
                           if "_PET_phonon_output_" in d]
    if not vibroml_output_dirs:
        print(f"  No PET phonon output dir found in {out_dir}")
        return None

    vibroml_dir = os.path.join(out_dir, vibroml_output_dirs[0])

    # Parse UQ summary
    summaries = [f for f in os.listdir(vibroml_dir) if f.endswith("_phonon_uq_summary.txt")]
    if not summaries:
        print(f"  No summary file found in {vibroml_dir}")
        return None

    summary_path = os.path.join(vibroml_dir, summaries[0])
    summary = parse_summary(summary_path)

    # Parse npz for soft_frac
    npz_files = [f for f in os.listdir(vibroml_dir) if f.endswith("_phonon_uq_data.npz")]
    pet_soft_frac = None
    pet_soft_frac_15 = None
    if npz_files:
        try:
            npz = np.load(os.path.join(vibroml_dir, npz_files[0]), allow_pickle=True)
            mean_freqs = np.concatenate([seg.ravel() for seg in npz['mean_freqs']])
            pet_soft_frac = float((mean_freqs < -0.1).sum() / len(mean_freqs))
            pet_soft_frac_15 = float((mean_freqs < -1.5).sum() / len(mean_freqs))
        except Exception as exc:
            print(f"  Warning: could not compute soft frac from npz ({exc})")

    # --- Step 3: Compute I_curv on uniform q-mesh ---
    print(f"  Step 3: Computing I_curv on uniform {I_CURV_Q_MESH} mesh...")
    try:
        if isinstance(sc, int):
            sc_tuple = (sc, sc, sc)
        else:
            sc_tuple = tuple(sc)
        i_curv_result = compute_i_curv(relaxed_path, pet_calc, sc_tuple)
    except Exception as e:
        print(f"  I_CURV COMPUTATION FAILED: {e}")
        import traceback
        traceback.print_exc()
        i_curv_result = {"i_curv": None, "omega_inst": None, "omega_min": None,
                         "soft_frac": None, "n_modes": 0, "n_imag": 0, "n_excluded": 0}

    elapsed = time.time() - t_start

    # Save I_curv to a JSON in the vibroml output dir
    i_curv_path = os.path.join(vibroml_dir, "i_curv_metric.json")
    with open(i_curv_path, "w") as f:
        json.dump(i_curv_result, f, indent=2)

    print(f"  DONE in {elapsed:.1f}s | min_freq={summary.get('min_freq', '?')} THz, "
          f"I_curv={i_curv_result.get('i_curv', '?')} THz², "
          f"max_std={summary.get('max_std', '?')} THz")

    return {
        "idx": idx,
        "category": category,
        "compound": compound,
        "space_group": entry.get("space_group", ""),
        "space_group_no": sg_no,
        "space_group_sym": sg_sym,
        "space_group_no_relaxed": sg_no_r,
        "space_group_sym_relaxed": sg_sym_r,
        "natoms": natoms,
        "supercell": sc,
        # MACE stats (from manifest)
        "mace_omega_min": entry.get("mace_omega_min"),
        "mace_soft_frac": entry.get("mace_soft_frac"),
        "mace_soft_frac_below_1p5": entry.get("mace_soft_frac_below_1p5"),
        # MACE I_curv (from pre-computed uniform q-mesh)
        "mace_i_curv": mace_icurv.get("mace_i_curv"),
        "mace_omega_inst": mace_icurv.get("mace_omega_inst"),
        "mace_i_curv_omega_min": mace_icurv.get("mace_i_curv_omega_min"),
        "mace_i_curv_soft_frac": mace_icurv.get("mace_i_curv_soft_frac"),
        "mace_i_curv_n_imag": mace_icurv.get("mace_i_curv_n_imag"),
        "mace_i_curv_n_modes": mace_icurv.get("mace_i_curv_n_modes"),
        # PET UQ stats (from vibroml summary)
        "pet_omega_min": summary.get("min_freq"),
        "pet_soft_frac": pet_soft_frac,
        "pet_soft_frac_below_1p5": pet_soft_frac_15,
        "pet_max_std_thz": summary.get("max_std"),
        "pet_std_at_min_thz": summary.get("std_at_min"),
        "pet_min_minus_3std": summary.get("min_minus_3std"),
        "pet_n_ensemble": summary.get("n_ensemble"),
        # I_curv (new metric, from uniform q-mesh)
        "pet_i_curv": i_curv_result.get("i_curv"),
        "pet_omega_inst": i_curv_result.get("omega_inst"),
        "pet_i_curv_omega_min": i_curv_result.get("omega_min"),
        "pet_i_curv_soft_frac": i_curv_result.get("soft_frac"),
        "pet_i_curv_n_imag": i_curv_result.get("n_imag"),
        "pet_i_curv_n_modes": i_curv_result.get("n_modes"),
        # Metadata
        "time_s": elapsed,
        "output_dir": vibroml_dir,
        "summary_file": summary_path,
        "npz_file": os.path.join(vibroml_dir, npz_files[0]) if npz_files else None,
        "i_curv_file": i_curv_path,
        "relaxed_cif": relaxed_path,
    }


def load_completed_results():
    """Load previously completed results from summary CSV."""
    completed = {}
    if os.path.exists(SUMMARY_CSV):
        with open(SUMMARY_CSV) as f:
            reader = csv.DictReader(f)
            for row in reader:
                try:
                    completed[int(row["idx"])] = row
                except (ValueError, KeyError):
                    pass
    return completed


def load_mace_icurv():
    """Load MACE I_curv results from the pre-computed JSON.
    Returns a dict keyed by idx (1-based)."""
    mace_icurv = {}
    if os.path.exists(MACE_ICURV_JSON):
        with open(MACE_ICURV_JSON) as f:
            data = json.load(f)
        for r in data:
            # idx in the JSON is 0-based, convert to 1-based to match manifest
            idx_1based = r.get("idx", -1) + 1
            mace_icurv[idx_1based] = {
                "mace_i_curv": r.get("i_curv"),
                "mace_omega_inst": r.get("omega_inst"),
                "mace_i_curv_omega_min": r.get("omega_min"),
                "mace_i_curv_soft_frac": r.get("soft_frac"),
                "mace_i_curv_n_imag": r.get("n_imag"),
                "mace_i_curv_n_modes": r.get("n_modes"),
            }
        print(f"  Loaded MACE I_curv for {len(mace_icurv)} structures from {MACE_ICURV_JSON}")
    else:
        print(f"  Warning: MACE I_curv JSON not found at {MACE_ICURV_JSON}")
    return mace_icurv


def save_summary_csv(results):
    """Save results to summary CSV (called incrementally and at end)."""
    fields = [
        "idx", "category", "compound",
        "space_group", "space_group_no", "space_group_sym",
        "space_group_no_relaxed", "space_group_sym_relaxed",
        "natoms", "supercell",
        # MACE stats (from manifest)
        "mace_omega_min", "mace_soft_frac", "mace_soft_frac_below_1p5",
        # MACE I_curv (from pre-computed uniform q-mesh)
        "mace_i_curv", "mace_omega_inst", "mace_i_curv_omega_min",
        "mace_i_curv_soft_frac", "mace_i_curv_n_imag", "mace_i_curv_n_modes",
        # PET UQ stats (from vibroml)
        "pet_omega_min", "pet_soft_frac", "pet_soft_frac_below_1p5",
        "pet_max_std_thz", "pet_std_at_min_thz", "pet_min_minus_3std",
        "pet_n_ensemble",
        # PET I_curv (from uniform q-mesh)
        "pet_i_curv", "pet_omega_inst", "pet_i_curv_omega_min",
        "pet_i_curv_soft_frac", "pet_i_curv_n_imag", "pet_i_curv_n_modes",
        # Metadata
        "time_s", "output_dir", "relaxed_cif",
    ]
    with open(SUMMARY_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for r in sorted(results, key=lambda x: int(x.get("idx", 0))):
            writer.writerow(r)
    print(f"  Summary saved: {SUMMARY_CSV} ({len(results)} entries)")


def main():
    with open(STRUCTURES_JSON) as f:
        structures = json.load(f)

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    dry_run = "--dry-run" in sys.argv
    if dry_run:
        print("DRY RUN - would execute:\n")
        for i, entry in enumerate(structures, 1):
            omega = entry.get("mace_omega_min")
            sf = entry.get("mace_soft_frac")
            omega_str = f"{omega:.3f}" if omega is not None else "?"
            sf_str = f"{100*sf:.1f}%" if sf is not None else "?"
            print(f"  [{i:3d}/{len(structures)}] {entry['category']:>20s} {entry['compound']:>15s} "
                  f"ω_min={omega_str:>8s} THz, sf={sf_str:>6s}, "
                  f"{entry.get('natoms','?'):>3} atoms, sc={entry['supercell']}")
        print(f"\nTotal: {len(structures)} structures")
        return

    # Load previously completed results for resume
    completed = load_completed_results()
    if completed:
        print(f"Resuming: {len(completed)} structures already completed, "
              f"{len(structures) - len(completed)} remaining")

    # Load MACE I_curv results (pre-computed on uniform q-mesh)
    print("\nLoading MACE I_curv results...")
    mace_icurv_all = load_mace_icurv()

    results = []
    total = len(structures)

    # Add completed results to the list
    for idx, row in completed.items():
        entry = structures[idx - 1]
        row.setdefault("category", entry.get("category", ""))
        row.setdefault("compound", entry["compound"])
        row.setdefault("space_group", entry.get("space_group", ""))
        row.setdefault("natoms", entry.get("natoms", 0))
        row.setdefault("supercell", entry["supercell"])
        row.setdefault("mace_omega_min", entry.get("mace_omega_min"))
        row.setdefault("mace_soft_frac", entry.get("mace_soft_frac"))
        row.setdefault("mace_soft_frac_below_1p5", entry.get("mace_soft_frac_below_1p5"))
        # Merge MACE I_curv for completed rows too
        mace_icurv = mace_icurv_all.get(idx, {})
        row.setdefault("mace_i_curv", mace_icurv.get("mace_i_curv"))
        row.setdefault("mace_omega_inst", mace_icurv.get("mace_omega_inst"))
        row.setdefault("mace_i_curv_omega_min", mace_icurv.get("mace_i_curv_omega_min"))
        row.setdefault("mace_i_curv_soft_frac", mace_icurv.get("mace_i_curv_soft_frac"))
        row.setdefault("mace_i_curv_n_imag", mace_icurv.get("mace_i_curv_n_imag"))
        row.setdefault("mace_i_curv_n_modes", mace_icurv.get("mace_i_curv_n_modes"))
        row.setdefault("time_s", 0)
        results.append(row)

    # Load PET-MAD calculator once
    print("\nLoading PET-MAD calculator...")
    from upet.calculator import UPETCalculator
    pet_calc = UPETCalculator(model="pet-mad-xs", device="cuda")

    for i, entry in enumerate(structures, 1):
        if i in completed:
            continue

        # Skip large structures (>150 atoms) on debug-gpu to avoid timeouts
        # These will be handled separately on the gpu partition
        natoms = entry.get('natoms', 0)
        if natoms > 150 and os.environ.get('VIBROML_SKIP_LARGE', '0') == '1':
            print(f"\n[{i}/{total}] {entry['compound']} ({entry.get('category','?')}) "
                  f"SKIPPED ({natoms} atoms > 150, VIBROML_SKIP_LARGE=1)")
            continue

        # Skip small structures (<=150 atoms) when running large-only batch
        if natoms <= 150 and os.environ.get('VIBROML_SKIP_SMALL', '0') == '1':
            continue

        mace_icurv = mace_icurv_all.get(i, {})
        result = run_single_structure(entry, i, total, pet_calc, mace_icurv)
        if result:
            results.append(result)
            save_summary_csv(results)

    # Write final summary CSV
    if results:
        save_summary_csv(results)

        # Print summary table
        print(f"\n{'idx':>3} {'category':>20} {'compound':>15} "
              f"{'MACE_Icurv':>10} {'PET_Icurv':>10} "
              f"{'MACE_ωmin':>9} {'PET_ωmin':>9} {'σ@min':>7} {'σ_max':>7}")
        for r in sorted(results, key=lambda x: int(x.get("idx", 0))):
            mace_ic = f"{r.get('mace_i_curv',0):.4f}" if r.get('mace_i_curv') is not None else "?"
            pet_ic = f"{r.get('pet_i_curv',0):.4f}" if r.get('pet_i_curv') is not None else "?"
            mace_w = f"{r.get('mace_omega_min',0):.3f}" if r.get('mace_omega_min') is not None else "?"
            pet_w = f"{r.get('pet_omega_min',0):.3f}" if r.get('pet_omega_min') is not None else "?"
            sig = f"{r.get('pet_std_at_min_thz',0):.3f}" if r.get('pet_std_at_min_thz') is not None else "?"
            smax = f"{r.get('pet_max_std_thz',0):.3f}" if r.get('pet_max_std_thz') is not None else "?"
            cat = r.get('category', '')[:20]
            print(f"{r.get('idx',0):>3} {cat:>20} {r['compound']:>15} "
                  f"{mace_ic:>10} {pet_ic:>10} "
                  f"{mace_w:>9} {pet_w:>9} {sig:>7} {smax:>7}")
    else:
        print("\nNo successful results!")


if __name__ == "__main__":
    main()
