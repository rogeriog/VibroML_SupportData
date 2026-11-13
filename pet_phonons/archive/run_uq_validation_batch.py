#!/usr/bin/env python3
"""
Batch PET-MAD UQ validation for VibroML stability thresholds.

Runs phonon UQ calculations on 105 structures:
  - 43 benchmark polymorphs from TableS1 (MACE-OMAT)
  - 2 ProtoCSP final validated structures
  - 50 randomly selected quaternary abandoned structures
  - 10 randomly selected quinary abandoned structures

Then collects results into a summary CSV for figure generation.

Resume logic: detects completed structures from existing output dirs and CSV,
skipping them on re-submission. This allows chunked execution across multiple
debug-gpu submissions.

Usage:
  python run_uq_validation_batch.py [--dry-run]
"""

import os, sys, json, subprocess, time, csv
import numpy as np
from pathlib import Path

STRUCTURES_JSON = os.path.join(os.path.dirname(__file__), "uq_validation_structures.json")
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "uq_validation_runs")
SUMMARY_CSV = os.path.join(os.path.dirname(__file__), "uq_validation_summary.csv")

# Use orb_env Python (has upet/metatrain); VibroML is on PYTHONPATH
VIBROML_PYTHON = os.environ.get("VIBROML_PYTHON", sys.executable)
VIBROML_MODULE = "vibroml.main"


def copy_original_phonon_data(cif_path, out_dir, entry):
    """Copy original MACE phonon band structure data and plots alongside the PET UQ results."""
    import shutil, re

    # For benchmark structures, the CIF is inside a phonon_output dir
    # For quaternary/quinary abandoned, the CIF is next to the phonon_output dir
    fname = os.path.basename(cif_path)
    m = re.search(r'energy_(m\d+p\d+)', fname)
    energy_tag = m.group(1) if m else ""

    # Determine the phonon output directory based on category
    category = entry.get("category", "")
    if category in ("quaternary_abandoned", "quinary_abandoned"):
        # CIF is in the base dir, phonon output is a sibling dir
        base_dir = os.path.dirname(cif_path)
        base_name = fname.replace(".cif", "")
        phonon_dirs = [d for d in os.listdir(base_dir)
                       if d.startswith(base_name) and "_MACE_GA_phonon_output_" in d]
        if not phonon_dirs:
            return
        phonon_dir = os.path.join(base_dir, phonon_dirs[0])
        src_dir = phonon_dir  # band structure files are directly in the phonon dir
    else:
        # Benchmark: CIF is in final_structures or final_structures_recovered
        cif_dir = os.path.dirname(cif_path)
        if "final_structures" in cif_dir:
            phonon_dir = os.path.dirname(cif_dir)
        else:
            phonon_dir = cif_dir

        # Find matching final_phonon_analysis dir
        if not os.path.isdir(phonon_dir):
            return
        analysis_dirs = [d for d in os.listdir(phonon_dir)
                         if d.startswith("final_phonon_analysis_") and energy_tag in d]
        if analysis_dirs:
            src_dir = os.path.join(phonon_dir, analysis_dirs[0])
        else:
            src_dir = phonon_dir  # fallback: band structure files directly in phonon dir

    if not os.path.isdir(src_dir):
        return

    dest_dir = os.path.join(out_dir, "original_mace_phonon_data")
    os.makedirs(dest_dir, exist_ok=True)

    # Copy band structure data, plots, and summary
    for f in os.listdir(src_dir):
        if any(f.endswith(ext) for ext in (".txt", ".png", ".svg", ".json")):
            if "band_structure" in f or "dos" in f or "phonon_run_summary" in f \
               or "special_point" in f or "k_point" in f or f.endswith(".png") \
               or f.endswith(".svg"):
                try:
                    shutil.copy2(os.path.join(src_dir, f), os.path.join(dest_dir, f))
                except:
                    pass


def run_single_structure(entry, idx, total):
    """Run PET-MAD UQ on a single structure."""
    cif_path = entry["filepath"]
    compound = entry["compound"]
    category = entry.get("category", "unknown")
    sc = entry["supercell"]
    name = entry["filename"].replace("_conventional", "").replace(".cif", "")
    # Short prefix for output
    prefix = f"{idx:03d}_{compound}_{category[:4]}_{name[:40]}"

    out_dir = os.path.join(OUTPUT_DIR, prefix)
    os.makedirs(out_dir, exist_ok=True)

    # Copy original MACE phonon data (band structure, plots, summary) into out_dir
    copy_original_phonon_data(cif_path, out_dir, entry)

    cmd = [
        VIBROML_PYTHON, "-m", VIBROML_MODULE,
        "--cif", cif_path,
        "--engine", "pet",
        "--model_name", "pet-mad-xs",
        "--supercell", str(sc),
        "--delta", "0.05",
        "--fmax", "1e-4",
        "--units", "THz",
        "--phonon-uq",
        "--relaxation-patience", "15",
    ]

    mace_omega = entry.get("mace_omega_min")
    mace_str = f"ω={mace_omega:.3f}" if mace_omega is not None else "ω=?"
    print(f"\n[{idx}/{total}] {compound} ({category}) {mace_str} THz, "
          f"{entry['natoms']} atoms, sc={sc}x{sc}x{sc}")
    print(f"  Output: {out_dir}")

    start = time.time()
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=3600,
            cwd=out_dir,  # run from the output dir so vibroml creates its folder here
            env={**os.environ, "PYTHONUNBUFFERED": "1"}
        )
        elapsed = time.time() - start

        if result.returncode != 0:
            print(f"  FAILED (exit {result.returncode}) in {elapsed:.1f}s")
            print(f"  stderr: {result.stderr[:500]}")
            return None

        # vibroml creates a *_PET_phonon_output_* subdirectory in CWD (out_dir)
        # Find it and locate the summary file inside
        vibroml_output_dirs = [d for d in os.listdir(out_dir)
                               if d.endswith("_PET_phonon_output_") or "_PET_phonon_output_" in d]
        if not vibroml_output_dirs:
            print(f"  No PET phonon output dir found in {out_dir}")
            return None

        vibroml_dir = os.path.join(out_dir, vibroml_output_dirs[0])
        summaries = [f for f in os.listdir(vibroml_dir) if f.endswith("_phonon_uq_summary.txt")]
        if not summaries:
            print(f"  No summary file found in {vibroml_dir}")
            return None

        summary_path = os.path.join(vibroml_dir, summaries[0])
        summary = parse_summary(summary_path)

        # Also check for npz data and compute soft mode fraction from PET mean
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

        if pet_soft_frac is not None:
            print(f"  DONE in {elapsed:.1f}s | min_freq={summary.get('min_freq', '?')} THz, "
                  f"soft_frac={100*pet_soft_frac:.1f}%")
        else:
            print(f"  DONE in {elapsed:.1f}s | min_freq={summary.get('min_freq', '?')} THz")

        return {
            "idx": idx,
            "category": category,
            "compound": compound,
            "space_group": entry.get("space_group", ""),
            "natoms": entry["natoms"],
            "supercell": sc,
            # MACE stats (from corrected manifest)
            "mace_omega_min": entry.get("mace_omega_min"),
            "mace_soft_frac": entry.get("mace_soft_frac"),
            "mace_soft_frac_below_1p5": entry.get("mace_soft_frac_below_1p5"),
            # PET stats (from UQ ensemble mean)
            "pet_omega_min": summary.get("min_freq"),
            "pet_soft_frac": pet_soft_frac,
            "pet_soft_frac_below_1p5": pet_soft_frac_15,
            "pet_max_std_thz": summary.get("max_std"),
            "pet_std_at_min_thz": summary.get("std_at_min"),
            "pet_min_minus_3std": summary.get("min_minus_3std"),
            "time_s": elapsed,
            "output_dir": vibroml_dir,
            "summary_file": summary_path,
            "npz_file": os.path.join(vibroml_dir, npz_files[0]) if npz_files else None,
        }

    except subprocess.TimeoutExpired:
        print(f"  TIMEOUT after {time.time()-start:.1f}s")
        return None
    except Exception as e:
        print(f"  ERROR: {e}")
        return None


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
    return data


def load_completed_results():
    """Load previously completed results from summary CSV and output dirs."""
    completed = {}
    if os.path.exists(SUMMARY_CSV):
        with open(SUMMARY_CSV) as f:
            reader = csv.DictReader(f)
            for row in reader:
                try:
                    completed[int(row["idx"])] = row
                except (ValueError, KeyError):
                    pass

    # Also check output dirs for summary files (in case CSV wasn't written)
    if os.path.exists(OUTPUT_DIR):
        for d in os.listdir(OUTPUT_DIR):
            dpath = os.path.join(OUTPUT_DIR, d)
            if not os.path.isdir(dpath):
                continue
            # vibroml creates a *_PET_phonon_output_* subdirectory; summary is inside
            vibroml_subdirs = [sd for sd in os.listdir(dpath)
                               if "_PET_phonon_output_" in sd]
            search_dir = os.path.join(dpath, vibroml_subdirs[0]) if vibroml_subdirs else dpath
            if not os.path.isdir(search_dir):
                continue
            summaries = [f for f in os.listdir(search_dir) if f.endswith("_phonon_uq_summary.txt")]
            if summaries:
                try:
                    idx = int(d.split("_")[0])
                    if idx not in completed:
                        summary = parse_summary(os.path.join(search_dir, summaries[0]))
                        # Also extract soft frac from npz
                        pet_sf = None
                        pet_sf_15 = None
                        npz_files = [f for f in os.listdir(search_dir) if f.endswith("_phonon_uq_data.npz")]
                        if npz_files:
                            try:
                                npz = np.load(os.path.join(search_dir, npz_files[0]), allow_pickle=True)
                                mf = np.concatenate([seg.ravel() for seg in npz['mean_freqs']])
                                pet_sf = float((mf < -0.1).sum() / len(mf))
                                pet_sf_15 = float((mf < -1.5).sum() / len(mf))
                            except:
                                pass
                        completed[idx] = {
                            "idx": idx,
                            "compound": d.split("_")[1] if "_" in d else "",
                            "pet_omega_min": summary.get("min_freq"),
                            "pet_soft_frac": pet_sf,
                            "pet_soft_frac_below_1p5": pet_sf_15,
                            "pet_max_std_thz": summary.get("max_std"),
                            "pet_std_at_min_thz": summary.get("std_at_min"),
                            "pet_min_minus_3std": summary.get("min_minus_3std"),
                        }
                except (ValueError, IndexError):
                    pass
    return completed


def save_summary_csv(results):
    """Save results to summary CSV (called incrementally and at end)."""
    fields = ["idx", "category", "compound", "space_group", "natoms", "supercell",
              "mace_omega_min", "mace_soft_frac", "mace_soft_frac_below_1p5",
              "pet_omega_min", "pet_soft_frac", "pet_soft_frac_below_1p5",
              "pet_max_std_thz", "pet_std_at_min_thz", "pet_min_minus_3std",
              "time_s"]
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
                  f"{entry['natoms']:>3d} atoms, sc={entry['supercell']}x{entry['supercell']}x{entry['supercell']}")
        print(f"\nTotal: {len(structures)} structures")
        return

    # Load previously completed results for resume
    completed = load_completed_results()
    if completed:
        print(f"Resuming: {len(completed)} structures already completed, "
              f"{len(structures) - len(completed)} remaining")

    results = []
    total = len(structures)

    # Add completed results to the list
    for idx, row in completed.items():
        entry = structures[idx - 1]
        row.setdefault("category", entry.get("category", ""))
        row.setdefault("compound", entry["compound"])
        row.setdefault("space_group", entry.get("space_group", ""))
        row.setdefault("natoms", entry["natoms"])
        row.setdefault("supercell", entry["supercell"])
        row.setdefault("mace_omega_min", entry.get("mace_omega_min"))
        row.setdefault("mace_soft_frac", entry.get("mace_soft_frac"))
        row.setdefault("mace_soft_frac_below_1p5", entry.get("mace_soft_frac_below_1p5"))
        row.setdefault("time_s", 0)
        results.append(row)

    for i, entry in enumerate(structures, 1):
        if i in completed:
            continue

        result = run_single_structure(entry, i, total)
        if result:
            results.append(result)
            save_summary_csv(results)

    # Write final summary CSV
    if results:
        save_summary_csv(results)

        # Print summary table
        print(f"\n{'idx':>3} {'category':>20} {'compound':>15} {'MACE_ωmin':>9} {'MACE_sf%':>8} {'PET_ωmin':>9} {'PET_sf%':>8} {'σ@min':>7}")
        for r in sorted(results, key=lambda x: x.get("idx", 0)):
            mace_w = f"{r['mace_omega_min']:.3f}" if r.get('mace_omega_min') is not None else "?"
            mace_sf = f"{100*r['mace_soft_frac']:.1f}" if r.get('mace_soft_frac') is not None else "?"
            pet_w = f"{r['pet_omega_min']:.3f}" if r.get('pet_omega_min') is not None else "?"
            pet_sf = f"{100*r['pet_soft_frac']:.1f}" if r.get('pet_soft_frac') is not None else "?"
            sig = f"{r['pet_std_at_min_thz']:.3f}" if r.get('pet_std_at_min_thz') is not None else "?"
            cat = r.get('category', '')[:20]
            print(f"{r.get('idx',0):>3} {cat:>20} {r['compound']:>15} {mace_w:>9} {mace_sf:>8} {pet_w:>9} {pet_sf:>8} {sig:>7}")
    else:
        print("\nNo successful results!")


if __name__ == "__main__":
    main()
