"""Parse final-step DFT forces from the white-space verification OUTCARs.

Reads `dft_whitespace_results.csv` to enumerate sampled systems, locates
each run's OUTCAR under a VASP inputs tree, extracts the final TOTAL-FORCE
block and POSCAR species, and writes:

- ``dft_whitespace_forces.csv.gz`` : one row per atom
  (material, side, atom_index, element, fx, fy, fz in eV/Angstrom)
- ``dft_whitespace_force_stats.csv`` : one row per structure
  (material, side, n_atoms, f_max, f_rms, e0)

Usage:
    python parse_dft_forces.py --results dft_whitespace_results.csv \
        --vasp-tree /path/to/dft_inputs/whitespace_winners \
        --out-dir .
"""
import argparse
import glob
import os
import re

import numpy as np
import pandas as pd


def parse_outcar_forces(path):
    """Return (forces Nx3 ndarray, e0 float or None) from the last ionic step."""
    forces, e0 = None, None
    block = []
    in_block = False
    with open(path, errors="replace") as fh:
        for line in fh:
            if "TOTAL-FORCE" in line:
                in_block, block = True, []
                continue
            if in_block:
                if line.strip().startswith("---") and len(block) == 0:
                    continue  # separator right after the header
                if re.match(r"^\s*-{10,}\s*$", line):
                    if block:
                        in_block = False
                        forces = np.array(block, dtype=float)
                    continue
                cols = line.split()
                if len(cols) >= 6:
                    try:
                        block.append([float(c) for c in cols[3:6]])
                    except ValueError:
                        pass
            if "free  energy   TOTEN" in line or re.search(r"E0=\s*[-\d.E+]+", line):
                m = re.search(r"E0=\s*([-\d.E+]+)", line)
                if m:
                    e0 = float(m.group(1))
    return forces, e0


def parse_poscar_species(path):
    """Return list of element symbols per site (VASP5 POSCAR)."""
    with open(path) as fh:
        lines = [l.strip() for l in fh.readlines() if l.strip()]
    symbols = lines[5].split()
    counts = [int(c) for c in lines[6].split()]
    return [s for s, c in zip(symbols, counts) for _ in range(c)]


def find_run_dirs(vasp_tree, run_dir):
    """Return {'initial': path, 'winner': path} by prefix-matching dir names.

    Some run directory names were truncated at filesystem limits, so match
    on prefix rather than exact equality.
    """
    hits = {}
    for side in ("initial", "winner"):
        cand = sorted(glob.glob(os.path.join(vasp_tree, f"{run_dir}*_{side}*"),
                               recursive=False))
        cand += sorted(glob.glob(os.path.join(vasp_tree, f"{run_dir}_{side[:3]}*")))
        for d in cand:
            base = os.path.basename(d)
            if os.path.isdir(d) and (base.endswith(side) or base.endswith(side[:3])):
                hits[side] = d
                break
    return hits


_SIDE_SUFFIX = {"initial": ("initial", "init", "ini", "in"),
                "winner": ("winner", "winn", "win", "wi")}


def find_spilled_outcar(vasp_tree, run_dir, side):
    """VASP 5.4.4 truncates long output paths: for the ~110-char run dirs
    the OUTCAR is written one level up under a truncated name. Return the
    best spilled-file candidate whose name is a prefix of run_dir_<side>.
    """
    best = None
    for cand in os.listdir(vasp_tree):
        p = os.path.join(vasp_tree, cand)
        if not os.path.isfile(p):
            continue
        base = f"{run_dir}_{side}"
        if not base.startswith(cand):
            continue
        suffix = cand[len(run_dir) + 1:] if cand.startswith(run_dir + "_") else cand
        if not any(suffix.startswith(s) for s in _SIDE_SUFFIX[side]):
            continue
        if best is None or len(cand) > len(os.path.basename(best)):
            best = p
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="dft_whitespace_results.csv")
    ap.add_argument("--vasp-tree", required=True,
                    help="Directory holding <run_dir>_{initial,winner}/ OUTCAR+POSCAR")
    ap.add_argument("--out-dir", default=".")
    args = ap.parse_args()

    res = pd.read_csv(args.results)
    atoms_rows, stats_rows = [], []
    missing = []
    for _, row in res.iterrows():
        rd, mat = row["run_dir"], row["material"]
        dirs = find_run_dirs(args.vasp_tree, rd)
        for side in ("initial", "winner"):
            d = dirs.get(side)
            oc = os.path.join(d, "OUTCAR") if d else None
            forces = e0 = None
            if oc and os.path.exists(oc):
                forces, e0 = parse_outcar_forces(oc)
            if forces is None or len(forces) == 0:
                spilled = find_spilled_outcar(args.vasp_tree, rd, side)
                if spilled:
                    forces, e0 = parse_outcar_forces(spilled)
            if (forces is None or len(forces) == 0) or not d:
                missing.append((rd, side))
                continue
            species = parse_poscar_species(os.path.join(d, "POSCAR"))
            if forces is None or len(forces) == 0:
                missing.append((rd, side))
                continue
            norms = np.linalg.norm(forces, axis=1)
            stats_rows.append(dict(material=mat, run_dir=rd, side=side,
                                   n_atoms=len(forces), e0=e0,
                                   f_max=float(norms.max()),
                                   f_rms=float(np.sqrt((norms ** 2).mean()))))
            for i, (el, f) in enumerate(zip(species, forces)):
                atoms_rows.append(dict(material=mat, side=side, atom_index=i,
                                       element=el, fx=f[0], fy=f[1], fz=f[2]))
    pd.DataFrame(atoms_rows).to_csv(
        os.path.join(args.out_dir, "dft_whitespace_forces.csv.gz"),
        index=False, compression="gzip")
    stats = pd.DataFrame(stats_rows)
    stats.to_csv(os.path.join(args.out_dir, "dft_whitespace_force_stats.csv"),
                 index=False)
    print(f"wrote {len(atoms_rows)} atom rows, {len(stats)} structure rows; "
          f"{len(missing)} missing/incomplete")


if __name__ == "__main__":
    main()
