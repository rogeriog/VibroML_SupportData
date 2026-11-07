#!/usr/bin/env python3
"""Scan every GA structure in the explorer DB for steric clashes.

A structure is flagged when its minimum interatomic distance (periodic
images included) falls below CLASH_RATIO x (sum of covalent radii) or below
CLASH_ABS_ANGSTROM. Flags are written to steric_clash_flags.csv with one row
per scanned structure, keyed by ga_structures.id.

Usage: python scan_steric_clashes.py [--workers N] [--limit N]
"""

import argparse
import os
import sqlite3
import sys
import time
from multiprocessing import Pool

import numpy as np
from ase.data import atomic_numbers, covalent_radii
from ase.io import read

BASE_DIR = "/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/02_runs_and_data/candidates_lowcoverage_cif"
DB_PATH = "/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/02_runs_and_data/VibroML_SupportData/whitespace_data/vibroml_database.db"
OUT_CSV = "/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/02_runs_and_data/VibroML_SupportData/whitespace_data/steric_clash_flags.csv"

CLASH_RATIO = 0.5
CLASH_ABS_ANGSTROM = 1.0


def check_structure(row):
    """row = (id, cif_path, material_name, struct_type, rank, energy)."""
    sid, rel_path, name, stype, rank, energy = row
    full = os.path.join(BASE_DIR, rel_path) if rel_path else None
    if not full or not os.path.exists(full):
        return (sid, name, stype, rank, energy, None, None, None, "missing")
    try:
        atoms = read(full)
        d = atoms.get_all_distances(mic=True)
        np.fill_diagonal(d, 1e9)
        mind = float(d.min())
        i, j = np.unravel_index(d.argmin(), d.shape)
        s1, s2 = atoms[i].symbol, atoms[j].symbol
        ratio = mind / (covalent_radii[atomic_numbers[s1]] + covalent_radii[atomic_numbers[s2]])
        clashed = (ratio < CLASH_RATIO) or (mind < CLASH_ABS_ANGSTROM)
        return (sid, name, stype, rank, energy, round(mind, 4), round(ratio, 4),
                f"{s1}-{s2}", "clashed" if clashed else "ok")
    except Exception:
        return (sid, name, stype, rank, energy, None, None, None, "error")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    con = sqlite3.connect(DB_PATH)
    rows = con.execute(
        """SELECT g.id, g.cif_path, m.material_name, g.struct_type, g.rank,
                  g.energy_per_atom_eV
           FROM ga_structures g JOIN materials m ON m.id = g.material_id
           WHERE g.cif_path IS NOT NULL"""
    ).fetchall()
    con.close()
    if args.limit:
        rows = rows[: args.limit]
    print(f"scanning {len(rows)} structures with {args.workers} workers")

    t0 = time.time()
    n_written = 0
    with open(OUT_CSV, "w") as out, Pool(args.workers) as pool:
        out.write("ga_structure_id,material_name,struct_type,rank,energy_per_atom_eV,"
                  "min_distance_angstrom,min_distance_over_covalent,closest_pair,status\n")
        for res in pool.imap_unordered(check_structure, rows, chunksize=200):
            out.write(",".join(str(x) if x is not None else "" for x in res) + "\n")
            n_written += 1
            if n_written % 20000 == 0:
                rate = n_written / (time.time() - t0)
                eta = (len(rows) - n_written) / rate / 60
                print(f"  {n_written}/{len(rows)}  ({rate:.0f}/s, ETA {eta:.1f} min)", flush=True)
    print(f"done in {(time.time()-t0)/60:.1f} min -> {OUT_CSV}")


if __name__ == "__main__":
    main()
