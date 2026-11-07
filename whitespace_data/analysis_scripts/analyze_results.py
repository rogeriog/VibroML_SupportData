"""Analyze the white-space DFT verification results (GMR item B).

Pairs each system's _initial and _winner single points, extracts E0 from
OSZICAR, and compares the DFT stabilization sign/magnitude against the
MACE prediction stored in sample_manifest.json.

Outputs dft_whitespace_results.csv and prints aggregate statistics:
fraction of winners DFT-confirmed below the initial structure, MAE of the
DFT-vs-MACE energy reduction, and per-dataset breakdown.
"""
import os
import re
import json
import csv
import statistics

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "dft_inputs", "whitespace_winners")
MANIFEST = json.load(open(os.path.join(HERE, "sample_manifest.json")))


def e0_from_oszicar(path):
    # VASP writes Fortran-style floats like "-.14678295E+03"; only accept
    # an E0 if the SCF converged (last electronic step |dE| < EDIFF).
    e0 = None
    last_de = None
    with open(path) as f:
        for line in f:
            m = re.search(r"E0=\s*(-?\d*\.?\d+(?:[Ee][+-]?\d+)?)", line)
            if m:
                e0 = float(m.group(1))
            m2 = re.match(r"\s*(?:DAV|RMM):\s+\d+\s+\S+\s+(\S+)", line)
            if m2:
                last_de = float(m2.group(1))
    if last_de is None or abs(last_de) >= 1e-6:
        return None
    return e0


def natoms(path):
    with open(path) as f:
        return sum(int(x) for x in f.readlines()[6].split())


def e0_from_spilled_outcar(run_dir):
    """VASP 5.4.4 truncates long output paths: for the ~110-char run dirs
    the OUTCAR is written one level up under a truncated name. Recover
    E0 from 'energy(sigma->0)' of any complete spilled OUTCAR whose name
    is a prefix of run_dir."""
    d = os.path.basename(run_dir)
    best = None
    for cand in os.listdir(ROOT):
        p = os.path.join(ROOT, cand)
        if os.path.isfile(p) and d.startswith(cand):
            try:
                txt = open(p, errors="replace").read()
            except OSError:
                continue
            if "Voluntary context switches" not in txt:
                continue
            es = re.findall(
                r"energy\(sigma->0\)\s*=\s*(-?\d+\.?\d*(?:[Ee][+-]?\d+)?)",
                txt)
            if es:
                best = float(es[-1])
    return best


rows = []
for r in MANIFEST:
    base = r["run_dir"]
    rec = dict(r)
    for tag in ("initial", "winner"):
        d = os.path.join(ROOT, f"{base}_{tag}")
        osz = os.path.join(d, "OSZICAR")
        rec[f"E0_{tag}"] = e0_from_oszicar(osz) if os.path.exists(osz) else None
        if rec[f"E0_{tag}"] is None:
            rec[f"E0_{tag}"] = e0_from_spilled_outcar(d)
        rec[f"natoms_{tag}"] = (
            natoms(os.path.join(d, "POSCAR"))
            if os.path.exists(os.path.join(d, "POSCAR")) else None)
    rows.append(rec)

done, rows_out = 0, []
for rec in rows:
    ei, ew = rec.get("E0_initial"), rec.get("E0_winner")
    ni, nw = rec.get("natoms_initial"), rec.get("natoms_winner")
    if ei is not None and ew is not None and ni and nw:
        # compare per atom on the winner's cell size; also record raw totals
        rec["dE_dft_total"] = ew - ei
        rec["dE_dft_per_atom"] = ew / nw - ei / ni
        done += 1
    else:
        rec["dE_dft_total"] = rec["dE_dft_per_atom"] = None
    rows_out.append(rec)

fields = ["run_dir", "material", "dataset", "elements", "e_mace_init",
          "e_mace_winner", "dE_mace", "E0_initial", "E0_winner",
          "natoms_initial", "natoms_winner", "dE_dft_total",
          "dE_dft_per_atom"]
with open(os.path.join(HERE, "dft_whitespace_results.csv"), "w",
          newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    for rec in rows_out:
        w.writerow({k: rec.get(k) for k in fields})

have = [r for r in rows_out if r.get("dE_dft_per_atom") is not None
        and r.get("dE_mace") is not None]
print(f"paired results: {done}/{len(rows_out)}")
if have:
    conf = [r for r in have if r["dE_dft_per_atom"] < 0]
    conf_mace = [r for r in have if r["dE_mace"] < 0]
    both = [r for r in have if r["dE_mace"] < 0 and r["dE_dft_per_atom"] < 0]
    mae = statistics.mean(
        abs(r["dE_dft_per_atom"] - r["dE_mace"]) for r in have)
    print(f"DFT-confirmed stabilization: {len(conf)}/{len(have)}")
    print(f"MACE claimed stabilization:  {len(conf_mace)}/{len(have)}")
    print(f"both agree (winner below init): {len(both)}/{len(have)}")
    print(f"MAE |dE_dft - dE_mace|: {mae*1000:.1f} meV/atom")
    from collections import defaultdict
    by = defaultdict(list)
    for r in have:
        by[r["dataset"]].append(r)
    for ds, rs in by.items():
        c = sum(1 for r in rs if r["dE_dft_per_atom"] < 0)
        print(f"  {ds}: {c}/{len(rs)} DFT-confirmed")
