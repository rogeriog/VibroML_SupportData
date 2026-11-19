#!/usr/bin/env python3
"""Scan every certified GA winner in the explorer DB for loss of
structural coherence ("decomposition-trajectory" signatures).

For each (seed, winner) pair a periodic bond graph is built with bonds
defined as contacts below 1.30 x (sum of covalent radii). The scan writes
per-structure metrics to structural_coherence_flags.csv:

  frac_unbound   fraction of sites with zero bonded neighbours
  frac_lowcoord  fraction of sites with <2 bonded neighbours
  frac_finite    fraction of atoms in bond-graph components containing no
                 periodic edge (fragments that do not wrap the cell)
  frag_n         number of such finite components
  loccomp_std    std over sites of the local metal fraction in a 5 A
                 sphere (spatial homogeneity; phase separation raises it)
  d_loccomp      loccomp_std(winner) - loccomp_std(seed)
  hetbond        hetero-species fraction of bonded contacts
  cn_cv          max over species of CV of bonded coordination number
  void_frac      fraction of ~12k random cell points >2.8 A from the
                 nearest atom (27-image minimum)

Flagged (incoherent) = unbound site, non-percolating fragment, or
void_frac > 0.15. Seeds provide the physical baseline (~0 on all axes).

Usage: python scan_structural_coherence.py
"""

Design: the target is NOT family retention (valid polymorphs legitimately
leave the seed family). The target is physical incoherence — a structure
on a decomposition trajectory. Metrics:

AXIS A — connectivity integrity (periodic bond graph, bonds =
sub-covalent contacts d < 1.30*(rcov_i + rcov_j)):
  frac_unbound   sites with zero bonded neighbours / N
  frac_lowcoord  sites with <2 bonded neighbours / N
  frac_finite    sites in bond-graph components containing no periodic
                 edge (image != 0,0,0) i.e. fragments that do NOT wrap the
                 cell / N. A real crystal — even layered — has every atom
                 in a component that wraps at least one direction; a
                 finite fragment is a cluster floating in a box.
  frag_n         number of such finite components

AXIS B — spatial homogeneity (decomposition = phase separation):
  loccomp_std    std over sites of the local metal fraction in a 5 A
                 periodic sphere
  d_loccomp      loccomp_std_winner - loccomp_std_seed (segregation delta;
                 uniform cluster compounds stay near seed level)

AXIS C — bond-level chemistry:
  hetbond        hetero-species fraction of bonded contacts in winner
  d_hetbond      hetbond_winner - hetbond_seed (interpenetration loss)

AXIS D — coordination heterogeneity + voids:
  cn_cv          max over species of std(CN_bond)/mean(CN_bond)
  d_cn_cv        cn_cv_winner - cn_cv_seed (heterogeneous environments)
  void_frac      fraction of ~12k uniform random cell points whose
                 nearest atom is > 2.8 A away (27-image minimum)

Seeds are the physical baseline: unbound/finite/void should be ~0 for
every real structure. Deltas are seed-referenced so family-dependent
baselines cancel. Output frankenstein_panel.csv. Pool(8).
"""

import csv
import os
import sqlite3
import sys
import warnings
from multiprocessing import Pool

import numpy as np
from pymatgen.core import Element, Structure

warnings.filterwarnings("ignore")

BASE_DIR = "/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/02_runs_and_data/candidates_lowcoverage_cif"
OUT_CSV = "/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/02_runs_and_data/VibroML_SupportData/whitespace_data/structural_coherence_flags.csv"
BASE = BASE_DIR
DB = os.path.join(BASE_DIR, "database", "vibroml_database.db")
OUT = OUT_CSV
ANIONS = {"N", "O", "F", "P", "S", "Cl", "As", "Se", "Br", "Sb", "Te", "I", "H"}
BOND_K = 1.30
LSPHERE = 5.0
VOID_D = 2.8
RNG = np.random.default_rng(7)


def metrics(s):
    n = len(s.sites)
    syms = [st.specie.symbol for st in s.sites]
    rad = [Element(e).atomic_radius_calculated for e in syms]
    fc = s.frac_coords
    nbrs = s.get_all_neighbors(LSPHERE)
    deg = [0] * n
    parent = list(range(n))
    wrap_root = {}          # root -> has periodic edge
    het = tot = 0
    lf = []                  # local metal fraction
    cats = {i for i in range(n) if syms[i] not in ANIONS}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    wrap_edges = []          # (root-at-time-of-add is stale; store site)
    edges = []
    for i in range(n):
        nc = ntt = 0
        for nb in nbrs[i]:
            j = nb.index
            d = nb.nn_distance
            if d < LSPHERE:
                ntt += 1
                nc += j in cats
            if j <= i:
                continue
            if d < BOND_K * (rad[i] + rad[j]):
                deg[i] += 1
                deg[j] += 1
                tot += 1
                het += syms[i] != syms[j]
                img = getattr(nb, 'image', None)
                if img is None:
                    img = np.round(nb.frac_coords - fc[j])
                edges.append((i, j, np.asarray(img)))
                ri, rj = find(i), find(j)
                if ri != rj:
                    parent[ri] = rj
                    w = wrap_root.get(ri, False) or wrap_root.get(rj, False)
                    wrap_root.pop(ri, None)
                    wrap_root.pop(rj, None)
                    wrap_root[rj] = w
                if np.any(img != 0):
                    wrap_root[find(j)] = True
        if ntt:
            lf.append(nc / ntt)
    comps = {}
    for i in range(n):
        comps.setdefault(find(i), []).append(i)
    finite_atoms = frag_n = 0
    for root, mem in comps.items():
        w = wrap_root.get(root, False) or any(
            np.any(e[2] != 0) for e in edges
            if find(e[0]) == root and e[0] in mem)
        if not w:
            frag_n += 1
            finite_atoms += len(mem)
    cncv = 0.0
    for e in set(syms):
        v = np.array([deg[i] for i in range(n) if syms[i] == e])
        m = v.mean()
        if m > 0:
            cncv = max(cncv, float(v.std() / m))
    # void sampling
    lat = s.lattice
    pts = RNG.random((12000, 3))
    fdiff = pts[:, None, :] - fc[None, :, :]
    mind = np.linalg.norm(fdiff @ lat.matrix, axis=-1).min(axis=1)
    for img in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1),
                (0, 0, -1), (1, 1, 0), (-1, -1, 0), (1, -1, 0), (-1, 1, 0),
                (1, 0, 1), (-1, 0, -1), (1, 0, -1), (-1, 0, 1),
                (0, 1, 1), (0, -1, -1), (0, 1, -1), (0, -1, 1),
                (1, 1, 1), (-1, -1, -1), (1, 1, -1), (-1, -1, 1),
                (1, -1, 1), (-1, 1, -1), (1, -1, -1), (-1, 1, 1)):
        dd = np.linalg.norm((fdiff + np.asarray(img)) @ lat.matrix, axis=-1)
        mind = np.minimum(mind, dd.min(axis=1))
    void = float(np.mean(mind > VOID_D))
    return dict(frac_unbound=sum(1 for d in deg if d == 0) / n,
                frac_lowcoord=sum(1 for d in deg if d < 2) / n,
                frac_finite=finite_atoms / n, frag_n=frag_n,
                hetbond=het / tot if tot else np.nan, cn_cv=cncv,
                loccomp_std=float(np.std(lf)) if lf else np.nan,
                void_frac=void)


def work(t):
    stem, ds, ci, cf, de = t
    rec = dict(stem=stem, dataset=ds, dE=de)
    fi, ff = os.path.join(BASE, ci or ''), os.path.join(BASE, cf or '')
    if not (ci and cf and os.path.exists(fi) and os.path.exists(ff)):
        rec['err'] = 'missing'
        return rec
    try:
        si = Structure.from_file(fi)
        sf = Structure.from_file(ff)
        mi = metrics(si)
        mf = metrics(sf)
    except Exception as e:
        rec['err'] = str(e)[:60]
        return rec
    for k in ('frac_unbound', 'frac_lowcoord', 'frac_finite', 'frag_n',
              'hetbond', 'cn_cv', 'loccomp_std', 'void_frac'):
        rec[k + '_i'] = mi[k]
        rec[k + '_f'] = mf[k]
    rec['d_hetbond'] = mf['hetbond'] - mi['hetbond']
    rec['d_loccomp'] = mf['loccomp_std'] - mi['loccomp_std']
    rec['d_cn_cv'] = mf['cn_cv'] - mi['cn_cv']
    rec['d_void'] = mf['void_frac'] - mi['void_frac']
    return rec


def main():
    con = sqlite3.connect(DB)
    rows = con.execute(
        """SELECT cif_stem, dataset, cif_path, cif_optimized_path,
                  energy_reduction_vs_initial_eV_per_atom
           FROM materials WHERE run_status='ga_completed'""").fetchall()
    print(len(rows))
    with Pool(8) as p:
        out = []
        for n, rec in enumerate(p.imap_unordered(work, rows, chunksize=15), 1):
            out.append(rec)
            if n % 300 == 0:
                print(n, flush=True)
    keys = sorted({k for r in out for k in r})
    with open(OUT, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(out)
    print('DONE')


if __name__ == '__main__':
    sys.exit(main())
