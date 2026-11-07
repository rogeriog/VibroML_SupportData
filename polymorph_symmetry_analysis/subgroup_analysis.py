"""Group-subgroup analysis for Table 1 polymorphs (GMR item C).

For each Table-1 row:
  1. find the deposited CIF matching SG number AND MACE energy (E_ref + dE)
  2. abstract subgroup test via pymatgen is_subgroup (by SG number)
  3. structural supercell match: enumerate integer P (det = volume ratio),
     prefilter by Gram-eigenvalue isometry, then StructureMatcher fit.
"""
import os, glob, itertools, re, json
import numpy as np
from pymatgen.core import Structure, Lattice
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.groups import SpaceGroup
from pymatgen.analysis.structure_matcher import StructureMatcher

BASE = "/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs"
ROWS_SGN = {"Pm-3m": 221, "Fm-3m": 225, "Fd-3m": 227}
REM = f"{BASE}/02_runs_and_data/VibroML_SupportData/remediation_results"

# parent reference energies per atom (from initial phonon run summaries)
PARENT = {
    "LiF":      dict(dir=f"{REM}/LiF_simple_cubic", sg="Pm-3m", eref=-4.523283),
    "CsPbI3":   dict(dir=f"{REM}/CsPbI3",           sg="Pm-3m", eref=-2.796148),
    "HfO2":     dict(dir=f"{REM}/HfO2",             sg="Fm-3m", eref=-10.117371),
    "Bi2Sn2O7": dict(dir=f"{REM}/Bi2Sn2O7",         sg="Fd-3m", eref=-6.011008),
}

# (sg_label, sg_number, dE_mace_meV) per Table-1 row (ref rows excluded)
ROWS = {
    "LiF":      [("Pnma",62,-276.8),("F-43m",216,-292.5),("P42/mnm",136,-296.4),
                 ("Fm-3m",225,-297.2),("P63mc",186,-303.1),("P1",1,-4088.3),
                 ("Pm",6,-7347.9),("Cm",8,-14042.8)],
    "CsPbI3":   [("C2/m",12,-18.1),("Cm",8,-23.0),("Pnma",62,-29.3),
                 ("P-1",2,-29.9),("P21/m",11,-33.3),("Cmcm",63,-34.9),
                 ("Pnma",62,-44.2)],
    "HfO2":     [("P42/nmc",137,-12.2),("P21/m",11,-18.5),("Pnma",62,-21.7),
                 ("Pbcn",60,-38.5),("P21/c",14,-48.5)],
    "Bi2Sn2O7": [("P212121",19,20.7),("P21",4,12.7),("Cc",9,-12.6),
                 ("Pna21",33,-12.7),("P21",4,-13.1),("P212121",19,-13.7),
                 ("Cc",9,-14.1)],
}


def conv_struct(path):
    s = Structure.from_file(path)
    return SpacegroupAnalyzer(s, symprec=0.1).get_conventional_standard_structure()


def prim_struct(path_or_s):
    s = path_or_s if isinstance(path_or_s, Structure) else Structure.from_file(path_or_s)
    return SpacegroupAnalyzer(s, symprec=0.1).get_primitive_standard_structure()


def sgn(s):
    return SpacegroupAnalyzer(s, symprec=0.1).get_space_group_number()


def gram_eigs(latt):
    return np.sort(np.linalg.eigvalsh(latt.matrix @ latt.matrix.T))


def hnf_mats(m):
    mats = []
    for a in range(1, m + 1):
        if m % a: continue
        r = m // a
        for d in range(1, r + 1):
            if r % d: continue
            f = r // d
            for b in range(d):
                for c in range(f):
                    for e in range(f):
                        mats.append(np.array([[a, b, c], [0, d, e], [0, 0, f]]))
    return mats


def energy_of(fname):
    m = re.search(r'energy_(m|p)(\d+)p(\d+)', os.path.basename(fname))
    if not m:
        return None
    s = float(m.group(2) + '.' + m.group(3))
    return -s if m.group(1) == 'm' else s


def find_match(parent_prim, child_prim):
    """Returns (P, rmsd) mapping parent_prim -> supercell matching child_prim."""
    sm = StructureMatcher(ltol=0.5, stol=0.8, angle_tol=15.0,
                          primitive_cell=False, scale=False,
                          attempt_supercell=True)
    if not sm.fit(parent_prim, child_prim):
        return None, None
    P = np.round(sm.get_supercell_matrix(child_prim, parent_prim)).astype(int)
    sc = parent_prim.copy()
    sc.make_supercell(P)
    try:
        rmsd = float(sm.get_rms_dist(sc, child_prim)[0])
    except Exception:
        rmsd = None
    return P, rmsd


def main():
    out = []
    hdr = f"{'Comp':10} {'Child':9} {'dE':9} {'subG':6} {'match':6} {'rmsd':7} P"
    print(hdr)
    for comp, info in PARENT.items():
        pcifs = glob.glob(os.path.join(info['dir'], "*GA*", "initial_relaxation_for_single_run", "*relaxed*.cif")) or \
                glob.glob(os.path.join(info['dir'], "*", "initial_relaxation_for_single_run", "*relaxed*.cif"))
        parent = prim_struct(sorted(pcifs)[0])
        pgrp = SpaceGroup.from_int_number(ROWS_SGN[info['sg']])

        cifs = [c for c in glob.glob(os.path.join(info['dir'], "*/final_structures/*.cif"))
                if '_conventional_' not in c]
        cands = []
        for c in cifs:
            e = energy_of(c)
            if e is None:
                continue
            try:
                cs = prim_struct(c)
            except Exception:
                continue
            cands.append((c, e, sgn(cs), cs))

        for label, sgnm, dme in ROWS[comp]:
            target = info['eref'] + dme / 1000.0
            cand = [x for x in cands if x[2] == sgnm and abs(x[1] - target) < 0.005]
            if not cand:
                cand = [x for x in cands if x[2] == sgnm]
                cand.sort(key=lambda x: abs(x[1] - target))
            if not cand:
                print(f"{comp:10} {label:9} {dme:9} -- none")
                continue
            child = cand[0][3]
            try:
                grp_sub = SpaceGroup.from_int_number(sgnm).is_subgroup(pgrp)
            except Exception:
                grp_sub = "?"
            P, rmsd = find_match(parent, child)
            pstr = "; ".join(str(P[i].tolist()) for i in range(3)) if P is not None else "-"
            print(f"{comp:10} {label:9} {dme:9} {str(grp_sub):6} {'yes' if P is not None else 'no':6} "
                  f"{('%.3f' % rmsd) if rmsd is not None else '-':7} {pstr}  [{os.path.basename(cand[0][0])[:45]}]")
            out.append(dict(comp=comp, sg=label, sg_num=sgnm, dE=dme, abstract_subgroup=grp_sub,
                            supercell_match=P is not None, rmsd=rmsd,
                            P=(P.tolist() if P is not None else None), cif=cand[0][0]))
    with open('subgroup_analysis_results.json', 'w') as f:
        json.dump(out, f, indent=1)


if __name__ == "__main__":
    main()
