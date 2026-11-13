#!/usr/bin/env python3
"""Compare VibroML-discovered HfO2 polymorphs against literature phases.

Reproduces Supplementary Tables S17 (Pnma comparison with Zhao, Zhou & Rappe,
2026) and S18 (cross-match against Qi & Rabe, 2025).

Reference structures: Qi & Rabe, Phys. Rev. Lett. 135, 046101 (2025),
downloaded from https://github.com/yuboqiuab/unstableflatband
(commit a438e4ecf63cddfa68d4ae8ea83d4c19bd238969, directory figure3/).
"""

import glob
import os
from pymatgen.core import Structure
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

HERE = os.path.dirname(os.path.abspath(__file__))

# Qi & Rabe Fig. 3 panel -> (label, is "known" phase)
QR_LABELS = {
    "3a": "t", "3b": "oIII", "3c": "m", "3d": "oIV", "3e": "o-AP",
    "3f": "m-dom", "3g": "C2/c", "3h": "Pc", "3i": "Cc",
    "3j": "Pmna", "3k": "P2_1", "3l": "dm-dm",
}


def main() -> None:
    cifs = sorted(glob.glob(os.path.join(HERE, "HfO2_SG*.cif")))
    print(f"Found {len(cifs)} VibroML HfO2 structures:\n")
    print(f"{'Space Group':<16} {'SG#':<6} {'N_conv':<8} {'a (A)':<8} {'b (A)':<8} {'c (A)':<8} {'File'}")
    print("-" * 96)
    vibroml = {}
    for path in cifs:
        s = Structure.from_file(path)
        sga = SpacegroupAnalyzer(s, symprec=0.1)
        conv = sga.get_conventional_standard_structure()
        a, b, c = conv.lattice.abc
        print(
            f"{sga.get_space_group_symbol():<16} "
            f"{sga.get_space_group_number():<6} "
            f"{len(conv):<8} "
            f"{a:<8.3f} {b:<8.3f} {c:<8.3f} "
            f"{os.path.basename(path)}"
        )
        vibroml[os.path.basename(path)] = s

    refs = {}
    for path in sorted(glob.glob(os.path.join(HERE, "qi_rabe_reference", "*.vasp"))):
        refs[os.path.basename(path)] = Structure.from_file(path)

    if not refs:
        print("\nNo Qi-Rabe reference files found in qi_rabe_reference/.")
        return

    print(f"\n{'='*96}\nCross-match vs Qi & Rabe reference structures "
          f"(StructureMatcher, ltol=0.10, stol=0.20, angle_tol=3; "
          f"best match by RMS):\n")
    matcher = StructureMatcher(
        primitive_cell=True, attempt_supercell=True,
        ltol=0.10, stol=0.20, angle_tol=3,
    )
    print(f"{'VibroML structure':<72} {'Match':<10} {'RMS (A)':<8} {'QR label'}")
    print("-" * 96)
    for name, s in vibroml.items():
        best = None
        for rname, q in refs.items():
            if matcher.fit(q, s):
                rms = matcher.get_rms_dist(q, s)[0]
                if best is None or rms < best[0]:
                    best = (rms, rname)
        if best is None:
            print(f"{name:<72} {'---':<10} {'---':<8} (VibroML-only)")
        else:
            qr = best[1].replace(".vasp", "")
            print(f"{name:<72} {qr:<10} {best[0]:<8.3f} {QR_LABELS.get(qr, '?')}")


if __name__ == "__main__":
    main()
