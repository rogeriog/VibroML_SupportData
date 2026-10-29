#!/usr/bin/env python3
"""Compare VibroML-discovered HfO2 polymorphs against literature phases.

Reproduces Supplementary Tables S17 (Pnma comparison with Zhao, Zhou & Rappe,
2026) and S18 (cross-match against Qi & Rabe, 2025).
"""

import glob
import os
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

HERE = os.path.dirname(os.path.abspath(__file__))


def main() -> None:
    cifs = sorted(glob.glob(os.path.join(HERE, "HfO2_SG*.cif")))
    print(f"Found {len(cifs)} VibroML HfO2 structures:\n")
    print(f"{'Space Group':<16} {'SG#':<6} {'N_conv':<8} {'a (A)':<8} {'b (A)':<8} {'c (A)':<8} {'File'}")
    print("-" * 96)
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


if __name__ == "__main__":
    main()
