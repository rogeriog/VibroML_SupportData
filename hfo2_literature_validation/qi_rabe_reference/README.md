# Qi & Rabe Reference Structures (HfO2)

Relaxed reference structures published by Qi & Rabe for the 12 locally
stable phases derived from the unstable flat phonon bands of the
$Cmma$ HfO$_2$ reference structure.

## Source

- Publication: Y. Qi and K. M. Rabe, *Competing Phases of HfO$_2$ from
  Unstable Flat Phonon Bands of an Unconventional High-Symmetry
  Structure*, Phys. Rev. Lett. **135**, 046101 (2025).
  DOI: 10.1103/PhysRevLett.135.046101; arXiv:2412.16792.
- Repository: https://github.com/yuboqiuab/unstableflatband
  (linked from the paper's Supplemental Material note: "Data,
  structures, and pseudopotentials involved in this study can be found
  at https://github.com/yuboqiuab/unstableflatband").
- Snapshot taken: commit `a438e4ecf63cddfa68d4ae8ea83d4c19bd238969`
  (2025-06-18), `figure3/` directory, files `3a.vasp`–`3l.vasp`
  (VASP POSCAR format, DFT-LDA relaxed coordinates).

## File-to-phase mapping (Fig. 3 of the publication)

| File     | Space group  | Atoms | Label  | Notes |
|----------|--------------|-------|--------|-------|
| 3a.vasp  | P4_2/nmc (137) | 12  | t      | known phase |
| 3b.vasp  | Pca2_1 (29)    | 12  | oIII   | ferroelectric |
| 3c.vasp  | P2_1/c (14)    | 12  | m      | baddeleyite ground state |
| 3d.vasp  | Pmn2_1 (31)    | 12  | oIV    | known phase |
| 3e.vasp  | Pbca (61)      | 24  | o-AP   | antipolar oIII domain structure |
| 3f.vasp  | Pbca (61)      | 24  | m-dom  | unit-cell-scale domain structure of m |
| 3g.vasp  | C2/c (15)      | 24  | new    | |
| 3h.vasp  | Pc (7)         | 24  | new    | |
| 3i.vasp  | Cc (9)         | 24  | new    | |
| 3j.vasp  | Pmna (62)      | 24  | new    | period-2 oIV superlattice |
| 3k.vasp  | P2_1 (4)       | 24  | new    | |
| 3l.vasp  | P2_1/c (14)    | 24  | dm-dm  | two mirror-related distorted-m domains |

## Verification performed

`compare_hfo2_phases.py` (parent directory) loads these reference
structures and cross-matches them against the VibroML-recovered CIFs
with `pymatgen.analysis.structure_matcher.StructureMatcher`. Confirmed
matches (RMS displacement after symmetry alignment):

- VibroML `top_1_iter3_sample26` (P2_1/c, 12 atoms) == `3c.vasp` (m).
- VibroML `unique_1_iter2_sample63` (P2_1/c, 24 atoms) == `3l.vasp`
  (dm-dm), RMS = 0.019 Angstrom.
