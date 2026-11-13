# HfO2 Literature Validation

Relaxed CIF structures for the $\ce{HfO2}$ polymorphs recovered by VibroML
(MACE-OMAT-0 medium) starting from the cubic fluorite ($Fm\overline{3}m$)
aristotype, corresponding to Supplementary Section S16 (Tables S17 and S18) in
the manuscript.

## Contents

- `qi_rabe_reference/` — the 12 relaxed reference structures (POSCAR format)
  published by Qi & Rabe for their Fig. 3 phases, downloaded from their
  public repository https://github.com/yuboqiuab/unstableflatband
  (commit a438e4ecf63cddfa68d4ae8ea83d4c19bd238969). See the README inside
  that directory for the file-to-phase mapping.
- `HfO2_SG62_Pnma_*.cif` — nonpolar orthorhombic $Pnma$ polymorph matching the
  independently reported `Pnma(b)` phase of Zhao, Zhou & Rappe (2026) (Table S17).
- `HfO2_SG60_Pbcn_*.cif` — antiferroelectric $Pbcn$ parent structure of
  ferroelectric $Pca2_1$ discussed by Raeliarijaona & Cohen (2023) and Qi &
  Rabe (2025).
- `HfO2_SG29_Pca2_1_*.cif`, `HfO2_SG31_Pmn2_1_*.cif`,
  `HfO2_SG14_P2_1_c_top_1_iter3_sample26_*.cif`,
  `HfO2_SG137_P4_2_nmc_*.cif` — the four known locally stable phases (oIII,
  oIV, m, t) also recovered by Qi & Rabe (2025) (Table S18).
- `HfO2_SG14_P2_1_c_unique_1_iter2_sample63_*.cif` — 24-atom $P2_1/c$
  period-2 domain superlattice (two mirror-related distorted-monoclinic
  layers), matching Qi & Rabe's predicted dm-dm phase (their Fig. 3(l))
  with a symmetry-aligned RMS displacement of 0.02 Angstrom (Table S18).
  This is a distinct polymorph from the 12-atom baddeleyite $m$ phase;
  both share space group $P2_1/c$ (SG 14).
- `HfO2_SG141_I4_1_amd_*.cif`, `HfO2_SG136_P4_2_mnm_*.cif`,
  `HfO2_SG11_P2_1_m_*.cif`, `HfO2_SG43_Fdd2_*.cif` — additional locally stable
  basins uncovered by VibroML's unconstrained search outside the $Cmma$
  flat-band subspace (Table S18).
- `compare_hfo2_phases.py` — script that inspects space groups, conventional
  cell atom counts, and lattice parameters across all recovered $\ce{HfO2}$
  structures, then cross-matches them against the Qi & Rabe reference
  structures with `pymatgen.analysis.structure_matcher.StructureMatcher`
  (ltol=0.10, stol=0.20, angle_tol=3, best match by RMS).

## Verified cross-match (re-run `compare_hfo2_phases.py` to reproduce)

| VibroML structure                    | Qi & Rabe | Label | RMS (A) |
|--------------------------------------|-----------|-------|---------|
| top_1_iter3_sample26 (P2_1/c, 12)    | 3c        | m     | 0.062   |
| unique_1_iter2_sample63 (P2_1/c, 24) | 3l        | dm-dm | 0.019   |
| top_5_iter3_sample6 (Pca2_1, 12)     | 3b        | oIII  | 0.028   |
| unique_6_iter4_sample68 (Pmn2_1, 6)  | 3d        | oIV   | 0.032   |
| unique_5_iter1_sample23 (P4_2/nmc, 6)| 3a        | t     | 0.019   |

VibroML recovered 5 of Qi & Rabe's 12 locally stable phases. The remaining
VibroML structures (Pbcn, I4_1/amd, P4_2/mnm, P2_1/m, Fdd2, and the
Zhao-matched Pnma) have no counterpart among the 12.
