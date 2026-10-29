# HfO2 Literature Validation

Relaxed CIF structures for the $\ce{HfO2}$ polymorphs recovered by VibroML
(MACE-OMAT-0 medium) starting from the cubic fluorite ($Fm\overline{3}m$)
aristotype, corresponding to Supplementary Section S16 (Tables S17 and S18) in
the manuscript.

## Contents

- `HfO2_SG62_Pnma_*.cif` — nonpolar orthorhombic $Pnma$ polymorph matching the
  independently reported `Pnma(b)` phase of Zhao, Zhou & Rappe (2026) (Table S17).
- `HfO2_SG60_Pbcn_*.cif` — antiferroelectric $Pbcn$ parent structure of
  ferroelectric $Pca2_1$ discussed by Raeliarijaona & Cohen (2023) and Qi &
  Rabe (2025).
- `HfO2_SG29_Pca2_1_*.cif`, `HfO2_SG31_Pmn2_1_*.cif`, `HfO2_SG14_P2_1_c_*.cif`,
  `HfO2_SG137_P4_2_nmc_*.cif` — the four known locally stable phases also
  recovered by Qi & Rabe (2025) (Table S18).
- `HfO2_SG141_I4_1_amd_*.cif`, `HfO2_SG136_P4_2_mnm_*.cif`,
  `HfO2_SG11_P2_1_m_*.cif`, `HfO2_SG43_Fdd2_*.cif` — additional locally stable
  basins uncovered by VibroML's unconstrained search outside the $Cmma$
  flat-band subspace (Table S18).
- `compare_hfo2_phases.py` — script that inspects space groups, conventional
  cell atom counts, and lattice parameters across all recovered $\ce{HfO2}$
  structures.
