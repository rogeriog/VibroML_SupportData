# B2 Intermetallic Data

Phonon evaluation outputs and DFT relaxations for the three B2
(CsCl-type) benchmark systems: CuN, ZnCu, and ZrCu.

## Phonon metrics (Table S12)

Each `b2_metrics_*` directory contains the `phonon_only` evaluation of
the corresponding engine's GA rank-one structure (primitive cell) with
the indicated MLIP engine (MACE-OMAT, eSEN, or UMA). All nine entries
have $f_\mathrm{soft} = 0$ and $I_\mathrm{curv} = 0$ on the sampled mesh,
confirming dynamical stability.

| System | Engine | Rank-one structure | Directory |
|--------|--------|--------------------|-----------|
| CuN | MACE | $P6_3mc$ | `b2_metrics_CuN_MACE_*` |
| CuN | eSEN | Amm2 | `b2_metrics_CuN_ESEN_*` |
| CuN | UMA | $P1$ (64-atom cell) | `b2_metrics_CuN_UMA_guidance_check/` |
| ZnCu | MACE | Pmma | `b2_metrics_ZnCu_MACE_rerun3x3x3/` |
| ZnCu | eSEN | Pmma | `b2_metrics_ZnCu_ESEN_*` |
| ZnCu | UMA | Pnma (8-atom cell) | `b2_metrics_ZnCu_UMA_*` |
| ZrCu | MACE | $P\bar{1}$ | `b2_metrics_ZrCu_MACE_*` |
| ZrCu | eSEN | $P2_1/c$ | `b2_metrics_ZrCu_ESEN_*` |
| ZrCu | UMA | $P2_1/c$ | `b2_metrics_ZrCu_UMA_*` |

The ZnCu–MACE-OMAT entry comes from a GA run launched with a 3×3×3
initial screen; the automatic 2×2×2 sweep terminated at the cubic B2
parent before launching the GA (see `zncu_mace_screen_termination/` for the run
artifacts cited in the table caption). The CuN–UMA entry is the
iteration-3 certified structure from the August GA campaign, evaluated
in its guidance phonon check.

These data correspond to Supplementary Table S_b2_phonon_metrics in the
manuscript.

## DFT polymorph relaxations

`polymorph_relaxations.tar.gz` contains the VASP inputs and outputs of
full DFT relaxations of all GA-discovered polymorph sets (31 polymorphs,
62 job directories including custodian restarts), extracted in the same
layout as the source tree:

```
CuN/<polymorph_dir>/job_*/   ZnCu/<polymorph_dir>/job_*/   ZrCu/<polymorph_dir>/job_*/
```

Each job directory holds `POSCAR`, `POSCAR.orig`, `CONTCAR`, `INCAR`,
`INCAR.orig`, `KPOINTS`, `KPOINTS.orig`, `IBZKPT`, `OUTCAR`, `OSZICAR`,
`vasprun.xml`, and `custodian.json`. **POTCARs are excluded** (licensed
VASP PAWs); the required POTCAR species are recorded in each `OUTCAR`
header (`POTCAR:` TITEL lines).

The three loose `polymorph_summary_*.csv` files list
`material, run_name, formula, n_atoms, energy, energy_per_atom, status`
for every relaxed polymorph. Within the conventional nitride and
intermetallic polymorph families, the DFT-relaxed lowest-energy
structures are Amm2 (CuN, −5.0458 eV/atom), Pmma (ZnCu, −2.5115 eV/atom),
and $P2_1/c$ (ZrCu, −6.2834 eV/atom); the cubic B2 references relax to
−4.6310, −2.5086, and −6.2372 eV/atom respectively.

## CuN UMA $P1$ 64-atom basin

`CuN_UMA_P1_dft_relax/` contains the DFT relaxation of UMA's rank-1
CuN candidate, a 64-atom triclinic structure outside the nitride
polymorph family. It relaxes to −5.1129 eV/atom (stage-2, converged to
EDIFFG −0.01; max residual force 0.009 eV/Å), 67 meV/atom below Amm2,
and contains nearly molecular N2 units (two N–N pairs at 1.12–1.14 Å,
near the N≡N triple-bond distance) with partial Cu clustering
(minimum Cu–Cu 2.49 Å), consistent with incipient decomposition toward
Cu and N2. `CuN_UMA_P1_analysis.md` in that directory documents the
settings, energetics, and structural analysis; the DFT-relaxed
geometry is deposited as `CuN_UMA_P1_64atom_dft_relaxed.cif`
(`CuN_UMA_P1_64atom.cif` is the UMA input structure).
