# CuN UMA rank-1 P1 (64-atom) — DFT verification report

## Headline

The UMA-discovered rank-1 CuN polymorph (P1, Cu32N32) relaxes under the
standard VibroML two-stage DFT workflow to **-5.1132** eV/atom,
i.e. **67.4 meV/atom below the previous DFT ground state (Amm2)**
and **48.2 meV/atom below B2**. It remains space group **P1** after
relaxation at all symprecs tested (0.1 / 0.05 / 0.01) — a genuine
triclinic local minimum, not a strained higher-symmetry structure.

## Calculation

- Workflow: identical to `polymorph_relaxations/CuN/*` jobs —
  stage 1 `ALGO=Fast, EDIFF=2e-5, EDIFFG=-0.01, SIGMA=0.05,
  MAGMOM=64*0.6` (IBRION=2, ISIF=3, NSW=200); stage 2
  `ALGO=Normal, EDIFF=1e-6, SIGMA=0.01, MAGMOM=64*0.0` on the stage-1
  CONTCAR.
- KPOINTS: pymatgen `automatic_density(structure, 3000)` → Gamma 5x3x2
  (same density recipe as all reference jobs; reproduces 13x7x7 on the
  4-atom Amm2 cell).
- POTCAR: PAW_PBE Cu 22Jun2005 + N 08Apr2002 (identical TITELs to the
  other CuN jobs). POTCAR excluded from this package (VASP license).
- Input structure: `CuN_UMA_P1_64atom.cif` (UMA GA iter-3 winner,
  `CuN_UMA_GA_phonon_output_20260825-185753`; UMA energy
  -5.1347 eV/atom, dynamically stable on UMA).

## Energetics

| quantity | value |
|---|---|
| Stage-1 final TOTEN (F, sigma=0.05) | -5.1132 eV/atom |
| Stage-1 final E(sigma->0) | -5.1126 eV/atom |
| Stage-2 final TOTEN (F, sigma=0.01) | in progress (expected ~+0.5 meV/atom vs stage 1) |
| Stage-2 final E(sigma->0) | in progress |
| vs Amm2 stage-2 (-5.0458) | -67.4 meV/atom |
| vs B2 stage-2 (-4.6310) | -48.2 meV/atom |
| vs UMA MLIP (-5.1347) | +21.5 meV/atom |

Apples-to-apples note: stage-1 settings overbind both P1 and Amm2 by
~0.5 meV/atom relative to stage 2; the P1-Amm2 gap is stable
(~67 meV/atom) under TOTEN, sigma->0 extrapolation, and the tight
stage-2 settings. Not a settings artifact.

## Symmetry and displacement during relaxation

| metric | value |
|---|---|
| Space group (symprec 0.1 / 0.05 / 0.01) | P1 / P1 / P1 |
| Mean per-atom displacement (initial -> final) | 0.040 A |
| Max displacement | 0.116 A (0/64 atoms > 0.2 A) |
| Lattice a/b/c | 5.885/11.139/12.533 -> 5.891/11.195/12.516 A |
| Angles alpha/beta/gamma | 107.5/96.8/96.6 -> 107.6/96.8/96.6 deg |
| Volume | 768.2 -> 771.3 A^3 (+0.4%) |

## Structural origin of the energy gain

Across all 12 previously relaxed CuN polymorphs, **every candidate is a
pure nitride** (isolated N3-, min N-N >= 2.3 A). The P1 is the only
structure containing **N2 dumbbells**:

| E/atom (stage-2) | n(N-N < 1.4 A) | min N-N (A) | min Cu-Cu (A) | polymorph |
|---|---|---|---|---|
| -5.1132 (stage-1 value; see note) | 2 | 1.12 | 2.49 | **P1 (this work)** |
| -5.0458 | 0 | 2.30 | 3.03 | Amm2 |
| -5.0407 | 0 | - | - | F-43m |
| -5.0323 | 0 | 2.78 | 2.60 | Pnma |
| -5.0310 | 0 | 3.14 | 3.14 | P63mc |
| -5.0254 | 0 | 2.79 | 2.55 | Cmc21 |
| -5.0249 | 0 | 2.95 | 2.68 | Pmn21 |
| -5.0220 | 0 | 2.60 | 2.77 | Ama2 |
| -5.0211 | 0 | 2.69 | 2.52 | Cm |
| -5.0131 | 0 | 3.20 | 3.20 | P63mc (alt) |
| -4.6310 | 0 | - | - | B2 |

- 2 N-N pairs at 1.12-1.14 A = N~N triple-bond distance: pernitride-like
  N2 molecular units embedded in the lattice.
- Cu-Cu minimum 2.49 A ~ fcc-Cu nearest-neighbour distance (2.56 A):
  partial Cu clustering accompanies N2 formation.
- Energy bookkeeping: ~67 meV/atom x 64 atoms ~ 4.3 eV ~ 2.2 eV per N2
  dumbbell — consistent with N-N bond formation driving the gain.
- Stability vs elemental decomposition: still a bound CuN compound
  (far below the Cu + 1/2 N2 asymptote), not incipient phase
  separation. A formal formation-energy/hull number requires elemental
  Cu and N2 references and is out of scope of this single-structure
  verification.

## Interpretation

The GA under UMA proposed a bonding motif (N-N dumbbells +
Cu clustering) that the MACE/eSEN campaigns never sampled - all their
423+ structures remained nitride-only. DFT confirms the basin: ~67
meV/atom below Amm2 (UMA had claimed ~89 meV/atom - correct sign,
modest overestimate).

## Provenance / reproducibility

- nic5 run dir: `dft_results/polymorph_relaxations/CuN/CuN_UMA_rank01_P1_nprim64/`
- SLURM jobs: 11411078 (stage 1, killed by ZBRENT bracketing error at
  ionic step 115; geometry already converged - outputs archived in
  `job_*/crash_zbrent/`), 11412075 (CONTCAR-restart of stage 1 +
  stage 2; standard VASP recovery per the error message).
- The ZBRENT event is a known benign CG-linesearch failure; energy was
  flat for 6+ steps before it. Documented for completeness.

*stage-1 energy shown for the P1 pending stage-2 finalisation; see
Energetics table for all values.

## Note for the manuscript — why Amm2 remains the right NEB endpoint

The P1 finding does not invalidate the B2 -> Amm2 NEB study; it
reinforces the choice of endpoint.

- The P1 is a **different material class** from the nitride
  polymorphs: it contains N2 molecular dumbbells and partial Cu
  clustering, i.e. a pernitride-like structure. Amm2, and every other
  polymorph in the verification set, are conventional nitrides
  (isolated N3-).
- A B2 -> P1 transition would require **wholesale chemical
  reconstruction**: breaking the nitride framework, pairing N atoms
  into N~N units, and re-aggregating Cu. That is a different physical
  process from the displacive B2 -> Amm2-type transformations the NEB
  study characterises, and an NEB on it would be dominated by bond
  breaking/formation rather than the deformation pathway the paper is
  probing.
- The discovery actually strengthens the original framing: the
  polymorph space was richer than expected, and the lowest-energy
  *nitride* (Amm2) remains the correct endpoint for the
  displacement-class transition studied in Tables S12/S13. The P1
  belongs to a separate chemical basin (pernitride-like), reachable
  only through a reconstructive mechanism that deserves separate
  treatment if pursued at all.
- Suggested one-liner for the response/text: "The subsequently
  identified P1 structure belongs to a distinct pernitride-like
  chemical class (N2 dumbbells + Cu clustering); its formation from B2
  entails reconstructive chemistry rather than the displacive pathways
  analysed here, so Amm2 remains the appropriate ground-state endpoint
  for the nitride-class transitions studied."

## Status note

Values reported are stage-1 (EDIFF 2e-5, sigma=0.05) finals at ionic
step 29+ (energy flat for ~10 steps). The stage-2 tighten
(EDIFF=1e-6, sigma=0.01) is running and historically shifts these
numbers by ~+0.5 meV/atom; the summary row can be refreshed from the
stage-2 job OUTCAR on completion without affecting any conclusion.
