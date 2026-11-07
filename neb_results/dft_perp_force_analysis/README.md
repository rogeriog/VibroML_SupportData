# NEB Perpendicular-Force Decomposition (GMR item F)

Reviewer request: "add the perpendicular force component, or at least the
endpoint forces, next to the maximum total force" in SI Table S12.

## Status

- `dft_inputs/neb_images_to_gs/` — 63 VASP input dirs (POSCAR/INCAR/KPOINTS/
  POTCAR), extracted from `vibroml_runs/dft_neb_statics_to_gs.tar.gz`.
  These are the same MLIP-NEB geometries used for Table S12.
- The original OUTCARs were not retained, so force vectors must be recomputed.
- No VASP installation exists on Lucia; these single points ran on an
  external machine previously.

## Workflow

1. Copy `dft_inputs/` + `run_all_vasp.sh` to a machine with `vasp_std`.
2. `bash run_all_vasp.sh 24` — runs/skips the 63 statics (resumable).
3. Return the tree (only `OUTCAR` needed per image dir).
4. `python decompose_forces.py` — writes:
   - `perp_force_per_image.csv` — per-image 3N norms of the total,
     path-parallel, and path-perpendicular DFT force plus max per-atom total.
   - `perp_force_summary.csv` — per-path `F_maxatom_tot`,
     `F_3N_tot_max_interior`, `F_perp_3N_max_interior`, `F_endpoint_max`
     (columns destined for SI Table S12).

## Method

Interior image i: tangent t_i = (R_{i+1} − R_{i−1}) normalized in the 3N
configuration space; F_par = (F·t̂)t̂, F_perp = F − F_par. Endpoint images
report the unprojected residual (the full force should vanish there).

## What is needed for the paper

**Manuscript locations affected:**
- `tables/tableS_neb_comprehensive.tex` (SI Table S12) — currently shows
  only `F_max^DFT` (total) and `F_max^MLIP`. Two columns will be added:
  `F_perp^DFT` (max perpendicular residual over interior images) and
  `F_end^DFT` (max residual at endpoint images).
- SI Section S9 text (~SI line 514) — the sentence asserting that the
  perpendicular component vanishes on an exact MEP will be backed by the
  numbers rather than asserted.
- `response_to_referees.tex` item F — point-by-point reply citing the new
  columns.

**Numbers consumed per path (7 paths: CuN/ZnCu/ZrCu × MACE/eSEN/UMA):**
- `perp_force_summary.csv` columns: `F_maxatom_tot` (existing table
  column, recomputed as a cross-check against the published values
  1.87/0.37/0.59/0.30/0.79/1.19/0.57/0.40/0.37 eV/Å),
  `F_perp_3N_max_interior` (new `F_perp^DFT` column),
  `F_endpoint_max` (new `F_end^DFT` column).

**What must come back:** `OUTCAR` from each of the 63 image dirs
(image_00..image_08 under each of the 7 path dirs). If only a partial run
is feasible, the 14 endpoint images (image_00 and image_08 per path)
satisfy the reviewer's minimum ask ("or at least the endpoint forces").
