# NEB force decomposition (SI Table S13)

Decomposition of the residual DFT forces on MLIP-generated NEB images into
path-parallel and path-perpendicular components, in the 3N configuration
space. This produces the `F_perp^DFT`, `F_par^DFT`, and `F_end^DFT` columns
of SI Table S13 for all nine B2 NEB paths (CuN / ZnCu / ZrCu x
MACE-OMAT / eSEN / UMA).

## Status

Complete. All nine paths (81 image statics) were decomposed; Table S13 has
no missing entries.

- Seven path groups (CuN x 3 engines, ZnCu/ZrCu x eSEN/UMA) came from the
  archived static single points in
  `vibroml_runs/dft_neb_statics_to_gs.tar.gz`.
- The two MACE-OMAT groups (`ZnCu_MACE_to_gs`, `ZrCu_MACE_to_gs`) were
  recovered from the original VASP run tree on nic5, where they were stored
  under the earlier naming `ZnCu_MACE_poly02_Pmma` and
  `ZrCu_MACE_poly01_P2_1`. Endpoint identity was verified by structure
  matching against the eSEN `to_gs` groups (RMS ~ 0.03 A, different cell
  setting) and by reproducing the Table S13 scalars exactly:
  ZnCu/MACE E_b^DFT = 448 meV/atom and F_max^DFT = 0.30 eV/A;
  ZrCu/MACE E_b^DFT = 59 meV/atom and F_max^DFT = 0.57 eV/A.

## Contents

- `decompose_forces.py` — parser + decomposition script.
- `neb_images_to_gs/` — the two recovered MACE-OMAT groups only
  (18 image dirs with `OUTCAR`, `POSCAR`, `INCAR`, `KPOINTS`).
  `vasprun.xml`/`CONTCAR`/`OSZICAR` are omitted for size; `OUTCAR` is all
  the analysis needs. The seven other groups' OUTCARs are retained in the
  working tree (`vibroml_runs/02_runs_and_data/dft_neb_perp_analysis/`),
  not duplicated here.
- `perp_force_per_image.csv` — per-image norms: `F_tot_3N`, `F_par_3N`,
  `F_perp_3N`, `F_max_per_atom`, `F_perp_max_per_atom`,
  `F_par_max_per_atom`.
- `perp_force_summary.csv` — per-path summary: `F_maxatom_tot`,
  `F_perp_max_per_atom`, `F_perp_3N_max_interior`, `F_par_max_per_atom`,
  `F_endpoint_max`.
- `run_all_vasp.sh` — the batch script originally used to run the statics.

## Method and table conventions

Interior image i: path tangent t_i = (R_{i+1} - R_{i-1}) normalized in the
3N configuration space; F_par_vec = (F.t_hat)t_hat, F_perp_vec = F - F_par_vec.
Endpoint images report the unprojected residual (the full force should
vanish there).

Mapping to Table S13 columns:

- `F_perp^DFT` = `F_perp_max_per_atom`: max over interior images of the
  largest per-atom norm of F_perp_vec.
- `F_par^DFT` = `F_par_max_per_atom`: max over interior images of the
  largest per-atom norm of F_par_vec.
- `F_end^DFT` = `F_endpoint_max`: max 3N norm of the residual at the two
  endpoint images.

All published values (F_max^DFT cross-check, F_perp, F_par, F_end)
reproduce from `perp_force_summary.csv` at two decimals for the full
9-path table.

## Usage

`python decompose_forces.py` inside the working-tree analysis directory
(`vibroml_runs/02_runs_and_data/dft_neb_perp_analysis/`), which contains
all nine path groups under `dft_inputs/neb_images_to_gs/`.
