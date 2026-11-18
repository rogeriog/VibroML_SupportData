# Evidence behind the ZnCu–MACE-OMAT Table S12 row

Table S12 (`tab:b2_phonon_metrics`) lists all nine engine–compound
combinations. The ZnCu–MACE-OMAT entry required a finer initial screen
than the automatic sweep used by default; the files below are the run
artifacts cited in the SI text and the table caption.

## `ZnCu` with MACE-OMAT — automatic sweep terminated at cubic B2

- `b2init_ZnCu_MACE.log` — driver log of the `--auto` supercell sweep.
- `ZnCu_b2_phonon_run_summary.txt` — run summary of the first screened
  supercell (`N2x2x2_D0.05_F0.001`). It reports
  `Most Negative Frequency (overall): 0.0539 THz` and `1296 total,
  0 imaginary` — no imaginary modes anywhere on that grid, because the
  B2 M-point instability requires fourth-nearest-neighbor interactions
  (on `3x3x3` MACE-OMAT gives omega_min = -0.80 THz at M). With no soft
  modes below threshold the automatic workflow exited at the cubic B2
  parent without launching the GA.

The Table S12 row for this combination was therefore produced by a GA
run launched with a 3×3×3 initial screen
(`b2_metrics_ZnCu_MACE_rerun3x3x3/` in this directory).

## `CuN` with UMA — earlier attempt log (superseded)

- `log_uma_ga_fresh.txt` — GA run log of a June attempt.
- `log_cun_uma_resume_debug.txt` — resume log ending with
  `JOB 15330609 ON cna049 CANCELLED AT 2026-08-26T14:28:36 DUE TO TIME
  LIMIT`, while relaxing a 772 Ang^3 distorted supercell (step 169).

The August campaign (`CuN_UMA_GA_phonon_output_20260825-185753`) completed
five iterations and produced the certified 64-atom P1 entry
(`b2_metrics_CuN_UMA_guidance_check/`); these logs are retained only as
run history.

For the Table S13 NEB force components: the ZnCu and ZrCu MACE-OMAT
OUTCARs were recovered from the original nic5 VASP run tree and all nine
paths are decomposed; see `neb_results/dft_perp_force_analysis/`
(OUTCARs, script, per-image and per-path CSVs) and the Table S13 caption.
