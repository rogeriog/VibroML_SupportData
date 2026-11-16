# Evidence for the two absent engine–compound rows of Table S12

Table S12 (`tab:b2_phonon_metrics`) lists the seven engine–compound GA
combinations that completed. Two are absent; the files below are the run
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

## `CuN` with UMA — GA run hit the SLURM walltime limit

- `log_uma_ga_fresh.txt` — GA run log.
- `log_cun_uma_resume_debug.txt` — resume log ending with
  `JOB 15330609 ON cna049 CANCELLED AT 2026-08-26T14:28:36 DUE TO TIME
  LIMIT`, while relaxing a 772 Ang^3 distorted supercell (step 169). The
  run never reached final-structure selection, so no UMA-GA lowest-energy
  candidate exists for `CuN`.

For the Table S13 NEB force components: the ZnCu and ZrCu MACE-OMAT
OUTCARs were later recovered from the original nic5 VASP run tree and
all nine paths are now decomposed; see
`neb_results/dft_perp_force_analysis/` (OUTCARs, script, per-image and
per-path CSVs) and the Table S13 caption.
