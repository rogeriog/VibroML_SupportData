# LiF Steric-Clash Anomaly Analysis

Linear interpolation pathways and DFT vs. MLIP energy/force evaluations for the
three unphysical short-distance $\ce{LiF}$ configurations ($Cm$, $P1$, and
$Pm$) uncovered by stochastic search (`--method opt_random`), corresponding to
Supplementary Section S5 (Table S6 and Figures S8--S11).

## Contents

- `LiF_Cm/`, `LiF_P1/`, `LiF_Pm/` — per-pathway interpolated image geometries
  and evaluation scripts (`eval_energy_forces.py`) comparing DFT (VASP PBE)
  against MACE-OMAT-0 (medium), MACE-MP, eSEN-30M-OMat, and UMA-m-1p1.
- `neb_image_details_DFT.csv` — image-by-image DFT reference energies and
  maximum forces.
- `final_path_statistics_summary.csv`, `final_path_statistics_detailed.csv`,
  `final_path_statistics_with_errors.csv`, `mlip_error_overview.csv` —
  aggregated barrier and force error statistics (Table S6).
- `get_statistics_path.py` — script that compiles the summary CSV tables.
- `plot_LiF_anomaly_comparison_broken_allMLIPs.py` and `comparison_plots/` —
  plotting script and generated broken-axis comparison figures (Figures S9--S11).
