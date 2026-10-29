# NEB Phase-Transition Pathway Benchmarks

Nudged elastic band (NEB) minimum-energy pathways driven by MACE-OMAT-0
(medium), eSEN-30M-OMat, and UMA-m-1p1 and verified by single-point DFT (VASP
PBE) calculations across the benchmark systems, corresponding to main-text
Section "Evaluating Phase Transition Pathways" and Supplementary Section S9
(Figures S13--S19).

## Contents

- `<Compound>_<ENGINE>_neb_summary.txt` — per-path summary metrics (reaction
  energy error, barrier error, path average energy, relative-energy MAE, force
  MAE, and endpoint residual forces) for $\ce{LiF}$, $\ce{CsPbI3}$,
  $\ce{HfO2}$, and $\ce{Bi2Sn2O7}$.
- `all_neb_images/` — relaxed NEB image structures (`.cif`) along each pathway.
- `vasp_plots_and_neb_data/` — per-image MLIP and DFT energies and maximum
  forces (`*_NEB_plot_data.csv`) and profile plots (`*_NEB_profile.png`).
- `comparison_plots/` — combined side-by-side energy and force comparison
  panels (`*_comparison.png`).
- `plot_comparison_DFT_MLIP_neb_paths.py` — script used to generate the
  comparison figures from the CSV data.
