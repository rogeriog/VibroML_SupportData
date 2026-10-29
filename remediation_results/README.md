# Remediation Algorithm Benchmarks

Outputs from the four dynamical-instability remediation workflows—Genetic
Algorithm (`GA`), Single Soft-Mode Following (`TRAD`), Coupled Soft-Mode
Following (`ALL`), and Stochastic Search (`RAND`)—evaluated with MACE-OMAT-0
(medium) on $\ce{LiF}$, $\ce{CsPbI3}$, $\ce{HfO2}$, and $\ce{Bi2Sn2O7}$.
Corresponds to Table 1 and Supplementary Tables S5, S14, and S15.

## Contents

- `LiF_simple_cubic/`, `CsPbI3/`, `HfO2/`, `Bi2Sn2O7/` — per-compound
  directories containing the run outputs for each remediation algorithm,
  including `final_structures_recovered/`, `overall_relaxation_summary.txt`,
  and `initial_settings.json`.
- `compile_unique_structures.py` and `unique_structures_summary.txt` — script
  and summary table counting unique polymorphs (>0.2 meV/atom energy threshold)
  recovered by each algorithm (Supplementary Table S5).
- `compile_timing.py` and `vibroml_timing_table.csv` — script and summary CSV
  reporting single-GPU execution times across algorithms and compounds
  (Supplementary Table S15).
