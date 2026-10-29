# VibroML Support Data — Repository Manifest

This manifest maps each manuscript display item and analysis section to the
corresponding directories, data files, and scripts in `VibroML_SupportData` and
the companion Figshare deposit (`https://doi.org/10.6084/m9.figshare.34790283`).

## Model Checkpoints Used

| Label in Manuscript | Checkpoint / Package |
|---------------------|----------------------|
| `MACE-OMAT` / `MACE-OMAT-0 (medium)` | `mace-omat-0-medium` (`mace-torch`) |
| `MACE-MP` | `mace-mp-0-medium` (`mace-torch`) |
| `eSEN` | `esen_30m_omat.pt` (`fairchem-core` 1.10) |
| `UMA` | `uma-m-1p1.pt` (`fairchem-core` 2.12) |
| `PET-MAD (xs)` | `pet-mad-xs` (128-member ensemble via `upet`) |

## Display Item & Section Mapping

| Manuscript / SI Item | Repository Path | Key Files / Scripts |
|----------------------|-----------------|---------------------|
| **Table 1** (Benchmark remediation) | `remediation_results/`, `phonon_and_md_benchmarking/` | `remediation_results/<Compound>/`, `vibroml_master_eval_summary.csv` |
| **Table 2, Fig. 4, SI Sec. S15, Fig. S24** (6,574 Alexandria systems) | `whitespace_data/` | `quaternary_ga_summary.csv`, `quintenary_ga_summary.csv`, `ehull_approx_data.csv`, `vibroml_database.db.gz` |
| **Table 3** (200 white-space elpasolites) | `whitespace_data/` | `elpasolite_sample_200.json`, `elpasolite_whitespaces_stats.json`, `vibroml_database.db.gz` |
| **SI Sec. S1, Tables S1--S3, S7--S9** (Phonon metrics & rankings) | `phonon_and_md_benchmarking/` | `MACE_evals/`, `MACEMP_evals/`, `ESEN_evals/`, `UMA_evals/`, `extract_evals_all.py` |
| **SI Sec. S2, Tables S4 & S10** (MD thermal stability) | `phonon_and_md_benchmarking/` | `MACE_MD_evals/`, `MACEMP_MD_evals/`, `ESEN_MD_evals/`, `UMA_MD_evals/`, `extract_md_evals_all.py` |
| **SI Sec. S4, Table S5** (Unique structure counts) | `remediation_results/` | `compile_unique_structures.py`, `unique_structures_summary.txt` |
| **SI Sec. S5, Table S6, Figs. S8--S11** (LiF steric-clash anomaly) | `LiF_anomaly_tests/` | `LiF_Cm/`, `LiF_P1/`, `LiF_Pm/`, `get_statistics_path.py`, `final_path_statistics_summary.csv` |
| **SI Sec. S8, Table S11** (B2 intermetallic phonon metrics) | `b2_intermetallics/` | `b2_metrics_<System>_<Engine>_*/` |
| **SI Sec. S9, Table S12, Figs. S13--S19** (NEB pathways) | `neb_results/` | `*_neb_summary.txt`, `vasp_plots_and_neb_data/`, `plot_comparison_DFT_MLIP_neb_paths.py` |
| **SI Sec. S11, Tables S14--S15** (Remediation settings & timing) | `remediation_results/` | `compile_timing.py`, `vibroml_timing_table.csv` |
| **SI Sec. S16, Tables S17--S18** (HfO2 literature validation) | `hfo2_literature_validation/` | `HfO2_SG*.cif`, `compare_hfo2_phases.py` |
| **83,562 GA local-minima XYZ dataset & Explorer blob store** | Figshare (`10.6084/m9.figshare.34790283`) | `vibroml_ga_local_minima_structures.xyz.gz`, `vibroml_explorer_bundle.tar.gz` |
