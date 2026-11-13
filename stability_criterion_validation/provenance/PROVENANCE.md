# Provenance Record

## Computational Settings

| Parameter | Value |
|-----------|-------|
| Q-mesh | 6×6×6 (unshifted Monkhorst–Pack, no Γ point) |
| Supercell | Per-structure (smallest commensurate with q-mesh) |
| Displacement delta | 0.05 Å |
| Acoustic exclusion | 3 Γ acoustic modes by index — no-op on the default even mesh (Γ is not sampled); retained for Γ-centred or odd meshes |
| Frequency threshold | None (only boundary is ω = 0) |
| Classification thresholds | Stable: f_soft < 2%, Quasi-Stable: 2% ≤ f_soft < 8%, Unstable: f_soft ≥ 8% |

## MLIP Engine Versions

| Engine | Model | Version | Checkpoint |
|--------|-------|---------|------------|
| MACE-OMAT | medium-omat-0 | MACE foundation model | Downloaded via `mace_mp()` |
| MACE-MP | medium | MACE-MP-0 medium | Downloaded via `mace_mp()` |
| PET-MAD | pet-mad-xs | PET-MAD ensemble | Downloaded via `UPETCalculator` |
| eSEN | esen-30m-omat | eSEN-30M-OMAT | `fairchem_models/esen_30m_omat.pt` |
| UMA | uma-m-1p1 | UMA-m-1.1 | `fairchem_models/uma-m-1p1.pt` |

## Geometry Provenance

- **MACE-OMAT**: evaluated the original MACE-OMAT-relaxed structures from the VibroML search.
- **PET-MAD, MACE-MP, eSEN, UMA**: evaluated PET-MAD symmetry-constrained relaxed versions of the same MACE-OMAT-discovered structures. This ensures all engines evaluate the same geometry rather than independently relaxing to different local minima. This is an intentional design choice.

## Benchmark Structures

- Total structures: 44
  - 4 `high_sym_unstable` (known-unstable aristotypes)
  - 20 `literature_stable` (lower-symmetry polymorphs whose space groups are reported for the same compound)
  - 4 `literature_quasi` (literature-discussed phases with borderline stability)
  - 13 `ga_discovered` (structures recovered by VibroML's search)
  - 3 `anomaly` (LiF steric-clash artifacts, excluded from ground-truth)
- Ground-truth subset: 24 (4 known unstable + 20 literature stable)
- Anomalies excluded from GT: idx 9, 10, 11 (LiF steric clash)
- Crosswalk: `structure_crosswalk.csv` maps ip_variants idx ↔ Table S1 row

### Label-basis caveat

The 20 `literature_stable` entries record that their *space group* is a reported polymorph of the compound; the atomic coordinates are MLIP-relaxed, not experimentally refined. No entry has had its dynamical stability confirmed by DFT phonon calculations in this work. The label-basis provenance for every reference structure is given in SI Table `tab:si_label_basis`. Index 6 (labelled LiF $Fm\bar{3}m$) has a metadata inconsistency (actual detected space group $P1$, anomalous energy) and is flagged as Tier X in that table.

## Git State

- VibroML commit: `9295052`
- Branch: `feature/parallel-ga-evaluation`
- Key module: `vibroml/utils/stability_metrics.py`

## SLURM Jobs

| Engine | Job ID | Partition | Status |
|--------|--------|-----------|--------|
| MACE-OMAT | 15399046 | gpu | completed |
| PET-MAD | 15399047 | gpu | completed |
| MACE-MP | 15399048 | gpu | completed |
| eSEN | 15399048 | gpu | completed |
| UMA | 15399048 | gpu | completed |

(Initial jobs 15392081–15392083 failed due to ASE cache-key errors; the fix in `compute_ip_variants.py` — per-structure `tempfile.mkdtemp` — allowed successful resume.)

## Environment Specification

| Environment | Path | Used for |
|-------------|------|----------|
| vibroml_env | `04_envs_and_archives/vibroml_env` | MACE-OMAT, MACE-MP |
| esen_env | `04_envs_and_archives/esen_env` | eSEN |
| uma_env | `04_envs_and_archives/uma_env` | UMA |
| orb_env | `conda_envs/orb_env` | PET-MAD (contains `upet`) |

All environments require CUDA 11.7 and cuDNN 8.4.1.50.

## Known Limitations

1. **No independent DFT phonon ground truth**: DFT verification covers relative formation energies only; no DFT phonon spectra were computed for the reference structures.
2. **Internal benchmark**: The reference labels encode which polymorphs are reported for each compound; the assessment demonstrates internal consistency and cross-architecture reproducibility, not calibration against an independent dynamical-stability ground truth.
3. **Bi₂Sn₂O₇ zero-dispersion issue**: Structures idx 34 and 39 previously used a [1,1,1] supercell fallback, producing q-independent spectra. Corrected reruns at [2,2,2] completed for MACE-OMAT (job 15401136), MACE-MP (15401137), PET-MAD (15401424), UMA (15401139), and eSEN idx=34 (15401425). eSEN idx=39 (1408 atoms) OOMs on GPU (both 40 GB and 80 GB A100s) and was not feasible on CPU within reasonable time; this entry retains the original [1,1,1] value for eSEN only. The corrected data has genuine q-dispersion (band std 0.66–0.94 THz vs 0.000 before). Merged files are `ip_variants_*_merged.json`; originals backed up as `ip_variants_*.json.bak_sc111`.
4. **PET-MAD LiF outlier**: PET-MAD predicts f_soft≈0.81 for the LiF aristotype vs ≈0.29–0.31 for the other four engines. All agree on the Unstable classification; the quantitative discrepancy is noted but not diagnosed.
5. **MACE-MP Bi₂Sn₂O₇ outlier**: MACE-MP assigns idx 40 to Unstable (ω_min = −74 THz) while the other four engines give ω_min ≈ −3.5 THz. This is described as a model-specific outlier, not a diagnosed artifact.

## Date

Generated: 2026-08-30 (last updated 2026-08-31)
