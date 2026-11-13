# Stability Criterion Validation — Referee-Facing Data Package

This directory contains all data, scripts, and outputs for the validation
of the dynamical stability criterion used in the VibroML manuscript.

## Directory Structure

```
stability_criterion_validation/
├── README.md                    ← this file
├── MANIFEST.md                  ← file inventory with checksums
├── REFEREE_GUIDE.md             ← quick-start guide for referees
├── raw_frequencies/             ← raw phonon frequencies (.npz) per engine
├── processed_metrics/           ← derived metrics (soft_frac, I_curv, etc.)
├── analysis_outputs/            ← ROC/AUC, confusion matrices, calibration
├── scripts/                     ← all scripts to reproduce results
├── figures/                     ← generated figures (PDF + PNG)
├── tables/                      ← generated LaTeX tables
└── provenance/                  ← model versions, settings, git commits
```

## Quick Start

1. **Reproduce all metrics from raw frequencies:**
   ```bash
   cd scripts/
   python analyze_stability_criterion.py --engines mace pet macemp esen uma
   ```

2. **Regenerate all LaTeX tables:**
   ```bash
   cd scripts/
   python generate_tables.py
   ```

3. **Regenerate the main-text figure:**
   ```bash
   cd scripts/
   python generate_main_figure.py --engines mace pet macemp esen uma
   ```

4. **Recompute raw frequencies from scratch (requires GPU):**
   ```bash
   cd scripts/
   sbatch submit_ip_variants_macemp.sh
   sbatch submit_ip_variants_esen.sh
   sbatch submit_ip_variants_uma.sh
   ```

## Key Results

- **Criterion:** soft_frac < 2% (Stable), 2-8% (Quasi-Stable), >= 8% (Unstable)
- **Benchmark:** 24 ground-truth structures (4 unstable, 20 stable)
- **Engines:** MACE-OMAT, MACE-MP, PET-MAD, eSEN, UMA (all 44/44 structures complete)
- **AUC:** 0.95-1.00 for soft_frac across engines
- **Dangerous errors:** 0 in 4 of 5 engines; MACE-MP has 1 false-unstable (Bi₂Sn₂O₇ idx 40, model-specific outlier not reproduced by other engines)
- **LOO-CV accuracy:** 83-92% across engines
- **Q-mesh:** uniform 6×6×6
- **Acoustic handling:** 3 Γ acoustic modes excluded by index when Γ is sampled; the default unshifted even 6×6×6 mesh does not contain Γ, so this exclusion is a no-op for the current runs. No frequency threshold is applied.

## Provenance

See `provenance/` for:
- MLIP model versions and checkpoint paths
- Git commit hashes
- SLURM job IDs
- q-mesh and supercell settings
- Date of computation

## Contact

For questions about this data package, contact the corresponding author.
