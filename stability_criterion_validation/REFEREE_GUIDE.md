# Referee Quick-Start Guide

## What to look at first

1. **`processed_metrics/stability_metrics_all_engines.csv`** — the master table
   with soft_frac, I_curv, omega_min, and classification for every structure
   and every engine.

2. **`analysis_outputs/roc_auc_per_engine.csv`** — AUC for each metric and engine.

3. **`analysis_outputs/three_tier_confusion.csv`** — confusion matrix showing
   zero dangerous errors in 4 of 5 engines (MACE-MP has 1 model-specific outlier).

4. **`figures/figure_criterion_validation.pdf`** — the main-text figure.

## How to verify the criterion

The three-tier rule is:
- Stable: soft_frac < 0.02
- Quasi-Stable: 0.02 ≤ soft_frac < 0.08
- Unstable: soft_frac ≥ 0.08

Check that:
1. All 4 known-unstable structures have soft_frac ≥ 0.08 (or are in Quasi-Stable)
2. All 20 literature-stable structures have soft_frac < 0.08 (or are in Quasi-Stable)
3. No stable structure is classified as Unstable
4. No unstable structure is classified as Stable

```bash
python -c "
import pandas as pd
df = pd.read_csv('processed_metrics/stability_metrics_all_engines.csv')
for eng in df['engine'].unique():
    sub = df[df['engine'] == eng]
    dangerous = ((sub['ground_truth']=='stable') & (sub['classification']=='unstable')).sum() + \
                ((sub['ground_truth']=='unstable') & (sub['classification']=='stable')).sum()
    print(f'{eng}: {len(sub)} structures, {dangerous} dangerous errors')
"
```

## How to check the I_curv computation

The I_curv formula is:
    I_curv = (1/M) * sum_{omega<0} omega^2   [THz^2]

where M = total modes (Γ acoustic excluded when Γ is sampled; the default unshifted even 6×6×6 mesh does not contain Γ, so no exclusion is applied for the current runs).

```bash
python -c "
import numpy as np
d = np.load('raw_frequencies/ip_variants_freqs_mace/idx_001_LiF_freqs.npz', allow_pickle=True)
freqs = d['freqs']
gamma = d['gamma_indices']
# Zero out acoustic at Gamma
for gi in gamma:
    freqs[gi, :3] = 0.0
flat = freqs.ravel()
M = len(flat)
neg = flat[flat < 0]
i_curv = (neg**2).sum() / M
print(f'LiF Pm-3m: M={M}, n_imag={len(neg)}, I_curv={i_curv:.6f} THz^2')
print(f'soft_frac={len(neg)/M:.6f}')
"
```

## How to reproduce from scratch

See README.md for full instructions. The key scripts are:
1. `scripts/compute_ip_variants.py` — computes raw frequencies
2. `scripts/analyze_stability_criterion.py` — computes all metrics and analysis
3. `scripts/generate_tables.py` — generates LaTeX tables
4. `scripts/generate_main_figure.py` — generates the main-text figure
