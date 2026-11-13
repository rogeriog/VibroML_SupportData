#!/usr/bin/env python3
"""
Analyze PET-MAD UQ results vs MACE.
Focus: justify the -1.5 THz threshold using ensemble uncertainty.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats
from collections import Counter

BASE = Path(__file__).parent
CSV = BASE / "uq_validation_summary.csv"
OUTPUT_DIR = BASE / "uq_analysis_plots"

df = pd.read_csv(CSV)
print(f"Loaded {len(df)} structures")

# Use pet_i_curv_soft_frac as the soft_frac (it's computed on uniform q-mesh, complete)
if 'pet_i_curv_soft_frac' in df.columns:
    df['pet_soft_frac'] = df['pet_i_curv_soft_frac']

# Also use mace_i_curv_soft_frac if mace_soft_frac is missing
if 'mace_i_curv_soft_frac' in df.columns and df['mace_soft_frac'].isna().any():
    df['mace_soft_frac'] = df['mace_soft_frac'].fillna(df['mace_i_curv_soft_frac'])

# Filter to rows with both MACE and PET omega_min
df_valid = df.dropna(subset=['pet_omega_min', 'mace_omega_min']).copy()
print(f"Valid rows (both MACE and PET ω_min): {len(df_valid)}")

# Compute soft_frac if still missing
df_valid = df_valid.dropna(subset=['pet_soft_frac', 'mace_soft_frac'])
print(f"Valid rows (with soft_frac): {len(df_valid)}")

OUTPUT_DIR.mkdir(exist_ok=True)

# ============================================================
# Classification
# ============================================================
mace_stable = (df_valid['mace_omega_min'] > -1.5) & (df_valid['mace_soft_frac'] < 0.02)
mace_unstable = (df_valid['mace_omega_min'] <= -1.5) & (df_valid['mace_soft_frac'] >= 0.02)
mace_quasi = ~(mace_stable | mace_unstable)

pet_stable = (df_valid['pet_omega_min'] > -1.5) & (df_valid['pet_soft_frac'] < 0.02)
pet_unstable = (df_valid['pet_omega_min'] <= -1.5) & (df_valid['pet_soft_frac'] >= 0.02)
pet_quasi = ~(pet_stable | pet_unstable)

mace_labels = np.where(mace_stable, 'Stable', np.where(mace_unstable, 'Unstable', 'Quasi'))
pet_labels = np.where(pet_stable, 'Stable', np.where(pet_unstable, 'Unstable', 'Quasi'))

# ============================================================
# Figure 1: Per-class uncertainty distributions (THE key figure)
# Shows that ensemble σ justifies the discrete criteria
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Panel (a): ω_min colored by MACE classification, with error bars
ax = axes[0, 0]
df_sorted = df_valid.sort_values('mace_omega_min').reset_index(drop=True)
x = range(len(df_sorted))
colors_pt = ['#27ae60' if s else '#f39c12' if q else '#c0392b' 
             for s, q in zip(mace_stable.sort_values(), mace_quasi.sort_values())]
# Actually need to sort the labels too
df_sorted['mace_label'] = mace_labels[df_valid.index.get_indexer(df_sorted.index)]
colors_pt = df_sorted['mace_label'].map({'Stable': '#27ae60', 'Quasi': '#f39c12', 'Unstable': '#c0392b'})

ax.bar(x, df_sorted['pet_omega_min'].values, color=colors_pt, alpha=0.7, edgecolor='none')
ax.errorbar(x, df_sorted['pet_omega_min'].values,
            yerr=3*df_sorted['pet_std_at_min_thz'].values,
            fmt='none', color='black', capsize=1, alpha=0.3, linewidth=0.5)
ax.axhline(y=-1.5, color='red', linestyle='--', linewidth=2, label='Threshold ($-$1.5 THz)')
ax.set_xlabel('Structure (sorted by MACE $\\omega_{min}$)', fontsize=12)
ax.set_ylabel('PET-MAD $\\omega_{min}$ (THz)', fontsize=12)
ax.set_title('(a) PET-MAD $\\omega_{min}$ ± 3σ, colored by MACE class', fontsize=13, fontweight='bold')
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color='#27ae60', label='Stable'),
                   Patch(color='#f39c12', label='Quasi-stable'),
                   Patch(color='#c0392b', label='Unstable'),
                   plt.Line2D([0], [0], color='red', linestyle='--', label='Threshold')], 
          fontsize=9)

# Panel (b): σ@ω_min distribution by class (box plot)
ax = axes[0, 1]
classes = ['Stable', 'Quasi', 'Unstable']
sigma_data = {c: df_valid.loc[df_valid['mace_label'] == c if 'mace_label' in df_valid.columns else 
                              (np.where(mace_stable, 'Stable', np.where(mace_unstable, 'Unstable', 'Quasi')) == c), 
                              'pet_std_at_min_thz'].dropna().values 
              for c in classes}
# Simpler approach
df_valid['mace_class'] = mace_labels
sigma_data = [df_valid.loc[df_valid['mace_class'] == c, 'pet_std_at_min_thz'].dropna().values for c in classes]
bp = ax.boxplot(sigma_data, labels=classes, patch_artist=True, widths=0.6,
                showmeans=True, meanprops={'marker': 'D', 'markerfacecolor': 'white', 'markersize': 6})
for patch, color in zip(bp['boxes'], ['#27ae60', '#f39c12', '#c0392b']):
    patch.set_facecolor(color)
    patch.set_alpha(0.6)
ax.axhline(y=0.5, color='gray', linestyle=':', linewidth=1, label='0.5 THz')
ax.set_ylabel('PET-MAD σ at $\\omega_{min}$ (THz)', fontsize=12)
ax.set_title('(b) Uncertainty by stability class', fontsize=13, fontweight='bold')
ax.legend(fontsize=10)

# Panel (c): ω_min vs σ@ω_min, colored by class
ax = axes[1, 0]
for c, color in zip(classes, ['#27ae60', '#f39c12', '#c0392b']):
    mask = df_valid['mace_class'] == c
    ax.scatter(df_valid.loc[mask, 'pet_omega_min'], df_valid.loc[mask, 'pet_std_at_min_thz'],
               c=color, alpha=0.7, s=60, edgecolors='black', linewidth=0.5, label=c)
ax.axvline(x=-1.5, color='red', linestyle='--', linewidth=1.5, label='Threshold')
ax.set_xlabel('PET-MAD $\\omega_{min}$ (THz)', fontsize=12)
ax.set_ylabel('PET-MAD σ at $\\omega_{min}$ (THz)', fontsize=12)
ax.set_title('(c) Uncertainty vs frequency by class', fontsize=13, fontweight='bold')
ax.legend(fontsize=10)

# Panel (d): Conservative estimate ω_min - 3σ vs threshold
ax = axes[1, 1]
df_3std = df_valid.dropna(subset=['pet_min_minus_3std']).sort_values('pet_min_minus_3std').reset_index(drop=True)
df_3std['mace_class'] = df_3std.index.map(
    {i: mace_labels[df_valid.index.get_indexer([idx])[0]] for i, idx in enumerate(df_3std['idx'])})
# Simpler: just use the class from df_valid
df_3std = df_valid.dropna(subset=['pet_min_minus_3std']).sort_values('pet_min_minus_3std').reset_index(drop=True)
df_3std['class'] = df_3std['mace_class'].values if 'mace_class' in df_3std.columns else 'Stable'
# Actually let's just recompute
df_3std = df_valid.dropna(subset=['pet_min_minus_3std']).copy()
df_3std = df_3std.sort_values('pet_min_minus_3std').reset_index(drop=True)
colors_bar = df_3std['mace_class'].map({'Stable': '#27ae60', 'Quasi': '#f39c12', 'Unstable': '#c0392b'})
x = range(len(df_3std))
ax.bar(x, df_3std['pet_min_minus_3std'].values, color=colors_bar, alpha=0.7, edgecolor='none')
ax.axhline(y=-1.5, color='red', linestyle='--', linewidth=2, label='Threshold ($-$1.5 THz)')
ax.set_xlabel('Structure (sorted by $\\omega_{min}$ - 3σ)', fontsize=12)
ax.set_ylabel('$\\omega_{min}$ - 3σ (THz)', fontsize=12)
ax.set_title('(d) Conservative stability estimate', fontsize=13, fontweight='bold')
ax.legend(handles=[Patch(color='#27ae60', label='Stable'),
                   Patch(color='#f39c12', label='Quasi-stable'),
                   Patch(color='#c0392b', label='Unstable'),
                   plt.Line2D([0], [0], color='red', linestyle='--', label='Threshold')],
          fontsize=9)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'uq_threshold_validation.png', dpi=300, bbox_inches='tight')
plt.savefig(OUTPUT_DIR / 'uq_threshold_validation.svg', bbox_inches='tight')
plt.close()
print(f"Saved: {OUTPUT_DIR / 'uq_threshold_validation.png'}")

# ============================================================
# Figure 2: I_curv as continuous proxy + correlation
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Panel (a): ω_min correlation (secondary, just for reference)
ax = axes[0]
df_plot = df_valid.dropna(subset=['mace_omega_min', 'pet_omega_min'])
colors = df_plot['mace_class'].map({'Stable': '#27ae60', 'Quasi': '#f39c12', 'Unstable': '#c0392b'})
ax.scatter(df_plot['mace_omega_min'], df_plot['pet_omega_min'],
           c=colors, alpha=0.6, s=50, edgecolors='black', linewidth=0.5)
lim = max(abs(df_plot['mace_omega_min'].max()), abs(df_plot['pet_omega_min'].max())) * 1.1
ax.plot([-lim, lim], [-lim, lim], 'k--', linewidth=1, label='y = x')
ax.axvline(x=-1.5, color='red', linestyle=':', linewidth=1.5)
ax.axhline(y=-1.5, color='red', linestyle=':', linewidth=1.5)
if len(df_plot) > 2:
    slope, intercept, r_value, p_value, _ = stats.linregress(
        df_plot['mace_omega_min'], df_plot['pet_omega_min'])
    ax.plot([-lim, lim], [slope*(-lim)+intercept, slope*lim+intercept], 'g-',
            linewidth=1.5, label=f'Fit (R²={r_value**2:.3f})')
ax.set_xlabel('MACE $\\omega_{min}$ (THz)', fontsize=13)
ax.set_ylabel('PET-MAD $\\omega_{min}$ (THz)', fontsize=13)
ax.set_title('$\\omega_{min}$: MACE vs PET-MAD', fontsize=14, fontweight='bold')
ax.legend(fontsize=9, handles=[plt.Line2D([0],[0],marker='o',color='w',markerfacecolor='#27ae60',markersize=8,label='Stable'),
                               plt.Line2D([0],[0],marker='o',color='w',markerfacecolor='#f39c12',markersize=8,label='Quasi'),
                               plt.Line2D([0],[0],marker='o',color='w',markerfacecolor='#c0392b',markersize=8,label='Unstable'),
                               plt.Line2D([0],[0],color='k',linestyle='--',label='y=x'),
                               plt.Line2D([0],[0],color='green',linestyle='-',label=f'Fit (R²={r_value**2:.3f})')])

# Panel (b): I_curv as continuous proxy — show it captures both ω_min and soft_frac
ax = axes[1]
df_curv = df_valid.dropna(subset=['mace_i_curv', 'pet_i_curv'])
if len(df_curv) > 0:
    # Plot I_curv vs ω_min for PET-MAD, colored by soft_frac
    sc = ax.scatter(df_curv['pet_i_curv'], df_curv['pet_omega_min'],
                    c=df_curv['pet_soft_frac'], cmap='YlOrRd', s=60,
                    edgecolors='black', linewidth=0.5, vmin=0, vmax=0.3)
    ax.axhline(y=-1.5, color='red', linestyle='--', linewidth=1.5, label='$\\omega_{min}$ threshold')
    ax.set_xlabel('PET-MAD $I_{curv}$ (THz²)', fontsize=13)
    ax.set_ylabel('PET-MAD $\\omega_{min}$ (THz)', fontsize=13)
    ax.set_title('$I_{curv}$ vs $\\omega_{min}$ (color = soft frac)', fontsize=14, fontweight='bold')
    plt.colorbar(sc, ax=ax, label='Soft mode fraction')
    ax.legend(fontsize=10)
    ax.set_xscale('log')
    ax.set_ylim(-20, 2)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'uq_correlation.png', dpi=300, bbox_inches='tight')
plt.savefig(OUTPUT_DIR / 'uq_correlation.svg', bbox_inches='tight')
plt.close()
print(f"Saved: {OUTPUT_DIR / 'uq_correlation.png'}")

# ============================================================
# Print summary statistics
# ============================================================
print(f"\n{'='*60}")
print(f"UQ Summary Statistics ({len(df_valid)} structures)")
print(f"{'='*60}")
print(f"\nClassification (MACE criteria):")
print(f"  Stable:     {mace_stable.sum()}")
print(f"  Quasi:      {mace_quasi.sum()}")
print(f"  Unstable:   {mace_unstable.sum()}")
print(f"\nClassification (PET-MAD criteria):")
print(f"  Stable:     {pet_stable.sum()}")
print(f"  Quasi:      {pet_quasi.sum()}")
print(f"  Unstable:   {pet_unstable.sum()}")

agreement = sum(mace_labels == pet_labels)
print(f"\nAgreement: {agreement}/{len(df_valid)} ({agreement/len(df_valid)*100:.1f}%)")

for cat in df_valid['category'].unique():
    mask = df_valid['category'] == cat
    cat_agree = sum(mace_labels[mask] == pet_labels[mask])
    print(f"  {cat}: {cat_agree}/{mask.sum()} ({cat_agree/mask.sum()*100:.1f}%)")

print(f"\nCorrelations:")
corr_omega = stats.pearsonr(df_valid['mace_omega_min'], df_valid['pet_omega_min'])
print(f"  ω_min Pearson:  r={corr_omega[0]:.4f}, p={corr_omega[1]:.4e}")
spearman_omega = stats.spearmanr(df_valid['mace_omega_min'], df_valid['pet_omega_min'])
print(f"  ω_min Spearman: ρ={spearman_omega[0]:.4f}, p={spearman_omega[1]:.4e}")

df_curv = df_valid.dropna(subset=['mace_i_curv', 'pet_i_curv'])
if len(df_curv) > 2:
    corr_curv = stats.pearsonr(df_curv['mace_i_curv'], df_curv['pet_i_curv'])
    print(f"  I_curv Pearson:  r={corr_curv[0]:.4f}, p={corr_curv[1]:.4e}")
    spearman_curv = stats.spearmanr(df_curv['mace_i_curv'], df_curv['pet_i_curv'])
    print(f"  I_curv Spearman: ρ={spearman_curv[0]:.4f}, p={spearman_curv[1]:.4e}")

print(f"\nPET-MAD Uncertainty by class:")
for c in classes:
    mask = df_valid['mace_class'] == c
    if mask.sum() > 0:
        sigmas = df_valid.loc[mask, 'pet_std_at_min_thz']
        print(f"  {c:10s}: σ@ω_min = {sigmas.mean():.3f} ± {sigmas.std():.3f} (mean ± std), "
              f"range [{sigmas.min():.3f}, {sigmas.max():.3f}], n={mask.sum()}")

print(f"\nThreshold Robustness:")
borderline = ((df_valid['pet_omega_min'] <= -1.5) & (df_valid['pet_min_minus_3std'] > -1.5)).sum()
reverse = ((df_valid['pet_omega_min'] > -1.5) & (df_valid['pet_min_minus_3std'] <= -1.5)).sum()
print(f"  Borderline (ω_min ≤ -1.5 but ω_min-3σ > -1.5): {borderline}")
print(f"  Reverse (ω_min > -1.5 but ω_min-3σ ≤ -1.5): {reverse}")

# Per-class analysis of the -1.5 THz threshold
print(f"\n{'='*60}")
print(f"THRESHOLD ANALYSIS: Is -1.5 THz justified?")
print(f"{'='*60}")
for c in classes:
    mask = df_valid['mace_class'] == c
    if mask.sum() == 0:
        continue
    sub = df_valid[mask]
    pet_omegas = sub['pet_omega_min']
    pet_sigmas = sub['pet_std_at_min_thz']
    pet_3std = sub['pet_min_minus_3std']
    
    print(f"\n  {c} (n={mask.sum()}):")
    print(f"    PET ω_min:     [{pet_omegas.min():.2f}, {pet_omegas.max():.2f}] THz")
    print(f"    PET σ@ω_min:   [{pet_sigmas.min():.2f}, {pet_sigmas.max():.2f}] THz")
    print(f"    PET ω_min-3σ:  [{pet_3std.min():.2f}, {pet_3std.max():.2f}] THz")
    
    # How many in this class have ω_min > -1.5?
    above = (pet_omegas > -1.5).sum()
    below = (pet_omegas <= -1.5).sum()
    print(f"    Above -1.5 THz: {above}/{mask.sum()}")
    print(f"    Below -1.5 THz: {below}/{mask.sum()}")
    
    # How many have ω_min-3σ > -1.5?
    above_3std = (pet_3std > -1.5).sum()
    below_3std = (pet_3std <= -1.5).sum()
    print(f"    ω_min-3σ > -1.5:  {above_3std}/{mask.sum()}")
    print(f"    ω_min-3σ ≤ -1.5:  {below_3std}/{mask.sum()}")

# Check what happens with different thresholds
print(f"\n{'='*60}")
print(f"THRESHOLD SENSITIVITY ANALYSIS")
print(f"{'='*60}")
print(f"\n  Threshold | Stable | Quasi | Unstable | Borderline(ω≤th, ω-3σ>th)")
print(f"  ----------|--------|-------|----------|---------------------------")
for thresh in [-0.5, -1.0, -1.5, -2.0, -3.0]:
    mace_s = (df_valid['mace_omega_min'] > thresh) & (df_valid['mace_soft_frac'] < 0.02)
    mace_u = (df_valid['mace_omega_min'] <= thresh) & (df_valid['mace_soft_frac'] >= 0.02)
    mace_q = ~(mace_s | mace_u)
    
    pet_s = (df_valid['pet_omega_min'] > thresh) & (df_valid['pet_soft_frac'] < 0.02)
    pet_u = (df_valid['pet_omega_min'] <= thresh) & (df_valid['pet_soft_frac'] >= 0.02)
    pet_q = ~(pet_s | pet_u)
    
    borderline = ((df_valid['pet_omega_min'] <= thresh) & (df_valid['pet_min_minus_3std'] > thresh)).sum()
    
    agree = sum(np.where(mace_s, 'S', np.where(mace_u, 'U', 'Q')) == 
                np.where(pet_s, 'S', np.where(pet_u, 'U', 'Q')))
    
    print(f"  {thresh:7.1f}   | {pet_s.sum():6d} | {pet_q.sum():5d} | {pet_u.sum():8d} | {borderline:25d}  agree={agree}/{len(df_valid)}")
