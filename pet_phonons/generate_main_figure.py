#!/usr/bin/env python3
"""
Generate the main-text criterion validation figure (4 panels).

Panel A: soft_frac distributions for stable vs unstable (violin/box)
Panel B: ROC curves for soft_frac, I_curv, omega_min (best engine)
Panel C: I_curv vs soft_frac scatter (colored by ground truth)
Panel D: Calibrated probability vs soft_frac

Usage:
  python generate_main_figure.py
  python generate_main_figure.py --engines mace pet --output figure_criterion_validation.pdf
"""
import argparse, os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneOut
from sklearn.preprocessing import StandardScaler

PET_PHONONS_DIR = "/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/02_runs_and_data/pet_phonons"
ANOMALY_IDXS = [9, 10, 11]


def load_engine_data(engine_label):
    import json
    path = os.path.join(PET_PHONONS_DIR, f"ip_variants_{engine_label}.json")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


def to_dataframe(data):
    rows = []
    for r in data:
        rows.append({
            'idx': r['idx'],
            'compound': r['compound'],
            'label': r.get('stability_label', 'unknown'),
            'omega_min': r.get('omega_min'),
            'soft_frac': r.get('soft_frac_no_thresh'),
            'i_curv': r.get('i_p_1.0'),
            'omega_inst': r.get('omega_p_1.0'),
        })
    return pd.DataFrame(rows)


def get_ground_truth(df):
    mask = df['label'].isin(['high_sym_unstable', 'literature_stable'])
    mask &= ~df['idx'].isin(ANOMALY_IDXS)
    result = df[mask].copy()
    result['unstable'] = (result['label'] == 'high_sym_unstable').astype(int)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--engines", nargs='+', default=['mace', 'pet'])
    parser.add_argument("--output", default=None)
    parser.add_argument("--fig-dir", default=None)
    args = parser.parse_args()

    # Load data
    engine_dfs = {}
    for eng in args.engines:
        data = load_engine_data(eng)
        if data is None:
            continue
        engine_dfs[eng] = to_dataframe(data)

    if not engine_dfs:
        print("No data found")
        sys.exit(1)

    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # ---- Panel A: soft_frac distributions ----
    ax = axes[0, 0]
    all_stable = []
    all_unstable = []
    for eng, df in engine_dfs.items():
        df_gt = get_ground_truth(df)
        stable = df_gt[df_gt['unstable'] == 0]['soft_frac'].values
        unstable = df_gt[df_gt['unstable'] == 1]['soft_frac'].values
        all_stable.extend(stable)
        all_unstable.extend(unstable)

    # Convert to percentages
    all_stable_pct = [s * 100 for s in all_stable]
    all_unstable_pct = [s * 100 for s in all_unstable]

    positions = [1, 2]
    bp = ax.boxplot([all_stable_pct, all_unstable_pct], positions=positions,
                    widths=0.6, patch_artist=True, showfliers=True)
    bp['boxes'][0].set_facecolor('lightblue')
    bp['boxes'][1].set_facecolor('lightcoral')
    ax.set_xticks(positions)
    ax.set_xticklabels(['Stable\n(literature)', 'Unstable\n(known)'], fontsize=13)
    ax.set_ylabel(r'$f_\mathrm{soft}$ (%)', fontsize=14)
    ax.set_title('(A) Soft-mode fraction distributions', fontsize=15)
    ax.axhline(4.5, color='gray', ls='--', lw=1.5, label=r'$f_\mathrm{soft} \approx 4.5\%$ (Youden $J$)')
    ax.legend(fontsize=12, loc='upper left')
    ax.tick_params(labelsize=12)
    ax.set_ylim(-1, max(max(all_unstable_pct) + 2, 35))

    # ---- Panel B: ROC curves ----
    ax = axes[0, 1]
    # Use the first engine for ROC curves
    eng = list(engine_dfs.keys())[0]
    df_gt = get_ground_truth(engine_dfs[eng])
    y = df_gt['unstable'].values

    metrics_to_plot = [
        ('soft_frac', r'$f_\mathrm{soft}$', 'tab:blue'),
        ('i_curv', r'$I_\mathrm{curv}$', 'tab:orange'),
        ('omega_min', r'$\omega_\mathrm{min}$ (negated)', 'tab:red'),
    ]

    for col, label, color in metrics_to_plot:
        if col == 'omega_min':
            scores = -df_gt[col].values
        else:
            scores = df_gt[col].values
        if len(np.unique(y)) < 2:
            continue
        fpr, tpr, _ = roc_curve(y, scores)
        auc = roc_auc_score(y, scores)
        ax.plot(fpr, tpr, label=f'{label} (AUC={auc:.3f})', color=color, lw=2.5)

    ax.plot([0, 1], [0, 1], 'k--', lw=0.5)
    ax.set_xlabel('False Positive Rate', fontsize=14)
    ax.set_ylabel('True Positive Rate', fontsize=14)
    ax.set_title(f'(B) ROC curves ({eng})', fontsize=15)
    ax.legend(fontsize=12)
    ax.set_xlim(-0.01, 1.01)
    ax.set_ylim(-0.01, 1.01)
    ax.tick_params(labelsize=12)

    # ---- Panel C: I_curv vs soft_frac scatter ----
    ax = axes[1, 0]
    colors = {0: 'tab:blue', 1: 'tab:red'}
    markers = {'mace': 'o', 'pet': 's', 'macemp': 'D', 'esen': '^', 'uma': 'v'}

    for eng, df in engine_dfs.items():
        df_gt = get_ground_truth(df)
        marker = markers.get(eng, 'o')
        for gt_val in [0, 1]:
            subset = df_gt[df_gt['unstable'] == gt_val]
            ax.scatter(subset['soft_frac'] * 100, subset['i_curv'],
                      c=colors[gt_val], marker=marker, alpha=0.6, s=60,
                      label=f'{eng}: {"unstable" if gt_val else "stable"}')

    ax.axvline(4.5, color='gray', ls='--', lw=1.5, alpha=0.7)
    ax.set_xlabel(r'$f_\mathrm{soft}$ (%)', fontsize=14)
    ax.set_ylabel(r'$I_\mathrm{curv}$ (THz$^2$)', fontsize=14)
    ax.set_title('(C) Continuous severity vs categorical metric', fontsize=15)
    ax.set_yscale('symlog', linthresh=0.01)
    ax.legend(fontsize=10, loc='upper left', ncol=2)
    ax.tick_params(labelsize=12)

    # ---- Panel D: Calibrated probability ----
    ax = axes[1, 1]
    eng = list(engine_dfs.keys())[0]
    df_gt = get_ground_truth(engine_dfs[eng])
    y = df_gt['unstable'].values
    X = df_gt[['soft_frac']].values
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    loo = LeaveOneOut()
    probs = np.zeros(len(y))
    for train_idx, test_idx in loo.split(X_scaled):
        X_train, y_train = X_scaled[train_idx], y[train_idx]
        X_test = X_scaled[test_idx]
        if len(np.unique(y_train)) < 2:
            probs[test_idx] = 0.5
            continue
        model = LogisticRegression(C=1.0, solver='lbfgs')
        model.fit(X_train, y_train)
        probs[test_idx] = model.predict_proba(X_test)[0, 1]

    ax.scatter(df_gt[df_gt['unstable'] == 0]['soft_frac'] * 100,
              probs[y == 0], c='tab:blue', marker='o', alpha=0.7, s=70,
              label='Stable (literature)')
    ax.scatter(df_gt[df_gt['unstable'] == 1]['soft_frac'] * 100,
              probs[y == 1], c='tab:red', marker='s', alpha=0.7, s=70,
              label='Unstable (known)')
    ax.axhline(0.5, color='k', ls='--', lw=0.5)
    ax.axvline(4.5, color='gray', ls='--', lw=1.5, alpha=0.7)
    ax.set_xlabel(r'$f_\mathrm{soft}$ (%)', fontsize=14)
    ax.set_ylabel(r'$P(\mathrm{unstable}|f_\mathrm{soft})$ (LOO-CV)', fontsize=14)
    ax.set_title(f'(D) Calibrated probability ({eng}, LOO-CV)', fontsize=15)
    ax.legend(fontsize=12)
    ax.set_ylim(-0.05, 1.05)
    ax.tick_params(labelsize=12)

    plt.tight_layout()

    # Save
    output_dir = args.fig_dir or os.path.join(PET_PHONONS_DIR, "figures")
    os.makedirs(output_dir, exist_ok=True)
    output = args.output or "figure_criterion_validation"
    if not output.endswith('.pdf'):
        output += '.pdf'
    out_path = os.path.join(output_dir, output)
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    fig.savefig(out_path.replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
    print(f"Saved: {out_path}")
    print(f"Saved: {out_path.replace('.pdf', '.png')}")
    plt.close(fig)


if __name__ == "__main__":
    main()
