#!/usr/bin/env python3
"""
Generate all SI figures for the stability criterion validation section.

Produces:
  1. fig_si_per_engine_roc.pdf/png  — Per-engine ROC curves for soft_frac
  2. fig_si_auc_vs_p.pdf/png        — AUC vs p for all 5 engines
  3. fig_si_marginal_mode.pdf/png   — Marginal-mode contribution bar chart
  4. fig_si_pet_uq_dispersion.pdf/png — PET-MAD UQ dispersion for 4 aristotypes

Usage:
  python generate_si_figures.py [--fig-dir figures_5eng] [--data-dir .]
"""

import argparse, json, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import csv

PET_PHONONS_DIR = os.path.dirname(os.path.abspath(__file__))

ENGINE_COLORS = {
    'mace':   '#1f77b4',
    'pet':    '#ff7f0e',
    'macemp': '#2ca02c',
    'esen':   '#d62728',
    'uma':    '#9467bd',
}
ENGINE_LABELS = {
    'mace':   'MACE-OMAT',
    'pet':    'PET-MAD',
    'macemp': 'MACE-MP',
    'esen':   'eSEN',
    'uma':    'UMA',
}


def load_engine_data(engine, data_dir):
    path = os.path.join(data_dir, f"ip_variants_{engine}.json")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


def get_ground_truth(data):
    """Return (labels, scores) for high_sym_unstable vs literature_stable."""
    labels = []
    scores = []
    for r in data:
        lab = r.get('stability_label', '')
        if lab in ('high_sym_unstable', 'literature_stable'):
            labels.append(1 if lab == 'high_sym_unstable' else 0)
            scores.append(r.get('soft_frac_no_thresh', 0))
    return np.array(labels), np.array(scores)


def compute_roc(labels, scores):
    """Compute ROC curve. Returns (fpr, tpr, auc)."""
    from sklearn.metrics import roc_curve, auc
    fpr, tpr, _ = roc_curve(labels, scores)
    a = auc(fpr, tpr)
    return fpr, tpr, a


def plot_per_engine_roc(engines, data_dir, fig_dir):
    """SI Figure: Per-engine ROC curves for soft_frac."""
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    # Panel A: All engines on one plot
    ax = axes[0]
    for eng in engines:
        data = load_engine_data(eng, data_dir)
        if data is None:
            continue
        labels, scores = get_ground_truth(data)
        if len(np.unique(labels)) < 2:
            continue
        fpr, tpr, a = compute_roc(labels, scores)
        ax.plot(fpr, tpr, color=ENGINE_COLORS[eng], label=f'{ENGINE_LABELS[eng]} (AUC={a:.3f})', linewidth=1.5)
    ax.plot([0, 1], [0, 1], 'k--', alpha=0.3)
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.set_title('(A) ROC: $f_\\mathrm{soft}$ (all engines)')
    ax.legend(fontsize=7, loc='lower right')
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)

    # Panel B: ROC for different metrics (MACE-OMAT as representative)
    ax = axes[1]
    eng = 'mace'
    data = load_engine_data(eng, data_dir)
    if data:
        for metric, label, color in [
            ('soft_frac_no_thresh', '$f_\\mathrm{soft}$', '#1f77b4'),
            ('i_p_1.0', '$I_\\mathrm{curv}$', '#ff7f0e'),
            ('omega_min', '$\\omega_\\mathrm{min}$', '#d62728'),
            ('omega_p_0.5', '$\\Omega_{0.5}$', '#2ca02c'),
        ]:
            labels_list = []
            scores_list = []
            for r in data:
                lab = r.get('stability_label', '')
                if lab in ('high_sym_unstable', 'literature_stable'):
                    labels_list.append(1 if lab == 'high_sym_unstable' else 0)
                    val = r.get(metric, 0)
                    if metric == 'omega_min':
                        val = -val  # negate so higher = more unstable
                    scores_list.append(val)
            if len(np.unique(labels_list)) < 2:
                continue
            fpr, tpr, a = compute_roc(np.array(labels_list), np.array(scores_list))
            ax.plot(fpr, tpr, color=color, label=f'{label} (AUC={a:.3f})', linewidth=1.5)
    ax.plot([0, 1], [0, 1], 'k--', alpha=0.3)
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.set_title(f'(B) ROC: metrics ({ENGINE_LABELS[eng]})')
    ax.legend(fontsize=7, loc='lower right')
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)

    fig.tight_layout()
    for ext in ['pdf', 'png']:
        path = os.path.join(fig_dir, f'fig_si_per_engine_roc.{ext}')
        fig.savefig(path, dpi=300, bbox_inches='tight')
        print(f"Saved: {path}")
    plt.close(fig)


def plot_auc_vs_p(engines, data_dir, fig_dir):
    """SI Figure: AUC vs p for all engines."""
    csv_path = os.path.join(data_dir, 'analysis_outputs_5eng', 'p_sweep.csv')
    if not os.path.exists(csv_path):
        print(f"WARNING: p_sweep.csv not found at {csv_path}")
        return

    # Read p-sweep data
    rows = []
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)

    # Group by engine
    fig, ax = plt.subplots(figsize=(6, 5))
    p_values = sorted(set(float(r['p']) for r in rows if r['p'] not in ('inf', 'soft_frac')))

    for eng in engines:
        eng_data = [(float(r['p']), float(r['auc'])) for r in rows if r['engine'] == eng and r['p'] not in ('inf', 'soft_frac')]
        eng_data.sort()
        if eng_data:
            ps = [d[0] for d in eng_data]
            aucs = [d[1] for d in eng_data]
            ax.plot(ps, aucs, 'o-', color=ENGINE_COLORS[eng], label=ENGINE_LABELS[eng],
                    markersize=4, linewidth=1.5)

    # Add soft_frac as horizontal reference
    for eng in engines:
        sf_rows = [r for r in rows if r['engine'] == eng and r['p'] == 'soft_frac']
        if sf_rows:
            sf_auc = float(sf_rows[0]['auc'])
            ax.axhline(y=sf_auc, color=ENGINE_COLORS[eng], linestyle='--', alpha=0.3)

    ax.set_xlabel('Exponent $p$ in $I_p = [\\frac{1}{M}\\sum|\\omega_i|^{2p}]^{1/p}$')
    ax.set_ylabel('AUC')
    ax.set_title('AUC vs. depth-weighting exponent $p$')
    ax.legend(fontsize=8)
    ax.set_ylim(0.5, 1.05)
    ax.grid(True, alpha=0.3)
    ax.text(0.95, 0.05, 'Dashed lines: $f_\\mathrm{soft}$ AUC', transform=ax.transAxes,
            ha='right', fontsize=7, color='gray')

    fig.tight_layout()
    for ext in ['pdf', 'png']:
        path = os.path.join(fig_dir, f'fig_si_auc_vs_p.{ext}')
        fig.savefig(path, dpi=300, bbox_inches='tight')
        print(f"Saved: {path}")
    plt.close(fig)


def plot_marginal_mode(engines, data_dir, fig_dir):
    """SI Figure: Marginal-mode contribution bar chart."""
    fig, ax = plt.subplots(figsize=(8, 5))

    x = np.arange(len(engines))
    width = 0.25

    pct_marginal = []
    pct_i05 = []
    pct_i10 = []

    for eng in engines:
        csv_path = os.path.join(data_dir, 'analysis_outputs_5eng', f'marginal_mode_sensitivity_{eng}.csv')
        if os.path.exists(csv_path):
            with open(csv_path) as f:
                reader = csv.DictReader(f)
                rows = list(reader)
            if rows:
                # Average across all structures (convert fractions to percentages)
                pm = np.mean([float(r.get('frac_marginal_modes', 0)) for r in rows]) * 100
                pi05 = np.mean([float(r.get('frac_contrib_p05', 0)) for r in rows]) * 100
                pi10 = np.mean([float(r.get('frac_contrib_p10', 0)) for r in rows]) * 100
                pct_marginal.append(pm)
                pct_i05.append(pi05)
                pct_i10.append(pi10)
            else:
                pct_marginal.append(0)
                pct_i05.append(0)
                pct_i10.append(0)
        else:
            pct_marginal.append(0)
            pct_i05.append(0)
            pct_i10.append(0)

    bars1 = ax.bar(x - width, pct_marginal, width, label='% marginal modes', color='#1f77b4', alpha=0.8)
    bars2 = ax.bar(x, pct_i05, width, label='% contrib to $I_{0.5}$', color='#ff7f0e', alpha=0.8)
    bars3 = ax.bar(x + width, pct_i10, width, label='% contrib to $I_\\mathrm{curv}$', color='#2ca02c', alpha=0.8)

    ax.set_xlabel('Engine')
    ax.set_ylabel('Percentage (%)')
    ax.set_title('Marginal-mode sensitivity ($|\\omega| < 0.5$ THz)')
    ax.set_xticks(x)
    ax.set_xticklabels([ENGINE_LABELS[e] for e in engines], rotation=20, ha='right')
    ax.legend(fontsize=8)
    ax.grid(axis='y', alpha=0.3)

    # Add value labels
    for bars in [bars1, bars2, bars3]:
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax.annotate(f'{height:.0f}%',
                            xy=(bar.get_x() + bar.get_width() / 2, height),
                            xytext=(0, 3), textcoords="offset points",
                            ha='center', va='bottom', fontsize=6)

    fig.tight_layout()
    for ext in ['pdf', 'png']:
        path = os.path.join(fig_dir, f'fig_si_marginal_mode.{ext}')
        fig.savefig(path, dpi=300, bbox_inches='tight')
        print(f"Saved: {path}")
    plt.close(fig)


def plot_pet_uq_dispersion(data_dir, fig_dir):
    """SI Figure: PET-MAD UQ dispersion for 4 aristotypes."""
    aristotypes = {
        'LiF': 'uq_validation_runs/001_LiF_benc_LiFsimplecubic/LiFsimplecubic_symrelaxed_PET_phonon_output_20260827-174404/LiFsimplecubic_symrelaxed_phonon_uq_data.npz',
        'CsPbI3': 'uq_validation_runs/106_CsPbI3_benc_CsPbI3_relaxed/CsPbI3_relaxed_symrelaxed_PET_phonon_output_20260829-223622/CsPbI3_relaxed_symrelaxed_phonon_uq_data.npz',
        'HfO2': 'uq_validation_runs/022_HfO2_benc_HfO2_Fm3m_relaxed/HfO2_Fm3m_relaxed_symrelaxed_PET_phonon_output_20260827-182403/HfO2_Fm3m_relaxed_symrelaxed_phonon_uq_data.npz',
        'Bi2Sn2O7': 'uq_validation_runs/034_Bi2Sn2O7_benc_Sn2Bi2O7_relaxed/Sn2Bi2O7_relaxed_symrelaxed_PET_phonon_output_20260827-184038/Sn2Bi2O7_relaxed_symrelaxed_phonon_uq_data.npz',
    }

    # Use LaTeX-formatted compound names with proper subscripts
    compound_labels = {
        'LiF': r'LiF',
        'CsPbI3': r'CsPbI$_3$',
        'HfO2': r'HfO$_2$',
        'Bi2Sn2O7': r'Bi$_2$Sn$_2$O$_7$',
    }

    fig = plt.figure(figsize=(14, 11))
    gs = GridSpec(2, 2, figure=fig, hspace=0.30, wspace=0.25)

    for i, (compound, npz_path) in enumerate(aristotypes.items()):
        ax = fig.add_subplot(gs[i // 2, i % 2])
        full_path = os.path.join(data_dir, npz_path)

        if not os.path.exists(full_path):
            ax.text(0.5, 0.5, f'{compound}\n(data not found)', transform=ax.transAxes,
                    ha='center', va='center', fontsize=14)
            ax.set_title(compound_labels.get(compound, compound), fontsize=14)
            continue

        data = np.load(full_path, allow_pickle=True)
        mean_freqs = data['mean_freqs']  # shape: (n_paths, n_qpoints, n_bands)
        std_freqs = data['std_freqs']

        # Use the first path (high-symmetry path)
        mean = mean_freqs[0]  # (n_qpoints, n_bands)
        std = std_freqs[0]

        n_q, n_bands = mean.shape
        colors = plt.cm.tab20(np.linspace(0, 1, min(n_bands, 20)))

        q_idx = np.arange(n_q)

        for b in range(min(n_bands, 20)):
            ax.fill_between(q_idx, mean[:, b] - 2*std[:, b], mean[:, b] + 2*std[:, b],
                            color=colors[b], alpha=0.15)
            ax.plot(q_idx, mean[:, b], color=colors[b], linewidth=0.6)

        ax.axhline(y=0, color='k', linestyle='-', linewidth=0.6)
        ax.set_ylabel('Frequency (THz)', fontsize=13)
        ax.set_xlabel('q-point index', fontsize=13)
        ax.set_title(f'{compound_labels.get(compound, compound)} ($n_\\mathrm{{bands}}={n_bands}$)',
                     fontsize=14)
        ax.set_ylim(bottom=min(mean.min() - 2*std.max(), -5))
        ax.tick_params(labelsize=11)

        # Add text with key stats — larger font
        omega_min = mean.min()
        max_std = std.max()
        std_at_min = std[np.unravel_index(np.argmin(mean), std.shape)]
        ax.text(0.02, 0.98, f'$\\omega_\\mathrm{{min}}$={omega_min:.2f} THz\n$\\sigma_\\mathrm{{max}}$={max_std:.2f} THz',
                transform=ax.transAxes, fontsize=10, va='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # No suptitle — caption in LaTeX handles the description

    for ext in ['pdf', 'png']:
        path = os.path.join(fig_dir, f'fig_si_pet_uq_dispersion.{ext}')
        fig.savefig(path, dpi=300, bbox_inches='tight')
        print(f"Saved: {path}")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description='Generate SI figures for criterion validation')
    parser.add_argument('--engines', nargs='+', default=['mace', 'pet', 'macemp', 'esen', 'uma'])
    parser.add_argument('--data-dir', default=PET_PHONONS_DIR)
    parser.add_argument('--fig-dir', default=os.path.join(PET_PHONONS_DIR, 'figures_5eng'))
    args = parser.parse_args()

    os.makedirs(args.fig_dir, exist_ok=True)

    print(f"Data dir: {args.data_dir}")
    print(f"Fig dir: {args.fig_dir}")
    print(f"Engines: {args.engines}")
    print()

    # Need sklearn for ROC
    try:
        from sklearn.metrics import roc_curve, auc
    except ImportError:
        print("WARNING: scikit-learn not available. Implementing ROC manually.")
        # Fallback ROC implementation
        def compute_roc(labels, scores):
            thresholds = np.unique(scores)
            thresholds = np.sort(thresholds)[::-1]
            fpr_list = [0.0]
            tpr_list = [0.0]
            for thresh in thresholds:
                pred = scores >= thresh
                tp = np.sum((pred == 1) & (labels == 1))
                fp = np.sum((pred == 1) & (labels == 0))
                tn = np.sum((pred == 0) & (labels == 0))
                fn = np.sum((pred == 0) & (labels == 1))
                fpr_list.append(fp / (fp + tn) if (fp + tn) > 0 else 0)
                tpr_list.append(tp / (tp + fn) if (tp + fn) > 0 else 0)
            fpr_list.append(1.0)
            tpr_list.append(1.0)
            fpr = np.array(fpr_list)
            tpr = np.array(tpr_list)
            # Compute AUC via trapezoidal rule
            idx = np.argsort(fpr)
            a = np.trapz(tpr[idx], fpr[idx])
            return fpr, tpr, a
        globals()['compute_roc'] = compute_roc

    print("1. Per-engine ROC curves...")
    plot_per_engine_roc(args.engines, args.data_dir, args.fig_dir)

    print("\n2. AUC vs p sweep...")
    plot_auc_vs_p(args.engines, args.data_dir, args.fig_dir)

    print("\n3. Marginal-mode bar chart...")
    plot_marginal_mode(args.engines, args.data_dir, args.fig_dir)

    print("\n4. PET-MAD UQ dispersion...")
    plot_pet_uq_dispersion(args.data_dir, args.fig_dir)

    print("\nDone.")


if __name__ == '__main__':
    main()
