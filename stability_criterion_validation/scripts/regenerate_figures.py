#!/usr/bin/env python3
"""
Regenerate the three stability-criterion figures with larger fonts (~1.5x)
and updated convergence plot (Omega_inst in panel B + Youden's J line).

Reads existing data files — does NOT recompute phonons.

Usage:
  conda activate /gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/04_envs_and_archives/vibroml_env
  python regenerate_figures.py
"""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score, roc_curve

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PET_DIR = os.path.join(SCRIPT_DIR, "..", "..", "pet_phonons")
CONV_JSON = os.path.join(SCRIPT_DIR, "..", "analysis_outputs", "convergence_sweep.json")
FIG_DIR = os.path.join(SCRIPT_DIR, "..", "..", "..", "03_paper_and_rebuttal",
                       "VIBROML_PAPER", "figures")

ANOMALY_IDXS = [9, 10, 11]

# Font sizes — ~1.5x the originals
FS_TITLE = 16   # was ~11
FS_LABEL = 14   # was ~9-10
FS_LEGEND = 12  # was ~7-8
FS_TICK = 12    # was ~8-9
FS_ANNOT = 11   # was ~7

SOFT_FRAC_REF = 0.045  # 4.5% Youden's J for f_soft

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


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def load_engine_data(engine_label):
    path = os.path.join(PET_DIR, f"ip_variants_{engine_label}.json")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


def to_dataframe(data):
    import pandas as pd
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


def compute_youdens_j(labels, scores):
    """Compute Youden's J = max(TPR - FPR) and the corresponding threshold."""
    fpr, tpr, thresholds = roc_curve(labels, scores)
    j_values = tpr - fpr
    best_idx = np.argmax(j_values)
    return j_values[best_idx], thresholds[best_idx]


# ---------------------------------------------------------------------------
# Figure 1: Main criterion validation (4 panels)
# ---------------------------------------------------------------------------
def generate_main_figure():
    print("Generating main criterion validation figure...")
    engine_dfs = {}
    for eng in ['mace', 'pet']:
        data = load_engine_data(eng)
        if data is not None:
            engine_dfs[eng] = to_dataframe(data)
    if not engine_dfs:
        print("  No data found!")
        return

    fig, axes = plt.subplots(2, 2, figsize=(14, 11))

    # Panel A: soft_frac distributions
    ax = axes[0, 0]
    all_stable, all_unstable = [], []
    for eng, df in engine_dfs.items():
        df_gt = get_ground_truth(df)
        all_stable.extend(df_gt[df_gt['unstable'] == 0]['soft_frac'].values * 100)
        all_unstable.extend(df_gt[df_gt['unstable'] == 1]['soft_frac'].values * 100)

    bp = ax.boxplot([all_stable, all_unstable], positions=[1, 2],
                    widths=0.6, patch_artist=True, showfliers=True)
    bp['boxes'][0].set_facecolor('lightblue')
    bp['boxes'][1].set_facecolor('lightcoral')
    ax.set_xticks([1, 2])
    ax.set_xticklabels(['Stable\n(literature)', 'Unstable\n(known)'], fontsize=FS_TICK)
    ax.set_ylabel(r'$f_\mathrm{soft}$ (%)', fontsize=FS_LABEL)
    ax.set_title('(A) Soft-mode fraction distributions', fontsize=FS_TITLE)
    ax.axhline(4.5, color='gray', ls='--', lw=1.2,
               label='optimal separation\n(Youden $J$, pooled)')
    ax.legend(fontsize=FS_LEGEND, loc='upper left')
    ax.set_ylim(-1, max(max(all_unstable) + 2, 35))
    ax.tick_params(axis='both', labelsize=FS_TICK)

    # Panel B: ROC curves
    ax = axes[0, 1]
    eng = list(engine_dfs.keys())[0]
    df_gt = get_ground_truth(engine_dfs[eng])
    y = df_gt['unstable'].values
    for col, label, color in [
        ('soft_frac', r'$f_\mathrm{soft}$', 'tab:blue'),
        ('i_curv', r'$I_\mathrm{curv}$', 'tab:orange'),
        ('omega_min', r'$\omega_\mathrm{min}$ (negated)', 'tab:red'),
    ]:
        scores = -df_gt[col].values if col == 'omega_min' else df_gt[col].values
        if len(np.unique(y)) < 2:
            continue
        fpr, tpr, _ = roc_curve(y, scores)
        auc = roc_auc_score(y, scores)
        ax.plot(fpr, tpr, label=f'{label} (AUC={auc:.3f})', color=color, lw=2.2)
    ax.plot([0, 1], [0, 1], 'k--', lw=0.5)
    ax.set_xlabel('False Positive Rate', fontsize=FS_LABEL)
    ax.set_ylabel('True Positive Rate', fontsize=FS_LABEL)
    ax.set_title(f'(B) ROC curves ({ENGINE_LABELS[eng]})', fontsize=FS_TITLE)
    ax.legend(fontsize=FS_LEGEND, loc='lower right')
    ax.set_xlim(-0.01, 1.01)
    ax.set_ylim(-0.01, 1.01)
    ax.tick_params(axis='both', labelsize=FS_TICK)

    # Panel C: I_curv vs soft_frac scatter
    ax = axes[1, 0]
    colors_gt = {0: 'tab:blue', 1: 'tab:red'}
    markers = {'mace': 'o', 'pet': 's', 'macemp': 'D', 'esen': '^', 'uma': 'v'}
    for eng, df in engine_dfs.items():
        df_gt = get_ground_truth(df)
        marker = markers.get(eng, 'o')
        for gt_val in [0, 1]:
            subset = df_gt[df_gt['unstable'] == gt_val]
            ax.scatter(subset['soft_frac'] * 100, subset['i_curv'],
                      c=colors_gt[gt_val], marker=marker, alpha=0.6, s=60,
                      label=f'{eng}: {"unstable" if gt_val else "stable"}')
    ax.axvline(4.5, color='gray', ls='--', lw=1.2, alpha=0.5)
    ax.set_xlabel(r'$f_\mathrm{soft}$ (%)', fontsize=FS_LABEL)
    ax.set_ylabel(r'$I_\mathrm{curv}$ (THz$^2$)', fontsize=FS_LABEL)
    ax.set_title(r'(C) Severity ($I_\mathrm{curv}$) vs breadth ($f_\mathrm{soft}$)',
                fontsize=FS_TITLE)
    ax.set_yscale('symlog', linthresh=0.01)
    ax.legend(fontsize=FS_LEGEND - 1, loc='upper left', ncol=2)
    ax.tick_params(axis='both', labelsize=FS_TICK)

    # Panel D: Calibrated probability
    ax = axes[1, 1]
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import LeaveOneOut
    from sklearn.preprocessing import StandardScaler
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
    ax.axvline(4.5, color='gray', ls='--', lw=1.2, alpha=0.5)
    ax.set_xlabel(r'$f_\mathrm{soft}$ (%)', fontsize=FS_LABEL)
    ax.set_ylabel(r'$P(\mathrm{unstable}|f_\mathrm{soft})$ (LOO-CV)', fontsize=FS_LABEL)
    ax.set_title(f'(D) Calibrated probability ({eng}, LOO-CV)', fontsize=FS_TITLE)
    ax.legend(fontsize=FS_LEGEND)
    ax.set_ylim(-0.05, 1.05)
    ax.tick_params(axis='both', labelsize=FS_TICK)

    plt.tight_layout()
    out = os.path.join(FIG_DIR, "figure_criterion_validation.pdf")
    fig.savefig(out, dpi=300, bbox_inches='tight')
    fig.savefig(out.replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
    print(f"  Saved: {out}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 2: SI per-engine ROC (2 panels)
# ---------------------------------------------------------------------------
def generate_si_roc_figure():
    print("Generating SI per-engine ROC figure...")
    engines = ['mace', 'pet', 'macemp', 'esen', 'uma']
    engine_dfs = {}
    for eng in engines:
        data = load_engine_data(eng)
        if data is not None:
            engine_dfs[eng] = to_dataframe(data)
    if not engine_dfs:
        print("  No data found!")
        return

    fig, axes = plt.subplots(1, 2, figsize=(13, 6))

    # Panel A: All engines
    ax = axes[0]
    for eng in engines:
        if eng not in engine_dfs:
            continue
        df_gt = get_ground_truth(engine_dfs[eng])
        y = df_gt['unstable'].values
        scores = df_gt['soft_frac'].values
        if len(np.unique(y)) < 2:
            continue
        fpr, tpr, _ = roc_curve(y, scores)
        auc = roc_auc_score(y, scores)
        ax.plot(fpr, tpr, color=ENGINE_COLORS[eng],
                label=f'{ENGINE_LABELS[eng]} (AUC={auc:.3f})', linewidth=2)
    ax.plot([0, 1], [0, 1], 'k--', alpha=0.3)
    ax.set_xlabel('False Positive Rate', fontsize=FS_LABEL)
    ax.set_ylabel('True Positive Rate', fontsize=FS_LABEL)
    ax.set_title('(A) ROC: $f_\\mathrm{soft}$ (all engines)', fontsize=FS_TITLE)
    ax.legend(fontsize=FS_LEGEND, loc='lower right')
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.tick_params(axis='both', labelsize=FS_TICK)

    # Panel B: Different metrics (MACE-OMAT)
    ax = axes[1]
    eng = 'mace'
    if eng in engine_dfs:
        data = load_engine_data(eng)
        for metric, label, color in [
            ('soft_frac_no_thresh', '$f_\\mathrm{soft}$', '#1f77b4'),
            ('i_p_1.0', '$I_\\mathrm{curv}$', '#ff7f0e'),
            ('omega_min', '$\\omega_\\mathrm{min}$', '#d62728'),
            ('omega_p_0.5', '$\\Omega_{0.5}$', '#2ca02c'),
        ]:
            labels_list, scores_list = [], []
            for r in data:
                lab = r.get('stability_label', '')
                if lab in ('high_sym_unstable', 'literature_stable'):
                    labels_list.append(1 if lab == 'high_sym_unstable' else 0)
                    val = r.get(metric, 0)
                    if metric == 'omega_min':
                        val = -val
                    scores_list.append(val)
            if len(np.unique(labels_list)) < 2:
                continue
            fpr, tpr, a = compute_roc(np.array(labels_list), np.array(scores_list))
            ax.plot(fpr, tpr, color=color, label=f'{label} (AUC={a:.3f})', linewidth=2)
    ax.plot([0, 1], [0, 1], 'k--', alpha=0.3)
    ax.set_xlabel('False Positive Rate', fontsize=FS_LABEL)
    ax.set_ylabel('True Positive Rate', fontsize=FS_LABEL)
    ax.set_title(f'(B) ROC: metrics ({ENGINE_LABELS[eng]})', fontsize=FS_TITLE)
    ax.legend(fontsize=FS_LEGEND, loc='lower right')
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.tick_params(axis='both', labelsize=FS_TICK)

    fig.tight_layout()
    out = os.path.join(FIG_DIR, "fig_si_per_engine_roc.pdf")
    fig.savefig(out, dpi=300, bbox_inches='tight')
    fig.savefig(out.replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
    print(f"  Saved: {out}")
    plt.close(fig)


def compute_roc(labels, scores):
    fpr, tpr, _ = roc_curve(labels, scores)
    auc = roc_auc_score(labels, scores)
    return fpr, tpr, auc


# ---------------------------------------------------------------------------
# Figure 3: Convergence sweep (2 panels, with Omega_inst + Youden's J)
# ---------------------------------------------------------------------------
def generate_convergence_figure():
    print("Generating convergence sweep figure...")
    if not os.path.exists(CONV_JSON):
        print(f"  Convergence data not found: {CONV_JSON}")
        return
    with open(CONV_JSON) as f:
        results = json.load(f)
    results = [r for r in results if r.get("soft_frac") is not None]
    if not results:
        print("  No valid convergence results!")
        return

    # Compute Youden's J for omega_inst using the full benchmark
    print("  Computing Youden's J for omega_inst...")
    mace_data = load_engine_data('mace')
    labels_list, omega_inst_list = [], []
    for r in mace_data:
        lab = r.get('stability_label', '')
        if lab in ('high_sym_unstable', 'literature_stable') and r['idx'] not in ANOMALY_IDXS:
            labels_list.append(1 if lab == 'high_sym_unstable' else 0)
            omega_inst_list.append(r.get('omega_p_1.0', 0))
    labels_arr = np.array(labels_list)
    omega_inst_arr = np.array(omega_inst_list)
    j_omega, thresh_omega = compute_youdens_j(labels_arr, omega_inst_arr)
    print(f"    Youden's J for omega_inst: J={j_omega:.3f}, threshold={thresh_omega:.3f} THz")

    # Organize convergence results
    by_compound = {}
    for r in results:
        comp = r["compound"]
        q = r["q_mesh_size"]
        by_compound.setdefault(comp, {})[q] = r
    compounds = sorted(by_compound.keys())
    cmap = plt.get_cmap("tab10")
    colors = {comp: cmap(i % 10) for i, comp in enumerate(compounds)}
    markers = {comp: m for comp, m in zip(compounds, ["o", "s", "^", "D", "v", "P"])}

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), constrained_layout=True)

    # Panel A: soft_frac
    axA = axes[0]
    for comp in compounds:
        qs = sorted(by_compound[comp].keys())
        ys = [by_compound[comp][q]["soft_frac"] * 100 for q in qs]
        axA.plot(qs, ys, marker=markers.get(comp, "o"), color=colors[comp],
                 label=comp, linewidth=2.2, markersize=9)
    axA.axhline(SOFT_FRAC_REF * 100, color="gray", linestyle="--", linewidth=1.2,
                alpha=0.7, label=f"optimal sep. ({SOFT_FRAC_REF*100:.1f}%)")
    axA.set_xlabel("q-mesh size (n × n × n)", fontsize=FS_LABEL)
    axA.set_ylabel(r"$f_\mathrm{soft}$ (%)", fontsize=FS_LABEL)
    axA.set_title("(A) Soft-mode fraction vs q-mesh", fontsize=FS_TITLE)
    axA.set_xticks([4, 6, 8])
    axA.legend(fontsize=FS_LEGEND, loc="best", framealpha=0.9)
    axA.grid(True, alpha=0.3)
    axA.tick_params(axis='both', labelsize=FS_TICK)

    # Panel B: Omega_inst (= sqrt(I_curv))
    axB = axes[1]
    for comp in compounds:
        qs = sorted(by_compound[comp].keys())
        ys = [np.sqrt(max(by_compound[comp][q]["I_curv"], 0)) for q in qs]
        axB.plot(qs, ys, marker=markers.get(comp, "o"), color=colors[comp],
                 label=comp, linewidth=2.2, markersize=9)
    axB.set_xlabel("q-mesh size (n × n × n)", fontsize=FS_LABEL)
    axB.set_ylabel(r"$\Omega_\mathrm{inst}$ (THz)", fontsize=FS_LABEL)
    axB.set_title(r"(B) Instability severity $\Omega_\mathrm{inst}=\sqrt{I_\mathrm{curv}}$ vs q-mesh",
                 fontsize=FS_TITLE)
    axB.set_xticks([4, 6, 8])
    axB.set_ylim(bottom=0)
    axB.legend(fontsize=FS_LEGEND, loc="best", framealpha=0.9)
    axB.grid(True, alpha=0.3)
    axB.tick_params(axis='both', labelsize=FS_TICK)

    out = os.path.join(FIG_DIR, "convergence_sweep.pdf")
    fig.savefig(out, dpi=300, bbox_inches='tight')
    fig.savefig(out.replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
    print(f"  Saved: {out}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    os.makedirs(FIG_DIR, exist_ok=True)
    print(f"Figure output dir: {FIG_DIR}")
    print()

    generate_main_figure()
    print()
    generate_si_roc_figure()
    print()
    generate_convergence_figure()

    print("\nDone.")


if __name__ == "__main__":
    main()
