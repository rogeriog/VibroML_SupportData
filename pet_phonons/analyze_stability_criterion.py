#!/usr/bin/env python3
"""
Analyze the stability criterion across multiple MLIP engines.

Consumes ip_variants_*.json files and uq_validation_structures.json labels,
produces:
  - Per-engine ROC/AUC for soft_frac, I_curv, Omega_p, omega_min
  - 3-tier confusion matrices + dangerous-error counts
  - p-sweep AUC table
  - Calibrated logistic regression with LOO-CV
  - Marginal-mode sensitivity analysis
  - All outputs as CSV + JSON for LaTeX generation

Usage:
  python analyze_stability_criterion.py
  python analyze_stability_criterion.py --engines mace pet macemp esen uma
"""
import argparse, json, os, sys
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneOut
from sklearn.preprocessing import StandardScaler
from scipy.special import expit  # sigmoid

# ============================================================================
# Constants
# ============================================================================

PET_PHONONS_DIR = "/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/02_runs_and_data/pet_phonons"

# Ground-truth labels for binary evaluation
GT_LABELS = ['high_sym_unstable', 'literature_stable']

# Structures to exclude
ANOMALY_IDXS = [9, 10, 11]  # LiF steric-clash anomalies (non-physical)

# p values for the sweep
P_VALUES = [0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0, 20.0, 50.0]
P_VALUES_INF = P_VALUES + ['inf']

# Criterion thresholds (must match stability_metrics.py)
SOFT_FRAC_STABLE_MAX = 0.02
SOFT_FRAC_UNSTABLE_MIN = 0.08


# ============================================================================
# Data loading
# ============================================================================

def load_engine_data(engine_label, data_dir=PET_PHONONS_DIR):
    """Load ip_variants JSON for a given engine."""
    path = os.path.join(data_dir, f"ip_variants_{engine_label}.json")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


def to_dataframe(data):
    """Convert ip_variants JSON list to a DataFrame."""
    rows = []
    for r in data:
        row = {
            'idx': r['idx'],
            'compound': r['compound'],
            'label': r.get('stability_label', 'unknown'),
            'space_group': r.get('space_group', ''),
            'omega_min': r.get('omega_min'),
            'soft_frac': r.get('soft_frac_no_thresh'),
            'i_curv': r.get('i_p_1.0'),
            'omega_inst': r.get('omega_p_1.0'),
        }
        for p in P_VALUES_INF:
            key_i = f"i_p_{p}"
            key_o = f"omega_p_{p}"
            row[key_i] = r.get(key_i)
            row[key_o] = r.get(key_o)
        rows.append(row)
    return pd.DataFrame(rows)


def get_ground_truth(df, exclude_anomalies=True, exclude_outliers=None):
    """Filter to ground-truth structures (high_sym_unstable + literature_stable)."""
    mask = df['label'].isin(GT_LABELS)
    if exclude_anomalies:
        mask &= ~df['idx'].isin(ANOMALY_IDXS)
    if exclude_outliers:
        mask &= ~df['idx'].isin(exclude_outliers)
    result = df[mask].copy()
    result['unstable'] = (result['label'] == 'high_sym_unstable').astype(int)
    return result


# ============================================================================
# ROC/AUC analysis
# ============================================================================

def compute_auc(y_true, scores):
    """Compute AUC, handling edge cases."""
    if len(np.unique(y_true)) < 2:
        return float('nan')
    # For metrics where higher = more unstable, we use scores directly
    # For omega_min, lower = more unstable, so negate
    try:
        return roc_auc_score(y_true, scores)
    except ValueError:
        return float('nan')


def per_engine_auc(df_gt):
    """Compute AUC for all metrics."""
    y = df_gt['unstable'].values
    results = {}

    # soft_frac: higher = more unstable
    results['soft_frac'] = compute_auc(y, df_gt['soft_frac'].values)

    # I_curv (p=1): higher = more unstable
    results['i_curv'] = compute_auc(y, df_gt['i_curv'].values)

    # omega_min: lower = more unstable → negate for AUC
    results['omega_min'] = compute_auc(y, -df_gt['omega_min'].values)

    # Omega_p for each p
    for p in P_VALUES_INF:
        col = f"omega_p_{p}"
        if col in df_gt.columns and df_gt[col].notna().any():
            results[f'omega_p_{p}'] = compute_auc(y, df_gt[col].values)

    return results


# ============================================================================
# 3-tier classification
# ============================================================================

def classify_3tier(soft_frac):
    """3-tier classification using soft_frac only."""
    if soft_frac >= SOFT_FRAC_UNSTABLE_MIN:
        return 'unstable'
    elif soft_frac >= SOFT_FRAC_STABLE_MAX:
        return 'quasi'
    else:
        return 'stable'


def evaluate_3tier(df_gt):
    """Evaluate 3-tier classification.

    Returns:
        dict with exact, quasi, dangerous, danger_rate
    """
    tiers = df_gt['soft_frac'].apply(classify_3tier)
    y = df_gt['unstable']

    danger = 0  # stable→unstable or unstable→stable
    buffered = 0  # landed in quasi
    exact = 0  # stable→stable or unstable→unstable

    for t, yy in zip(tiers, y):
        if yy == 0 and t == 'unstable':
            danger += 1
        elif yy == 1 and t == 'stable':
            danger += 1
        elif t == 'quasi':
            buffered += 1
        else:
            exact += 1

    n = len(df_gt)
    return {
        'n': n,
        'exact': exact,
        'quasi': buffered,
        'danger': danger,
        'danger_rate': danger / n if n > 0 else 0,
        'tiers': tiers.tolist(),
    }


def get_dangerous_errors(df_gt):
    """Return rows that are dangerously misclassified."""
    tiers = df_gt['soft_frac'].apply(classify_3tier)
    mask = ((df_gt['unstable'] == 0) & (tiers == 'unstable')) | \
           ((df_gt['unstable'] == 1) & (tiers == 'stable'))
    return df_gt[mask].copy()


# ============================================================================
# Calibrated logistic regression
# ============================================================================

def calibrated_logistic_loocv(df_gt, features='soft_frac'):
    """Leave-one-out CV logistic regression.

    Returns calibrated probabilities and LOO accuracy.
    """
    y = df_gt['unstable'].values
    X = df_gt[[features]].values
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    loo = LeaveOneOut()
    probs = np.zeros(len(y))
    correct = 0

    for train_idx, test_idx in loo.split(X_scaled):
        X_train, y_train = X_scaled[train_idx], y[train_idx]
        X_test = X_scaled[test_idx]

        if len(np.unique(y_train)) < 2:
            probs[test_idx] = 0.5
            continue

        model = LogisticRegression(C=1.0, solver='lbfgs')
        model.fit(X_train, y_train)
        prob = model.predict_proba(X_test)[0, 1]
        probs[test_idx] = prob

        pred = 1 if prob >= 0.5 else 0
        if pred == y[test_idx[0]]:
            correct += 1

    return {
        'probabilities': probs,
        'loo_accuracy': correct / len(y),
        'y_true': y,
    }


# ============================================================================
# Marginal-mode sensitivity
# ============================================================================

def marginal_mode_sensitivity(freqs_dir, engine_label, threshold=0.5):
    """Analyze how much marginal modes (|omega| < threshold) contribute to I_p.

    Returns per-structure: fraction of modes marginal, fraction of I_0.5
    and I_1.0 contributed by marginal modes.
    """
    import glob
    freqs_path = os.path.join(freqs_dir, f"ip_variants_freqs_{engine_label}")
    files = sorted(glob.glob(os.path.join(freqs_path, "*.npz")))

    results = []
    for f in files:
        d = np.load(f, allow_pickle=True)
        freqs = d['freqs']
        gamma_idx = d['gamma_indices']

        # Zero out acoustic at Gamma
        freqs = freqs.copy()
        for gi in gamma_idx:
            freqs[gi, :3] = 0.0

        flat = freqs.ravel()
        neg = flat[flat < 0]
        if len(neg) == 0:
            continue

        c05 = np.abs(neg)
        c10 = neg ** 2

        marginal_mask = np.abs(neg) < threshold
        frac_marginal_modes = marginal_mask.sum() / len(neg)
        frac_contrib_p05 = c05[marginal_mask].sum() / c05.sum() if c05.sum() > 0 else 0
        frac_contrib_p10 = c10[marginal_mask].sum() / c10.sum() if c10.sum() > 0 else 0

        results.append({
            'idx': int(d['idx']),
            'compound': str(d['compound']),
            'n_neg': len(neg),
            'frac_marginal_modes': frac_marginal_modes,
            'frac_contrib_p05': frac_contrib_p05,
            'frac_contrib_p10': frac_contrib_p10,
        })

    return pd.DataFrame(results)


# ============================================================================
# Main analysis
# ============================================================================

def main():
    parser = argparse.ArgumentParser(description="Analyze stability criterion across engines")
    parser.add_argument("--engines", nargs='+', default=['mace', 'pet'],
                        help="Engine labels to analyze (e.g. mace pet macemp esen uma)")
    parser.add_argument("--data-dir", default=PET_PHONONS_DIR)
    parser.add_argument("--output-dir", default="analysis_outputs")
    parser.add_argument("--exclude-outliers", nargs='*', type=int, default=None,
                        help="Structure idxs to exclude as known outliers (e.g. 6 for LiF Fm-3m)")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Load all engine data
    engine_dfs = {}
    for eng in args.engines:
        data = load_engine_data(eng, args.data_dir)
        if data is None:
            print(f"WARNING: No data for engine '{eng}' (file not found)")
            continue
        df = to_dataframe(data)
        engine_dfs[eng] = df
        print(f"Loaded {eng}: {len(df)} structures")

    if not engine_dfs:
        print("ERROR: No engine data found")
        sys.exit(1)

    # ================================================================
    # 1. Per-engine AUC
    # ================================================================
    print("\n" + "=" * 80)
    print("1. Per-engine AUC (ground truth: high_sym_unstable vs literature_stable)")
    print("=" * 80)

    auc_rows = []
    for eng, df in engine_dfs.items():
        df_gt = get_ground_truth(df, exclude_outliers=args.exclude_outliers)
        print(f"\n  {eng}: {len(df_gt)} GT structures ({df_gt['unstable'].sum()} unstable)")

        auc = per_engine_auc(df_gt)
        auc['engine'] = eng
        auc['n_gt'] = len(df_gt)
        auc['n_unstable'] = int(df_gt['unstable'].sum())
        auc_rows.append(auc)

        for metric, val in sorted(auc.items()):
            if metric not in ('engine', 'n_gt', 'n_unstable'):
                print(f"    {metric:>20}: AUC = {val:.4f}")

    auc_df = pd.DataFrame(auc_rows)
    auc_df.to_csv(os.path.join(args.output_dir, "roc_auc_per_engine.csv"), index=False)
    print(f"\n  Saved: {args.output_dir}/roc_auc_per_engine.csv")

    # ================================================================
    # 2. 3-tier confusion
    # ================================================================
    print("\n" + "=" * 80)
    print("2. 3-tier classification (soft_frac: <2% stable, 2-8% quasi, >=8% unstable)")
    print("=" * 80)

    tier_rows = []
    for eng, df in engine_dfs.items():
        df_gt = get_ground_truth(df, exclude_outliers=args.exclude_outliers)
        result = evaluate_3tier(df_gt)
        result['engine'] = eng
        tier_rows.append(result)

        print(f"\n  {eng}: exact={result['exact']}, quasi={result['quasi']}, "
              f"DANGEROUS={result['danger']}/{result['n']}")

        # Show dangerous errors
        errors = get_dangerous_errors(df_gt)
        if len(errors) > 0:
            print(f"    Dangerous errors:")
            for _, r in errors.iterrows():
                print(f"      idx={int(r['idx'])}: {r['compound']}, {r['label']}, "
                      f"sf={r['soft_frac']:.4f}, wmin={r['omega_min']:.2f}")

    tier_df = pd.DataFrame(tier_rows)
    tier_df.to_csv(os.path.join(args.output_dir, "three_tier_confusion.csv"), index=False)
    print(f"\n  Saved: {args.output_dir}/three_tier_confusion.csv")

    # ================================================================
    # 3. p-sweep AUC (per engine)
    # ================================================================
    print("\n" + "=" * 80)
    print("3. p-sweep AUC")
    print("=" * 80)

    p_sweep_rows = []
    for eng, df in engine_dfs.items():
        df_gt = get_ground_truth(df, exclude_outliers=args.exclude_outliers)
        y = df_gt['unstable'].values

        for p in P_VALUES_INF:
            col = f"omega_p_{p}"
            if col not in df_gt.columns or df_gt[col].isna().all():
                continue
            auc = compute_auc(y, df_gt[col].values)
            p_sweep_rows.append({'engine': eng, 'p': str(p), 'auc': auc})

    # Add soft_frac for comparison
    for eng, df in engine_dfs.items():
        df_gt = get_ground_truth(df, exclude_outliers=args.exclude_outliers)
        y = df_gt['unstable'].values
        auc = compute_auc(y, df_gt['soft_frac'].values)
        p_sweep_rows.append({'engine': eng, 'p': 'soft_frac', 'auc': auc})

    p_sweep_df = pd.DataFrame(p_sweep_rows)
    p_sweep_df.to_csv(os.path.join(args.output_dir, "p_sweep.csv"), index=False)

    # Print summary table
    print(f"\n  {'p':>10} | " + " | ".join(f"{eng:>10}" for eng in engine_dfs.keys()))
    print("  " + "-" * (12 + 13 * len(engine_dfs)))
    for p in P_VALUES_INF + ['soft_frac']:
        vals = []
        for eng in engine_dfs.keys():
            row = p_sweep_df[(p_sweep_df['engine'] == eng) & (p_sweep_df['p'] == str(p))]
            if len(row) > 0:
                vals.append(f"{row.iloc[0]['auc']:.4f}")
            else:
                vals.append("—")
        print(f"  {str(p):>10} | " + " | ".join(f"{v:>10}" for v in vals))

    print(f"\n  Saved: {args.output_dir}/p_sweep.csv")

    # ================================================================
    # 4. Calibrated logistic regression (LOO-CV)
    # ================================================================
    print("\n" + "=" * 80)
    print("4. Calibrated logistic regression (LOO-CV)")
    print("=" * 80)

    for eng, df in engine_dfs.items():
        df_gt = get_ground_truth(df, exclude_outliers=args.exclude_outliers)
        if len(df_gt) < 10 or df_gt['unstable'].nunique() < 2:
            print(f"\n  {eng}: insufficient data for LOO-CV")
            continue

        result = calibrated_logistic_loocv(df_gt, features='soft_frac')
        print(f"\n  {eng}: LOO accuracy = {result['loo_accuracy']:.4f} "
              f"(N={len(df_gt)}, {df_gt['unstable'].sum()} unstable)")

        # Save calibrated probabilities
        prob_df = pd.DataFrame({
            'idx': df_gt['idx'].values,
            'compound': df_gt['compound'].values,
            'label': df_gt['label'].values,
            'unstable': result['y_true'],
            'soft_frac': df_gt['soft_frac'].values,
            'calibrated_prob': result['probabilities'],
        })
        prob_df.to_csv(os.path.join(args.output_dir, f"calibration_loo_{eng}.csv"), index=False)

    # ================================================================
    # 5. Marginal-mode sensitivity
    # ================================================================
    print("\n" + "=" * 80)
    print("5. Marginal-mode sensitivity (|omega| < 0.5 THz)")
    print("=" * 80)

    for eng in engine_dfs.keys():
        freqs_dir = os.path.join(args.data_dir, f"ip_variants_freqs_{eng}")
        if not os.path.exists(freqs_dir):
            print(f"\n  {eng}: no frequency data directory")
            continue

        mm_df = marginal_mode_sensitivity(args.data_dir, eng, threshold=0.5)
        if len(mm_df) > 0:
            mm_df.to_csv(os.path.join(args.output_dir, f"marginal_mode_sensitivity_{eng}.csv"),
                        index=False)
            avg_marginal = mm_df['frac_marginal_modes'].mean()
            avg_contrib_p05 = mm_df['frac_contrib_p05'].mean()
            avg_contrib_p10 = mm_df['frac_contrib_p10'].mean()
            print(f"\n  {eng}: {len(mm_df)} structures with imaginary modes")
            print(f"    avg % marginal modes: {100*avg_marginal:.1f}%")
            print(f"    avg % contrib to I_0.5: {100*avg_contrib_p05:.1f}%")
            print(f"    avg % contrib to I_1.0: {100*avg_contrib_p10:.1f}%")

    # ================================================================
    # 6. Master metrics table (all engines)
    # ================================================================
    print("\n" + "=" * 80)
    print("6. Master metrics table")
    print("=" * 80)

    master_rows = []
    for eng, df in engine_dfs.items():
        df_gt = get_ground_truth(df, exclude_outliers=args.exclude_outliers)
        for _, r in df_gt.iterrows():
            master_rows.append({
                'engine': eng,
                'idx': int(r['idx']),
                'compound': r['compound'],
                'label': r['label'],
                'ground_truth': 'unstable' if r['unstable'] else 'stable',
                'soft_frac': r['soft_frac'],
                'omega_min': r['omega_min'],
                'i_curv': r['i_curv'],
                'omega_inst': r['omega_inst'],
                'classification': classify_3tier(r['soft_frac']),
            })

    master_df = pd.DataFrame(master_rows)
    master_df.to_csv(os.path.join(args.output_dir, "stability_metrics_all_engines.csv"), index=False)
    print(f"\n  Saved: {args.output_dir}/stability_metrics_all_engines.csv ({len(master_df)} rows)")

    print("\n" + "=" * 80)
    print("Analysis complete.")
    print(f"Outputs in: {args.output_dir}/")
    print("=" * 80)


if __name__ == "__main__":
    main()
