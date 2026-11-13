#!/usr/bin/env python3
"""
Analyze I_p variants: find the best p value and threshold
that works cross-model (MACE and PET-MAD).
"""
import json, numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score, accuracy_score, confusion_matrix, f1_score

P_VALUES = [0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0, 20.0, 50.0]
P_VALUES_INF = P_VALUES + ['inf']

def load_results(path):
    with open(path) as f:
        return json.load(f)

def main():
    mace_path = "/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/02_runs_and_data/pet_phonons/ip_variants_mace.json"
    pet_path = "/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/02_runs_and_data/pet_phonons/ip_variants_pet.json"

    mace_data = load_results(mace_path)
    pet_data = load_results(pet_path)

    print(f"MACE: {len(mace_data)} structures")
    print(f"PET:  {len(pet_data)} structures")

    # Build DataFrames
    def to_df(data):
        rows = []
        for r in data:
            row = {'idx': r['idx'], 'compound': r['compound'],
                   'label': r['stability_label'], 'natoms': r['natoms']}
            for p in P_VALUES_INF:
                key = f"i_p_{p}" if p != 'inf' else "i_p_inf"
                okey = f"omega_p_{p}" if p != 'inf' else "omega_p_inf"
                row[f"i_p_{p}"] = r.get(key)
                row[f"omega_p_{p}"] = r.get(okey)
            row['omega_min'] = r.get('omega_min')
            row['soft_frac_no_thresh'] = r.get('soft_frac_no_thresh')
            row['soft_frac_0p05'] = r.get('soft_frac_0p05')
            rows.append(row)
        return pd.DataFrame(rows)

    df_mace = to_df(mace_data)
    df_pet = to_df(pet_data)

    # Ground truth: high_sym_unstable=1, literature_stable=0
    # Exclude anomalies and ga_discovered
    gt_labels = ['high_sym_unstable', 'literature_stable']
    anomaly_idxs = [9, 10, 11]

    df_mace_gt = df_mace[(df_mace['label'].isin(gt_labels)) & (~df_mace['idx'].isin(anomaly_idxs))].copy()
    df_pet_gt = df_pet[(df_pet['label'].isin(gt_labels)) & (~df_pet['idx'].isin(anomaly_idxs))].copy()

    df_mace_gt['unstable'] = (df_mace_gt['label'] == 'high_sym_unstable').astype(int)
    df_pet_gt['unstable'] = (df_pet_gt['label'] == 'high_sym_unstable').astype(int)

    print(f"\nGround truth: {len(df_mace_gt)} MACE, {len(df_pet_gt)} PET")
    print(f"  Unstable: {df_mace_gt['unstable'].sum()}")
    print(f"  Stable: {(1-df_mace_gt['unstable']).sum()}")

    # ============================================================
    # AUC for each p, each model
    # ============================================================
    print(f"\n{'='*90}")
    print(f"{'p':>6} | {'MACE I_p AUC':>12} {'MACE Ω_p AUC':>12} {'MACE best_acc':>13} | {'PET I_p AUC':>12} {'PET Ω_p AUC':>12} {'PET best_acc':>13} | {'min_acc':>8}")
    print(f"{'='*90}")

    best_results = []

    for p in P_VALUES_INF:
        col_i = f"i_p_{p}"
        col_o = f"omega_p_{p}"

        # MACE
        mace_i = df_mace_gt[col_i].dropna()
        mace_o = df_mace_gt[col_o].dropna()
        y_m = df_mace_gt.loc[mace_i.index, 'unstable']

        if len(mace_i) < 5 or y_m.nunique() < 2:
            print(f"{p:>6} | {'N/A':>12} {'N/A':>12} {'N/A':>13} | ", end="")
            mace_auc_i = mace_auc_o = mace_best_acc = 0
        else:
            mace_auc_i = roc_auc_score(y_m, mace_i)
            mace_auc_o = roc_auc_score(y_m, mace_o)
            # Best threshold accuracy
            mace_best_acc = 0
            for th in np.linspace(mace_o.min(), mace_o.max(), 200):
                pred = (mace_o >= th).astype(int)
                acc = accuracy_score(y_m, pred)
                if acc > mace_best_acc:
                    mace_best_acc = acc
            print(f"{p:>6} | {mace_auc_i:>12.4f} {mace_auc_o:>12.4f} {mace_best_acc:>13.3f} | ", end="")

        # PET
        pet_i = df_pet_gt[col_i].dropna()
        pet_o = df_pet_gt[col_o].dropna()
        y_p = df_pet_gt.loc[pet_i.index, 'unstable']

        if len(pet_i) < 5 or y_p.nunique() < 2:
            print(f"{'N/A':>12} {'N/A':>12} {'N/A':>13} | {'N/A':>8}")
            pet_auc_i = pet_auc_o = pet_best_acc = 0
        else:
            pet_auc_i = roc_auc_score(y_p, pet_i)
            pet_auc_o = roc_auc_score(y_p, pet_o)
            pet_best_acc = 0
            for th in np.linspace(pet_o.min(), pet_o.max(), 200):
                pred = (pet_o >= th).astype(int)
                acc = accuracy_score(y_p, pred)
                if acc > pet_best_acc:
                    pet_best_acc = acc
            min_acc = min(mace_best_acc, pet_best_acc)
            print(f"{pet_auc_i:>12.4f} {pet_auc_o:>12.4f} {pet_best_acc:>13.3f} | {min_acc:>8.3f}")

        best_results.append({
            'p': p, 'mace_auc': mace_auc_o, 'pet_auc': pet_auc_o,
            'mace_acc': mace_best_acc, 'pet_acc': pet_best_acc,
            'min_acc': min(mace_best_acc, pet_best_acc)
        })

    # Also check soft_frac and omega_min
    for metric_name, col in [('soft_frac (no thresh)', 'soft_frac_no_thresh'),
                              ('soft_frac (0.05)', 'soft_frac_0p05'),
                              ('omega_min (neg)', 'omega_min')]:
        mace_v = df_mace_gt[col].dropna()
        pet_v = df_pet_gt[col].dropna()
        y_m = df_mace_gt.loc[mace_v.index, 'unstable']
        y_p = df_pet_gt.loc[pet_v.index, 'unstable']
        if col == 'omega_min':
            mace_v = -mace_v  # negate
            pet_v = -pet_v
        mace_auc = roc_auc_score(y_m, mace_v) if y_m.nunique() > 1 else 0
        pet_auc = roc_auc_score(y_p, pet_v) if y_p.nunique() > 1 else 0
        mace_acc = 0
        pet_acc = 0
        for th in np.linspace(mace_v.min(), mace_v.max(), 200):
            pred = (mace_v >= th).astype(int)
            acc = accuracy_score(y_m, pred)
            if acc > mace_acc: mace_acc = acc
        for th in np.linspace(pet_v.min(), pet_v.max(), 200):
            pred = (pet_v >= th).astype(int)
            acc = accuracy_score(y_p, pred)
            if acc > pet_acc: pet_acc = acc
        print(f"{metric_name:>20} | {mace_auc:>12.4f} {'—':>12} {mace_acc:>13.3f} | {pet_auc:>12.4f} {'—':>12} {pet_acc:>13.3f} | {min(mace_acc,pet_acc):>8.3f}")

    # ============================================================
    # Find best cross-model threshold for each p
    # ============================================================
    print(f"\n{'='*90}")
    print("=== Cross-model shared threshold analysis ===")
    print(f"{'='*90}")
    print(f"{'p':>6} | {'Best shared Ω_p thresh':>22} | {'MACE acc':>8} {'MACE TP':>7} {'MACE FP':>7} {'MACE FN':>7} {'MACE TN':>7} | {'PET acc':>8} {'PET TP':>7} {'PET FP':>7} {'PET FN':>7} {'PET TN':>7} | {'min_acc':>8}")
    print("-" * 150)

    for p in P_VALUES_INF:
        col_o = f"omega_p_{p}"
        mace_o = df_mace_gt[col_o].dropna()
        pet_o = df_pet_gt[col_o].dropna()
        y_m = df_mace_gt.loc[mace_o.index, 'unstable']
        y_p = df_pet_gt.loc[pet_o.index, 'unstable']

        # Combined threshold range
        all_vals = np.concatenate([mace_o.values, pet_o.values])
        thresholds = np.linspace(all_vals.min(), all_vals.max(), 300)

        best_min_acc = 0
        best_th = 0
        for th in thresholds:
            pred_m = (mace_o >= th).astype(int)
            pred_p = (pet_o >= th).astype(int)
            acc_m = accuracy_score(y_m, pred_m)
            acc_p = accuracy_score(y_p, pred_p)
            if min(acc_m, acc_p) > best_min_acc:
                best_min_acc = min(acc_m, acc_p)
                best_th = th
                best_acc_m = acc_m
                best_acc_p = acc_p

        pred_m = (mace_o >= best_th).astype(int)
        pred_p = (pet_o >= best_th).astype(int)
        cm_m = confusion_matrix(y_m, pred_m, labels=[0, 1])
        cm_p = confusion_matrix(y_p, pred_p, labels=[0, 1])

        print(f"{str(p):>6} | {best_th:>22.6f} | {best_acc_m:>8.3f} {cm_m[1,1]:>7d} {cm_m[0,1]:>7d} {cm_m[1,0]:>7d} {cm_m[0,0]:>7d} | {best_acc_p:>8.3f} {cm_p[1,1]:>7d} {cm_p[0,1]:>7d} {cm_p[1,0]:>7d} {cm_p[0,0]:>7d} | {best_min_acc:>8.3f}")

    # ============================================================
    # Show sorted values for the best p
    # ============================================================
    best_p_row = max(best_results, key=lambda x: x['min_acc'])
    best_p = best_p_row['p']
    print(f"\n{'='*90}")
    print(f"=== Best p = {best_p} (cross-model min_acc = {best_p_row['min_acc']:.3f}) ===")
    print(f"{'='*90}")

    col_o = f"omega_p_{best_p}"
    print(f"\n{'idx':>4} {'compound':>10} {'label':>22} {'MACE Ω_p':>10} {'PET Ω_p':>10} {'MACE ωmin':>10} {'PET ωmin':>10}")
    print("-" * 90)

    # Merge MACE and PET by idx
    merged = df_mace_gt[['idx', 'compound', 'label', col_o, 'omega_min']].merge(
        df_pet_gt[['idx', col_o, 'omega_min']], on='idx', suffixes=('_mace', '_pet'))

    for _, r in merged.sort_values(f'{col_o}_mace').iterrows():
        print(f"{int(r['idx']):>4} {r['compound']:>10} {r['label']:>22} "
              f"{r[f'{col_o}_mace']:>10.4f} {r[f'{col_o}_pet']:>10.4f} "
              f"{r['omega_min_mace']:>10.4f} {r['omega_min_pet']:>10.4f}")


if __name__ == "__main__":
    main()
