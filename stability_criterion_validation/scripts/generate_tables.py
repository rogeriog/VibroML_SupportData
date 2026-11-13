#!/usr/bin/env python3
"""
Generate LaTeX tables for the manuscript from analysis outputs.

Consumes the CSV files produced by analyze_stability_criterion.py and
writes LaTeX table files into the manuscript tables/ directory.

Usage:
  python generate_tables.py
  python generate_tables.py --analysis-dir analysis_outputs --tables-dir ../../03_paper_and_rebuttal/VIBROML_PAPER/tables
"""
import argparse, os, sys
import pandas as pd
import numpy as np

# ============================================================================
# Paths
# ============================================================================

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_ANALYSIS_DIR = os.environ.get(
    "VIBROML_ANALYSIS_DIR",
    os.path.join(_SCRIPT_DIR, "..", "..", "pet_phonons", "analysis_outputs_5eng"))
DEFAULT_TABLES_DIR = os.environ.get(
    "VIBROML_TABLES_DIR",
    os.path.join(_SCRIPT_DIR, "..", "..", "..", "03_paper_and_rebuttal",
                 "VIBROML_PAPER", "tables"))

# Pooled optimal separation (descriptive, not an operational threshold)
SOFT_FRAC_OPTIMAL_SEP = 0.045


# ============================================================================
# Table S_roc_auc: ROC/AUC per engine per metric
# ============================================================================

def generate_roc_auc_table(analysis_dir, tables_dir):
    """Generate TableS_roc_auc.tex from roc_auc_per_engine.csv."""
    csv_path = os.path.join(analysis_dir, "roc_auc_per_engine.csv")
    if not os.path.exists(csv_path):
        print(f"  SKIP: {csv_path} not found")
        return

    df = pd.read_csv(csv_path)
    engines = df['engine'].tolist()

    # Select metrics to display
    display_metrics = ['soft_frac', 'omega_p_0.5', 'i_curv', 'omega_p_2.0', 'omega_min']

    latex = []
    latex.append(r"\begin{table}[htbp]")
    latex.append(r"\centering")
    latex.append(r"\caption{ROC AUC for distinguishing stable from unstable structures "
                 r"across MLIP engines and scalar metrics. The soft-mode fraction "
                 r"($f_\mathrm{soft}$) achieves the highest AUC across all engines. "
                 r"Depth-weighted metrics ($\Omega_p$ with larger $p$, $\omega_\mathrm{min}$) "
                 r"are increasingly vulnerable to MLIP interpolation artifacts.}")
    latex.append(r"\label{tab:si_roc_auc}")
    latex.append(r"\begin{tabular}{l" + "c" * len(engines) + "}")
    latex.append(r"\toprule")
    latex.append("Metric & " + " & ".join(engines) + r" \\")
    latex.append(r"\midrule")

    metric_labels = {
        'soft_frac': r'$f_\mathrm{soft}$',
        'omega_p_0.5': r'$\Omega_{0.5}$',
        'i_curv': r'$I_\mathrm{curv}$ ($p=1$)',
        'omega_p_2.0': r'$\Omega_{2}$',
        'omega_min': r'$\omega_\mathrm{min}$',
    }

    for metric in display_metrics:
        if metric not in df.columns:
            continue
        label = metric_labels.get(metric, metric)
        vals = []
        for eng in engines:
            row = df[df['engine'] == eng]
            if len(row) > 0 and not pd.isna(row[metric].iloc[0]):
                vals.append(f"{row[metric].iloc[0]:.3f}")
            else:
                vals.append("—")
        latex.append(f"{label} & " + " & ".join(vals) + r" \\")

    latex.append(r"\bottomrule")
    latex.append(r"\end{tabular}")
    latex.append(r"\end{table}")

    out_path = os.path.join(tables_dir, "TableS_roc_auc.tex")
    with open(out_path, 'w') as f:
        f.write("\n".join(latex))
    print(f"  Wrote: {out_path}")


# ============================================================================
# Table S_marginal_mode: marginal mode sensitivity
# ============================================================================

def generate_marginal_mode_table(analysis_dir, tables_dir):
    """Generate TableS_marginal_mode.tex from marginal_mode_sensitivity_*.csv."""
    eng_files = []
    for f in os.listdir(analysis_dir):
        if f.startswith("marginal_mode_sensitivity_") and f.endswith(".csv"):
            eng_files.append(f)

    if not eng_files:
        print(f"  SKIP: no marginal_mode_sensitivity files found")
        return

    latex = []
    latex.append(r"\begin{table}[htbp]")
    latex.append(r"\centering")
    latex.append(r"\caption{Marginal-mode sensitivity analysis. "
                 r"``\% marginal'' = fraction of imaginary modes with $|\omega|<0.5$~THz. "
                 r"``\% contrib $I_{0.5}$'' and ``\% contrib $I_\mathrm{curv}$'' = "
                 r"fraction of the metric contributed by these marginal modes. "
                 r"The $p=1$ formulation is less sensitive to near-zero noise than $p=0.5$.}")
    latex.append(r"\label{tab:si_marginal}")
    latex.append(r"\begin{tabular}{lccc}")
    latex.append(r"\toprule")
    latex.append(r"Engine & \% marginal modes & \% contrib $I_{0.5}$ & \% contrib $I_\mathrm{curv}$ \\")
    latex.append(r"\midrule")

    for f in sorted(eng_files):
        eng = f.replace("marginal_mode_sensitivity_", "").replace(".csv", "")
        df = pd.read_csv(os.path.join(analysis_dir, f))
        if len(df) == 0:
            continue
        avg_marg = 100 * df['frac_marginal_modes'].mean()
        avg_p05 = 100 * df['frac_contrib_p05'].mean()
        avg_p10 = 100 * df['frac_contrib_p10'].mean()
        latex.append(f"{eng} & {avg_marg:.1f}\\% & {avg_p05:.1f}\\% & {avg_p10:.1f}\\% \\\\")

    latex.append(r"\bottomrule")
    latex.append(r"\end{tabular}")
    latex.append(r"\end{table}")

    out_path = os.path.join(tables_dir, "TableS_marginal_mode.tex")
    with open(out_path, 'w') as f:
        f.write("\n".join(latex))
    print(f"  Wrote: {out_path}")


# ============================================================================
# Main
# ============================================================================

def main():
    parser = argparse.ArgumentParser(description="Generate LaTeX tables from analysis outputs")
    parser.add_argument("--analysis-dir", default=DEFAULT_ANALYSIS_DIR)
    parser.add_argument("--tables-dir", default=DEFAULT_TABLES_DIR)
    args = parser.parse_args()

    os.makedirs(args.tables_dir, exist_ok=True)

    print(f"Analysis dir: {args.analysis_dir}")
    print(f"Tables dir: {args.tables_dir}")
    print()

    generate_roc_auc_table(args.analysis_dir, args.tables_dir)
    generate_marginal_mode_table(args.analysis_dir, args.tables_dir)

    print("\nDone.")


if __name__ == "__main__":
    main()
