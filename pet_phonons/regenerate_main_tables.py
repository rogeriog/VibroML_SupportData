#!/usr/bin/env python3
"""
Regenerate Table 1 and Table S1 from the uniform 6x6x6 mesh data,
reporting continuous metrics (f_soft, I_curv, omega_min) instead of
categorical stability labels.

Preserves the DFT energy columns and remediation-algorithm checkmarks in Table 1
(these are from the original VibroML search and are not affected by the metric change).
"""
import json, os, re, shutil

HERE = os.path.dirname(os.path.abspath(__file__))
TABLES_DIR = os.path.join(HERE, "..", "..", "03_paper_and_rebuttal",
                         "VIBROML_PAPER", "tables")
TABLES_DIR = os.path.abspath(TABLES_DIR)


def load_ip_variants(engine):
    fn = os.path.join(HERE, f"ip_variants_{engine}.json")
    with open(fn) as f:
        data = json.load(f)
    return {r["idx"]: r for r in data}


def load_structures():
    fn = os.path.join(HERE, "uq_validation_structures.json")
    with open(fn) as f:
        data = json.load(f)
    return {s["order"]: s for s in data}


# ---------------------------------------------------------------------------
# Table 1: Stability Remediation Benchmark (MACE-OMAT only)
# ---------------------------------------------------------------------------

# Hand-curated DFT energies and remediation algorithm assignments from the
# original VibroML search.  The f_soft and I_curv columns are regenerated from
# the uniform-mesh data.
#
# Each entry: (compound, natoms_label, space_group_latex, dE_mace, dE_init,
#              dE_relax, idx_in_44set_or_None, ga, trad, all_, rand)
#
# idx_in_44set maps to ip_variants_mace.json.  Where the Table 1 entry does
# not correspond to a 44-set structure (e.g. Fm-3m LiF which is the rocksalt
# ground state, not the anomalous idx=6), idx=None and we fall back to the
# original mace_soft_frac from uq_validation_structures.json.

TABLE1_ROWS = [
    # LiF
    ("LiF", 2, r"Pm\bar{3}m", 0.0, 0.0, 0.0, 1, "ref", "ref", "ref", "ref"),
    ("", 2, r"Pnma", -276.8, -275.4, -275.4, 2, r"\checkmark", "", "", ""),
    ("", 2, r"F\bar{4}3m", -292.5, -289.9, -289.9, 3, r"\checkmark", "", "", ""),
    ("", 2, r"P4_2/mnm", -296.4, -294.5, -294.5, 4, r"\checkmark", "", "", ""),
    ("", 2, r"Fm\bar{3}m", -297.2, -296.1, -296.1, None, r"\checkmark", r"\checkmark", r"\checkmark", r"\checkmark"),
    ("", 2, r"P6_3mc", -303.1, -300.7, -300.6, 8, r"\checkmark", "", r"\checkmark", r"\checkmark"),
    ("", 2, r"P1", -4088.3, -3121.4, "--", 9, "", "", "", r"\checkmark"),
    ("", 2, r"Pm", -7347.9, -7691.9, "--", 10, "", "", "", r"\checkmark"),
    ("", 2, r"Cm", -14042.8, -15170.5, "--", 11, "", "", "", r"\checkmark"),
    # CsPbI3
    ("CsPbI$_3$", 5, r"Pm\bar{3}m", 0.0, 0.0, -5.4, 44, "ref", "ref", "ref", "ref"),
    ("", 5, r"C2/m", -18.1, -21.5, -21.5, 13, "", r"\checkmark", "", ""),
    ("", 5, r"Cm", -23.0, -23.0, -23.5, 14, "", r"\checkmark", r"\checkmark", ""),
    ("", 5, r"Pnma", -29.3, -29.1, -27.9, 16, r"\checkmark", "", r"\checkmark", ""),
    ("", 5, r"P\bar{1}", -29.9, -28.1, -28.5, 15, "", "", "", r"\checkmark"),
    ("", 5, r"P2_1/m", -33.3, -26.2, -46.7, 19, r"\checkmark", "", "", ""),
    ("", 5, r"Cmcm", -34.9, -26.7, -26.8, 20, r"\checkmark", "", "", ""),
    ("", 5, r"Pnma", -44.2, -37.9, -38.0, 21, r"\checkmark", "", "", ""),
    # HfO2
    ("HfO$_2$", 12, r"Fm\bar{3}m", 0.0, 0.0, -5.4, 22, "ref", "ref", "ref", "ref"),
    ("", 12, r"P4_2/nmc", -12.2, -38.0, -39.7, 24, "", r"\checkmark", r"\checkmark", r"\checkmark"),
    ("", 12, r"P2_1/m", -18.5, -68.0, -68.8, 26, r"\checkmark", "", "", r"\checkmark"),
    ("", 12, r"Pnma", -21.7, -66.0, -66.4, 27, r"\checkmark", "", "", r"\checkmark"),
    ("", 12, r"Pbcn", -38.5, -79.4, -79.7, 31, r"\checkmark", r"\checkmark", "", ""),
    ("", 12, r"P2_1/c", -48.5, -94.4, -95.1, 28, r"\checkmark", "", "", r"\checkmark"),
    # Bi2Sn2O7
    ("Bi$_2$Sn$_2$O$_7$", 88, r"Fd\bar{3}m", 0.0, 0.0, -10.2, 34, "ref", "ref", "ref", "ref"),
    ("", 88, r"P2_12_12_1", 20.7, 6.1, 5.4, 35, r"\checkmark", "", "", ""),
    ("", 88, r"P2_1", 12.7, 1.0, 0.5, 36, r"\checkmark", "", "", r"\checkmark"),
    ("", 88, r"Cc", -12.6, -23.4, -24.1, 39, "", r"\checkmark", r"\checkmark", ""),
    ("", 88, r"Pna2_1", -12.7, -23.7, -24.0, 40, r"\checkmark", "", r"\checkmark", ""),
    ("", 88, r"P2_1", -13.1, -24.1, -24.4, 41, r"\checkmark", "", r"\checkmark", ""),
    ("", 88, r"P2_12_12_1", -13.7, -24.5, -24.8, 42, r"\checkmark", "", r"\checkmark", r"\checkmark"),
    ("", 88, r"Cc", -14.1, -25.7, -25.8, 43, r"\checkmark", r"\checkmark", r"\checkmark", ""),
]


def generate_table1(mace_data, structures):
    """Regenerate Table 1 with continuous f_soft and I_curv metrics."""
    lines = []
    lines.append("")
    lines.append("")
    lines.append("% Table 1: Stability Remediation Benchmark "
                 "(regenerated from uniform 6x6x6 mesh, continuous metrics)")
    lines.append(r"\begin{table}[htbp]")
    lines.append(r"    \centering")
    lines.append(r"    \caption{\textbf{Benchmark of VibroML instability remediation "
                 r"algorithms across four test systems using MACE-OMAT-0 (medium) "
                 r"model as MLIP.} The table compares the relative energies "
                 r"($\Delta E$) of the most relevant polymorphs identified starting "
                 r"from high-symmetry unstable phases. Energies are reported in "
                 r"meV/atom relative to the initial unstable structure (set to 0.0). "
                 r"$\Delta E_{f,\text{MACE}}$ denotes the energy predicted by the "
                 r"MLIP, while $\Delta E_{f,\text{init}}$ and "
                 r"$\Delta E_{f,\text{relax}}$ correspond to DFT verification values "
                 r"for the initial and relaxed structures, respectively. Dynamical "
                 r"stability is quantified by two continuous descriptors: the "
                 r"soft-mode fraction $f_\mathrm{soft}$ (fraction of imaginary modes "
                 r"on a uniform $6\times6\times6$ q-point mesh) and the mean negative "
                 r"curvature $I_\mathrm{curv}$ (in THz$^2$). The columns under "
                 r"``Remediation Algorithm'' indicate which method successfully "
                 r"recovered the specific polymorph: Genetic Algorithm (GA), Single "
                 r"Mode Following (TRAD), Coupled Mode Following (ALL), and "
                 r"Stochastic Search (RAND).}")
    lines.append(r"    \label{tab:stability_remediation}")
    lines.append("")
    lines.append(r"\setlength{\tabcolsep}{0pt}")
    lines.append(r"\renewcommand{\arraystretch}{1.3}")
    lines.append(r"\footnotesize")
    lines.append(r"\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}} l c c c c c c c c c c }")
    lines.append(r"    \toprule")
    lines.append(r"    \textbf{Compound} & ")
    lines.append(r"    \multirow{2}{*}{\textbf{Space Group}} & ")
    lines.append(r"    \textbf{$\Delta E_{f,\text{MACE}}$} & ")
    lines.append(r"    \textbf{$\Delta E_{f,\text{init}}$} & ")
    lines.append(r"    \textbf{$\Delta E_{f,\text{relax}}$} & ")
    lines.append(r"    \textbf{$f_\mathrm{soft}$} & ")
    lines.append(r"    \textbf{$I_\mathrm{curv}$} & ")
    lines.append(r"    \multicolumn{4}{c}{\textbf{Remediation Algorithm}} \\")
    lines.append(r"    \cmidrule(lr){8-11}")
    lines.append(r"    \textbf{(\# atoms)} & & (meV/at) & (meV/at) & (meV/at) & (\%) & (THz$^2$) & \textbf{GA} & \textbf{TRAD} & \textbf{ALL} & \textbf{RAND} \\")
    lines.append(r"    \midrule")

    _pending_natoms = None
    for row in TABLE1_ROWS:
        (compound, natoms, sg, dE_mace, dE_init, dE_relax,
         idx, ga, trad, all_, rand) = row

        # Determine soft_frac and i_curv
        if idx is not None and idx in mace_data:
            sf = mace_data[idx]["soft_frac_no_thresh"]
            icurv = mace_data[idx].get("i_p_1.0", 0.0)
        elif idx is not None and idx in structures:
            sf = structures[idx].get("mace_soft_frac", 0.0)
            icurv = 0.0
        else:
            sf = 0.0
            icurv = 0.0

        sf_str = f"{sf*100:.1f}"
        icurv_str = f"{icurv:.2f}"

        # Format: first row of a compound has the compound name;
        # second row has the natoms label; subsequent rows are blank in col 1.
        if compound:
            comp_cell = f"\\textbf{{{compound}}}"
            natoms_cell = f"\\textbf{{({natoms})}}"
            lines.append(f"    {comp_cell} & ${sg}$ & {dE_mace} & {dE_init} & {dE_relax} & {sf_str} & {icurv_str} & {ga} & {trad} & {all_} & {rand} \\\\")
            _pending_natoms = natoms_cell
        else:
            prefix = _pending_natoms if _pending_natoms else " "
            lines.append(f"    {prefix} & ${sg}$ & {dE_mace} & {dE_init} & {dE_relax} & {sf_str} & {icurv_str} & {ga} & {trad} & {all_} & {rand} \\\\")
            _pending_natoms = None

        # Check if next row is a new compound
        next_idx = TABLE1_ROWS.index(row) + 1
        if next_idx < len(TABLE1_ROWS) and TABLE1_ROWS[next_idx][0]:
            lines.append(r"    \midrule")

    lines.append(r"    \bottomrule")
    lines.append(r"\end{tabular*}")
    lines.append(r"\end{table}")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Table S1: Full comparison across all 5 engines
# ---------------------------------------------------------------------------

ENGINE_ORDER = ["mace", "macemp", "esen", "uma", "pet"]
ENGINE_LABELS = {
    "mace": "MACE-OMAT",
    "macemp": "MACE-MP",
    "esen": "eSEN",
    "uma": "UMA",
    "pet": "PET-MAD",
}

# Compound ordering and structure selection for Table S1
# (same structures as Table 1, but showing all 5 engines)
TABLE_S1_GROUPS = [
    ("LiF", [(1, r"Pm\bar{3}m"), (2, r"Pnma"), (3, r"F\bar{4}3m"),
             (4, r"P4_2/mnm"), (5, r"P6_3mc"), (6, r"Fm\bar{3}m"),
             (9, r"P1"), (10, r"Pm"), (11, r"Cm")]),
    ("CsPbI$_3$", [(44, r"Pm\bar{3}m"), (13, r"C2/m"), (14, r"Cm"),
                   (16, r"Pnma"), (15, r"P\bar{1}"), (19, r"P2_1/m"),
                   (20, r"Cmcm"), (21, r"Pnma")]),
    ("HfO$_2$", [(22, r"Fm\bar{3}m"), (24, r"P4_2/nmc"), (26, r"P2_1/m"),
                 (27, r"Pnma"), (31, r"Pbcn"), (28, r"P2_1/c")]),
    ("Bi$_2$Sn$_2$O$_7$", [(34, r"Fd\bar{3}m"), (35, r"P2_12_12_1"),
                           (36, r"P2_1"), (39, r"Cc"), (40, r"Pna2_1"),
                           (41, r"P2_1"), (42, r"P2_12_12_1"), (43, r"Cc")]),
]


def generate_table_s1(all_engine_data):
    """Regenerate Table S1 with continuous metrics, all 5 engines."""
    lines = []
    lines.append(r"\begin{landscape}")
    lines.append(r"\small")
    lines.append(r"\setlength{\tabcolsep}{5pt}")
    lines.append(r"\renewcommand{\arraystretch}{1.2}")
    lines.append(r"\setlength{\LTcapwidth}{\linewidth}")
    lines.append(r"\begin{longtable}{ll l c c c}")
    lines.append(r"\caption{\textbf{Comprehensive benchmarking of dynamical "
                 r"stability across five foundational MLIPs.} The data is grouped "
                 r"by compound, starting with the reference high-symmetry "
                 r"aristotype, followed by the remaining polymorphs sorted by "
                 r"MACE-OMAT relative energy from highest to lowest. "
                 r"(1) \textbf{Engine}: the MLIP used for force-constant "
                 r"evaluation; (2) \textbf{Soft Frac} ($f_\mathrm{soft}$): "
                 r"the proportion of imaginary phonon frequencies integrated "
                 r"across the sampled Brillouin zone on a uniform "
                 r"$6\times6\times6$ q-point mesh; (3) \textbf{$I_\mathrm{curv}$}: "
                 r"the mean negative curvature in THz$^2$, quantifying the average "
                 r"depth of the imaginary modes; (4) \textbf{Soft Freq} "
                 r"($\omega_\mathrm{min}$): the lowest (most negative) phonon "
                 r"frequency in THz. "
                 r"\label{tab:fullcomp_MLIPs_phonon_evals}} \\")
    lines.append(r"\toprule")
    lines.append(r"\textbf{Compound} & \textbf{Space Group} & \textbf{Engine} & "
                 r"\textbf{Soft Frac} & \textbf{$I_\mathrm{curv}$} & "
                 r"\textbf{Soft Freq} \\")
    lines.append(r" & & & (\%) & (THz$^2$) & (THz) \\")
    lines.append(r"\midrule")
    lines.append(r"\endfirsthead")
    lines.append(r"\multicolumn{6}{c}{{\bfseries \tablename\ \thetable{} "
                 r"-- continued from previous page}} \\")
    lines.append(r"\toprule")
    lines.append(r"\textbf{Comp} & \textbf{SG} & \textbf{Engine} & "
                 r"\textbf{Soft Frac} & \textbf{$I_\mathrm{curv}$} & "
                 r"\textbf{Soft Freq} \\")
    lines.append(r"\midrule")
    lines.append(r"\endhead")
    lines.append(r"\bottomrule")
    lines.append(r"\multicolumn{6}{r}{{Continued on next page}} \\")
    lines.append(r"\endfoot")
    lines.append(r"\bottomrule")
    lines.append(r"\endlastfoot")

    for compound, entries in TABLE_S1_GROUPS:
        n_entries = len(entries)
        lines.append(f"\\multirow{{{n_entries*5}}}{{*}}{{\\large \\textbf{{{compound}}}}}")

        for i, (idx, sg) in enumerate(entries):
            if i > 0:
                lines.append("")
            lines.append(f"& \\multirow{{5}}{{*}}{{${sg}$}}")

            for j, eng in enumerate(ENGINE_ORDER):
                data = all_engine_data[eng]
                if idx in data:
                    r = data[idx]
                    sf = r["soft_frac_no_thresh"]
                    om = r["omega_min"]
                    icurv = r.get("i_p_1.0", 0.0)
                    sf_str = f"{sf*100:.1f}"
                    icurv_str = f"{icurv:.2f}"
                    om_str = f"{om:.2f}"
                else:
                    sf_str = "—"
                    icurv_str = "—"
                    om_str = "—"

                eng_label = ENGINE_LABELS[eng]
                if j < len(ENGINE_ORDER) - 1:
                    lines.append(f"  & {eng_label} & {sf_str} & {icurv_str} & {om_str} \\\\*")
                else:
                    lines.append(f"  & {eng_label} & {sf_str} & {icurv_str} & {om_str} \\\\")

            if i < len(entries) - 1:
                lines.append(r"\cmidrule{2-6}")

        lines.append(r"\midrule")

    # Remove last midrule and replace with bottomrule
    if lines[-1] == r"\midrule":
        lines.pop()

    lines.append(r"\end{longtable}")
    lines.append(r"\end{landscape}")

    return "\n".join(lines)


def main():
    print("Loading data...")
    mace_data = load_ip_variants("mace")
    structures = load_structures()
    all_engine_data = {eng: load_ip_variants(eng) for eng in ENGINE_ORDER}

    # Backup originals
    for name in ["table1_stability_remediation.tex",
                 "TableS1_fullcomp_MLIPs_phonon_evals.tex"]:
        src = os.path.join(TABLES_DIR, name)
        bak = src + ".bak_legacy"
        if os.path.exists(src) and not os.path.exists(bak):
            shutil.copy2(src, bak)
            print(f"  Backed up: {bak}")

    # Generate Table 1
    print("\nGenerating Table 1...")
    t1 = generate_table1(mace_data, structures)
    t1_path = os.path.join(TABLES_DIR, "table1_stability_remediation.tex")
    with open(t1_path, "w") as f:
        f.write(t1 + "\n")
    print(f"  Wrote: {t1_path}")

    # Generate Table S1
    print("\nGenerating Table S1...")
    ts1 = generate_table_s1(all_engine_data)
    ts1_path = os.path.join(TABLES_DIR, "TableS1_fullcomp_MLIPs_phonon_evals.tex")
    with open(ts1_path, "w") as f:
        f.write(ts1 + "\n")
    print(f"  Wrote: {ts1_path}")

    print("\nDone.")


if __name__ == "__main__":
    main()
