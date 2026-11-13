#!/usr/bin/env python3
"""
Generate the reference-label provenance table (SI Table: tab:si_label_basis).

For every structure in the ground-truth reference set, record:
  - the labelled space group (the basis of the "stable"/"unstable" claim)
  - the space group actually detected in the evaluated geometry
  - the basis tier of the label

Tiers
-----
A : the space group is an experimentally observed polymorph of that compound
    -> defensible external anchor for the "dynamically stable" label
B : the space group is discussed in the literature for that compound but its
    experimental status is less clear -> flagged for author confirmation
C : the space group is, to our knowledge, a hypothetical polymorph recovered by
    the search and not experimentally reported -> NOT a valid external anchor
X : metadata inconsistency; the entry should not be used as ground truth

NOTE: the A/B/C assignment is a crystallographic judgement and is recorded here
explicitly so that it can be reviewed and overridden. Edit TIERS below.

Usage:
  python generate_label_basis_table.py [--tables-dir DIR]
"""

import argparse, json, os

HERE = os.path.dirname(os.path.abspath(__file__))

# --- label-basis assignment (review and edit as needed) ----------------------
TIERS = {
    # Tier A: experimentally observed polymorph
    16: ("A", "$\\delta$/$\\gamma$-\\ce{CsPbI3}, experimentally reported"),
    21: ("A", "$\\delta$/$\\gamma$-\\ce{CsPbI3}, experimentally reported"),
    28: ("A", "monoclinic baddeleyite, ambient phase of \\ce{HfO2}"),
    33: ("A", "monoclinic baddeleyite, ambient phase of \\ce{HfO2}"),
    29: ("A", "polar orthorhombic (o-III) ferroelectric \\ce{HfO2}"),
    39: ("A", "$\\alpha$-\\ce{Bi2Sn2O7}, low-temperature phase"),
    43: ("A", "$\\alpha$-\\ce{Bi2Sn2O7}, low-temperature phase"),
    # Tier B: literature-discussed, experimental status to be confirmed
    23: ("B", "polar orthorhombic \\ce{HfO2}; confirm experimental status"),
    25: ("B", "polar \\ce{HfO2} phase; confirm experimental status"),
    32: ("B", "not a standard \\ce{HfO2} polymorph; confirm"),
    35: ("B", "$\\beta$-\\ce{Bi2Sn2O7} symmetry is debated; confirm"),
    42: ("B", "$\\beta$-\\ce{Bi2Sn2O7} symmetry is debated; confirm"),
    40: ("B", "confirm experimental status"),
    # Tier C: hypothetical, recovered by search only
    2:  ("C", "\\ce{LiF} retains B1 to very high pressure; not reported"),
    3:  ("C", "zincblende \\ce{LiF}; hypothetical"),
    4:  ("C", "hypothetical"),
    5:  ("C", "wurtzite \\ce{LiF}; hypothetical"),
    7:  ("C", "hypothetical"),
    8:  ("C", "wurtzite \\ce{LiF}; hypothetical"),
    # Tier X: metadata inconsistency -- see note in the caption
    6:  ("X", "labelled $Fm\\overline{3}m$ but evaluated geometry is $P1$; "
              "recorded energy ($-8.61$~eV) is far below the physical "
              "polymorphs ($\\approx-4.82$~eV) and clusters with the "
              "steric-clash anomalies"),
}

TIER_DESC = {
    "A": "experimentally observed polymorph",
    "B": "literature-discussed, status to confirm",
    "C": "hypothetical, recovered by search only",
    "X": "metadata inconsistency -- excluded",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--structures",
                    default=os.path.join(HERE, "uq_validation_structures.json"))
    ap.add_argument("--tables-dir", default=os.path.join(HERE, "tables_5eng"))
    args = ap.parse_args()
    os.makedirs(args.tables_dir, exist_ok=True)

    with open(args.structures) as f:
        S = {s["order"]: s for s in json.load(f)
             if s.get("category") == "benchmark"}

    pos = [o for o, s in S.items()
           if s.get("stability_label") == "high_sym_unstable"]
    neg = [o for o, s in S.items()
           if s.get("stability_label") == "literature_stable"]

    L = []
    L.append(r"\begin{table}[htbp]")
    L.append(r"\centering\footnotesize")
    L.append(r"\caption{\textbf{Provenance of every reference-set label.} "
             r"``SG (labelled)'' is the space group on which the stability "
             r"label is based; ``SG (evaluated)'' is the space group detected "
             r"in the geometry actually used for the phonon calculation. "
             r"Basis tiers: A, the space group is an experimentally observed "
             r"polymorph of that compound; B, the space group is discussed in "
             r"the literature but its experimental status is less certain; C, "
             r"the space group is a hypothetical polymorph recovered by the "
             r"search and is therefore not an independent anchor for a "
             r"``dynamically stable'' label; X, the entry contains a metadata "
             r"inconsistency and is excluded from the reference set. No entry "
             r"has had its dynamical stability confirmed by DFT phonons in "
             r"this work (Section~\ref{sec:si_label_provenance}).}")
    L.append(r"\label{tab:si_label_basis}")
    L.append(r"\begin{tabular}{r l l l c p{5.4cm}}")
    L.append(r"\toprule")
    L.append(r"Idx & Compound & SG (labelled) & SG (evaluated) & Basis & Note \\")
    L.append(r"\midrule")
    L.append(r"\multicolumn{6}{l}{\textit{Known-unstable aristotypes "
             r"(positives)}} \\")
    for o in sorted(pos):
        s = S[o]
        L.append(f"{o} & \\ce{{{s['compound']}}} & ${s.get('space_group_labeled','?')}$ "
                 f"& {_fmt_sg(s.get('space_group_actual'))} & A & "
                 f"high-symmetry aristotype, unstable at ambient \\\\")
    L.append(r"\midrule")
    L.append(r"\multicolumn{6}{l}{\textit{Reported lower-symmetry polymorphs "
             r"(negatives)}} \\")
    for o in sorted(neg):
        s = S[o]
        tier, note = TIERS.get(o, ("?", "unassigned"))
        L.append(f"{o} & \\ce{{{s['compound']}}} & ${s.get('space_group_labeled','?')}$ "
                 f"& {_fmt_sg(s.get('space_group_actual'))} & {tier} & {note} \\\\")
    L.append(r"\bottomrule")
    L.append(r"\end{tabular}")
    L.append(r"\end{table}")

    out = os.path.join(args.tables_dir, "TableS_label_basis.tex")
    with open(out, "w") as f:
        f.write("\n".join(L) + "\n")
    print(f"Wrote: {out}")

    # console summary
    from collections import Counter
    c = Counter(TIERS.get(o, ("?", ""))[0] for o in neg)
    print("\nNegative-label basis breakdown:")
    for t in "ABCX?":
        if c.get(t):
            print(f"  Tier {t} ({TIER_DESC.get(t,'?')}): {c[t]}")
    usable = c.get("A", 0)
    print(f"\nReference set sizes:")
    print(f"  as published      : {len(pos)} positives + {len(neg)} negatives")
    print(f"  Tier A only       : {len(pos)} + {usable}")
    print(f"  Tier A+B          : {len(pos)} + {usable + c.get('B',0)}")
    print(f"  excluding Tier X  : {len(pos)} + {len(neg) - c.get('X',0)}")


def _fmt_sg(v):
    if v is None:
        return "--"
    return str(v).replace("_", r"\_")


if __name__ == "__main__":
    main()
