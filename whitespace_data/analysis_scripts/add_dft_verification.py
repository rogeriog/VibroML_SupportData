"""Add the `dft_verification` table to vibroml_database.db.

Populates one row per system in the 200-system DFT single-point
verification sample (seed 42) from:

- ``dft_whitespace_results.csv``   (energies; converged rows only)
- ``dft_whitespace_failures.csv``  (the 13 non-converged systems)
- ``dft_whitespace_force_stats.csv`` (per-structure force magnitudes,
  produced by ``parse_dft_forces.py``)
- ``steric_clash_flags.csv``       (clash flag of the top-ranked winner)

Columns written: material_id, run_dir, dataset, converged, clash_flag,
confirmed, e0_initial, e0_winner, natoms_initial, natoms_winner,
de_dft_per_atom, de_mlip_per_atom, fmax_initial, frms_initial,
fmax_winner, frms_winner, failure_note.

Usage:
    python add_dft_verification.py --db ../vibroml_database.db --data-dir ..
"""
import argparse
import sqlite3

import pandas as pd

SCHEMA = """
CREATE TABLE IF NOT EXISTS dft_verification (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    material_id INTEGER REFERENCES materials(id),
    run_dir TEXT,
    dataset TEXT,
    converged INTEGER,
    clash_flag INTEGER,
    confirmed INTEGER,
    e0_initial REAL,
    e0_winner REAL,
    natoms_initial INTEGER,
    natoms_winner INTEGER,
    de_dft_per_atom REAL,
    de_mlip_per_atom REAL,
    fmax_initial REAL,
    frms_initial REAL,
    fmax_winner REAL,
    frms_winner REAL,
    failure_note TEXT
);
CREATE INDEX IF NOT EXISTS idx_dft_verification_material
    ON dft_verification(material_id);
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True)
    ap.add_argument("--data-dir", default=".")
    args = ap.parse_args()

    res = pd.read_csv(f"{args.data_dir}/dft_whitespace_results.csv")
    try:
        fail = pd.read_csv(f"{args.data_dir}/dft_whitespace_failures.csv")
    except FileNotFoundError:
        fail = pd.DataFrame(columns=["run_dir", "material", "reason"])
    stats = pd.read_csv(f"{args.data_dir}/dft_whitespace_force_stats.csv")
    flags = pd.read_csv(f"{args.data_dir}/steric_clash_flags.csv")
    flagged = set(
        flags[(flags.struct_type == "top") & (flags["rank"] == 1)
              & (flags.status == "clashed")].material_name)

    conn = sqlite3.connect(args.db)
    conn.executescript("DROP TABLE IF EXISTS dft_verification;" + SCHEMA)
    mats = pd.read_sql_query("SELECT id, material_name FROM materials", conn)
    name2id = dict(zip(mats.material_name, mats.id))

    note_col = "failure_mode" if "failure_mode" in fail.columns else "reason"
    notes = (fail.dropna(subset=[note_col])
                 .groupby("run_dir")[note_col]
                 .apply(lambda s: "; ".join(sorted(set(s.astype(str)))))
                 .to_dict()) if note_col in fail.columns else {}

    rows = []
    for _, r in res.iterrows():
        conv = pd.notna(r.get("dE_dft_per_atom"))
        st = stats[stats.run_dir == r.run_dir]
        si = st[st.side == "initial"].iloc[0] if (st.side == "initial").any() and len(st[st.side == "initial"]) else None
        sw = st[st.side == "winner"].iloc[0] if (st.side == "winner").any() and len(st[st.side == "winner"]) else None
        clash = int(r.material in flagged)
        confirmed = int(r.dE_dft_per_atom < 0) if conv else None
        rows.append(dict(
            material_id=int(name2id.get(r.material, -1)),
            run_dir=r.run_dir, dataset=r.dataset, converged=int(conv),
            clash_flag=clash, confirmed=confirmed,
            e0_initial=r.get("E0_initial"), e0_winner=r.get("E0_winner"),
            natoms_initial=r.get("natoms_initial"),
            natoms_winner=r.get("natoms_winner"),
            de_dft_per_atom=r.get("dE_dft_per_atom"),
            de_mlip_per_atom=r.get("dE_mace"),
            fmax_initial=None if si is None else si.f_max,
            frms_initial=None if si is None else si.f_rms,
            fmax_winner=None if sw is None else sw.f_max,
            frms_winner=None if sw is None else sw.f_rms,
            failure_note=notes.get(r.run_dir)))

    df = pd.DataFrame(rows)
    df.to_sql("dft_verification", conn, if_exists="append", index=False)
    n = conn.execute("SELECT COUNT(*) FROM dft_verification").fetchone()[0]
    nc = conn.execute(
        "SELECT COUNT(*) FROM dft_verification WHERE converged=1").fetchone()[0]
    conn.commit()
    print(f"dft_verification: {n} rows ({nc} converged)")


if __name__ == "__main__":
    main()
