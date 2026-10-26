#!/usr/bin/env python3
"""
VibroML GA Phonon Results Explorer - Streamlit Web Interface

Reads external CIF, band-structure, DOS, settings, and phonon-mode trajectory
files either directly from the filesystem (when running on the raw directory
tree) or seamlessly from the compact `database/vibroml_blobs.db` SQLite bundle
downloaded from Figshare.
"""

import os
import json
import zlib
import random
import sqlite3
from pathlib import Path

import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from ase.data import chemical_symbols

BASE_DIR = Path(os.environ.get("VIBROML_DATA_DIR", Path(__file__).resolve().parent))
DATABASE_PATH = BASE_DIR / "database" / "vibroml_database.db"
if not DATABASE_PATH.exists() and (BASE_DIR / "vibroml_database.db").exists():
    DATABASE_PATH = BASE_DIR / "vibroml_database.db"

BLOBS_DB_PATH = BASE_DIR / "database" / "vibroml_blobs.db"
if not BLOBS_DB_PATH.exists() and (BASE_DIR / "vibroml_blobs.db").exists():
    BLOBS_DB_PATH = BASE_DIR / "vibroml_blobs.db"


def to_absolute_path(rel_path):
    """Convert relative path to absolute path (prepend BASE_DIR)."""
    if rel_path is None:
        return None
    return str(BASE_DIR / rel_path)


def _blob_conn():
    if BLOBS_DB_PATH.exists():
        return sqlite3.connect(f"file:{BLOBS_DB_PATH}?mode=ro", uri=True)
    return None


def file_available(rel_path, table="text_files"):
    """Return True if rel_path exists on disk or inside vibroml_blobs.db."""
    if not rel_path:
        return False
    p = Path(to_absolute_path(rel_path))
    if p.exists():
        return True
    conn = _blob_conn()
    if conn is None:
        return False
    try:
        c = conn.cursor()
        c.execute(f"SELECT 1 FROM {table} WHERE rel_path = ? LIMIT 1", (str(rel_path),))
        return c.fetchone() is not None
    except Exception:
        return False
    finally:
        conn.close()


def read_text_resource(rel_path):
    """Read UTF-8 text from disk or from vibroml_blobs.db."""
    if not rel_path:
        return None
    p = Path(to_absolute_path(rel_path))
    if p.exists():
        try:
            return p.read_text()
        except Exception:
            pass
    conn = _blob_conn()
    if conn is None:
        return None
    try:
        c = conn.cursor()
        c.execute("SELECT data_z FROM text_files WHERE rel_path = ?", (str(rel_path),))
        row = c.fetchone()
        if row and row[0]:
            return zlib.decompress(row[0]).decode("utf-8", errors="replace")
    except Exception:
        pass
    finally:
        conn.close()
    return None


def read_array_resource(rel_path):
    """Read 2D numeric array from disk or from vibroml_blobs.db."""
    if not rel_path:
        return None
    p = Path(to_absolute_path(rel_path))
    if p.exists():
        try:
            return np.loadtxt(p)
        except Exception:
            pass
    conn = _blob_conn()
    if conn is None:
        return None
    try:
        c = conn.cursor()
        c.execute("SELECT shape_0, shape_1, data_z FROM array_files WHERE rel_path = ?", (str(rel_path),))
        row = c.fetchone()
        if row and row[2]:
            s0, s1, data_z = row
            arr = np.frombuffer(zlib.decompress(data_z), dtype=np.float32)
            return arr.reshape(s0, s1) if s1 > 1 else arr.reshape(s0)
    except Exception:
        pass
    finally:
        conn.close()
    return None


def find_kpoints_text(bs_rel_path):
    """Find special_k_points_*.txt content in the same folder as bs_rel_path."""
    if not bs_rel_path:
        return None
    p = Path(to_absolute_path(bs_rel_path))
    if p.parent.exists():
        kpts = list(p.parent.glob("special_k_points_*.txt"))
        if kpts:
            try:
                return kpts[0].read_text()
            except Exception:
                pass
    conn = _blob_conn()
    if conn is None:
        return None
    try:
        parent_rel = str(Path(bs_rel_path).parent)
        c = conn.cursor()
        c.execute(
            "SELECT data_z FROM text_files WHERE rel_path LIKE ? LIMIT 1",
            (f"{parent_rel}/special_k_points_%.txt",),
        )
        row = c.fetchone()
        if row and row[0]:
            return zlib.decompress(row[0]).decode("utf-8", errors="replace")
    except Exception:
        pass
    finally:
        conn.close()
    return None


def load_traj_frames_xyz(traj_rel_path):
    """Return multi-frame XYZ string for a phonon mode .traj from disk or vibroml_blobs.db."""
    if not traj_rel_path:
        return None
    p = Path(to_absolute_path(traj_rel_path))
    if p.exists():
        try:
            from ase.io import read
            traj = read(str(p), index=":")
            frames_xyz = ""
            for i, atoms in enumerate(traj):
                syms = atoms.get_chemical_symbols()
                pos_str = "\n".join(
                    f"{sym} {x:.6f} {y:.6f} {z:.6f}"
                    for sym, (x, y, z) in zip(syms, atoms.get_positions())
                )
                frames_xyz += f"{len(atoms)}\nFrame {i}\n{pos_str}\n"
            return frames_xyz
        except Exception:
            pass

    conn = _blob_conn()
    if conn is None:
        return None
    try:
        c = conn.cursor()
        c.execute(
            "SELECT n_frames, n_atoms, numbers, r0, uc, us FROM traj_files WHERE rel_path = ?",
            (str(traj_rel_path),),
        )
        row = c.fetchone()
        if not row:
            return None
        n_frames, n_atoms, num_z, r0_z, uc_z, us_z = row
        numbers = np.frombuffer(zlib.decompress(num_z), dtype=np.uint8)
        syms = [chemical_symbols[int(z)] for z in numbers]
        if uc_z is not None and us_z is not None:
            r0 = np.frombuffer(zlib.decompress(r0_z), dtype=np.float32).reshape(n_atoms, 3)
            uc = np.frombuffer(zlib.decompress(uc_z), dtype=np.float32).reshape(n_atoms, 3)
            us = np.frombuffer(zlib.decompress(us_z), dtype=np.float32).reshape(n_atoms, 3)
            angles = 2.0 * np.pi * np.arange(n_frames) / n_frames
            all_pos = (
                r0[None, :, :]
                + np.cos(angles)[:, None, None] * uc[None, :, :]
                + np.sin(angles)[:, None, None] * us[None, :, :]
            )
        else:
            all_pos = np.frombuffer(zlib.decompress(r0_z), dtype=np.float32).reshape(n_frames, n_atoms, 3)

        frames_xyz = []
        for i in range(n_frames):
            pos_i = all_pos[i]
            pos_str = "\n".join(
                f"{sym} {x:.6f} {y:.6f} {z:.6f}"
                for sym, (x, y, z) in zip(syms, pos_i)
            )
            frames_xyz.append(f"{n_atoms}\nFrame {i}\n{pos_str}\n")
        return "".join(frames_xyz)
    except Exception:
        return None
    finally:
        conn.close()


JMOL_COLORS = {
    'H':'#FFFFFF','He':'#D9FFFF','Li':'#CC80FF','Be':'#C2FF00','B':'#FFB5B5',
    'C':'#909090','N':'#3050F8','O':'#FF0D0D','F':'#90E050','Ne':'#B3E3F5',
    'Na':'#AB5CF2','Mg':'#8AFF00','Al':'#BFA6A6','Si':'#F0C8A0','P':'#FF8000',
    'S':'#FFFF30','Cl':'#1FF01F','Ar':'#80D1E3','K':'#8F40D4','Ca':'#3DFF00',
    'Ti':'#BFC2C7','V':'#A6A6AB','Cr':'#8A99C7','Mn':'#9C7AC7','Fe':'#E06633',
    'Co':'#F090A0','Ni':'#50D050','Cu':'#C88033','Zn':'#7D80B0','Ga':'#C28F8F',
    'Ge':'#668F8F','As':'#BD80E3','Se':'#FFA100','Br':'#A62929','Kr':'#5CB8D1',
    'Rb':'#702EB0','Sr':'#00FF00','Y':'#94FFFF','Zr':'#94E0E0','Nb':'#73C2C9',
    'Mo':'#54B5B5','Tc':'#3B9E9E','Ru':'#248F8F','Rh':'#0A7D8C','Pd':'#006985',
    'Ag':'#C0C0C0','Cd':'#FFD98F','In':'#A67573','Sn':'#668080','Sb':'#9E63B5',
    'Te':'#D47A00','I':'#940094','Cs':'#57178F','Ba':'#00C900','La':'#70D4FF',
    'Hf':'#4DC2FF','Ta':'#4DA6FF','W':'#2194D6','Re':'#267DAB','Os':'#266696',
    'Ir':'#175487','Pt':'#D0D0E0','Au':'#FFD123','Hg':'#B8B8D0','Pb':'#575961',
    'Bi':'#9E4FB5',
}


@st.cache_data
def load_data():
    conn = sqlite3.connect(DATABASE_PATH)
    df = pd.read_sql_query("SELECT * FROM materials", conn)
    conn.close()
    return df


def load_ga_structures(material_id):
    conn = sqlite3.connect(DATABASE_PATH)
    df = pd.read_sql_query(
        "SELECT * FROM ga_structures WHERE material_id = ? ORDER BY CASE struct_type WHEN 'top' THEN 0 ELSE 1 END, rank",
        conn, params=(material_id,)
    )
    conn.close()
    return df


def load_phonon_data(material_id):
    conn = sqlite3.connect(DATABASE_PATH)
    df = pd.read_sql_query(
        "SELECT * FROM phonon_data WHERE material_id = ? ORDER BY CASE data_type WHEN 'initial' THEN 0 WHEN 'final' THEN 1 WHEN 'unique' THEN 2 ELSE 3 END, iteration",
        conn, params=(material_id,)
    )
    conn.close()
    return df


def load_phonon_modes(phonon_data_id):
    conn = sqlite3.connect(DATABASE_PATH)
    df = pd.read_sql_query(
        "SELECT * FROM phonon_modes WHERE phonon_data_id = ?",
        conn, params=(phonon_data_id,)
    )
    conn.close()
    return df


def show_3d_structure(cif_path, height=380):
    """Display 3D structure from a CIF file path using py3Dmol."""
    cif_content = read_text_resource(cif_path)
    if not cif_content:
        st.warning(f"CIF file not found: {cif_path}")
        return
    try:
        import py3Dmol
        view = py3Dmol.view(width=700, height=height)
        view.addModel(cif_content, "cif")
        view.setStyle({"sphere": {"colorscheme": "Jmol", "scale": 0.28},
                       "stick": {"colorscheme": "Jmol", "radius": 0.12}})
        view.addUnitCell()
        view.zoomTo()
        raw_html = view._make_html()
        uid = f"v{random.randint(0, 999999)}"
        html = f"""
        <div id="{uid}_wrap" style="position:relative;border:1px solid #ddd;border-radius:4px">
          {raw_html}
          <div id="{uid}_ov" style="position:absolute;top:0;left:0;right:0;bottom:0;
               z-index:9999;cursor:pointer;display:flex;align-items:flex-end;
               justify-content:center;padding-bottom:6px">
            <span style="background:rgba(0,0,0,.5);color:#fff;padding:2px 10px;
                         border-radius:3px;font-size:11px">🖱 Click to interact</span>
          </div>
        </div>
        <script>
          (function(){{
            var ov=document.getElementById('{uid}_ov');
            var wrap=document.getElementById('{uid}_wrap');
            ov.addEventListener('click',function(){{ov.style.display='none';wrap.focus();}});
            document.addEventListener('click',function(e){{
              if(wrap&&!wrap.contains(e.target)){{ov.style.display='flex';}}
            }});
          }})();
        </script>
        """
        st.components.v1.html(html, height=height + 30, scrolling=False)
    except Exception as e:
        st.error(f"Error displaying 3D structure: {e}")


def show_phonon_mode(traj_path, height=380):
    """Display phonon mode animation from a .traj file or blob bundle using py3Dmol."""
    frames_xyz = load_traj_frames_xyz(traj_path)
    if not frames_xyz:
        st.warning(f"Trajectory not found: {traj_path}")
        return
    try:
        import py3Dmol
        view = py3Dmol.view(width=700, height=height)
        view.addModelsAsFrames(frames_xyz, "xyz")
        view.setStyle({"sphere": {"colorscheme": "Jmol", "scale": 0.28},
                       "stick": {"colorscheme": "Jmol", "radius": 0.12}})
        view.zoomTo()
        view.animate({'loop': 'backAndForth', 'interval': 100})
        raw_html = view._make_html()
        uid = f"v{random.randint(0, 999999)}"
        html = f"""
        <div id="{uid}_wrap" style="position:relative;border:1px solid #ddd;border-radius:4px">
          {raw_html}
          <div id="{uid}_ov" style="position:absolute;top:0;left:0;right:0;bottom:0;
               z-index:9999;cursor:pointer;display:flex;align-items:flex-end;
               justify-content:center;padding-bottom:6px">
            <span style="background:rgba(0,0,0,.5);color:#fff;padding:2px 10px;
                         border-radius:3px;font-size:11px">🖱 Click to interact</span>
          </div>
        </div>
        <script>
          (function(){{
            var ov=document.getElementById('{uid}_ov');
            var wrap=document.getElementById('{uid}_wrap');
            ov.addEventListener('click',function(){{ov.style.display='none';wrap.focus();}});
            document.addEventListener('click',function(e){{
              if(wrap&&!wrap.contains(e.target)){{ov.style.display='flex';}}
            }});
          }})();
        </script>
        """
        st.components.v1.html(html, height=height + 30, scrolling=False)
    except Exception as e:
        st.error(f"Error displaying phonon mode: {e}")


def parse_kpoints_text(text):
    """Parse special_k_points_*.txt text into (distances, labels)."""
    if not text:
        return None, None
    distances, labels = [], []
    in_dist = in_labels = False
    for line in text.splitlines():
        line = line.strip()
        if line == "Special K-point Distances:":
            in_dist, in_labels = True, False
        elif line == "Special K-point Labels:":
            in_dist, in_labels = False, True
        elif in_dist and line:
            try:
                distances.append(float(line))
            except ValueError:
                pass
        elif in_labels and line:
            labels.append("Γ" if line in ("G", "Gamma") else line)
    return (distances, labels) if distances else (None, None)


def plot_band_structure(bs_path, title="Phonon Band Structure"):
    try:
        data = read_array_resource(bs_path)
        if data is None or data.size == 0:
            return None
        n_kpts = len(data) if data.ndim == 1 else data.shape[0]
        fig, ax = plt.subplots(figsize=(7, 3.5))
        if data.ndim == 1:
            ax.plot(data, lw=0.8, color="steelblue")
        else:
            for col in range(data.shape[1]):
                ax.plot(data[:, col], lw=0.7, color="steelblue", alpha=0.75)
        ax.axhline(0, color="red", lw=1.2, ls="--", alpha=0.8)
        kpts_text = find_kpoints_text(bs_path)
        dists, kpt_labels = parse_kpoints_text(kpts_text)
        if dists and kpt_labels and len(dists) == len(kpt_labels):
            max_d = dists[-1]
            if max_d > 0:
                positions = [d / max_d * (n_kpts - 1) for d in dists]
                for pos in positions:
                    ax.axvline(pos, color="gray", lw=0.8, alpha=0.5)
                ax.set_xticks(positions)
                ax.set_xticklabels(kpt_labels, fontsize=9)
                ax.set_xlabel("k-path")
        else:
            ax.set_xlabel("k-point index")
        ax.set_ylabel("Frequency (THz)")
        ax.set_title(title)
        ax.grid(True, alpha=0.2)
        plt.tight_layout()
        return fig
    except Exception:
        return None


def load_initial_settings(output_dir_path):
    """Load initial_settings.json from output directory or blob bundle."""
    if not output_dir_path:
        return None
    text = read_text_resource(f"{output_dir_path}/initial_settings.json")
    if not text:
        return None
    try:
        return json.loads(text)
    except Exception:
        return None


STABILITY_ICONS = {
    "Stable": "🟢 Stable",
    "Quasi-Stable": "🟡 Quasi-Stable",
    "Unstable": "🔴 Unstable",
    "Unknown": "⚪ Unknown",
}


def stability_badge(s):
    return STABILITY_ICONS.get(s, f"⚪ {s}")


def fmt(v, decimals=4):
    try:
        return f"{float(v):.{decimals}f}"
    except Exception:
        return "—"


def main():
    st.set_page_config(page_title="VibroML Explorer", page_icon="🔬", layout="wide")
    st.title("🔬 VibroML GA Phonon Results Explorer")

    if not DATABASE_PATH.exists():
        st.error(f"Database not found: {DATABASE_PATH}")
        return

    if "view" not in st.session_state:
        st.session_state["view"] = "table"

    df = load_data()

    # ── Sidebar filters ──────────────────────────────────────────────
    st.sidebar.header("Filters")
    st.sidebar.markdown(f"**{len(df)} total materials**")

    dataset_opts = sorted(df["dataset"].dropna().unique().tolist())
    dataset = st.sidebar.multiselect("Dataset", options=dataset_opts, default=dataset_opts)

    status_opts = sorted(df["run_status"].dropna().unique().tolist())
    run_status = st.sidebar.multiselect("Run Status", options=status_opts, default=status_opts)

    init_opts = sorted(df["initial_stability"].dropna().unique().tolist())
    if "Unknown" not in init_opts:
        init_opts = ["Unknown"] + init_opts
    initial_stability = st.sidebar.multiselect("Initial Stability", options=init_opts, default=init_opts)

    final_opts = sorted(df["final_stability"].dropna().unique().tolist())
    if "Unknown" not in final_opts:
        final_opts = ["Unknown"] + final_opts
    final_stability = st.sidebar.multiselect("Final Stability", options=final_opts, default=final_opts)

    sg_min, sg_max = int(df["space_group"].dropna().min()), int(df["space_group"].dropna().max())
    sg_range = st.sidebar.slider("Space Group", sg_min, sg_max, (sg_min, sg_max))

    e_min = float(df["e_above_hull_alexandria_eV_per_atom"].dropna().min())
    e_max = float(df["e_above_hull_alexandria_eV_per_atom"].dropna().max())
    e_range = st.sidebar.slider("E_hull (eV/atom)", e_min, e_max, (e_min, e_max), step=0.001)

    search = st.sidebar.text_input("Search by name / AGM ID")

    # ── Apply filters ────────────────────────────────────────────────
    mask_init = df["initial_stability"].isin(initial_stability) | (
        ("Unknown" in initial_stability) & df["initial_stability"].isna()
    )
    mask_final = df["final_stability"].isin(final_stability) | (
        ("Unknown" in final_stability) & df["final_stability"].isna()
    )
    mask_sg = (
        ((df["space_group"] >= sg_range[0]) & (df["space_group"] <= sg_range[1]))
        | (df["space_group"].isna() if (sg_range[0] == sg_min and sg_range[1] == sg_max) else False)
    )
    mask_ehull = (
        ((df["e_above_hull_alexandria_eV_per_atom"] >= e_range[0]) & (df["e_above_hull_alexandria_eV_per_atom"] <= e_range[1]))
        | (df["e_above_hull_alexandria_eV_per_atom"].isna() if (e_range[0] == e_min and e_range[1] == e_max) else False)
    )

    df_filtered = df[
        df["dataset"].isin(dataset) &
        df["run_status"].isin(run_status) &
        mask_init & mask_final &
        mask_sg & mask_ehull
    ]

    if search:
        m = (
            df_filtered["material_name"].str.contains(search, case=False, na=False) |
            df_filtered["agm_id"].str.contains(search, case=False, na=False)
        )
        df_filtered = df_filtered[m]

    st.sidebar.markdown(f"**{len(df_filtered)} shown**")

    # ── View: Table ──────────────────────────────────────────────────
    if st.session_state["view"] == "table":
        st.subheader(f"Filtered: {len(df_filtered)} materials")
        st.caption("Use the selector below the table to jump to a material's detail view.")

        df_filtered = df_filtered.copy()
        df_filtered["predicted_final_ehull_eV_per_atom"] = (
            df_filtered["e_above_hull_alexandria_eV_per_atom"] + df_filtered["energy_reduction_vs_initial_eV_per_atom"]
        )
        df_filtered["energy_reduction_meV_per_atom"] = df_filtered["energy_reduction_vs_initial_eV_per_atom"] * 1000

        display_cols = [
            "cif_stem", "material_name", "dataset", "space_group",
            "e_above_hull_alexandria_eV_per_atom",
            "initial_stability", "final_stability", "run_status",
            "initial_phonon_min_freq_THz", "final_phonon_min_freq_THz",
            "energy_reduction_meV_per_atom", "predicted_final_ehull_eV_per_atom",
            "ga_n_main_iterations", "ga_total_relaxations", "ga_unique_structures",
            "runtime_hours",
        ]
        col_rename = {
            "material_name": "Name", "dataset": "Dataset", "space_group": "SG",
            "e_above_hull_alexandria_eV_per_atom": "E_hull",
            "initial_stability": "Init. Stability", "final_stability": "Final Stability",
            "run_status": "Status", "ga_n_main_iterations": "GA Iters",
            "ga_total_relaxations": "Relaxations", "ga_unique_structures": "Unique Structs",
            "runtime_hours": "Runtime (h)",
            "initial_phonon_min_freq_THz": "Init MinFreq (THz)",
            "final_phonon_min_freq_THz": "Final MinFreq (THz)",
            "energy_reduction_meV_per_atom": "Energy Reduction (meV/atom)",
            "predicted_final_ehull_eV_per_atom": "Predicted Final Ehull (eV/atom)",
        }
        df_show = df_filtered[display_cols].rename(columns=col_rename).copy()
        for col in ["E_hull", "Runtime (h)", "Init MinFreq (THz)", "Final MinFreq (THz)",
                   "Energy Reduction (meV/atom)", "Predicted Final Ehull (eV/atom)"]:
            if col in df_show.columns:
                df_show[col] = pd.to_numeric(df_show[col], errors="coerce").round(4)

        df_show = df_show.drop(columns=["cif_stem"])

        st.dataframe(df_show, width="stretch", height=500, hide_index=True)

        st.caption("Jump to a material's detail view:")
        jc1, jc2 = st.columns([4, 1])
        with jc1:
            material_options = {
                row["cif_stem"]: f"🔍 {row['material_name']} ({row['dataset']})"
                for _, row in df_filtered.iterrows()
            }
            selected_stem = st.selectbox(
                "Material",
                options=list(material_options.keys()),
                format_func=lambda x: material_options[x],
                key="material_quickjump",
                label_visibility="collapsed",
            )
        with jc2:
            if st.button("View Details →", key="quickjump_btn", use_container_width=True):
                st.session_state["selected_material"] = selected_stem
                st.session_state["view"] = "detail"
                st.rerun()

        # Statistics section
        st.divider()
        with st.expander("📈 Statistics for Filtered Materials", expanded=True):
            n_total = len(df_filtered)
            n_stable_no_soft = len(df_filtered[df_filtered["initial_stability"].str.lower() == "stable"])
            n_final_stable = len(df_filtered[df_filtered["final_stability"].str.lower() == "stable"])
            denom_vibroml = max(n_total - n_stable_no_soft, 0)

            col1, col2, col3, col4, col5 = st.columns(5)

            with col1:
                avg_ehull = df_filtered["e_above_hull_alexandria_eV_per_atom"].mean()
                st.metric("Average E_hull", f"{avg_ehull:.4f} eV/atom" if pd.notna(avg_ehull) else "—")

            with col2:
                pct_no_soft = (n_stable_no_soft / n_total * 100) if n_total else 0.0
                st.metric("Stable (no soft modes)", f"{n_stable_no_soft} / {n_total}")
                st.markdown(f"<span style='color:gray;font-size:0.8em'>{pct_no_soft:.1f}%</span>", unsafe_allow_html=True)

            with col3:
                pct_vibroml = (n_final_stable / denom_vibroml * 100) if denom_vibroml else 0.0
                st.metric("VibroML stabilized", f"{n_final_stable} / {denom_vibroml}")
                st.markdown(f"<span style='color:gray;font-size:0.8em'>{pct_vibroml:.1f}%</span>", unsafe_allow_html=True)

            with col4:
                avg_energy_reduction = df_filtered["energy_reduction_meV_per_atom"].mean()
                st.metric("Avg Energy Reduction", f"{avg_energy_reduction:.2f} meV/atom" if pd.notna(avg_energy_reduction) else "—")

            with col5:
                avg_runtime = df_filtered["runtime_hours"].mean()
                st.metric("Avg Runtime", f"{avg_runtime:.2f} h" if pd.notna(avg_runtime) else "—")

            fig_ehull, ax_ehull = plt.subplots(figsize=(10, 3))
            ehull_vals = df_filtered["e_above_hull_alexandria_eV_per_atom"].dropna()
            if not ehull_vals.empty:
                ax_ehull.hist(ehull_vals, bins=50, edgecolor="black", alpha=0.7, color="steelblue")
                ax_ehull.set_xlabel("E_hull (eV/atom)")
                ax_ehull.set_ylabel("Count")
                ax_ehull.set_title("E_hull Distribution (Initial)")
                xlo = min(0, float(ehull_vals.min()))
                xhi = max(0, float(ehull_vals.max()))
                ax_ehull.set_xlim(xlo, xhi if xhi > xlo else xlo + 1)
                ax_ehull.axvline(0, color="black", linewidth=1, linestyle="--")
                plt.tight_layout()
                st.pyplot(fig_ehull)
            plt.close(fig_ehull)

            fig_cmp, ax_cmp = plt.subplots(figsize=(10, 3))
            initial_vals = df_filtered["e_above_hull_alexandria_eV_per_atom"].dropna()
            predicted_vals = df_filtered["predicted_final_ehull_eV_per_atom"].dropna()
            all_vals = pd.concat([initial_vals, predicted_vals])
            if not all_vals.empty:
                xlo = min(0, float(all_vals.min()))
                xhi = max(0, float(all_vals.max()))
                bins = np.linspace(xlo, xhi if xhi > xlo else xlo + 1, 51)
                ax_cmp.hist(initial_vals, bins=bins, edgecolor="black", alpha=0.5, color="red", label="Initial E_hull")
                ax_cmp.hist(predicted_vals, bins=bins, edgecolor="black", alpha=0.5, color="steelblue", label="Predicted Final E_hull")
                ax_cmp.set_xlabel("E_hull (eV/atom)")
                ax_cmp.set_ylabel("Count")
                ax_cmp.set_title("Initial vs Predicted Final E_hull Distribution")
                ax_cmp.set_xlim(xlo, xhi if xhi > xlo else xlo + 1)
                ax_cmp.axvline(0, color="black", linewidth=1, linestyle="--")
                ax_cmp.legend()
                plt.tight_layout()
                st.pyplot(fig_cmp)
            plt.close(fig_cmp)

            fig_init_stab, ax_init_stab = plt.subplots(figsize=(10, 3))
            init_stab_counts = df_filtered["initial_stability"].value_counts(dropna=False)
            init_stab_counts = init_stab_counts[init_stab_counts.index.notna()]
            init_stab_counts.index = init_stab_counts.index.astype(str)
            ax_init_stab.bar(init_stab_counts.index, init_stab_counts.values, color="steelblue", alpha=0.7)
            ax_init_stab.set_xlabel("Initial Stability")
            ax_init_stab.set_ylabel("Count")
            ax_init_stab.set_title("Initial Stability Distribution")
            plt.xticks(rotation=45)
            plt.tight_layout()
            st.pyplot(fig_init_stab)
            plt.close(fig_init_stab)

            fig_omega, ax_omega = plt.subplots(figsize=(10, 3))
            ax_omega.hist(df_filtered["initial_phonon_min_freq_THz"].dropna(),
                         bins=50, edgecolor="black", alpha=0.7, color="steelblue")
            ax_omega.set_xlabel("Initial MinFreq (THz)")
            ax_omega.set_ylabel("Count")
            ax_omega.set_title("Initial Omega Min Distribution")
            plt.tight_layout()
            st.pyplot(fig_omega)
            plt.close(fig_omega)

            fig_soft, ax_soft = plt.subplots(figsize=(10, 3))
            ax_soft.hist(df_filtered["initial_softness_fraction"].dropna(),
                        bins=50, edgecolor="black", alpha=0.7, color="steelblue")
            ax_soft.set_xlabel("Softness Fraction")
            ax_soft.set_ylabel("Count")
            ax_soft.set_title("Softness Fraction Distribution")
            plt.tight_layout()
            st.pyplot(fig_soft)
            plt.close(fig_soft)

            fig_final_stab, ax_final_stab = plt.subplots(figsize=(10, 3))
            final_stab_counts = df_filtered["final_stability"].value_counts(dropna=False)
            final_stab_counts = final_stab_counts[final_stab_counts.index.notna()]
            final_stab_counts.index = final_stab_counts.index.astype(str)
            ax_final_stab.bar(final_stab_counts.index, final_stab_counts.values, color="steelblue", alpha=0.7)
            ax_final_stab.set_xlabel("Final Stability")
            ax_final_stab.set_ylabel("Count")
            ax_final_stab.set_title("Final Stability Distribution")
            plt.xticks(rotation=45)
            plt.tight_layout()
            st.pyplot(fig_final_stab)
            plt.close(fig_final_stab)

            fig_runtime, ax_runtime = plt.subplots(figsize=(10, 3))
            ax_runtime.hist(df_filtered["runtime_hours"].dropna(),
                          bins=50, edgecolor="black", alpha=0.7, color="steelblue")
            ax_runtime.set_xlabel("Runtime (hours)")
            ax_runtime.set_ylabel("Count")
            ax_runtime.set_title("Runtime Distribution")
            plt.tight_layout()
            st.pyplot(fig_runtime)
            plt.close(fig_runtime)

    # ── View: Detail ─────────────────────────────────────────────────
    elif st.session_state["view"] == "detail":
        if st.button("⬅ Back to table", key="back_to_table"):
            st.session_state["view"] = "table"
            st.rerun()

        if df_filtered.empty:
            st.info("No materials match current filters.")
            return

        stems = df_filtered["cif_stem"].tolist()
        default_idx = 0
        if "selected_material" in st.session_state and st.session_state["selected_material"] in stems:
            default_idx = stems.index(st.session_state["selected_material"])
        selected = st.selectbox("Select material", options=stems, index=default_idx)
        if not selected:
            return
        st.session_state["selected_material"] = selected

        with st.spinner("Loading material details..."):
            m = df_filtered[df_filtered["cif_stem"] == selected].iloc[0]
            mid = int(m["id"])

            c1, c2, c3, c4, c5, c6 = st.columns(6)
            c1.metric("Material", m["material_name"])
            c2.metric("AGM ID", m["agm_id"] or "—")
            c3.metric("Space Group", int(m["space_group"]) if pd.notna(m["space_group"]) else "—")
            c4.metric("E_hull (eV/atom)", fmt(m["e_above_hull_alexandria_eV_per_atom"]))
            c5.metric("Atoms", int(m["n_atoms"]) if pd.notna(m["n_atoms"]) else "—")
            c6.metric("Dataset", m["dataset"])

            st.divider()
            c1, c2, c3, c4, c5 = st.columns(5)
            c1.metric("Run Status", m["run_status"])
            c2.metric("GA Iterations", int(m["ga_n_main_iterations"]) if pd.notna(m["ga_n_main_iterations"]) else "—")
            c3.metric("Relaxations", int(m["ga_total_relaxations"]) if pd.notna(m["ga_total_relaxations"]) else "—")
            c4.metric("Unique Structs", int(m["ga_unique_structures"]) if pd.notna(m["ga_unique_structures"]) else "—")
            c5.metric("Runtime (h)", fmt(m["runtime_hours"], 2))

            st.divider()
            c1, c2 = st.columns(2)
            with c1:
                st.subheader("Initial Phonon")
                st.markdown(f"**{stability_badge(m['initial_stability'])}**")
                c1a, c1b = st.columns(2)
                c1a.metric("Min Freq (THz)", fmt(m["initial_phonon_min_freq_THz"]))
                c1b.metric("Softness", f"{float(m['initial_softness_fraction'])*100:.2f}%" if pd.notna(m["initial_softness_fraction"]) else "—")
            with c2:
                st.subheader("Final / Best Phonon")
                st.markdown(f"**{stability_badge(m['final_stability'])}**")
                c2a, c2b = st.columns(2)
                c2a.metric("Min Freq (THz)", fmt(m["final_phonon_min_freq_THz"]))
                c2b.metric("Softness", f"{float(m['final_softness_fraction'])*100:.2f}%" if pd.notna(m["final_softness_fraction"]) else "—")
                if m["last_phonon_source"] and m["last_phonon_source"] != "none":
                    st.caption(f"Source: `{m['last_phonon_source']}`")

            settings = load_initial_settings(m["output_dir_path"])
            if settings:
                with st.expander("⚙️ Computation Settings (initial_settings.json)", expanded=False):
                    keys = ["engine", "model_name", "supercell_n", "delta", "fmax", "no_relax", "units"]
                    cs = st.columns(len([k for k in keys if k in settings]))
                    ci = 0
                    for k in keys:
                        if k in settings:
                            cs[ci].metric(k, str(settings[k]))
                            ci += 1

            st.divider()
            st.subheader("3D Structure Viewer")
            struct_labels = ["📄 Original (Input CIF)"]
            if m["cif_optimized_path"] and file_available(m["cif_optimized_path"], "text_files"):
                struct_labels.append("✨ Best GA Structure")

            struct_tabs = st.tabs(struct_labels)
            with struct_tabs[0]:
                show_3d_structure(m["cif_path"])
                if m["cif_path"]:
                    st.caption(f"`{m['cif_path']}`")

            if len(struct_tabs) > 1:
                with struct_tabs[1]:
                    show_3d_structure(m["cif_optimized_path"])
                    st.caption(f"`{m['cif_optimized_path']}`")

            phonon_df = load_phonon_data(mid)
            ga_df = load_ga_structures(mid)

            st.divider()
            if not ga_df.empty:
                st.subheader(f"GA Structures — {len(ga_df)} unique")
                ga_show = ga_df[["struct_type", "rank", "iteration", "energy_per_atom_eV",
                                  "energy_reduction_vs_initial_eV_per_atom"]].copy()
                ga_show["Label"] = ga_show.apply(
                    lambda r: f"Top #{int(r['rank'])}" if r["struct_type"] == "top" else f"Unique #{int(r['rank'])}", axis=1
                )
                ga_show = ga_show[["Label", "struct_type", "iteration", "energy_per_atom_eV",
                                   "energy_reduction_vs_initial_eV_per_atom"]].copy()
                ga_show.columns = ["Label", "Type", "GA Iter", "Energy (eV/atom)", "ΔE vs Initial (eV/atom)"]
                for col in ["Energy (eV/atom)", "ΔE vs Initial (eV/atom)"]:
                    ga_show[col] = pd.to_numeric(ga_show[col], errors="coerce").round(4)
                st.dataframe(ga_show, width="stretch", height=280)

                def _rank_label(idx):
                    row_ = ga_df.loc[idx]
                    prefix = "Top" if row_.get("struct_type", "top") == "top" else "Unique"
                    return f"{prefix} #{int(row_['rank'])}  –  {fmt(row_['energy_per_atom_eV'])} eV/atom"

                selected_ga_idx = st.selectbox(
                    "View 3D structure",
                    options=ga_df.index.tolist(),
                    format_func=_rank_label
                )
                if selected_ga_idx is not None:
                    row_info = ga_df.loc[selected_ga_idx]
                    ga_cif_path = row_info["cif_path"]
                    selected_rank = int(row_info["rank"])
                    show_3d_structure(ga_cif_path)
                    st.caption(f"`{ga_cif_path}`")
                    if not phonon_df.empty:
                        stype = row_info.get("struct_type", "top")
                        phonon_type = "final" if stype == "top" else "unique"
                        rank_phonon = phonon_df[
                            (phonon_df["data_type"] == phonon_type) &
                            (phonon_df["iteration"] == selected_rank)
                        ]
                        type_label = f"Top #{selected_rank}" if phonon_type == "final" else f"Unique #{selected_rank}"
                        if not rank_phonon.empty:
                            row_p = rank_phonon.iloc[0]
                            st.markdown(f"**Phonon ({type_label}):** {stability_badge(row_p['stability'])}  |  {fmt(row_p['min_freq_THz'])} THz")
                            bs_path = row_p["band_structure_path"]
                            if bs_path and file_available(bs_path, "array_files"):
                                fig = plot_band_structure(bs_path, title=f"{type_label} \u2013 Phonon Band Structure")
                                if fig:
                                    st.pyplot(fig)
                                    plt.close(fig)
                        else:
                            st.info(f"No phonon analysis stored for {type_label}.")
            else:
                st.info("No GA structures in database for this material.")

            st.divider()
            st.subheader("Phonon Analyses")

            if not phonon_df.empty:
                for _, row in phonon_df.iterrows():
                    dt = row["data_type"]
                    it = int(row["iteration"]) if pd.notna(row["iteration"]) else None
                    if dt == "initial":
                        label = "Initial Phonon"
                    elif dt == "guidance":
                        label = f"Guidance – Iter {it}"
                    elif dt == "final":
                        label = f"Top #{it}"
                    elif dt == "unique":
                        label = f"Unique #{it}"
                    else:
                        label = f"{dt.capitalize()}" + (f" #{it}" if it else "")
                    min_f = fmt(row["min_freq_THz"]) if pd.notna(row["min_freq_THz"]) else "—"
                    with st.expander(f"{stability_badge(row['stability'])}  |  {label}  |  {min_f} THz",
                                     expanded=(row["data_type"] == "initial")):
                        ca, cb, cc = st.columns(3)
                        ca.metric("Min Freq (THz)", fmt(row["min_freq_THz"]))
                        cb.metric("Softness", f"{float(row['softness_fraction'])*100:.2f}%" if pd.notna(row["softness_fraction"]) else "—")
                        cc.metric("Stability", row["stability"])

                        bs_path = row["band_structure_path"]
                        if bs_path and file_available(bs_path, "array_files"):
                            fig = plot_band_structure(bs_path, title=f"{label} — Band Structure")
                            if fig:
                                st.pyplot(fig)
                                plt.close(fig)
                            st.caption(f"📄 `{bs_path}`")

                        dos_path = row["dos_path"]
                        if dos_path and file_available(dos_path, "array_files"):
                            st.caption(f"📄 DOS: `{dos_path}`")

                        modes_df = load_phonon_modes(row["id"])
                        if not modes_df.empty:
                            st.divider()
                            st.markdown("**Phonon Mode Visualization**")

                            softest_modes = modes_df[modes_df["mode_type"] == "softest"]
                            optical_modes = modes_df[modes_df["mode_type"] == "optical"]

                            mode_options = []
                            if not softest_modes.empty:
                                mode_options.extend([("softest", idx) for idx in softest_modes.index])
                            if not optical_modes.empty:
                                mode_options.extend([("optical", idx) for idx in optical_modes.index])

                            if mode_options:
                                def mode_label(idx):
                                    m_ = modes_df.loc[idx]
                                    name = m_["mode_name"]
                                    if m_["mode_type"] == "softest":
                                        return f"🔴 Softest: {name}"
                                    else:
                                        return f"🟡 Optical: {name}"

                                selected_mode_idx = st.selectbox(
                                    "Select phonon mode to animate",
                                    options=[idx for _, idx in mode_options],
                                    format_func=mode_label,
                                    key=f"mode_{row['id']}"
                                )

                                if selected_mode_idx is not None:
                                    selected_mode = modes_df.loc[selected_mode_idx]
                                    traj_path = selected_mode["traj_path"]
                                    show_phonon_mode(traj_path)
                                    st.caption(f"📄 Trajectory: `{traj_path}`")
                            else:
                                st.info("No phonon mode trajectories found.")
            else:
                st.info("No phonon data in database. Re-run ingestion to populate.")

            if m["output_dir_path"]:
                st.divider()
                st.caption(f"📁 Output directory: `{m['output_dir_path']}`")


if __name__ == "__main__":
    main()
