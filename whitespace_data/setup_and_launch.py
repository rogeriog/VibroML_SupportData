#!/usr/bin/env python3
"""
Seamless Download, Setup, and Launch Script for the VibroML Database Explorer.

Usage:
  1. Clone the repository:
       git clone https://github.com/rogeriog/VibroML_SupportData.git
       cd VibroML_SupportData

  2. Download the full phonon/structure bundle from Figshare and launch the app:
       python launch_explorer.py

  Optional flags:
    --url <FIGSHARE_DOWNLOAD_URL>   Override the default Figshare download URL
    --bundle <PATH_TO_TAR_GZ>       Use a locally downloaded vibroml_explorer_bundle.tar.gz
    --no-launch                     Download and set up the database without starting Streamlit
    --reconstruct-files             Materialize physical .cif, .txt, and .traj files on disk
                                    from vibroml_blobs.db (optional; the explorer reads
                                    directly from vibroml_blobs.db without needing extraction)
    --material <NAME_OR_STEM>       When used with --reconstruct-files, only materialize files
                                    for materials matching this substring
"""

import argparse
import gzip
import os
import shutil
import sqlite3
import subprocess
import sys
import tarfile
import urllib.request
import zlib
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
DB_DIR = HERE / "database"
DB_PATH = DB_DIR / "vibroml_database.db"
DB_GZ_IN_REPO = HERE / "vibroml_database.db.gz"
BLOBS_DB_PATH = DB_DIR / "vibroml_blobs.db"
EXPLORER_SCRIPT = HERE / "vibroml_explorer.py"

# Update this URL/DOI after publishing the Figshare item:
FIGSHARE_BUNDLE_URL = os.environ.get(
    "VIBROML_FIGSHARE_BUNDLE_URL",
    "https://ndownloader.figshare.com/files/YOUR_FIGSHARE_FILE_ID",
)


def download_with_progress(url: str, dest: Path):
    print(f"Downloading phonon & structure bundle from:\n  {url}")
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (VibroML-Setup)"})
    with urllib.request.urlopen(req) as resp, open(dest, "wb") as out:
        total = int(resp.headers.get("Content-Length", 0))
        downloaded = 0
        chunk_size = 1024 * 1024  # 1 MB
        while True:
            chunk = resp.read(chunk_size)
            if not chunk:
                break
            out.write(chunk)
            downloaded += len(chunk)
            if total > 0:
                pct = downloaded * 100.0 / total
                mb_done = downloaded / (1024 ** 2)
                mb_tot = total / (1024 ** 2)
                print(f"\r  {mb_done:7.1f} / {mb_tot:.1f} MB ({pct:5.1f}%)", end="", flush=True)
            else:
                mb_done = downloaded / (1024 ** 2)
                print(f"\r  {mb_done:7.1f} MB downloaded", end="", flush=True)
    print("\nDownload complete.")


def ensure_metadata_db():
    DB_DIR.mkdir(parents=True, exist_ok=True)
    if DB_PATH.exists():
        return
    if DB_GZ_IN_REPO.exists():
        print(f"Decompressing {DB_GZ_IN_REPO.name} -> {DB_PATH} ...")
        with gzip.open(DB_GZ_IN_REPO, "rb") as fin, open(DB_PATH, "wb") as fout:
            shutil.copyfileobj(fin, fout)
        print("Metadata database ready.")


def ensure_blobs_db(bundle_path: Path = None, url: str = None):
    DB_DIR.mkdir(parents=True, exist_ok=True)
    if BLOBS_DB_PATH.exists():
        print(f"Found phonon/structure blob database: {BLOBS_DB_PATH}")
        return

    # Check if a local archive exists in HERE or parent
    if bundle_path is None:
        for candidate in [
            HERE / "vibroml_explorer_bundle.tar.gz",
            HERE.parent / "vibroml_explorer_bundle.tar.gz",
            HERE / "vibroml_blobs.db.gz",
            HERE.parent / "vibroml_blobs.db.gz",
        ]:
            if candidate.exists():
                bundle_path = candidate
                break

    if bundle_path is None:
        target_url = url or FIGSHARE_BUNDLE_URL
        if "YOUR_FIGSHARE_FILE_ID" in target_url:
            print(
                "\n[NOTE] Full phonon/trajectory bundle (vibroml_blobs.db) is not yet in this folder.\n"
                "       Download `vibroml_explorer_bundle.tar.gz` from Figshare and place it in\n"
                f"       `{HERE}/` (or pass `--url <FIGSHARE_URL>` / `--bundle <PATH>`).\n"
                "       Launching with metadata database (summary tables & statistics)...\n"
            )
            return
        bundle_path = HERE / "vibroml_explorer_bundle.tar.gz"
        download_with_progress(target_url, bundle_path)

    print(f"Extracting {bundle_path.name} into {HERE} ...")
    if bundle_path.name.endswith(".tar.gz") or bundle_path.name.endswith(".tgz"):
        with tarfile.open(bundle_path, "r:gz") as tar:
            tar.extractall(path=HERE)
    elif bundle_path.name.endswith(".db.gz"):
        with gzip.open(bundle_path, "rb") as fin, open(BLOBS_DB_PATH, "wb") as fout:
            shutil.copyfileobj(fin, fout)
    print("Extraction complete.")


def reconstruct_files_to_disk(material_filter: str = None):
    """Optionally materialize physical .cif, .txt, and .traj files on disk from vibroml_blobs.db."""
    if not BLOBS_DB_PATH.exists():
        print(f"ERROR: {BLOBS_DB_PATH} not found. Cannot reconstruct files.")
        sys.exit(1)

    from ase import Atoms
    from ase.io.trajectory import TrajectoryWriter

    conn = sqlite3.connect(BLOBS_DB_PATH)
    c = conn.cursor()

    like_clause = f"%{material_filter}%" if material_filter else "%"

    print("Reconstructing text files (.cif, initial_settings.json, special_k_points_*.txt)...")
    c.execute("SELECT rel_path, data_z FROM text_files WHERE rel_path LIKE ?", (like_clause,))
    n_txt = 0
    for rel_path, data_z in c:
        out_p = HERE / rel_path
        if not out_p.exists():
            out_p.parent.mkdir(parents=True, exist_ok=True)
            out_p.write_bytes(zlib.decompress(data_z))
            n_txt += 1
    print(f"  Wrote {n_txt:,} text files.")

    print("Reconstructing array files (band_structure_energies_*.txt, dos_energies_*.txt)...")
    c.execute("SELECT rel_path, shape_0, shape_1, data_z FROM array_files WHERE rel_path LIKE ?", (like_clause,))
    n_arr = 0
    for rel_path, s0, s1, data_z in c:
        out_p = HERE / rel_path
        if not out_p.exists():
            out_p.parent.mkdir(parents=True, exist_ok=True)
            arr = np.frombuffer(zlib.decompress(data_z), dtype=np.float32)
            arr = arr.reshape(s0, s1) if s1 > 1 else arr.reshape(s0)
            np.savetxt(out_p, arr, fmt="%.6e")
            n_arr += 1
    print(f"  Wrote {n_arr:,} array files.")

    print("Reconstructing ASE .traj phonon mode animations...")
    c.execute(
        "SELECT rel_path, n_frames, n_atoms, numbers, cell, r0, uc, us FROM traj_files WHERE rel_path LIKE ?",
        (like_clause,),
    )
    n_tr = 0
    for rel_path, n_frames, n_atoms, num_z, cell_z, r0_z, uc_z, us_z in c:
        out_p = HERE / rel_path
        if not out_p.exists():
            out_p.parent.mkdir(parents=True, exist_ok=True)
            numbers = np.frombuffer(zlib.decompress(num_z), dtype=np.uint8)
            cell = np.frombuffer(zlib.decompress(cell_z), dtype=np.float32).reshape(3, 3)
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

            writer = TrajectoryWriter(str(out_p), mode="w")
            for i in range(n_frames):
                atoms = Atoms(numbers=numbers, positions=all_pos[i], cell=cell, pbc=True)
                writer.write(atoms)
            writer.close()
            n_tr += 1
    print(f"  Wrote {n_tr:,} .traj files.")
    conn.close()


def main():
    parser = argparse.ArgumentParser(description="Download, set up, and launch the VibroML Database Explorer.")
    parser.add_argument("--url", type=str, default=None, help="Figshare download URL for vibroml_explorer_bundle.tar.gz")
    parser.add_argument("--bundle", type=Path, default=None, help="Local path to vibroml_explorer_bundle.tar.gz")
    parser.add_argument("--no-launch", action="store_true", help="Set up database without launching Streamlit")
    parser.add_argument("--reconstruct-files", action="store_true", help="Materialize physical .cif/.txt/.traj files on disk")
    parser.add_argument("--material", type=str, default=None, help="Substring filter when using --reconstruct-files")
    parser.add_argument("--port", type=int, default=8501, help="Streamlit server port (default: 8501)")
    args = parser.parse_args()

    ensure_metadata_db()
    ensure_blobs_db(bundle_path=args.bundle, url=args.url)

    if args.reconstruct_files:
        reconstruct_files_to_disk(material_filter=args.material)

    if not args.no_launch:
        env = os.environ.copy()
        env["VIBROML_DATA_DIR"] = str(HERE)
        print(f"\nLaunching VibroML Explorer on port {args.port}...")
        subprocess.run(
            [sys.executable, "-m", "streamlit", "run", str(EXPLORER_SCRIPT), "--server.port", str(args.port)],
            env=env,
            check=True,
        )


if __name__ == "__main__":
    main()
