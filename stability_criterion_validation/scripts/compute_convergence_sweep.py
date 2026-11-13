#!/usr/bin/env python3
"""
Convergence sweep: phonon stability metrics vs q-mesh size.

For each of 4 representative high-symmetry aristotype structures (LiF, CsPbI3,
HfO2, Bi2Sn2O7) this script computes phonon frequencies on three uniform
q-meshes (4x4x4, 6x6x6, 8x8x8) using a single MLIP engine (MACE-OMAT by
default) and evaluates three stability metrics:

  * soft_frac  - fraction of modes with omega < 0 (3 Gamma acoustic modes
                 excluded by index)
  * I_curv     - mean negative curvature at p=1, i.e. mean(max(0, -omega^2))
  * omega_min  - minimum frequency (THz)

Results are written to a JSON file and a CSV file, and a 2-panel convergence
figure (PDF + PNG) is generated showing soft_frac and I_curv vs q-mesh size
for each structure.

Usage:
  python compute_convergence_sweep.py --engine mace --model-name medium-omat-0 \
      --device cuda --output convergence_sweep.json
"""
import argparse
import csv
import json
import os
import shutil
import sys
import time

import numpy as np

EV_TO_THZ = 4.135667696e-3

# q-mesh sizes to sweep
Q_MESH_SIZES = [4, 6, 8]

# soft_frac reference line drawn on the convergence plot (pooled optimal separation)
SOFT_FRAC_REFERENCE = 0.045

# Add VibroML to path for initialize_calculator
VIBROML_PATH = os.environ.get(
    "VIBROML_CODEBASE",
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..",
                 "01_codebase", "VibroML"))
if VIBROML_PATH not in sys.path:
    sys.path.insert(0, VIBROML_PATH)


# -----------------------------------------------------------------------------
# Calculator loading (delegates to VibroML's canonical initializer)
# -----------------------------------------------------------------------------
def load_calculator(engine, model_name=None, device="cuda"):
    """Load a calculator via VibroML's canonical initialize_calculator.

    This delegates to vibroml.utils.structure_utils.initialize_calculator
    which handles all engine-specific setup (checkpoint paths, model
    downloads, device placement).
    """
    from vibroml.utils.structure_utils import initialize_calculator

    defaults = {
        "mace": "medium-omat-0",
        "macemp": "medium",
        "esen": None,
        "uma": None,
        "pet": "pet-mad-xs",
    }

    actual_engine = engine
    if engine == "macemp":
        actual_engine = "mace"
        if model_name is None:
            model_name = "medium"
    elif engine == "mace" and model_name is None:
        model_name = "medium-omat-0"
    elif model_name is None:
        model_name = defaults.get(engine)

    return initialize_calculator(
        engine=actual_engine,
        model_name=model_name or "medium-omat-0",
    )


# -----------------------------------------------------------------------------
# Phonon computation (same pattern as compute_ip_variants.py)
# -----------------------------------------------------------------------------
def compute_phonon_freqs_uniform(atoms, calc, supercell=(2, 2, 2),
                                 delta=0.05, q_mesh=(6, 6, 6)):
    """Compute phonon frequencies on a uniform q-mesh.

    Uses ASE's finite-displacement Phonons with a temporary cache directory
    that is always cleaned up in the ``finally`` block.

    Parameters
    ----------
    atoms : ase.Atoms
        Structure (calculator will be attached internally).
    calc : Calculator
        MLIP calculator to use for the force evaluations.
    supercell : tuple of int
        Supercell for finite displacements.
    delta : float
        Displacement amplitude (Angstrom).
    q_mesh : tuple of int
        Monkhorst-Pack q-mesh dimensions.

    Returns
    -------
    freqs_THz : np.ndarray
        Phonon frequencies in THz, shape (n_qpoints, n_modes).
    gamma_indices : np.ndarray
        Indices of Gamma q-points within the mesh.
    """
    from ase.phonons import Phonons
    from ase.dft.kpoints import monkhorst_pack
    import tempfile

    atoms.calc = calc
    ph_dir = tempfile.mkdtemp(prefix="phonon_conv_", dir=".")

    try:
        ph = Phonons(atoms, calc, supercell=supercell, delta=delta, name=ph_dir)
        ph.run()
        ph.read()

        qpoints = monkhorst_pack(q_mesh)
        energies_eV = ph.band_structure(qpoints, modes=False, verbose=False)
        freqs_THz = energies_eV / EV_TO_THZ

        gamma_mask = np.all(np.abs(qpoints) < 1e-10, axis=1)
        gamma_indices = np.where(gamma_mask)[0]

        ph.clean()
    finally:
        shutil.rmtree(ph_dir, ignore_errors=True)

    return freqs_THz, gamma_indices


# -----------------------------------------------------------------------------
# Metric computation
# -----------------------------------------------------------------------------
def compute_metrics(freqs_THz, gamma_indices):
    """Compute soft_frac, I_curv (p=1), and omega_min from phonon frequencies.

    The 3 acoustic modes at each Gamma q-point are excluded by index before
    computing the metrics.

    Parameters
    ----------
    freqs_THz : np.ndarray
        Phonon frequencies in THz, shape (n_qpoints, n_modes).
    gamma_indices : np.ndarray
        Indices of Gamma q-points within the mesh.

    Returns
    -------
    dict with keys: soft_frac, I_curv, omega_min, n_modes_total,
                    n_modes_used, n_neg
    """
    freqs = np.asarray(freqs_THz)
    n_modes_per_q = freqs.shape[1] if freqs.ndim > 1 else len(freqs)

    all_freqs = freqs.ravel().copy()

    # Exclude 3 acoustic modes at each Gamma q-point by flattend index
    if len(gamma_indices) > 0 and n_modes_per_q >= 3:
        mask = np.ones(len(all_freqs), dtype=bool)
        for gi in gamma_indices:
            for j in range(3):
                idx = int(gi) * n_modes_per_q + j
                if idx < len(mask):
                    mask[idx] = False
        all_freqs = all_freqs[mask]

    n_used = len(all_freqs)
    neg_freqs = all_freqs[all_freqs < 0]
    n_neg = int(len(neg_freqs))

    soft_frac = float(n_neg / n_used) if n_used > 0 else 0.0
    # I_curv: mean negative curvature at p=1.
    # For stored negative frequencies (imaginary modes), curvature = omega^2.
    # For positive frequencies (stable), curvature = 0.
    # Matches vibroml.utils.stability_metrics.compute_i_curv(p=1).
    neg_curvature = np.where(all_freqs < 0, all_freqs ** 2, 0.0)
    i_curv = float(neg_curvature.mean()) if n_used > 0 else 0.0
    omega_min = float(np.min(all_freqs)) if n_used > 0 else 0.0

    return {
        "soft_frac": soft_frac,
        "I_curv": i_curv,
        "omega_min": omega_min,
        "n_modes_total": int(len(freqs.ravel())),
        "n_modes_used": int(n_used),
        "n_neg": n_neg,
    }


# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------
def make_convergence_plot(results, fig_dir):
    """Generate the 2-panel convergence figure (PDF + PNG).

    Panel A: soft_frac vs q-mesh size with per-structure lines and
             reference line at the pooled optimal separation (0.045).
    Panel B: I_curv vs q-mesh size with per-structure lines.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(fig_dir, exist_ok=True)

    # Organize results: compound -> {q: metrics}
    by_compound = {}
    for r in results:
        comp = r["compound"]
        q = r["q_mesh_size"]
        by_compound.setdefault(comp, {})[q] = r

    compounds = sorted(by_compound.keys())
    # Stable color cycle across structures
    cmap = plt.get_cmap("tab10")
    colors = {comp: cmap(i % 10) for i, comp in enumerate(compounds)}
    markers = {comp: m for comp, m in zip(compounds, ["o", "s", "^", "D", "v", "P"])}

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)

    # --- Panel A: soft_frac ---
    axA = axes[0]
    for comp in compounds:
        qs = sorted(by_compound[comp].keys())
        ys = [by_compound[comp][q]["soft_frac"] * 100 for q in qs]
        axA.plot(qs, ys, marker=markers.get(comp, "o"), color=colors[comp],
                 label=comp, linewidth=1.8, markersize=7)
    axA.axhline(SOFT_FRAC_REFERENCE * 100, color="gray", linestyle="--", linewidth=1.0,
                alpha=0.7, label=f"optimal sep. ({SOFT_FRAC_REFERENCE*100:.1f}%)")
    axA.set_xlabel("q-mesh size (n x n x n)")
    axA.set_ylabel(r"$f_\mathrm{soft}$ (%)")
    axA.set_title("Panel A: soft fraction vs q-mesh")
    axA.set_xticks(Q_MESH_SIZES)
    axA.legend(fontsize=8, loc="best", framealpha=0.9)
    axA.grid(True, alpha=0.3)

    # --- Panel B: I_curv ---
    axB = axes[1]
    for comp in compounds:
        qs = sorted(by_compound[comp].keys())
        ys = [by_compound[comp][q]["I_curv"] for q in qs]
        axB.plot(qs, ys, marker=markers.get(comp, "o"), color=colors[comp],
                 label=comp, linewidth=1.8, markersize=7)
    axB.set_xlabel("q-mesh size (n x n x n)")
    axB.set_ylabel(r"$I_{\mathrm{curv}}$ (THz$^2$)")
    axB.set_title("Panel B: mean negative curvature vs q-mesh")
    axB.set_xticks(Q_MESH_SIZES)
    axB.legend(fontsize=8, loc="best", framealpha=0.9)
    axB.grid(True, alpha=0.3)

    pdf_path = os.path.join(fig_dir, "convergence_sweep.pdf")
    png_path = os.path.join(fig_dir, "convergence_sweep.png")
    fig.savefig(pdf_path)
    fig.savefig(png_path, dpi=200)
    plt.close(fig)
    print(f"Convergence figure saved: {pdf_path}")
    print(f"Convergence figure saved: {png_path}")
    return pdf_path, png_path


# -----------------------------------------------------------------------------
# Main driver
# -----------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Convergence sweep: phonon stability metrics vs q-mesh size")
    parser.add_argument("--engine", default="mace",
                        help="MLIP engine (default: mace)")
    parser.add_argument("--model-name", default="medium-omat-0",
                        help="Model name override (default: medium-omat-0)")
    parser.add_argument("--device", default="cuda",
                        help="Device placement (default: cuda)")
    parser.add_argument("--output", default="convergence_sweep.json",
                        help="Output JSON file (default: convergence_sweep.json)")
    parser.add_argument("--structures", default="uq_validation_structures.json",
                        help="Structures JSON manifest")
    parser.add_argument("--delta", type=float, default=0.05,
                        help="Finite-displacement amplitude (Angstrom)")
    parser.add_argument("--use-symrelaxed", action="store_true", default=True,
                        help="Prefer _symrelaxed.cif if available (default: True)")
    parser.add_argument("--no-symrelaxed", dest="use_symrelaxed",
                        action="store_false",
                        help="Use original CIFs instead of symrelaxed")
    parser.add_argument("--fig-dir", default="figures_5eng",
                        help="Directory for convergence figure output")
    args = parser.parse_args()

    csv_output = args.output.replace(".json", ".csv")

    # Load structures and filter to high-symmetry aristotypes
    with open(args.structures) as f:
        structures = json.load(f)

    targets = [s for s in structures
               if s.get("stability_label") == "high_sym_unstable"]
    # Keep the 4 representative compounds in a stable order
    desired = ["LiF", "CsPbI3", "HfO2", "Bi2Sn2O7"]
    targets_by_comp = {t.get("compound"): t for t in targets}
    targets = [targets_by_comp[c] for c in desired if c in targets_by_comp]

    print("=" * 70)
    print("Convergence sweep: phonon stability metrics vs q-mesh size")
    print("=" * 70)
    print(f"Engine:      {args.engine}")
    print(f"Model:       {args.model_name}")
    print(f"Device:      {args.device}")
    print(f"Structures:  {len(targets)} ({[t['compound'] for t in targets]})")
    print(f"Q-meshes:    {Q_MESH_SIZES}")
    print(f"Output JSON: {args.output}")
    print(f"Output CSV:  {csv_output}")
    print(f"Figure dir:  {args.fig_dir}")
    print(f"Symrelaxed:  {args.use_symrelaxed}")
    print("=" * 70)

    # Load calculator once (reused across all structures / q-meshes)
    print(f"\nLoading calculator: engine={args.engine}, model={args.model_name}")
    calc = load_calculator(args.engine, model_name=args.model_name,
                           device=args.device)
    print(f"Calculator loaded: {type(calc).__name__}")

    from ase.io import read

    results = []

    for entry in targets:
        compound = entry.get("compound", "?")
        label = entry.get("stability_label", "unknown")
        cif_path = entry.get("filepath", "")

        if args.use_symrelaxed:
            symrelaxed = cif_path.replace(".cif", "_symrelaxed.cif")
            if os.path.exists(symrelaxed):
                cif_path = symrelaxed

        sc = entry.get("supercell", 2)
        if isinstance(sc, int):
            sc = (sc, sc, sc)
        else:
            sc = tuple(sc)

        print(f"\n--- {compound} ({label}) ---")
        print(f"  CIF: {cif_path}")
        print(f"  Supercell: {sc}")

        atoms = read(cif_path)
        print(f"  Natoms (primitive): {len(atoms)}")

        for qn in Q_MESH_SIZES:
            q_mesh = (qn, qn, qn)
            t0 = time.time()
            print(f"  q-mesh {q_mesh} ... ", end="", flush=True)
            try:
                freqs, gamma_indices = compute_phonon_freqs_uniform(
                    atoms, calc, supercell=sc, delta=args.delta,
                    q_mesh=q_mesh)
                metrics = compute_metrics(freqs, gamma_indices)
                dt = time.time() - t0

                record = {
                    "compound": compound,
                    "stability_label": label,
                    "cif_path": cif_path,
                    "supercell": list(sc),
                    "q_mesh": list(q_mesh),
                    "q_mesh_size": qn,
                    "natoms": len(atoms),
                    "engine": args.engine,
                    "model_name": args.model_name,
                    "delta": args.delta,
                    "time_s": dt,
                }
                record.update(metrics)
                results.append(record)

                print(f"soft_frac={metrics['soft_frac']:.4f} "
                      f"I_curv={metrics['I_curv']:.6f} "
                      f"omega_min={metrics['omega_min']:.4f} "
                      f"({dt:.1f}s)")

                # Save incrementally so partial results survive a crash
                with open(args.output, "w") as f:
                    json.dump(results, f, indent=2)

            except Exception as e:
                dt = time.time() - t0
                print(f"FAILED ({dt:.1f}s): {e}")
                import traceback
                traceback.print_exc()
                # Record the failure so the sweep continues
                results.append({
                    "compound": compound,
                    "stability_label": label,
                    "cif_path": cif_path,
                    "supercell": list(sc),
                    "q_mesh": list(q_mesh),
                    "q_mesh_size": qn,
                    "natoms": len(atoms),
                    "engine": args.engine,
                    "model_name": args.model_name,
                    "delta": args.delta,
                    "time_s": dt,
                    "error": str(e),
                    "soft_frac": None,
                    "I_curv": None,
                    "omega_min": None,
                })
                with open(args.output, "w") as f:
                    json.dump(results, f, indent=2)

    # ---------------- Write CSV ----------------
    fieldnames = [
        "compound", "stability_label", "q_mesh_size", "q_mesh",
        "supercell", "natoms", "soft_frac", "I_curv", "omega_min",
        "n_modes_total", "n_modes_used", "n_neg", "time_s",
        "engine", "model_name", "delta", "cif_path",
    ]
    with open(csv_output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for r in results:
            writer.writerow(r)
    print(f"\nCSV written: {csv_output}")

    # ---------------- Generate convergence plot ----------------
    # Only plot records that actually succeeded
    plot_results = [r for r in results if r.get("soft_frac") is not None]
    if plot_results:
        try:
            make_convergence_plot(plot_results, args.fig_dir)
        except Exception as e:
            print(f"WARNING: plot generation failed: {e}")
            import traceback
            traceback.print_exc()
    else:
        print("WARNING: no successful records to plot")

    # ---------------- Summary ----------------
    print("\n" + "=" * 70)
    print("Convergence sweep complete")
    print(f"  Records: {len(results)}")
    print(f"  JSON:    {args.output}")
    print(f"  CSV:     {csv_output}")
    print("=" * 70)


if __name__ == "__main__":
    main()
