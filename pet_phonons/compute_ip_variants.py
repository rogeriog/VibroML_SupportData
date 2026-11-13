#!/usr/bin/env python3
"""
Post-processing script: recompute uniform-mesh phonon frequencies for
benchmark structures and compute I_p for multiple p values.

This script:
  1. Loads the relaxed CIF for each benchmark structure
  2. Recomputes phonons on a uniform q-mesh (read-only on calc files)
  3. Saves raw frequencies to a NEW npz file (does not modify calc outputs)
  4. Computes I_p for p = 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0, inf
  5. Saves results to a new JSON/CSV

Usage:
  # MACE-OMAT (default model)
  python compute_ip_variants.py --engine mace --model-name medium-omat-0 \
      --structures uq_validation_structures.json --output ip_variants_mace.json

  # MACE-MP
  python compute_ip_variants.py --engine mace --model-name medium \
      --structures uq_validation_structures.json --output ip_variants_macemp.json

  # eSEN
  python compute_ip_variants.py --engine esen \
      --structures uq_validation_structures.json --output ip_variants_esen.json

  # UMA
  python compute_ip_variants.py --engine uma \
      --structures uq_validation_structures.json --output ip_variants_uma.json

  # PET-MAD
  python compute_ip_variants.py --engine pet \
      --structures uq_validation_structures.json --output ip_variants_pet.json
"""
import argparse, csv, json, os, shutil, sys, time
import numpy as np
from pathlib import Path

EV_TO_THZ = 4.135667696e-3

# p values to test
P_VALUES = [0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0, 20.0, 50.0]

# Add VibroML to path for initialize_calculator
VIBROML_PATH = "/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/01_codebase/VibroML"
if VIBROML_PATH not in sys.path:
    sys.path.insert(0, VIBROML_PATH)


def load_calculator(engine, model_name=None, device="cuda"):
    """Load a calculator via VibroML's canonical initialize_calculator.

    This delegates to vibroml.utils.structure_utils.initialize_calculator
    which handles all engine-specific setup (checkpoint paths, model
    downloads, device placement).
    """
    from vibroml.utils.structure_utils import initialize_calculator

    # Default model names per engine
    defaults = {
        "mace": "medium-omat-0",   # MACE-OMAT-0 medium
        "macemp": "medium",         # MACE-MP-0 medium
        "esen": None,               # Uses checkpoint path
        "uma": None,                # Uses checkpoint path
        "pet": "pet-mad-xs",
    }

    # Normalize engine name: "macemp" → "mace" with model "medium"
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


def get_space_group(atoms, symprec=0.01):
    try:
        import spglib
        cell = (atoms.cell.array, atoms.get_scaled_positions(), atoms.numbers)
        ds = spglib.get_symmetry_dataset(cell, symprec=symprec)
        return ds["number"], str(spglib.get_spacegroup(cell, symprec=symprec))
    except:
        return -1, "unknown"


def compute_phonon_freqs_uniform(atoms, calc, supercell=(2,2,2), delta=0.05,
                                  q_mesh=(6,6,6)):
    """Compute phonon frequencies on a uniform q-mesh. Returns raw freq array."""
    from ase.phonons import Phonons
    from ase.dft.kpoints import monkhorst_pack
    import tempfile

    atoms.calc = calc
    # Use a unique temp directory per call to avoid stale cache from prior structures
    ph_dir = tempfile.mkdtemp(prefix="phonon_ip_", dir=".")

    try:
        ph = Phonons(atoms, calc, supercell=supercell, delta=delta, name=ph_dir)
        ph.run()
        ph.read()

        qpoints = monkhorst_pack(q_mesh)
        energies_eV = ph.band_structure(qpoints, modes=False, verbose=False)
        freqs_THz = energies_eV / EV_TO_THZ

        # Gamma indices
        gamma_mask = np.all(np.abs(qpoints) < 1e-10, axis=1)
        gamma_indices = np.where(gamma_mask)[0]

        ph.clean()
    finally:
        # Ensure cleanup even on error
        shutil.rmtree(ph_dir, ignore_errors=True)

    return freqs_THz, gamma_indices


def compute_ip(freqs_THz, gamma_indices=None, p_values=P_VALUES):
    """
    Compute I_p for multiple p values.

    I_p = [ (1/M) * sum_i (max(0, -omega_i^2))^p ]^(1/p)
        = [ (1/M) * sum_{omega<0} |omega_i|^{2p} ]^(1/p)

    Omega_p = sqrt(I_p) = [ (1/M) * sum |omega|^{2p} ]^(1/(2p))

    For p->inf: I_inf = max(|omega|^2) = |omega_min|^2

    Acoustic modes at Gamma are excluded (3 translational modes).
    No frequency threshold is applied — the only boundary is lambda=0.
    """
    freqs = freqs_THz.copy()

    # Exclude 3 acoustic modes at Gamma by index
    if gamma_indices is not None and len(gamma_indices) > 0:
        for gi in gamma_indices:
            freqs[gi, :3] = 0.0  # set to exactly zero (not unstable)

    freqs_flat = freqs.ravel()
    M = len(freqs_flat)

    # Negative curvature: max(0, -omega^2) = |omega|^2 for omega < 0
    neg_curv = np.where(freqs_flat < 0, freqs_flat**2, 0.0)

    # omega_min from ALL modes (including acoustic, but they're 0 now)
    omega_min = float(min(freqs_flat.min(), 0.0))

    # soft_frac: fraction of modes with omega < 0 (no threshold)
    n_neg = int((freqs_flat < 0).sum())
    soft_frac_no_thresh = float(n_neg) / M

    # soft_frac with 0.05 THz threshold (for comparison)
    n_neg_thresh = int((freqs_flat < -0.05).sum())
    soft_frac_thresh = float(n_neg_thresh) / M

    results = {
        "omega_min": omega_min,
        "soft_frac_no_thresh": soft_frac_no_thresh,
        "soft_frac_0p05": soft_frac_thresh,
        "n_modes": M,
        "n_neg": n_neg,
        "n_neg_thresh": n_neg_thresh,
    }

    for p in p_values:
        if p >= 50:
            # p -> inf: max negative curvature
            i_p = float(neg_curv.max())
        else:
            i_p = float((neg_curv**p).mean()**(1.0/p))
        omega_p = float(np.sqrt(i_p))
        results[f"i_p_{p}"] = i_p
        results[f"omega_p_{p}"] = omega_p

    # p = infinity
    i_inf = float(neg_curv.max())
    results["i_p_inf"] = i_inf
    results["omega_p_inf"] = float(np.sqrt(i_inf))

    return results


def main():
    parser = argparse.ArgumentParser(description="Compute I_p variants for benchmark structures")
    parser.add_argument("--engine", choices=["mace", "macemp", "esen", "uma", "pet"],
                        required=True,
                        help="MLIP engine: mace (MACE-OMAT), macemp (MACE-MP), "
                             "esen, uma, pet (PET-MAD)")
    parser.add_argument("--model-name", default=None,
                        help="Model name override (e.g. medium-omat-0, medium, "
                             "pet-mad-xs). Defaults per engine.")
    parser.add_argument("--structures", required=True)
    parser.add_argument("--output", default=None)
    parser.add_argument("--q_mesh", type=int, nargs=3, default=[6, 6, 6])
    parser.add_argument("--delta", type=float, default=0.05)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--max_structures", type=int, default=None)
    parser.add_argument("--use-symrelaxed", action="store_true", default=True,
                        help="Use _symrelaxed.cif if available (default: True)")
    parser.add_argument("--only-idx", type=int, nargs='+', default=None,
                        help="Restrict to these structure 'order' indices "
                             "(e.g. --only-idx 34 39). Used to recompute a subset.")
    parser.add_argument("--supercell-override", type=int, default=None,
                        help="Force a cubic (N,N,N) supercell for every selected "
                             "structure, ignoring the value in the structures JSON. "
                             "Needed for large cells that defaulted to 1x1x1 and "
                             "therefore produced no phonon dispersion.")
    parser.add_argument("--freqs-dir", default=None,
                        help="Directory for raw-frequency NPZ output. Defaults to "
                             "ip_variants_freqs_<engine>. Override to avoid "
                             "overwriting a previously validated dataset.")
    args = parser.parse_args()

    # Output naming: use engine name for files
    engine_label = args.engine
    output = args.output or f"ip_variants_{engine_label}.json"
    csv_output = output.replace(".json", ".csv")
    freqs_dir = args.freqs_dir or f"ip_variants_freqs_{engine_label}"
    os.makedirs(freqs_dir, exist_ok=True)

    # Load structures
    with open(args.structures) as f:
        structures = json.load(f)

    # Filter to benchmark only
    bench = [s for s in structures if s.get('category') == 'benchmark']
    if args.only_idx:
        wanted = set(args.only_idx)
        bench = [s for s in bench if s.get('order') in wanted]
        print(f"Restricted to idx {sorted(wanted)}: {len(bench)} structure(s) selected")
    if args.max_structures:
        bench = bench[:args.max_structures]

    print(f"Engine: {args.engine}")
    print(f"Structures: {len(bench)} benchmark")
    print(f"Q-mesh: {args.q_mesh}")
    print(f"Output: {output}")
    print(f"Freqs dir: {freqs_dir}")

    # Resume logic
    results = []
    completed = set()
    if os.path.exists(output):
        with open(output) as f:
            results = json.load(f)
        completed = {r.get('idx', -1) for r in results}
        print(f"Resuming: {len(completed)} completed, {len(bench) - len(completed)} remaining")

    # Load calculator
    print(f"\nLoading calculator: engine={args.engine}, model={args.model_name or 'default'}")
    calc = load_calculator(args.engine, model_name=args.model_name, device=args.device)
    print(f"Calculator loaded: {type(calc).__name__}")

    # Process each structure
    for i, entry in enumerate(bench):
        idx = entry.get('order', i + 1)

        # Use the relaxed CIF if available (from previous UQ runs)
        # Otherwise use the original filepath
        cif_path = entry.get('filepath', '')

        # Use the symrelaxed CIF if available (for all engines)
        if args.use_symrelaxed:
            symrelaxed = cif_path.replace('.cif', '_symrelaxed.cif')
            if os.path.exists(symrelaxed):
                cif_path = symrelaxed

        if idx in completed:
            continue

        compound = entry.get('compound', '?')
        label = entry.get('stability_label', 'unknown')
        if args.supercell_override is not None:
            n = args.supercell_override
            sc = (n, n, n)
        else:
            sc = entry.get('supercell', 2)
            if isinstance(sc, int):
                sc = (sc, sc, sc)
            else:
                sc = tuple(sc)

        print(f"\n[{i+1}/{len(bench)}] idx={idx}: {compound} ({label})")

        from ase.io import read
        atoms = read(cif_path)
        sg_no, sg_sym = get_space_group(atoms)

        t_start = time.time()

        try:
            freqs, gamma_indices = compute_phonon_freqs_uniform(
                atoms, calc, supercell=sc, delta=args.delta,
                q_mesh=tuple(args.q_mesh))

            # Save raw frequencies
            freqs_path = os.path.join(freqs_dir, f"idx_{idx:03d}_{compound}_freqs.npz")
            np.savez(freqs_path, freqs=freqs, gamma_indices=gamma_indices,
                     compound=compound, idx=idx, q_mesh=args.q_mesh)

            # Compute I_p
            ip_results = compute_ip(freqs, gamma_indices=gamma_indices)
            ip_results['idx'] = idx
            ip_results['compound'] = compound
            ip_results['stability_label'] = label
            ip_results['space_group'] = sg_sym
            ip_results['natoms'] = len(atoms)
            ip_results['supercell'] = list(sc)
            ip_results['q_mesh'] = list(args.q_mesh)
            ip_results['cif_path'] = cif_path
            ip_results['time_s'] = time.time() - t_start
            ip_results['freqs_file'] = freqs_path

            results.append(ip_results)

            # Save incrementally
            with open(output, 'w') as f:
                json.dump(results, f, indent=2)

            print(f"  omega_min={ip_results['omega_min']:.4f}, "
                  f"I_p_1={ip_results['i_p_1.0']:.6f}, "
                  f"Omega_p_1={ip_results['omega_p_1.0']:.6f}, "
                  f"I_p_2={ip_results['i_p_2.0']:.6f}, "
                  f"soft_frac={ip_results['soft_frac_no_thresh']:.4f}, "
                  f"time={ip_results['time_s']:.1f}s")

        except Exception as e:
            print(f"  FAILED: {e}")
            import traceback
            traceback.print_exc()

    print(f"\n{'='*60}")
    print(f"Completed: {len(results)}/{len(bench)}")
    print(f"Results: {output}")
    print(f"Raw freqs: {freqs_dir}/")


if __name__ == "__main__":
    main()
