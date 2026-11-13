#!/usr/bin/env python3
"""
Compute the continuous instability metric I_curv for a set of structures
using both MACE-OMAT and PET-MAD.

I_curv = (1/M) * sum_i max(0, -omega_i^2)   [units: THz^2]

where the sum runs over all non-acoustic modes on a uniform q-point mesh.
This is the mean negative curvature of the mass-weighted dynamical matrix.

For PET-MAD, the script also computes ensemble UQ statistics:
  - max_std: maximum frequency std across the ensemble (THz)
  - std_at_min: std at the q-point where the mean freq is most negative (THz)
  - min_minus_3std: mean min freq - 3*std_at_min (conservative lower bound)

Results are saved incrementally to both JSON and CSV.

Usage:
  # MACE-OMAT (vibroml_env):
  vibroml_env_python compute_omega_inst.py --engine mace --structures list.json

  # PET-MAD with UQ + symmetry-constrained relaxation (orb_env):
  orb_env_python compute_omega_inst.py --engine pet --structures list.json \
      --symmetry_constrained --phonon_uq
"""

import argparse
import csv
import json
import os
import shutil
import time
import numpy as np
from pathlib import Path

# Constants
EV_TO_THZ = 4.135667696e-3  # eV per THz


def load_calculator(engine, model_name="medium-omat-0", device="cuda"):
    """Load the ASE calculator for the given engine."""
    if engine == "mace":
        from mace.calculators import mace_mp
        print(f"  Loading MACE calculator (model={model_name}, device={device})...")
        calc = mace_mp(model=model_name, dispersion=False, default_dtype="float64",
                       device=device, stress=True)
        return calc
    elif engine == "pet":
        from upet.calculator import UPETCalculator
        pet_model = model_name if model_name != "medium-omat-0" else "pet-mad-xs"
        print(f"  Loading PET-MAD calculator (model={pet_model}, device={device})...")
        calc = UPETCalculator(model=pet_model, device=device)
        return calc
    else:
        raise ValueError(f"Unknown engine: {engine}")


def get_space_group(atoms, symprec=0.01):
    """Get the space group number and symbol."""
    try:
        import spglib
        # spglib needs (cell, scaled_positions, numbers) tuple, not Atoms
        cell = (atoms.cell.array, atoms.get_scaled_positions(), atoms.numbers)
        sg = spglib.get_spacegroup(cell, symprec=symprec)
        ds = spglib.get_symmetry_dataset(cell, symprec=symprec)
        sg_no = ds["number"]
        return sg_no, str(sg)
    except Exception:
        return -1, "unknown"


def relax_structure(atoms, calc, fmax=1e-3, max_steps=200,
                    symmetry_constrained=False, symprec=0.01):
    """Relax the structure using FIRE optimizer.

    If symmetry_constrained=True, apply ASE FixSymmetry constraint to
    preserve the initial space group during relaxation. Both atomic
    positions and cell parameters are allowed to relax, but only within
    the symmetry operations of the initial structure.
    """
    from ase.optimize import FIRE

    atoms.calc = calc

    if symmetry_constrained:
        from ase.constraints import FixSymmetry
        constraint = FixSymmetry(atoms, symprec=symprec,
                                 adjust_positions=True,
                                 adjust_cell=True)
        atoms.set_constraint(constraint)
        print(f"    Symmetry-constrained relaxation (symprec={symprec})...")

    try:
        opt = FIRE(atoms, logfile=None)
        opt.run(fmax=fmax, steps=max_steps)
        print(f"    Relaxed: fmax={atoms.get_forces().max():.6f} eV/A")
    except Exception as e:
        print(f"    Relaxation failed: {e}, using original structure")
    return atoms


def compute_phonon_eigenvalues(atoms, calc, supercell=(2, 2, 2), delta=0.05,
                               q_mesh=(8, 8, 8)):
    """
    Compute phonon frequencies on a uniform q-point mesh using ASE Phonons.

    Returns:
        freqs_THz: ndarray of shape (N_q, 3*N_atoms) - frequencies in THz
                   (negative = imaginary)
        gamma_indices: indices of q-points at Gamma (q=0)
    """
    from ase.phonons import Phonons
    from ase.dft.kpoints import monkhorst_pack

    atoms.calc = calc

    ph_dir = "phonon_temp_omega_inst"
    if os.path.exists(ph_dir):
        shutil.rmtree(ph_dir)

    ph = Phonons(atoms, calc, supercell=supercell, delta=delta, name=ph_dir)

    print(f"    Computing force constants (supercell={supercell}, delta={delta})...")
    ph.run()
    ph.read()

    qpoints = monkhorst_pack(q_mesh)
    N_q = len(qpoints)
    print(f"    Computing eigenvalues on {q_mesh} mesh ({N_q} q-points)...")

    energies_eV = ph.band_structure(qpoints, modes=False, verbose=False)
    freqs_THz = energies_eV / EV_TO_THZ

    gamma_mask_q = np.all(np.abs(qpoints) < 1e-10, axis=1)
    gamma_indices = np.where(gamma_mask_q)[0]

    ph.clean()

    return freqs_THz, gamma_indices


def compute_pet_uq_phonons(atoms, model_name, supercell, delta, q_mesh, device,
                           workdir):
    """
    Compute PET-MAD ensemble phonon UQ using uqphonon's PhononEnsemble.

    Returns:
        mean_freqs: ndarray (N_q, 3N) - mean frequencies in THz
        std_freqs: ndarray (N_q, 3N) - std frequencies in THz
        gamma_indices: indices of q-points at Gamma
        n_ensemble: int, number of ensemble members
    """
    from uqphonon import PhononEnsemble
    from ase.dft.kpoints import monkhorst_pack
    import upet

    # Export model to TorchScript if not cached
    model_base, model_size = model_name.rsplit("-", 1)
    model_cache = os.path.join(os.path.expanduser("~"), ".cache", "vibroml_pet_models")
    os.makedirs(model_cache, exist_ok=True)
    model_pt = os.path.join(model_cache, f"{model_base.replace('-','_')}_{model_size}.pt")
    if not os.path.exists(model_pt):
        print(f"    Exporting PET model to TorchScript: {model_pt}")
        upet.save_upet(model=model_base, size=model_size, output=model_pt)
    else:
        print(f"    Reusing cached TorchScript model: {model_pt}")

    sc_matrix = np.diag(supercell)

    print(f"    Setting up PhononEnsemble (supercell={supercell}, delta={delta})...")
    ensemble = PhononEnsemble(
        atoms=atoms,
        supercell_matrix=sc_matrix,
        model=model_pt,
        device=device,
        primitive_matrix=np.eye(3),
    )
    ensemble.compute_displacements(distance=delta)

    uq_workdir = os.path.join(workdir, f"_uqphonon_omega_inst")
    print(f"    Running i-PI ensemble forces (workdir={uq_workdir})...")
    ensemble.run_forces(workdir=uq_workdir)

    # Compute bands on uniform q-point mesh
    qpoints = monkhorst_pack(q_mesh)
    # Pass as a list with a single segment (raw q-point array)
    print(f"    Computing ensemble phonon bands on {q_mesh} mesh ({len(qpoints)} q-points)...")
    ensemble.compute_bands(bandpath=[qpoints], labels=["G"])

    # Extract per-member frequencies
    member_freqs = []
    for ph_member in ensemble.phonons:
        # ph_member.band_structure.frequencies is a list of segments
        # Each segment is (N_q, 3N) in THz (phonopy native units)
        seg_freqs = np.concatenate([seg.ravel() for seg in ph_member.band_structure.frequencies])
        member_freqs.append(seg_freqs)
    member_freqs = np.array(member_freqs)  # (n_ensemble, N_q * 3N)

    # Mean phonon frequencies
    mean_ph = ensemble.mean_phonon
    mean_freqs_segs = mean_ph.band_structure.frequencies
    mean_freqs = np.concatenate([seg.ravel() for seg in mean_freqs_segs])
    # Reshape to (N_q, 3N)
    n_bands = len(atoms) * 3
    mean_freqs = mean_freqs.reshape(len(qpoints), n_bands)

    # Std across ensemble members
    all_std = member_freqs.std(axis=0)  # (N_q * 3N,)
    std_freqs = all_std.reshape(len(qpoints), n_bands)

    # Gamma indices
    gamma_mask_q = np.all(np.abs(qpoints) < 1e-10, axis=1)
    gamma_indices = np.where(gamma_mask_q)[0]

    n_ensemble = member_freqs.shape[0]

    return mean_freqs, std_freqs, gamma_indices, n_ensemble


def compute_instability_metric(freqs_THz, gamma_indices=None):
    """
    Compute the instability metric I_curv from phonon frequencies on a
    uniform q-point mesh.

    I_curv = (1/M) * sum_i max(0, -omega_i^2)   [units: THz^2]

    The 3 acoustic modes at Gamma are excluded (rigid translations).
    No frequency threshold is applied — the only boundary is lambda=0.

    omega_min is computed from ALL modes, clamped to 0 if positive.
    """
    freqs = freqs_THz.ravel()
    N_q, N_bands = freqs_THz.shape

    # omega_min: from ALL modes, clamped to 0 if positive
    omega_min = float(min(freqs.min(), 0.0))

    # Build exclusion mask: only the 3 acoustic modes at Gamma by index
    exclude_mask = np.zeros_like(freqs_THz, dtype=bool)
    if gamma_indices is not None and len(gamma_indices) > 0:
        for gi in gamma_indices:
            exclude_mask[gi, :3] = True

    kept = freqs_THz[~exclude_mask]
    n_modes = len(kept)
    n_excluded = int(exclude_mask.sum())

    if n_modes == 0:
        return {
            "i_curv": 0.0,
            "omega_inst": 0.0,
            "omega_min": omega_min,
            "soft_frac": 0.0,
            "n_modes": 0,
            "n_imag": 0,
            "n_excluded": n_excluded,
        }

    n_imag = int((kept < 0).sum())
    soft_frac = float(n_imag) / n_modes

    neg_curvature = np.where(kept < 0, kept ** 2, 0.0)
    i_curv = float(neg_curvature.mean())
    omega_inst = float(np.sqrt(i_curv))

    return {
        "i_curv": i_curv,
        "omega_inst": omega_inst,
        "omega_min": omega_min,
        "soft_frac": soft_frac,
        "n_modes": n_modes,
        "n_imag": n_imag,
        "n_excluded": n_excluded,
    }


def process_structure(cif_path, engine, calc, supercell, q_mesh, delta,
                      relax=True, fmax=1e-3, symmetry_constrained=False,
                      phonon_uq=False, model_name="medium-omat-0",
                      device="cuda", workdir="."):
    """Process a single structure: load, relax, compute phonons, compute metric."""
    from ase.io import read

    print(f"\n  Loading: {os.path.basename(cif_path)}")
    atoms = read(cif_path)

    # Get space group before relaxation
    sg_no, sg_sym = get_space_group(atoms)

    t_start = time.time()

    if relax:
        print(f"  Relaxing with {engine}...")
        atoms = relax_structure(atoms, calc, fmax=fmax,
                                symmetry_constrained=symmetry_constrained)
        if atoms.constraints:
            atoms.set_constraint()
    else:
        atoms.calc = calc

    # Check if structure changed space group
    sg_no_relaxed, sg_sym_relaxed = get_space_group(atoms)

    if phonon_uq and engine == "pet":
        # PET-MAD ensemble UQ path
        print(f"  Computing PET-MAD ensemble UQ phonons...")
        try:
            mean_freqs, std_freqs, gamma_indices, n_ens = compute_pet_uq_phonons(
                atoms, model_name, supercell, delta, q_mesh, device, workdir
            )
        except Exception as e:
            print(f"  PET UQ PHONON COMPUTATION FAILED: {e}")
            import traceback
            traceback.print_exc()
            return None

        print(f"  Computing instability metric from ensemble mean...")
        result = compute_instability_metric(mean_freqs, gamma_indices=gamma_indices)

        # UQ statistics
        # max_std: maximum std across all q-points and bands
        max_std = float(std_freqs.max())

        # std_at_min: std at the q-point/band where mean freq is most negative
        min_idx = np.unravel_index(np.argmin(mean_freqs), mean_freqs.shape)
        std_at_min = float(std_freqs[min_idx])

        # min_minus_3std: conservative lower bound on the minimum frequency
        min_minus_3std = float(mean_freqs[min_idx] - 3 * std_at_min)

        result["max_std"] = max_std
        result["std_at_min"] = std_at_min
        result["min_minus_3std"] = min_minus_3std
        result["n_ensemble"] = n_ens

        print(f"  UQ: max_std={max_std:.4f}, std_at_min={std_at_min:.4f}, "
              f"min_minus_3std={min_minus_3std:.4f} (n_ens={n_ens})")

    else:
        # Standard path (MACE or PET without UQ)
        print(f"  Computing phonons (supercell={supercell}, q_mesh={q_mesh})...")
        try:
            freqs, gamma_indices = compute_phonon_eigenvalues(
                atoms, calc, supercell=supercell, delta=delta, q_mesh=q_mesh
            )
        except Exception as e:
            print(f"  PHONON COMPUTATION FAILED: {e}")
            return None

        print(f"  Computing instability metric...")
        result = compute_instability_metric(freqs, gamma_indices=gamma_indices)

    t_elapsed = time.time() - t_start

    # Add metadata
    result["engine"] = engine
    result["cif_path"] = cif_path
    result["natoms"] = len(atoms)
    result["supercell"] = list(supercell) if hasattr(supercell, '__iter__') else [supercell]*3
    result["q_mesh"] = list(q_mesh)
    result["delta"] = delta
    result["relaxed"] = relax
    result["symmetry_constrained"] = symmetry_constrained
    result["space_group_no"] = sg_no
    result["space_group_sym"] = sg_sym
    result["space_group_no_relaxed"] = sg_no_relaxed
    result["space_group_sym_relaxed"] = sg_sym_relaxed
    result["time_s"] = t_elapsed

    print(f"  RESULT: I_curv={result['i_curv']:.4f} THz^2, "
          f"omega_min={result['omega_min']:.4f} THz, "
          f"soft_frac={result['soft_frac']:.4f}, "
          f"n_imag={result['n_imag']}/{result['n_modes']}, "
          f"time={t_elapsed:.1f}s")

    return result


# CSV field order — all possible columns
CSV_FIELDS = [
    "idx", "category", "compound",
    "space_group_no", "space_group_sym",
    "space_group_no_relaxed", "space_group_sym_relaxed",
    "natoms", "supercell", "q_mesh", "delta",
    "engine", "relaxed", "symmetry_constrained",
    "i_curv", "omega_inst", "omega_min", "soft_frac",
    "n_modes", "n_imag", "n_excluded",
    "max_std", "std_at_min", "min_minus_3std", "n_ensemble",
    "time_s", "cif_path",
]


def save_csv(results, csv_path):
    """Save results to CSV."""
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_FIELDS, extrasaction="ignore")
        writer.writeheader()
        for r in results:
            # Flatten list fields for CSV
            row = dict(r)
            if isinstance(row.get("supercell"), list):
                row["supercell"] = "x".join(str(x) for x in row["supercell"])
            if isinstance(row.get("q_mesh"), list):
                row["q_mesh"] = "x".join(str(x) for x in row["q_mesh"])
            writer.writerow(row)


def main():
    parser = argparse.ArgumentParser(
        description="Compute I_curv instability metric for structures."
    )
    parser.add_argument("--engine", choices=["mace", "pet"], required=True,
                        help="Calculator engine to use")
    parser.add_argument("--model_name", default="medium-omat-0",
                        help="Model name (default: medium-omat-0 for MACE, "
                             "pet-mad-xs for PET)")
    parser.add_argument("--structures", required=True,
                        help="JSON file with list of structures")
    parser.add_argument("--output", default="omega_inst_results.json",
                        help="Output JSON file for results")
    parser.add_argument("--csv_output", default=None,
                        help="Output CSV file (default: same as output with .csv)")
    parser.add_argument("--q_mesh", type=int, nargs=3, default=[8, 8, 8],
                        help="Uniform q-point mesh (default: 8 8 8)")
    parser.add_argument("--delta", type=float, default=0.05,
                        help="Finite displacement amplitude (default: 0.05 A)")
    parser.add_argument("--fmax", type=float, default=1e-3,
                        help="Relaxation force criterion (default: 1e-3 eV/A)")
    parser.add_argument("--no_relax", action="store_true",
                        help="Skip relaxation (use structure as-is)")
    parser.add_argument("--symmetry_constrained", action="store_true",
                        help="Preserve initial space group during relaxation")
    parser.add_argument("--phonon_uq", action="store_true",
                        help="Compute PET-MAD ensemble UQ (requires uqphonon + i-PI)")
    parser.add_argument("--device", default="cuda",
                        help="Device for calculator (default: cuda)")
    parser.add_argument("--max_structures", type=int, default=None,
                        help="Limit number of structures (for testing)")
    args = parser.parse_args()

    # Determine CSV output path
    csv_path = args.csv_output or args.output.replace(".json", ".csv")

    # Load structures
    with open(args.structures) as f:
        structures = json.load(f)

    if args.max_structures:
        structures = structures[:args.max_structures]

    # Resume logic
    results = []
    completed_keys = set()
    if os.path.exists(args.output):
        try:
            with open(args.output) as f:
                results = json.load(f)
            for r in results:
                key = f"{r.get('compound','')}_{r.get('idx',-1)}"
                completed_keys.add(key)
            print(f"Resuming: {len(completed_keys)} structures already completed, "
                  f"{len(structures) - len(completed_keys)} remaining")
        except Exception as e:
            print(f"Warning: could not load existing results ({e}), starting fresh")
            results = []

    print(f"Engine: {args.engine}")
    print(f"Model: {args.model_name}")
    print(f"Structures: {len(structures)}")
    print(f"Q-mesh: {args.q_mesh}")
    print(f"Delta: {args.delta}")
    print(f"Relax: {not args.no_relax}")
    print(f"Symmetry constrained: {args.symmetry_constrained}")
    print(f"Phonon UQ: {args.phonon_uq}")
    print(f"Output: {args.output}")
    print(f"CSV output: {csv_path}")

    # Load calculator
    calc = load_calculator(args.engine, args.model_name, args.device)

    workdir = os.path.dirname(os.path.abspath(args.output))

    # Process each structure
    for i, entry in enumerate(structures):
        cif_path = entry["filepath"]
        sc = entry.get("supercell", 2)
        if isinstance(sc, int):
            sc = (sc, sc, sc)
        else:
            sc = tuple(sc)

        compound = entry.get("compound", os.path.basename(cif_path))
        category = entry.get("category", "unknown")

        # Skip if already completed
        key = f"{compound}_{i}"
        if key in completed_keys:
            continue

        print(f"\n[{i+1}/{len(structures)}] {compound} ({category})")

        result = process_structure(
            cif_path, args.engine, calc, sc, tuple(args.q_mesh),
            args.delta, relax=not args.no_relax, fmax=args.fmax,
            symmetry_constrained=args.symmetry_constrained,
            phonon_uq=args.phonon_uq, model_name=args.model_name,
            device=args.device, workdir=workdir
        )

        if result is not None:
            result["compound"] = compound
            result["category"] = category
            result["idx"] = i
            results.append(result)

            # Save incrementally (both JSON and CSV)
            with open(args.output, "w") as f:
                json.dump(results, f, indent=2)
            save_csv(results, csv_path)

    # Final summary
    print(f"\n{'='*60}")
    print(f"Completed: {len(results)}/{len(structures)}")
    print(f"Results saved to: {args.output}")
    print(f"CSV saved to: {csv_path}")
    print(f"\nSummary:")
    print(f"{'compound':25s} {'cat':15s} {'I_curv':>10s} {'wmin':>8s} {'soft%':>7s} {'time':>7s}")
    for r in results:
        print(f"{r['compound']:25s} {r.get('category',''):15s} "
              f"{r['i_curv']:10.4f} {r['omega_min']:8.4f} "
              f"{100*r['soft_frac']:7.1f} {r.get('time_s',0):7.1f}s")


if __name__ == "__main__":
    main()
