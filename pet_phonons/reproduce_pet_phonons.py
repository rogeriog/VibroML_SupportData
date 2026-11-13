"""
Minimal reproduction of the PET-MAD phonon dispersion recipe.

Covers:
  1. Al (FCC) — deterministic phonons with FixSymmetry relaxation
  2. BaTiO3 cubic — phonon UQ with uqphonon ensemble (imaginary soft mode)

Run with:
    PYTHONPATH=/path/to/VibroML \
    /gpfs/scratch/acad/htforft/rgouvea/conda_envs/orb_env/bin/python3.11 \
    reproduce_pet_phonons.py

Plots are saved as PNG files in the current directory.
"""

import warnings
import tempfile
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from ase import Atoms
from ase.build import bulk
from ase.constraints import FixSymmetry
from ase.filters import FrechetCellFilter
from ase.optimize import LBFGS

import spglib
import upet
from upet.calculator import UPETCalculator
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms

from uqphonon import PhononEnsemble
from uqphonon._core import _phonopy_to_ase, _resolve_bandpath
from uqphonon._plot import plot_bands

warnings.filterwarnings("ignore", category=RuntimeWarning,
                        message="logm result may be inaccurate")

# ─────────────────────────────────────────────────────────────────────────────
# Global settings — kept small for speed; increase for production
# ─────────────────────────────────────────────────────────────────────────────
FMAX    = 1e-4
STEPS   = 500
DELTA   = 0.05

try:
    import torch as _torch
    DEVICE = "cuda" if _torch.cuda.is_available() else "cpu"
except ImportError:
    DEVICE = "cpu"
print(f"Using device: {DEVICE}")

MODEL_BASE    = "pet-mad"
MODEL_SIZE    = "xs"
MODEL_VERSION = "1.5.0"
MODEL_NAME    = f"{MODEL_BASE}-{MODEL_SIZE}"

SUPERCELL_AL  = (3, 3, 3)   # recipe uses (4,4,4); 3x3x3 is faster
SUPERCELL_BTO = (2, 2, 2)


# ─────────────────────────────────────────────────────────────────────────────
# Initialise calculator and export TorchScript model
# ─────────────────────────────────────────────────────────────────────────────
print("Initialising PET-MAD calculator …")
calc = UPETCalculator(
    model=MODEL_NAME,
    device=DEVICE,
    dtype="float32",
    version=MODEL_VERSION,
)

MODEL_PT = "model_xs.pt"
print(f"Exporting TorchScript model to {MODEL_PT} …")
upet.save_upet(model=MODEL_BASE, size=MODEL_SIZE, version=MODEL_VERSION,
               output=MODEL_PT)


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def report_symmetry(atoms, label=""):
    cell = (atoms.get_cell(), atoms.get_scaled_positions(), atoms.get_atomic_numbers())
    sg = spglib.get_spacegroup(cell, symprec=1e-2)
    print(f"  {label:35s}  {sg}")


def compute_phonons_ase(atoms, supercell, npoints=50, label=""):
    """Phonon dispersion using ASE + phonopy (single-model, no UQ)."""
    print(f"\n--- {label} ---")
    ph_atoms = PhonopyAtoms(
        symbols=atoms.get_chemical_symbols(),
        cell=atoms.cell[:],
        scaled_positions=atoms.get_scaled_positions(),
    )
    phonon = Phonopy(ph_atoms, supercell_matrix=np.diag(supercell),
                     primitive_matrix=np.eye(3))
    phonon.generate_displacements(distance=DELTA)

    forces = []
    for sc in phonon.supercells_with_displacements:
        sc_ase = Atoms(symbols=sc.symbols, cell=sc.cell,
                       scaled_positions=sc.scaled_positions, pbc=True)
        sc_ase.calc = atoms.calc
        forces.append(sc_ase.get_forces())

    phonon.forces = np.array(forces)
    phonon.produce_force_constants()
    phonon.symmetrize_force_constants()

    prim = _phonopy_to_ase(phonon.primitive)
    qpts, conns, lbls = _resolve_bandpath(prim, None, npoints)
    phonon.run_band_structure(qpts, path_connections=conns, labels=lbls)
    print(f"  {len(forces)} displacements | {len(forces[0])} atoms/supercell")
    return phonon


def compute_phonons_uq(atoms, supercell, npoints=50, label=""):
    """Phonon dispersion with UQ using uqphonon + i-PI.

    Each call runs the i-PI force evaluation in an isolated subprocess so that
    i-PI's global softexit singleton does not kill the parent process after the
    first ensemble run — which would prevent a second call from completing.
    """
    print(f"\n--- {label} (UQ ensemble) ---")

    # ── Step 1: compute displaced supercells in-process (no i-PI involved) ────
    ensemble = PhononEnsemble(
        atoms=atoms,
        supercell_matrix=np.diag(supercell),
        model=str(Path(MODEL_PT).resolve()),
        device=DEVICE,
        primitive_matrix=np.eye(3),
    )
    ensemble.compute_displacements(distance=DELTA)

    # ── Step 2: run i-PI forces in an isolated subprocess ─────────────────────
    # i-PI uses a global singleton 'softexit' that sends SIGTERM after the
    # simulation ends.  A second call in the same process is killed immediately.
    # Running in a subprocess gives each call a fresh Python interpreter state.
    import sys, subprocess, pickle, ase.io

    workdir = Path(tempfile.gettempdir()) / f"uqphonon_{label.replace(' ', '_')}"
    workdir.mkdir(parents=True, exist_ok=True)

    # Serialise displaced supercells to disk so the subprocess can read them
    displacements_file = workdir / "displaced_supercells.traj"
    from uqphonon._core import _phonopy_to_ase
    displacements_ase = [_phonopy_to_ase(sc) for sc in ensemble._phonon.supercells_with_displacements]
    ase.io.write(str(displacements_file), displacements_ase)

    supercell_ase = _phonopy_to_ase(ensemble._phonon.supercell)
    supercell_file = workdir / "supercell.xyz"
    ase.io.write(str(supercell_file), supercell_ase)

    forces_file = workdir / "forces.npy"

    # Small script that runs inside the subprocess
    worker_code = f"""
import sys, numpy as np, ase.io
from pathlib import Path
from uqphonon._ipi import run_ipi_forces

displacements = ase.io.read("{displacements_file}", index=":")
supercell     = ase.io.read("{supercell_file}")
forces = run_ipi_forces(
    supercell=supercell,
    displacements=displacements,
    model="{str(Path(MODEL_PT).resolve())}",
    device="{DEVICE}",
    workdir="{workdir / 'ipi_run'}",
)
np.save("{forces_file}", forces)
"""
    result = subprocess.run(
        [sys.executable, "-c", worker_code],
        capture_output=False,   # let i-PI output print normally
    )
    if result.returncode != 0:
        raise RuntimeError(f"i-PI subprocess failed with exit code {result.returncode}")

    forces = np.load(str(forces_file))

    # ── Step 3: inject forces back and compute bands ───────────────────────────
    # Bypass run_forces() and set the forces array directly
    ensemble._forces = forces
    n_disp, n_ens, n_at, _ = forces.shape
    print(f"  {n_ens} ensemble members | {n_disp} displacements | {n_at} atoms/supercell")

    ensemble.compute_bands(npoints=npoints)
    return ensemble


# ─────────────────────────────────────────────────────────────────────────────
# 1.  Al (FCC) — constrained relaxation + deterministic phonons
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("1. Al (FCC) — constrained relaxation + deterministic phonons")
print("=" * 60)

atoms_al = bulk("Al", "fcc", a=4.05)
atoms_al.set_constraint(FixSymmetry(atoms_al))
atoms_al.calc = calc

opt = LBFGS(FrechetCellFilter(atoms_al, mask=[True]*3 + [False]*3), logfile=None)
opt.run(fmax=FMAX, steps=STEPS)
report_symmetry(atoms_al, "Al FCC (relaxed, constrained)")
atoms_al.set_constraint(None)

phonon_al = compute_phonons_ase(atoms_al, SUPERCELL_AL, label="Al FCC")

fig, ax = plt.subplots(figsize=(8, 5))
plot_bands([phonon_al], phonon_al, mode="mean", unit="THz", ax=ax,
           color="tab:blue", mean_linewidth=1.5)
ax.set_title("Al (FCC) — PET-MAD XS, deterministic")
plt.tight_layout()
plt.savefig("al_phonons_deterministic.png", dpi=150)
plt.close(fig)
print("  Saved → al_phonons_deterministic.png")


# ─────────────────────────────────────────────────────────────────────────────
# 2.  Al (FCC) — UQ ensemble phonons
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("2. Al (FCC) — UQ ensemble phonons")
print("=" * 60)

ensemble_al = compute_phonons_uq(atoms_al, SUPERCELL_AL, label="Al FCC UQ")

fig, ax = plt.subplots(figsize=(8, 5))
ensemble_al.plot(ax=ax, mode="ensemble", unit="THz",
                 color="tab:blue", std_alpha=0.25)
ax.set_title("Al (FCC) — PET-MAD XS ensemble (UQ)")
plt.tight_layout()
plt.savefig("al_phonons_uq_ensemble.png", dpi=150)
plt.close(fig)
print("  Saved → al_phonons_uq_ensemble.png")

fig, ax = plt.subplots(figsize=(8, 5))
ensemble_al.plot(ax=ax, mode="mean+std", unit="THz",
                 color="tab:blue", std_alpha=0.35)
ax.set_title("Al (FCC) — PET-MAD XS mean ± std (UQ)")
plt.tight_layout()
plt.savefig("al_phonons_uq_mean_std.png", dpi=150)
plt.close(fig)
print("  Saved → al_phonons_uq_mean_std.png")


# ─────────────────────────────────────────────────────────────────────────────
# 3.  BaTiO3 cubic — constrained relaxation + UQ (imaginary modes expected)
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("3. BaTiO3 cubic — constrained relaxation + UQ phonons")
print("=" * 60)

a_bto = 4.00
bto_cubic = Atoms(
    symbols="BaTiO3",
    scaled_positions=[[0,0,0],[.5,.5,.5],[.5,.5,0],[.5,0,.5],[0,.5,.5]],
    cell=[a_bto, a_bto, a_bto],
    pbc=True,
)
bto_cubic.set_constraint(FixSymmetry(bto_cubic))
bto_cubic.calc = calc

opt = LBFGS(FrechetCellFilter(bto_cubic, mask=[True]*3 + [False]*3), logfile=None)
opt.run(fmax=FMAX, steps=STEPS)
report_symmetry(bto_cubic, "BTO cubic (relaxed, constrained)")
bto_cubic.set_constraint(None)

ensemble_bto = compute_phonons_uq(bto_cubic, SUPERCELL_BTO, label="BTO cubic UQ")

fig, ax = plt.subplots(figsize=(8, 5))
ensemble_bto.plot(ax=ax, mode="ensemble", unit="THz",
                  color="tab:red", std_alpha=0.2)
ax.axhline(0, color="k", lw=0.8, ls="--")
ax.set_ylim(-12, 26)
ax.set_title(r"BaTiO$_3$ cubic $Pm\bar{3}m$ — PET-MAD XS ensemble (UQ)")
plt.tight_layout()
plt.savefig("bto_cubic_phonons_uq_ensemble.png", dpi=150)
plt.close(fig)
print("  Saved → bto_cubic_phonons_uq_ensemble.png")

fig, ax = plt.subplots(figsize=(8, 5))
ensemble_bto.plot(ax=ax, mode="mean+std", unit="THz",
                  color="tab:red", std_alpha=0.35)
ax.axhline(0, color="k", lw=0.8, ls="--")
ax.set_ylim(-12, 26)
ax.set_title(r"BaTiO$_3$ cubic $Pm\bar{3}m$ — PET-MAD XS mean ± std (UQ)")
plt.tight_layout()
plt.savefig("bto_cubic_phonons_uq_mean_std.png", dpi=150)
plt.close(fig)
print("  Saved → bto_cubic_phonons_uq_mean_std.png")

print("\n" + "=" * 60)
print("All done. Saved plots:")
for f in ["al_phonons_deterministic.png",
          "al_phonons_uq_ensemble.png",
          "al_phonons_uq_mean_std.png",
          "bto_cubic_phonons_uq_ensemble.png",
          "bto_cubic_phonons_uq_mean_std.png"]:
    print(f"  {f}")
print("=" * 60)
