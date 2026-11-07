"""Decompose residual DFT forces on MLIP NEB images into path-parallel and
path-perpendicular components (3N configuration space), for GMR item F.

Per image i the path tangent is t_i = (R_{i+1} - R_{i-1}) / |...| (endpoints:
one-sided). The 3N force vector is split as F_par = (F . t_hat) t_hat,
F_perp = F - F_par. Endpoint images have no projection applied (the full
force should vanish there).

Outputs per-image CSV and per-path summary with F_tot_max, F_perp_max, and
the endpoint residual forces for Table S12.
"""
import os
import re
import numpy as np
from ase.io import read

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "dft_inputs", "neb_images_to_gs")
OUT = os.path.dirname(os.path.abspath(__file__))


def parse_outcar_forces(path):
    """Return (N,3) force array from the last TOTAL-FORCE block of OUTCAR."""
    with open(path, "r", errors="replace") as f:
        lines = f.readlines()
    blocks = [i for i, l in enumerate(lines) if "TOTAL-FORCE" in l]
    if not blocks:
        return None
    i = blocks[-1]
    forces = []
    j = i + 2  # skip header and dashed line
    while j < len(lines) and "----- " not in lines[j] and "---" not in lines[j]:
        parts = lines[j].split()
        if len(parts) >= 6:
            forces.append([float(parts[3]), float(parts[4]), float(parts[5])])
        j += 1
    return np.array(forces)


def image_dirs(path_dir):
    imgs = sorted(d for d in os.listdir(path_dir)
                  if re.search(r"image_\d+$", d))
    imgs.sort(key=lambda d: int(d.rsplit("_", 1)[1]))
    return [os.path.join(path_dir, d) for d in imgs]


def tangent_vector(pos_im1, pos_ip1):
    t = (pos_ip1 - pos_im1).ravel()
    n = np.linalg.norm(t)
    return t / n if n > 0 else None


rows = []
summary = []
for path_name in sorted(os.listdir(ROOT)):
    pdir = os.path.join(ROOT, path_name)
    if not os.path.isdir(pdir):
        continue
    dirs = image_dirs(pdir)
    if not dirs:
        continue
    pos = [read(os.path.join(d, "POSCAR")).get_positions() for d in dirs]
    frc = []
    for d in dirs:
        oc = os.path.join(d, "OUTCAR")
        f = parse_outcar_forces(oc) if os.path.exists(oc) else None
        frc.append(f)
    n_img = len(dirs)

    for i, (d, F) in enumerate(zip(dirs, frc)):
        if F is None:
            rows.append([path_name, i, "missing"] + [""] * 6)
            continue
        Fv = F.ravel()
        is_endpoint = (i == 0 or i == n_img - 1)
        if is_endpoint:
            Fpar, Fperp = 0.0, np.linalg.norm(Fv)
        else:
            t = tangent_vector(pos[i - 1], pos[i + 1])
            fpar_vec = np.dot(Fv, t) * t
            Fperp = np.linalg.norm(Fv - fpar_vec)
            Fpar = abs(np.dot(Fv, t))
        ftot_3n = np.linalg.norm(Fv)
        fa_max = np.linalg.norm(F, axis=1).max()
        rows.append([path_name, i, "ok",
                     f"{ftot_3n:.4f}", f"{Fpar:.4f}", f"{Fperp:.4f}",
                     f"{fa_max:.4f}"])

    done = [(i, r) for i, r in enumerate(rows[-n_img:])
            if r[2] == "ok"]
    if done:
        ftot_max = max(float(r[3]) for _, r in done[1:-1] or done)
        fperp_max = max(float(r[5]) for _, r in done[1:-1] or done)
        fatom_max = max(float(r[6]) for _, r in done)
        end_f = max(float(rows[-n_img][3]), float(rows[-1][3]))
        summary.append([path_name, f"{fatom_max:.3f}", f"{ftot_max:.3f}",
                        f"{fperp_max:.3f}", f"{end_f:.3f}"])

with open(os.path.join(OUT, "perp_force_per_image.csv"), "w") as f:
    f.write("path,image,status,F_tot_3N,F_par_3N,F_perp_3N,F_max_per_atom\n")
    for r in rows:
        f.write(",".join(map(str, r)) + "\n")

with open(os.path.join(OUT, "perp_force_summary.csv"), "w") as f:
    f.write("path,F_maxatom_tot,F_3N_tot_max_interior,F_perp_3N_max_interior,"
            "F_endpoint_max\n")
    for r in summary:
        f.write(",".join(map(str, r)) + "\n")

print(f"Parsed {sum(1 for r in rows if r[2]=='ok')}/{len(rows)} images across "
      f"{len(summary)} paths")
for r in summary:
    print(f"  {r[0]:28s} Ftot={r[1]:>6s}  Fpar3Nmax~{r[2]:>6s}  "
          f"Fperp3Nmax={r[3]:>6s}  Fend={r[4]:>6s}")
