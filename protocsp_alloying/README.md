# ProtoCSP solid-solution alloying data

Source data and figure pipeline for the targeted compositional-engineering
case studies in the VibroML manuscript: the
Cs$_2$(K$_{1-x}$Na$_x$)In(I$_{1-y}$Br$_y$)$_6$ and
(K$_{1-x}$Ba$_x$)(Ta$_{1-x}$Zr$_x$)(Se$_{1-y}$O$_y$)$_3$ solid-solution
maps (main-text alloying figures and their Supplementary counterparts).

## Contents

- `cs2kini6_phase_diagram_data.csv`, `ktase3_phase_diagram_data.csv` —
  aggregated formation-energy grids. One row per
  (composition `x`, halide/anion fraction `y`, polymorph) with columns
  `base_compound, formula, x, y, polymorph,
  formation_energy_per_atom_eV, total_energy_per_atom_eV`.
  `polymorph = base` is the target 3D framework; `top*` and `2d_top*`
  are the competing lower-symmetry polymorphs.
- `protocsp_runs/` — raw ProtoCSP output tree. Each
  `*/x_*/y_*/<formula>_mlip_energies.json` holds the ranked,
  MLIP-relaxed candidate list for that composition-polymorph pair
  (space group, volume, formation and total energy per atom).
- `combined_phonon_data.csv`, `cs2kini6_phonon_data.csv`,
  `ktase3_phonon_data.csv` — per-composition dynamical-stability
  metrics (soft-mode fraction, minimum phonon frequency) behind the
  phonon heatmap figures.
- `extract_cs2kini6_phase_diagram_data.py`,
  `extract_ktase3_phase_diagram_data.py` — scripts that aggregate the
  raw `*_mlip_energies.json` tree into the CSVs above.
- `plot_cs2kini6_phase_diagram_complete.py`,
  `plot_ktase3_phase_diagram_complete.py` — the exact scripts used for
  the published combined panels. Each takes the corresponding CSV in
  the working directory and writes `*.jpg` and `*.svg`.
- `plot_ktase3_phase_diagram_3d_2d.py`, `plot_phase_diagram_3d.py` —
  generators of the interactive Plotly phase diagrams.
- `*_3d.html` — interactive 3D phase-diagram views.
- `Cs2GaNbF6_protoCSP.txt` — example ProtoCSP retrieval log for a
  white-space elpasolite composition.

## Note on the difference metric

The contour panels plot, per composition, the energy of the target
framework **minus the minimum over the competing polymorphs**
(3D-vs-1D for KTaSe$_3$, and lowest-2D vs lowest-1D in the second
KTaSe$_3$ panel) — i.e., a fixed-composition polymorph comparison, not
a convex hull over compositions. The plotting scripts are the
corrected versions: earlier working copies compared against the mean
of the polymorphs and carried a "Convex Hull" title label that did not
match the plotted quantity.

## Regeneration

```bash
python plot_cs2kini6_phase_diagram_complete.py   # writes cs2kini6_combined_formation_energy_3panel.{jpg,svg}
python plot_ktase3_phase_diagram_complete.py     # writes ktase3_combined_formation_energy_5panel.{jpg,svg}
```

Requires `pandas`, `matplotlib`, `scipy`, `numpy` (the VibroML
`vibroml_env` environment is sufficient).
