# Whitespace Data: Abandoned Alexandria Systems

Per-system energetic and vibrational data for the full-scale remediation study
of 6,574 abandoned systems from the Alexandria database (6,345 quaternary +
229 quinary), plus 200 elpasolite compositions evaluated with two
structure-generation strategies (forced-cubic and transmuted).

## Files

| File | Description |
|------|-------------|
| `quaternary_ga_summary.csv` | Per-system GA remediation summary for 6,345 quaternary systems. Columns: material name, Alexandria ID, space group, initial Ehull, initial/final energies, phonon metrics, run status, GA iterations, runtime. |
| `quintenary_ga_summary.csv` | Per-system GA remediation summary for 229 quinary systems. Same columns as quaternary. |
| `ehull_approx_data.csv` | Approximate post-remediation Ehull for all 6,574 systems, computed as `Ehull_init + dE_MLIP`. Columns: material, dataset, ehull_init, dE_mlip, ehull_post_approx, run_status. |
| `elpasolite_sample_200.json` | The 200-composition elpasolite sample (seed=42) with composition, alkali, metal, and halogen labels. |
| `elpasolite_whitespaces_stats.json` | Aggregate statistics for the elpasolite whitespace runs (forced-cubic and transmuted). |
| `vibroml_database.db.gz` | Gzip-compressed SQLite database (19 MB, ~190 MB uncompressed) containing the full structured data: materials, GA structures, phonon data, phonon modes, and the DFT-verification table. Decompress with `gunzip vibroml_database.db.gz`, or let `setup_and_launch.py` handle it automatically. |
| `steric_clash_flags.csv` | Per-structure steric-clash screen of every GA structure in the database (see note below). Columns: ga_structure_id, material_name, struct_type, rank, energy_per_atom_eV, min_distance_angstrom, min_distance_over_covalent, closest_pair, status (`ok`/`clashed`/`missing`/`error`). |
| `corrected_clash_systems.csv` | Per-system corrections for the 134 systems whose top-ranked GA structure is flagged: reported best energy, best clash-free energy, energy inflation, and corrected energy reduction relative to the initial structure. |
| `dft_whitespace_results.csv` | DFT single-point verification of a 200-system uniform random subsample (seed 42): static VASP energies for each system's initial structure and top-ranked GA structure at the stored geometries, alongside the MACE energies. Columns: run_dir, material, dataset, elements, e_mace_init, e_mace_winner, dE_mace, E0_initial, E0_winner, natoms_initial, natoms_winner, dE_dft_total, dE_dft_per_atom. 187 of 200 systems have converged DFT on both sides; restricting to the 174 clash-free systems, 166 (95.4%) retain a lower DFT energy than the initial structure (median DFT stabilization -96 meV/atom vs -89 meV/atom predicted by MACE; 175/187 = 93.6% counting all converged systems). |
| `dft_whitespace_failures.csv` | Per-run-directory failure modes for the 13 sampled systems whose DFT did not converge on at least one side: SCF non-convergence across multiple solver choices (predominantly magnetic metallic systems), oversized supercells exceeding the 48-hour wall limit, or never-launched directories. |
| `dft_whitespace_forces.csv.gz` | Per-atom DFT force vectors (eV/Angstrom) on the stored geometries for every converged verification run: material, side (initial/winner), atom_index, element, fx, fy, fz. 391 converged sides, 11,138 atom rows. |
| `dft_whitespace_force_stats.csv` | Per-structure force magnitudes derived from the same runs: n_atoms, E0 (where recoverable), max |F|, RMS |F|. The large force outliers are exactly the steric-clash-flagged structures. |
| `analysis_scripts/` | Reproduction scripts: `scan_steric_clashes.py` (writes `steric_clash_flags.csv` and `corrected_clash_systems.csv`), `prepare_dft_sample.py` + `sample_manifest.json` (the 200-system DFT sample, seed 42), `analyze_results.py` (writes `dft_whitespace_results.csv`/`dft_whitespace_failures.csv`, requiring genuine SCF convergence), `parse_dft_forces.py` (writes `dft_whitespace_forces.csv.gz` and `dft_whitespace_force_stats.csv` from the VASP outputs), and `add_dft_verification.py` (migrates the `dft_verification` table into the database). Local paths inside the scripts need adjusting. |
| `vibroml_explorer.py` | Streamlit web application for interactive exploration of the database: material summaries, structure viewers, phonon band structures, DOS, and phonon-mode animations. |
| `setup_and_launch.py` | Setup script that decompresses the metadata database, optionally downloads the Figshare phonon/structure bundle (`vibroml_explorer_bundle.tar.gz`), and launches the explorer. |

## SQLite Database Schema

The database contains four tables:

- `materials`: one row per system (6,774 distinct systems: 6,574 Alexandria + 200 elpasolite compositions across forced-cubic and transmuted strategies)
- `ga_structures`: one row per GA-discovered structure (264,542 rows)
- `phonon_data`: one row per phonon evaluation (123,865 rows)
- `phonon_modes`: one row per q-point mode (305,998 rows)
- `dft_verification`: one row per system in the 200-system DFT verification
  subsample (187 converged): static VASP energies for the stored initial and
  top-ranked geometries, per-atom-derived max/RMS force magnitudes, the
  steric-clash flag of the winner, and the failure mode for non-converged
  runs

Query example:

```sql
SELECT material_name, final_stability, final_softness_fraction
FROM materials
WHERE final_stability = 'stable'
ORDER BY final_softness_fraction;
```

## Launching the Database Explorer

Dependencies: `streamlit`, `pandas`, `numpy`, `matplotlib`, `ase`, and
`py3Dmol` (for 3D structure and phonon-mode viewers).

From the repository root:

```bash
python launch_explorer.py
```

or directly:

```bash
python whitespace_data/setup_and_launch.py
```

The script decompresses `vibroml_database.db.gz` into
`whitespace_data/database/`, then looks for the Figshare bundle
`vibroml_explorer_bundle.tar.gz` (either already downloaded next to the
script, or fetched from the URL embedded in the script or passed via
`--url`). After extraction, the Streamlit explorer starts on port 8501.

The bundle contains `database/vibroml_blobs.db`, a compact SQLite store
holding every external file referenced by the database (input/optimized/GA
CIFs, run settings, band-structure and DOS arrays, and phonon-mode
trajectories). Phonon-mode trajectories are stored in a harmonic form
`r(t) = r0 + uc*cos(2*pi*t/M) + us*sin(2*pi*t/M)` that reconstructs all
frames with sub-micro-angstrom error; the explorer decodes them on the fly.
Pass `--reconstruct-files` to `setup_and_launch.py` to materialize the
original `.cif`, `.txt`, and `.traj` files on disk instead.

The explorer also runs against the metadata database alone (summary tables
and statistics) when the bundle is not present.

### DFT-verified entries

The sidebar includes a **DFT data** filter. Selecting "DFT verified" shows
the 187 systems whose initial and top-ranked GA structures were both
evaluated with static VASP single points; "DFT attempted" also includes the
13 sampled systems whose runs did not converge. Verified entries display
the MLIP and DFT energy reductions side by side, the sign-confirmation
flag, and the max per-atom residual force on each stored geometry.
Per-atom force vectors are provided in `dft_whitespace_forces.csv.gz`.
Work toward the full DFT dataset for all deposited structures is ongoing.

## Note on the steric-clash screen

A post-hoc geometric audit scanned every GA structure in the database for
unphysical short contacts, defined as any atom pair closer than half the sum
of the covalent radii (`min_distance_over_covalent < 0.5`). Approximately
8.5% of the ~155,000 structures that could be re-read from disk are flagged
(~74,000 database rows point to files no longer on disk and are marked
`missing`; the flagged count is therefore a lower bound). The artifact is
strongly enriched in light-cation compounds (predominantly H-bearing
hydrides and lithium compounds), where a small mobile cation can relax
onto a neighbouring site and produce a spurious deep MLIP minimum that
passes the energy- and volume-based termination criteria.

Consequences for downstream use:

- The top-ranked structure is flagged in 115 of the 3,225 GA-completed
  Alexandria systems (and a further 19 elpasolite systems). Energies and
  phonon metrics reported for these entries are unreliable; corrected
  energies based on the best clash-free structure per system are given in
  `corrected_clash_systems.csv` (101 Alexandria systems have a clash-free
  substitute within the same run; 14 retain their reported values).
- Corrected aggregate statistics (best clash-free substitution): mean
  stabilization -104 meV/atom (median -83), 4.2% of systems reaching the
  hull, 20.9% within 50 meV/atom. The f_soft < 2% completion rate is 96%.
- In the DFT-verification subsample, flagged systems separate into
  MLIP-only extrapolation errors (repulsive in DFT) and inherited DFT
  artifacts (the DFT well lying deeper than the MLIP prediction); the
  class breakdown is given in the manuscript SI.
- The deposited Figshare XYZ dataset is not filtered; join it against
  `steric_clash_flags.csv` to exclude flagged frames before using the
  structures as DFT/MLIP training data.
- A steric-clash termination criterion (`--steric_clash_ratio`, default
  0.5) has been added to the VibroML workflow so future runs reject such
  candidates automatically.

## Note on the assembled XYZ dataset

The full collection of 83,562 GA-discovered local-minimum structures
(including stable guidance-phonon structures for runs without a final GA
population) in extended XYZ format, together with the
`vibroml_explorer_bundle.tar.gz` database bundle used by the explorer, is
deposited on Figshare under DOI `10.6084/m9.figshare.34790283` (see Data
Availability statement in the manuscript).
