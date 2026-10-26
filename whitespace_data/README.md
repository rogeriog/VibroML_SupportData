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
| `vibroml_database.db.gz` | Gzip-compressed SQLite database (22 MB, 218 MB uncompressed) containing the full structured data: materials, GA structures, phonon data, and phonon modes. Decompress with `gunzip vibroml_database.db.gz`, or let `setup_and_launch.py` handle it automatically. |
| `vibroml_explorer.py` | Streamlit web application for interactive exploration of the database: material summaries, structure viewers, phonon band structures, DOS, and phonon-mode animations. |
| `setup_and_launch.py` | Setup script that decompresses the metadata database, optionally downloads the Figshare phonon/structure bundle (`vibroml_explorer_bundle.tar.gz`), and launches the explorer. |

## SQLite Database Schema

The database contains four tables:

- `materials`: one row per composition (6,990 total including elpasolite)
- `ga_structures`: one row per GA-discovered structure (264,542 rows)
- `phonon_data`: one row per phonon evaluation (123,865 rows)
- `phonon_modes`: one row per q-point mode (305,998 rows)

Query example:

```sql
SELECT material_name, final_stability, final_softness_fraction
FROM materials
WHERE final_stability = 'stable'
ORDER BY final_softness_fraction;
```

## Launching the Database Explorer

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

## Note on the assembled XYZ dataset

The full collection of 83,562 GA-discovered local-minimum structures
(including stable guidance-phonon structures for runs without a final GA
population) in extended XYZ format is deposited separately on Figshare (see
Data Availability statement in the manuscript).
