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
| `vibroml_database.db.gz` | Gzip-compressed SQLite database (18 MB, 190 MB uncompressed) containing the full structured data: materials, GA structures, phonon entries, and phonon modes. Decompress with `gunzip vibroml_database.db.gz`. |

## SQLite Database Schema

The database contains four tables:

- `materials`: one row per composition (6,990 total including elpasolite)
- `ga_structures`: one row per GA-discovered structure (229,659 rows)
- `phonon_entries`: one row per phonon evaluation (107,324 rows)
- `phonon_modes`: one row per q-point mode (235,671 rows)

Query example:

```sql
SELECT material_name, final_stability, final_softness_fraction
FROM materials
WHERE final_stability = 'stable'
ORDER BY final_softness_fraction;
```

## Note on the assembled XYZ dataset

The full collection of 83,441 GA-discovered local-minimum structures in
extended XYZ format is deposited separately on Figshare (see Data Availability
statement in the manuscript).
