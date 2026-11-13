# Headline-number provenance

Every numerical claim in the manuscript, SI, and response letter can be
reproduced from files in this repository. For each headline number below we
list the source file(s) and the **exact command** that confirms it.

Base paths used below:

- `WS/`   = `whitespace_data/` in this repository
- `DB`    = `WS/vibroml_database.db` (decompress `vibroml_database.db.gz` first)
- `BASE/` = the deposited GA output tree (`candidates_lowcoverage_cif/` in the
  Figshare/explorer bundle); relative `cif_path` values in the database resolve
  against it

Python snippets assume `sqlite3`, `csv`, `os`, `numpy` from a standard install
(pymatgen/ASE only needed where a structure file is opened).

> **Recovery note (2026-11-12).** 1,232 GA-completed Alexandria run directories
> were truncated during the Sept-15 scratch re-archive (`vibroml_runs.tar.zst`
> premature EOF). They were recovered from the manneback backup in two batches
> (975 + 328 directories) and re-ingested into `ga_structures` (276,336 rows,
> including 99,597 `_conventional`/`_primitive` cell-variant rows). Coverage is
> now 3,476 / 3,479 GA-completed systems with an auditable top-ranked
> structure; all earlier "1,973 retained / 1,252 unretained" values in this
> document's history are superseded. The complete incident record — root
> cause, recovery procedure, per-number before/after table, and a re-runnable
> verification register — is in **`MISSING_RUN_RECOVERY_REPORT.md`**.

---

## 1. Campaign counts

| Claim | Value | How to confirm |
|---|---|---|
| Total abandoned Alexandria systems | **6,574** | `sqlite3 $DB "SELECT COUNT(*) FROM materials WHERE dataset IN ('quaternary','quintenary')"` — or `wc -l` on `BASE/quaternary_ga_summary.csv` + `quintenary_ga_summary.csv` minus 2 header lines |
| Met cutoff without GA (`stable_no_ga`) | **3,137** (47.7%) | `sqlite3 $DB "SELECT COUNT(*) FROM materials WHERE dataset IN ('quaternary','quintenary') AND run_status='stable_no_ga'"` |
| Entered GA | **3,437** | `sqlite3 $DB "SELECT COUNT(*) FROM materials WHERE dataset IN ('quaternary','quintenary') AND run_status IN ('ga_completed','ga_incomplete','failed_early','never_ran')"` |
| GA-completed | **3,225** | `sqlite3 $DB "SELECT COUNT(*) FROM materials WHERE dataset IN ('quaternary','quintenary') AND run_status='ga_completed'"` |
| Elpasolite GA-completed | **254** (136 forced-cubic + 118 transmuted) | `sqlite3 $DB "SELECT dataset, COUNT(*) FROM materials WHERE dataset LIKE 'elpasolite%' AND run_status='ga_completed' GROUP BY dataset"` |
| Unique GA polymorphs | **156,671** (mean 49/system, max 190) | `sqlite3 $DB "SELECT SUM(ga_unique_structures), AVG(ga_unique_structures), MAX(ga_unique_structures) FROM materials WHERE dataset IN ('quaternary','quintenary') AND run_status='ga_completed'"` |
| Total GPU cost | **~4,080 GPU-h** (mean 1.3 h/system) | `sqlite3 $DB "SELECT SUM(runtime_hours), AVG(runtime_hours) FROM materials WHERE dataset IN ('quaternary','quintenary') AND run_status='ga_completed'"` |
| Auditable Alexandria runs | **3,222** (3,045 quat + 177 quin) | GA-completed Alexandria materials whose rank-1 `top` `cif_path` resolves under `BASE/` — see snippet in §2 below |
| Non-auditable GA-completed runs | **3** | `Ba2Fe3MoH12` (no top-ranked rows), `LiCdSbRu`, `LiNbSbRu` (all retained candidates clash-flagged; excluded in the audit table) |
| DFT sampling pool | **1,778** | runs whose initial-relaxed **and** top-ranked CIFs existed at DFT dispatch time; snapshot recorded in `WS/analysis_scripts/sample_manifest.json` provenance (1,352 quat + 166 quin + 260 elpasolite). The subsample was drawn before the manneback recovery. |

## 2. Steric-clash audit (manuscript: "267 Alexandria + 19 elpasolite flagged")

Per-structure flags live in `WS/steric_clash_flags.csv` (276,336 rows; one row
per `ga_structures.id`, statuses `clashed`/`ok`/`missing`/`error`). Clash
criterion (see `analysis_scripts/scan_steric_clashes.py`): flag iff
`d_min/(r_cov_i + r_cov_j) < 0.5` **or** `d_min < 1.0 Å`, periodic images via
ASE `get_all_distances(mic=True)`.

```bash
# Rebuild the flags CSV (full-DB variant used for the campaign DB):
python candidates_lowcoverage_cif/scan_steric_clashes_fulldb.py --workers 16

# Count flagged GA-completed runs per dataset (retained top-1 structures):
python - <<'EOF'
import csv, sqlite3
flags = {int(r['ga_structure_id']): r['status'] for r in
         csv.DictReader(open('whitespace_data/steric_clash_flags.csv'))}
db = sqlite3.connect('whitespace_data/vibroml_database.db')
q = """SELECT m.dataset, g.id FROM materials m JOIN ga_structures g
       ON g.material_id=m.id AND g.struct_type='top' AND g.rank=1
       WHERE m.run_status='ga_completed'"""
from collections import Counter
print(Counter(ds for ds, gid in db.execute(q) if flags.get(gid)=='clashed'))
# -> {'quaternary': 264, 'quintenary': 3,
#     'elpasolite_forced_cubic': 12, 'elpasolite_transmuted': 7}
EOF
```

| Claim | Value |
|---|---|
| Alexandria flagged | **267 / 3,222 = 8.3%** (264 quat + 3 quin) |
| Elpasolite flagged | **19 / 254 = 7.5%** (12 + 7) |
| Substituted with best clash-free candidate | **243** (224 Alex + 19 elp) |
| Excluded (no clash-free candidate retained) | **43** |

`WS/corrected_clash_systems.csv` (286 rows) lists every flagged run:
`reported_best_eV`, `best_clean_eV`, `n_clashed_structures`,
`energy_inflation_eV` (= best_clean − reported_best), and corrected reductions.
Regenerate with `candidates_lowcoverage_cif/generate_corrected_clash_systems.py`.

## 3. Audited energy & hull statistics (Figure 4A; N = 3,126)

Definition: audited (3,222) − all-clash (43) = 3,179 clash-free runs; restrict
to finite `e_above_hull` + `ΔE` + `final_softness_fraction` → **N = 3,126**.
For the 224 substituted runs ΔE is `corrected_reduction_eV`.

```bash
python - <<'EOF'
import csv, sqlite3, numpy as np
db = sqlite3.connect('whitespace_data/vibroml_database.db')
corr, excl = {}, set()
for r in csv.DictReader(open('whitespace_data/corrected_clash_systems.csv')):
    (corr.__setitem__(r['material_name'], float(r['corrected_reduction_eV']))
     if r['best_clean_eV'] else excl.add(r['material_name']))
# retained = top1 cif on disk; substitute corrected ΔE; drop excl + NaNs
# (full logic in generate_main_ga_remediation_figure.py — it prints all stats)
EOF
```

`generate_main_ga_remediation_figure.py` prints the full block:

| Claim | Value |
|---|---|
| mean ΔE | **−92 meV/atom** (−91.5) |
| median ΔE | **−80 meV/atom** (−80.3) |
| >100 meV/atom improvement | **1,160 (37.1%)** |
| >200 meV/atom improvement | **208 (6.7%)** |
| mean E_hull | **0.185 → 0.093 eV/atom** |
| median E_hull | **180 → 89 meV/atom** (180.5 → 89.1) |
| on hull (ΔE_hull ≤ 0) | **1.3%** (42 systems) |
| within 50 meV/atom of hull | **18.5%** (578) |

The unaudited values across all 3,160 GA-completed Alexandria systems with
energy records (kept in SI S15 for provenance): mean −112 (median −86)
meV/atom; 41% improving by >100 meV/atom; median hull 180→84 meV/atom; mean
hull →0.073 eV/atom; 6.5% at or below hull; 23.3% within 50 meV/atom.

The non-Γ soft-mode check (manuscript/response "64% zero imaginary modes at
non-Γ points"): `rescan_nongamma_modes.py` loads the top-1
`final_phonon_analysis_*/band_structure_energies_*.txt` of every quaternary
GA-completed run and counts runs where all negative frequencies sit at Γ rows
(acoustic |ω| < 0.05 THz) — **1,760 / 2,754 = 63.9%**.

## 4. f_soft statistics (Figure 4B — all 3,225 GA-completed)

```bash
sqlite3 $DB "SELECT
  SUM(initial_softness_fraction<0.02), SUM(final_softness_fraction<0.02),
  SUM(final_softness_fraction=0), COUNT(*) FROM materials
  WHERE dataset IN ('quaternary','quintenary') AND run_status='ga_completed'"
# -> 999 | 3096 | 56 | 3225
```

| Claim | Value |
|---|---|
| below 2% before GA | **999 / 3,225 = 31.0%** |
| below 2% after GA | **3,096 / 3,225 = 96.0%** (65 NULL finals count as not-passing) |
| f_soft = 0 exactly | **56** (54 × 4-atom + 2 × 6-atom primitive cells) |
| stable-no-GA systems at f_soft = 0 | **1 / 3,137** |

## 5. DFT verification subsample

```bash
sqlite3 $DB "SELECT COUNT(*), SUM(converged), SUM(clash_flag), SUM(confirmed)
             FROM dft_verification"   # -> 200 | 187 | 18 | 175
```

| Claim | Value | Source |
|---|---|---|
| subsample size / converged | 200 / 187 | `dft_verification` table; per-system: `analysis_scripts/dft_whitespace_results.csv` |
| clash-free confirmed | **166/174 = 95.4%** | `dft_whitespace_results.csv` (`converged=1 AND clash_flag=0` → `confirmed`); all-187 rate 175/187 = 93.6% |
| DFT vs MLIP median stabilization | **−96 vs −89 meV/atom** | `de_dft_per_atom`, `de_mlip_per_atom` columns |
| median absolute deviation | **14 meV/atom** | |`de_dft_per_atom − de_mlip_per_atom`| |
| largest sign reversal | CaMgGaSn +19.5 eV/atom | same file |

## 6. Other sections

| Claim | Source file(s) |
|---|---|
| Elpasolite: 254 runs, mean ΔE −240 meV/atom | `DB`: `SELECT AVG(energy_reduction_vs_initial_eV_per_atom) FROM materials WHERE dataset LIKE 'elpasolite%' AND run_status='ga_completed'` |
| B2 phonon metrics (7 deposited engine×compound pairs) | `b2_intermetallics/` (see its README; Table S12 rows) |
| NEB perpendicular/parallel force decomposition (7 paths) | `neb_results/dft_perp_force_analysis/perp_force_summary.csv` + `decompose_forces.py` |
| HfO₂ MD RMSD (recomputed) | `HFO2_RMSD_ARTIFACT.md` — full forensic record |
| HfO₂ Qi–Rabe crossmatch: 5/12 recovered incl. dm–dm (Table S18) | `hfo2_literature_validation/` — `compare_hfo2_phases.py` matches VibroML CIFs against the published reference POSCARs in `qi_rabe_reference/` (`3a`–`3l.vasp` from github.com/yuboqiuab/unstableflatband, commit `a438e4e`). `unique_1_iter2_sample63` == `3l` (dm–dm, RMS 0.019 Å); `top_1_iter3_sample26` == `3c` (m, RMS 0.062 Å) |
| **186,653-structure Figshare XYZ dataset** | `candidates_lowcoverage_cif/assemble_ga_structures_xyz_v2.py` rebuilds `ga_local_minima_structures_v2.xyz` from every `*/final_structures/` dir (top+unique, excl. `_conventional`/`_primitive` variants; per-dataset: 151,464 quat + 5,097 quin + 18,307 forced-cubic + 11,785 transmuted) — supersedes the 83,562-structure v1 deposit assembled before the manneback recovery |
| Table S8 five-engine phonon metrics (43 structure blocks; uniform 6^3 mesh, 2x2x2 supercell, 0.05 A displacement) | `pet_phonons/ip_variants_{mace,macemp,esen,uma}_merged.json` (deterministic engines); PET-MAD column = ensemble mean ± 1sigma from `pet_phonons/uq_validation_summary.csv` (per-run ensemble npz on Figshare: `pet_phonons_uq_validation_runs.tar.gz`). `pet_phonons/structure_crosswalk.csv` maps every table row to run directory and source CIF. Two documented patches: LiF Pm-3m uses the `convergence_sweep.json` 8^3 entry; LiF Fm-3m uses the true-rocksalt eval (see `QMESH_PROTOCOL_RESOLUTION.md`) |
| GA vs. RAND budget figure (Fig. S8) | `remediation_results/plot_ga_vs_rand_budget.py`; GA `overall_*_summary*.txt` in `LiFsimplecubic_GA5it_phonon_output_20250831-152432`, `CsPbI3_GA5it_phonon_output_20251003-002011`, `HfO2_Fm3m_GA_phonon_output_20250826-161333`, `Sn2Bi2O7_MACE_GA_phonon_output_20260107-205131`; RAND per-sample `relaxation_summary.txt` / run summaries in the corresponding `*_OPT_RANDOM_*` directories (LiF excludes 14 clash-artifact rows below −5 eV/atom) |

*Register maintained alongside manuscript revisions; last update covers the
manneback recovery + full-database re-audit (2026-11-12).*
