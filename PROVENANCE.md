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
| Audited GA-completed runs | **3,225** (all of them) | every `ga_completed` Alexandria run is covered by the geometric-integrity audit below; `Ba2Fe3MoH12` was recovered into the audited set after the Part-C database ingest |
| DFT verification subsample | **200** | uniform random draw (seed 42); sampled systems and sampling pool recorded in `WS/analysis_scripts/sample_manifest.json` |

## 2. Geometric-integrity audit

(manuscript: "269 Alexandria + 19 elpasolite flagged for steric clashes;
64 + 23 for loss of structural coherence; 29 carry both flags")

Two independent geometric checks are applied to the top-ranked structure of
every audited run:

- **Steric clash** — unphysically short interatomic contacts. Per-structure
  flags live in `WS/steric_clash_flags.csv` (276,448 rows; one row per
  `ga_structures.id`, statuses `clashed`/`ok`/`missing`/`error`). Criterion
  (see `analysis_scripts/scan_steric_clashes.py`): flag iff
  `d_min/(r_cov_i + r_cov_j) < 0.5` **or** `d_min < 1.0 Å`, periodic images
  via ASE `get_all_distances(mic=True)`.
- **Structural coherence** — unbound atoms, non-percolating bond-network
  fragments, or collapsed voids. Flags live in
  `WS/structural_coherence_flags.csv` (3,154 evaluated winner structures;
  statuses `ok`/`incoherent`). Criterion:
  `analysis_scripts/scan_structural_coherence.py`.

```bash
# Rebuild the clash flags CSV (full-DB variant used for the campaign DB):
python candidates_lowcoverage_cif/scan_steric_clashes_fulldb.py --workers 16

# Count flagged runs per dataset from the substitution registries:
python - <<'EOF'
import csv
from collections import Counter
clash = list(csv.DictReader(open('whitespace_data/corrected_clash_systems.csv')))
coh   = list(csv.DictReader(open('whitespace_data/coherence_substitutions.csv')))
print('clash flagged :', Counter(r['dataset'] for r in clash))
print('coher. flagged:', Counter(r['dataset'] for r in coh))
cn = {r['material_name'] for r in clash}
kn = {r['material_name'] for r in coh}
print('both flags    :', len(cn & kn))
EOF
```

| Claim | Value |
|---|---|
| Alexandria clash-flagged | **269 / 3,225 = 8.3%** (266 quat + 3 quin) |
| Elpasolite clash-flagged | **19 / 254 = 7.5%** (12 + 7) |
| Alexandria coherence-flagged | **64** (60 quat + 4 quin) |
| Elpasolite coherence-flagged | **23** (15 + 8) |
| Both flags | **29** |
| Substituted, clash audit | **245** (226 Alex + 19 elp; best clash-free candidate from the same run) |
| Substituted, coherence audit | **84** (61 Alex + 23 elp) |
| Excluded, no integrity-passing substitute | **43** (all Alexandria; includes `LiCdSbRu`, `LiNbSbRu`, `Ca2ZnH8S`, `Cs2Li3ZrH12`, `NaCu4(SbH3)2`) |

`WS/corrected_clash_systems.csv` (288 rows) lists every clash-flagged run:
`reported_best_eV`, `best_clean_eV`, `n_clashed_structures`,
`energy_inflation_eV` (= best_clean − reported_best), and corrected
reductions. `WS/coherence_substitutions.csv` (87 rows) lists every
coherence-flagged winner with its substitute `ga_structure_id` and path.
Regenerate with `candidates_lowcoverage_cif/generate_corrected_clash_systems.py`.

## 3. Audited energy & hull statistics (Figure 4A; N = 3,127)

Definition: 3,225 GA-completed − 43 without an integrity-passing substitute
= 3,182 integrity-passing runs; restrict to finite `e_above_hull` + `ΔE` +
`final_softness_fraction` → **N = 3,127**. For substituted runs ΔE uses the
corrected energy of the substitute from the same run.

`generate_main_ga_remediation_figure.py` prints the full block:

| Claim | Value |
|---|---|
| mean ΔE | **−91 meV/atom** |
| median ΔE | **−80 meV/atom** |
| >100 meV/atom improvement | **1,156 (37%)** |
| >200 meV/atom improvement | **207** |
| mean E_hull | **0.185 → 0.093 eV/atom** |
| median E_hull | **180 → 89 meV/atom** |
| within 50 meV/atom of hull | **18.4%** (~575) |
| at or below hull (post-GA `e_above_hull` ≤ 0) | **39 (1.2%)** |

The uncorrected reference values across all GA-completed Alexandria systems
with energy records (N = 3,163, kept in SI S15 for provenance): mean **−112**,
median **−86** meV/atom; **206 systems (6.5%)** at or below the hull — most of
these apparent crossings are flagged structures and disappear once
integrity-passing substitutes are used.

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
# -> 999 | 3097 | 56 | 3225
```

| Claim | Value |
|---|---|
| below 2% before GA | **999 / 3,225 = 31.0%** |
| below 2% after GA | **3,097 / 3,225 = 96.0%** (64 NULL finals count as not-passing) |
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
| audit-flagged converged systems | **16** (13 clash + 3 coherence-only) | intersect `dft_whitespace_results.csv` materials with `corrected_clash_systems.csv` + `coherence_substitutions.csv`; the flagged winner is what DFT evaluated (verifiable via `reported_energy_eV_per_atom` = `e_mace_winner`) |
| integrity-passing confirmed | **163/171 = 95.3%** | 187 converged − 16 flagged = 171; 175 confirmed − 12 flagged confirmations = 163; the 4 unconfirmed flagged systems are the sign reversals |
| 8 unconfirmed integrity-passing | marginal MLIP reductions <25 meV/atom, DFT residual <10 meV/atom | same file |
| DFT vs MLIP median stabilization | **−96 vs −89 meV/atom** | `de_dft_per_atom`, `de_mlip_per_atom` columns |
| median absolute deviation | **14 meV/atom** | |`de_dft_per_atom − de_mlip_per_atom`| |
| largest sign reversal | CaMgGaSn +19.5 eV/atom | same file |

## 6. Other sections

| Claim | Source file(s) |
|---|---|
| Elpasolite: 254 runs, mean ΔE −240 meV/atom | `DB`: `SELECT AVG(energy_reduction_vs_initial_eV_per_atom) FROM materials WHERE dataset LIKE 'elpasolite%' AND run_status='ga_completed'` |
| B2 phonon metrics (9 deposited engine×compound pairs) | `b2_intermetallics/` (see its README; Table S12 rows) |
| NEB perpendicular/parallel force decomposition (9 paths) | `neb_results/dft_perp_force_analysis/perp_force_summary.csv` + `decompose_forces.py` |
| HfO₂ Qi–Rabe crossmatch: 5/12 recovered incl. dm–dm (Table S18) | `hfo2_literature_validation/` — `compare_hfo2_phases.py` matches VibroML CIFs against the published reference POSCARs in `qi_rabe_reference/` (`3a`–`3l.vasp` from github.com/yuboqiuab/unstableflatband, commit `a438e4e`). `unique_1_iter2_sample63` == `3l` (dm–dm, RMS 0.019 Å); `top_1_iter3_sample26` == `3c` (m, RMS 0.062 Å) |
| **186,653-structure Figshare XYZ dataset** | `candidates_lowcoverage_cif/assemble_ga_structures_xyz_v2.py` rebuilds `ga_local_minima_structures_v2.xyz` from every `*/final_structures/` dir (top+unique, excl. `_conventional`/`_primitive` variants; per-dataset: 151,464 quat + 5,097 quin + 18,307 forced-cubic + 11,785 transmuted) |
| Table S8 five-engine phonon metrics (43 structure blocks; uniform 6^3 mesh, 2x2x2 supercell, 0.05 A displacement) | `pet_phonons/ip_variants_{mace,macemp,esen,uma}_merged.json` (deterministic engines); PET-MAD column = ensemble mean ± 1sigma from `pet_phonons/uq_validation_summary.csv` (per-run ensemble npz on Figshare: `pet_phonons_uq_validation_runs.tar.gz`). `pet_phonons/structure_crosswalk.csv` maps every table row to run directory and source CIF. Two documented patches: LiF Pm-3m uses the `convergence_sweep.json` 8^3 entry; LiF Fm-3m uses a dedicated evaluation on the true rocksalt structure |
| GA vs. RAND budget figure (Fig. S8) | `remediation_results/plot_ga_vs_rand_budget.py`; GA `overall_*_summary*.txt` in `LiFsimplecubic_GA5it_phonon_output_20250831-152432`, `CsPbI3_GA5it_phonon_output_20251003-002011`, `HfO2_Fm3m_GA_phonon_output_20250826-161333`, `Sn2Bi2O7_MACE_GA_phonon_output_20260107-205131`; RAND per-sample `relaxation_summary.txt` / run summaries in the corresponding `*_OPT_RANDOM_*` directories (LiF excludes 14 clash-artifact rows below −5 eV/atom) |

*Register maintained alongside manuscript revisions.*
