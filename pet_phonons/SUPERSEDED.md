# Superseded Outputs in this Directory

The following outputs were produced with the legacy `compute_omega_inst.py`
which applied a 0.05 THz all-mode cutoff (`near_zero_mask`) in addition to
the Gamma-acoustic exclusion. This contradicts the I_curv definition
(I_curv = Tr(D_-)/3N, no frequency threshold) and has been corrected.

The corrected implementation lives in:
  - `compute_ip_variants.py` (already correct: Gamma-acoustic-only exclusion)
  - `compute_omega_inst.py`   (fixed on 2026-08-30: removed near_zero_mask)

## Superseded files (do NOT use for published numbers)

| File | Issue | Replacement |
|------|-------|-------------|
| `omega_inst_full_mace.json` | Used 0.05 THz all-mode cutoff | `ip_variants_mace.json` |
| `omega_inst_full_pet_symconstrained.json` | Used 0.05 THz all-mode cutoff | `ip_variants_pet.json` |
| `omega_inst_full_pet_symconstrained.csv` | Same | (derived from ip_variants_pet.json) |
| `omega_inst_cspbi3_cubic_mace.json` | Same | `ip_variants_mace.json` (idx=44) |
| `omega_inst_cspbi3_cubic_mace.csv` | Same | (derived from ip_variants_mace.json) |
| `uq_validation_summary.csv` (i_curv columns) | Mixed: some entries used the old cutoff | `ip_variants_{mace,pet}.json` |
| `archive/test_omega_inst_*.json` | Early test outputs with the old cutoff | Not needed |

**These files are NOT deleted** — they are retained for provenance and
to allow comparison between the old and corrected methodologies.

## Correct outputs (use these)

| File | Status |
|------|--------|
| `ip_variants_mace.json` | Correct: Gamma-acoustic-only exclusion, p-sweep, uniform 6x6x6 mesh |
| `ip_variants_pet.json` | Correct: same |
| `ip_variants_freqs_mace/` | Raw frequencies (npz) for MACE |
| `ip_variants_freqs_pet/` | Raw frequencies (npz) for PET-MAD |
