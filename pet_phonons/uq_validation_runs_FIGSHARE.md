# `uq_validation_runs/` — deposited on Figshare

The per-structure PET-MAD ensemble runs (143 directories, each containing the
`*_phonon_uq_data.npz` ensemble-phonon payload plus the reference deterministic
phonon outputs) are deposited as a single tarball in the project Figshare
deposit:

- File: `pet_phonons_uq_validation_runs.tar.gz` (~268 MB, 2,337 entries)
- Deposit: `10.6084/m9.figshare.34790283` (same deposit as the GA local-minima
  structure dataset and explorer bundle)

The ensemble summary statistics used for the Table S8 PET-MAD column
(`pet_omega_min`, `pet_std_at_min_thz`) are in this directory in
`uq_validation_summary.csv`. The i-PI per-step scratch (`_uqphonon*`,
`ipi_forces.committee_*`) is excluded from the tarball; every run's
`phonon_uq_data.npz` payload is retained.
