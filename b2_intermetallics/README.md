# B2 Intermetallic Phonon Metrics

Phonon evaluation outputs for the three B2 (CsCl-type) intermetallic benchmark
systems: CuN, ZnCu, and ZrCu.

Each subdirectory contains the `phonon_only` evaluation of the GA-discovered
ground-state primitive cell with the indicated MLIP engine (MACE-OMAT, eSEN,
or UMA). All seven ground states have $f_\mathrm{soft} = 0$ and
$I_\mathrm{curv} = 0$, confirming dynamical stability.

## Contents

| System | Engine | Directory |
|--------|--------|-----------|
| CuN | MACE | `b2_metrics_CuN_MACE_*` |
| CuN | eSEN | `b2_metrics_CuN_ESEN_*` |
| ZnCu | eSEN | `b2_metrics_ZnCu_ESEN_*` |
| ZnCu | UMA | `b2_metrics_ZnCu_UMA_*` |
| ZrCu | MACE | `b2_metrics_ZrCu_MACE_*` |
| ZrCu | eSEN | `b2_metrics_ZrCu_ESEN_*` |
| ZrCu | UMA | `b2_metrics_ZrCu_UMA_*` |

These data correspond to Supplementary Table S_b2_phonon_metrics in the
manuscript.
