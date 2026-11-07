# Polymorph Symmetry Analysis (Table S6)

Group-subgroup analysis of the Table-1 polymorphs found by each VibroML
remediation algorithm, behind SI Table S6.

For each Table-1 row the script (`subgroup_analysis.py`):

1. Locates the deposited CIF matching the space-group number and the MACE
   energy (reference + reduction).
2. Applies an abstract subgroup test via `pymatgen` `is_subgroup` (by
   space-group number).
3. Enumerates integer transformation matrices `P` (det = volume ratio),
   prefilters by Gram-matrix eigenvalue isometry, and fits with
   `StructureMatcher` to report the supercell mapping, RMSD, active
   distortion, and relation classification (proper subgroup, order-order,
   distortive, or reconstructed).

Output: `subgroup_analysis_results.json` (the values typeset in Table S6).
The $Pna2_1$ Bi2Sn2O7 entry is regenerated from its recorded displacement
recipe before the RMSD is computed, as noted in the table footnote.

Dependencies: `pymatgen`, `numpy`. Local paths inside the script need
adjusting.
