# Vendored inputs

Copies of files my scripts read that were produced OUTSIDE kingman2021/, in a directory
since retired by the consolidation. Vendored so this folder is self-contained.

| File | Read by | Original producer | Status |
|---|---|---|---|
| `null_popperm_3sp.rds` | `R/12_overlap_emmax17.R` | `module_sticklebacks_LDscnR/permutation_null_3sp.R` | producer **deleted** in commit `3cb47ba`; recoverable via `git show 3cb47ba^:module_sticklebacks_LDscnR/permutation_null_3sp.R` |

Also note `data/regions_tau0.05_lmin10_rho0.60.csv` (read by `R/08`, `R/09`, `R/16`) was
produced by `module_sticklebacks_LDscnR/manhattan_regions.R`, likewise deleted in `3cb47ba`
and recoverable the same way. The CSV itself survives in `data/`.

Neither file is regenerable from anything now on disk. If either needs re-deriving, restore
the producer from git first.
