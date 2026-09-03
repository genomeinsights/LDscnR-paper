# Superseded — retained as audit trail, do not run

These two scripts produced figures that were in `manuscript/` until 2026-09-02.
Both have known errors; corrected replacements live in `kingman2021/R/`
(owned by the Kingman/EcoPeak session), which is where current results come from.

| superseded here                | corrected replacement              | what was wrong |
|--------------------------------|------------------------------------|----------------|
| `kingman_overlap_sweep.R`      | `kingman2021/R/17_overlap_sweep.R` | left precision/fold unmasked in cells holding a single region, producing a spurious bright band in the high-tau/high-l_min corner and inverting the apparent optimum |
| `kingman_cscore_enrichment.R`  | `kingman2021/R/18_cscore_enrichment.R` | figure title claimed a rising trend in C that the bins do not show (the data show a C=0 -> C>0 step) |

Kept so that a referee asking why a figure changed can be answered from the record.
See `kingman2021/doc/00_PROVENANCE.md` and §10 of the Kingman report.

The rest of this module's scripts were deleted on 2026-09-02; they remain in git
history. `results/` and `figures/` are preserved here because `kingman2021/R/`
and `LD-pruning-paper/R_3sp/` read from `results/` directly.

---

## The null producers, restored 2026-09-03

The module's scripts were deleted in 3cb47ba. Eight that WRITE into `results/`
are restored here, because `results/` is consumed live by `kingman2021/R` and
`LD-pruning-paper/R_3sp` and its `.rds` had no producer in the working tree.

Two audits (mine and ldscnr-d4's) independently concluded that
`null_popperm_3sp.rds` had **no producer anywhere**. Both were wrong, for the
same reason: every one of these scripts builds its output name with `sprintf()`
from an argument, so a grep for a producer of a literal filename returns nothing.
`permutation_null_3sp.R:35` is
`sprintf(".../results/null_%s_3sp.rds", TAG)` — invisible to that search.

**Output to invocation.** All take `commandArgs()`, run from the repo root:

| output | script | argument |
|---|---|---|
| `null_regionperm_3sp.rds` | `permutation_null_3sp.R` | *(default)* `pop_locality` |
| `null_popperm_3sp.rds` | `permutation_null_3sp.R` | **`none`** — not the default |
| `null_latent_3sp.rds` | `emmax_latent_null_3sp.R` | basis `3sp_latent_basis.rds` |
| `null_latent_thin250_3sp.rds` | `emmax_latent_null_3sp.R` | basis `3sp_latent_basis_thin250.rds` |
| `null_spatial_3sp.rds` | `spatial_null_3sp.R` | fixed output path |
| `null_lfmm_perm_3sp.rds` | `lfmm_permutation_null_3sp.R` | `B` (default 100) |
| `lfmm_b1_fournulls_Cs.rds` | `lfmm_b1_fournulls_3sp.R` | fixed output path |
| `null_cache_3sp.rds` | `manhattan_regions.R` / `run_3sp_LDscnR.R` | cache path |

Note the trap in row 2: the *default* run of `permutation_null_3sp.R` produces
`regionperm`, so regenerating `null_popperm_3sp.rds` requires passing `none`
explicitly. Nothing recorded that until now, and it is the single fact that made
the file look unregenerable.

The `null_emmax_*` outputs are not covered above; their arguments were not
recovered, which is why `results/` is now tracked in full rather than relied on
being reproducible. The basis files these scripts read moved to `3sp_data/` in
3cb47ba, so paths need checking before any re-run.

**These are superseded.** They are here as provenance and as a regeneration path
of last resort, not as current code.
