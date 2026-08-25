# Handoff — the input contract for `analyse_one_dataset.R`

For the session generating simulation data. Written 2026-08-25.

Division of labour: you produce the data, I analyse it. The interface is two
files per (cell x engine x null type). **Validate before handing over:**

```bash
Rscript module_sim_LDscnR/check_inputs.R <panel.rds> <pvals.rds>
```

Exit 0 = pass. Every check in it is something that would otherwise fail deep
inside the analysis, or not fail at all and quietly produce wrong regions.

## You may already be done

`run_sim_nulls.R` writes `cell_<id>.rds` + `pnull_<eng>_<ty>_<id>.rds`, which
carry everything except genotypes. `module_sim_LDscnR/adapt_pnull.R` converts a
pnull pair into the two files below, re-pooling `GTs` from the bundles with the
identical `R<i>_` rule your `pool_cell()` uses and verifying the result against
`context$markers`. So:

```bash
Rscript module_sim_LDscnR/run_sim_nulls.R 2 1 2
Rscript module_sim_LDscnR/adapt_pnull.R <out>/pnull_emmax_env_orth_V2_c1_env2.rds
Rscript module_sim_LDscnR/check_inputs.R <out>/panel_V2_c1_env2.rds <out>/pvals_V2_c1_env2_emmax_env_orth.rds
```

is the whole path. If that works, ignore the rest of this file — it is the spec
the adapter already satisfies.

**As of writing, no `pnull_*` or `cell_*` files exist anywhere**, so the nulls
stage has not yet been run against the regenerated bundles.

## FILE 1 — `panel_<cell>.rds`: the LD description of the study

Depends on the panel, never on the association method. One per cell, shared by
every engine x null type.

| field | requirement |
|---|---|
| `GTs` | numeric matrix, individuals x markers. `colnames(GTs)` **identical to `map$marker`, in order** |
| `map` | data.table with `marker`, `Chr`, `Pos`; one row per marker, same order as `GTs` columns |
| `ld_ws` | numeric matrix, markers x rho windows. `rownames` identical to `map$marker`, in order; `colnames` = the rho windows |
| `decay_sum` | per-chromosome LD-decay summary; needs `Chr` and `b` (and `c` if fitted) |
| `map$true_QTN` | OPTIONAL. Present -> the evaluation section runs. Absent -> treated as real data and skipped |

**`decay_sum` must cover every chromosome in `map`.** A missing chromosome, or one
with an `NA` fit, does not error: `ld_edges()` silently falls back to `r2 = 0.5`,
against a typical fitted value of ~0.27, so that chromosome quietly clusters
almost nothing. The validator checks this explicitly because nothing downstream
will.

`GTs` cannot be omitted. `decay_sum` sets the r^2 *threshold*; `GTs` is what the
actual pairwise r^2 is measured from.

## FILE 2 — `pvals_<cell>_<engine>_<basis>.rds`: one method's output

| field | requirement |
|---|---|
| `p_obs` | numeric, one p-value per marker, `names` identical to `map$marker` **in order** |
| `p_perm` | markers-by-B matrix (rownames = markers) **or** a list of B vectors |
| `basis` | label for how surrogates were built, e.g. `"env_orth"` |
| `engine` | label for the scanner, e.g. `"emmax"` |

Hard requirements:

- **p-values, not F-statistics or -log10p.** The validator range-checks [0, 1].
- **Observed and every surrogate must come from the SAME engine with the SAME
  settings.** Nothing downstream can detect a mismatch.
- **Order is the one thing that cannot be recovered.** A p-vector of the right
  length in the wrong order gives no error and wrong regions. Names are checked,
  never used to reorder — a mismatch stops the run. Send names.
- **B >= 100.** The smallest attainable region p-value is `1/(1+B)`, so B = 100
  floors it at 0.0099. Below 20 is rejected.
- `NA` p-values are tolerated (they simply never become hits).

## What I would like, in priority order

1. **`env_orth`, EMMAX, B = 100, all ten env cells.** Your own header calls
   env_orth the arbiter and the analogue of 3sp's regional permutation. Ten cells
   because env is the replicate axis and a single cell flukes — env2 alone is a
   draw, not a result.
2. **`genetic`, EMMAX, same cells.** EMMAX's home field: the specificity check.
   Silence here is informative on its own.
3. `latent`, LFMM, if the LFMM budget allows — lets the two engines be compared
   on identical footing.

`global_perm` and `spatial` are lower value: your current `SIM_NULL_TYPES`
default already drops them.

## One thing worth checking on your side

The env2 panel I built has **30,628 NA values in `ld_ws`** (~0.5% of cells). Not
fatal — `ld_cscore()` uses `na.rm = TRUE` for the quantile and an NA marker is
simply never a candidate in that window — but it means those markers are
invisible to the C-score at that rho. Worth knowing whether that is expected from
the decay/ld_w step or a symptom.

Also note `emx_gif` records lambda BEFORE genomic control, and
`regen_sim_data.R:81` applies GC only when `gif > 1.1`. See
`HANDOFF_lambda_sweep.md` — deflation below 1 is neither corrected nor flagged,
and 35 of 100 files currently have it.
