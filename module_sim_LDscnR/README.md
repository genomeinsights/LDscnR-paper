# module_sim_LDscnR

Clean, package-only benchmark of the LD-aware C-score outlier-region method on the
Nemo simulations — the sim counterpart of `module_sticklebacks_LDscnR/`, following
the same two-step principle (regenerate → analyse) and using only exported
`LDscnR` functions.

## The point

Unlike 3sp, the sims carry **truth** (`true_pos_QTN`), so this is a **TP/FP
benchmark**, not region discovery:

- the **C-score** (integrated over `rho × q* × alpha`) beats **single-SNP**
  association on PR-AUC, and the **≥2-SNP region filter** is the key precision lever;
- the **structure-aware null** returns a data-driven `(tau_C, l_min)` whose PR sits
  on the sweep's high-PR ridge — calibration lands in the right place;
- the truth-PR ridge over `tau_C × l_min` coincides with the null-quiet zone.

**Always replicate-average env1–5** (mean ± SE) — env1 alone repeatedly flukes.

## Run (two steps)

From the `LDscnR-paper` root (needs the `LDscnR` package installed):

```bash
# 1. regenerate the per-file bundles (heavy; one subprocess per file recommended)
Rscript module_sim_LDscnR/regen_sim_data.R 2 1 all all      # V2_c1, all env, all chr
# 2. benchmark + figures
Rscript module_sim_LDscnR/run_sim_LDscnR.R 2 1              # V2_c1
```

## Data regeneration (`regen_sim_data.R`)

Rebuilds each chromosome bundle from the parsed genotypes with the current LDscnR
machinery: `compute_LD_decay` (**corr**) → `compute_ld_w` → GRM → `emmax` (+GC) →
LFMM (K=5, GC). Saves the **GRM + its marker set** so the structured null's
surrogate EMMAX runs on the identical kinship (engine coherence). **nobgs only**
(bgs is broken — co-author to reimplement).

**GRM basis** (`SIM_GRM` env / `GRM_METHOD`): `"complexity_chain"` is the validated
sim default (`ld_complexity_reduction` rho=0.5 → `ld_prune_and_eMLG`, ~45% kept on
these weak-LD sims). `"ld_w_threshold"` is the 3sp-style direct local-LD filter and
is **under evaluation** — on 3sp the chain over-corrected, so we are checking
whether an `ld_w` threshold is the consistent choice for both datasets (see the sim
GRM PR-AUC test).

Writes per-file bundles to `$SIM_OUT` (default `/Volumes/Nemo/Nemo_sim/regen_sim_data`).

## Analysis (`run_sim_LDscnR.R`)

Pools 10 chromosome files per `(V, c, env)` (each keeps its own saved GRM), builds
the C-score, runs the pooled structured null (one spatial surrogate → per-file
EMMAX → pooled C), calibrates `(tau_C, l_min)`, and scores against truth with
`flag_true_qtns` / `qtn_ld_table` / `evaluate_ors` / `pr_auc`.

Outputs (`figures/`, `results/`):
- `fig1_pr_auc.png` — env-averaged PR-AUC, C-score vs single-SNP, by `l_min`.
- `fig2_heatmaps.png` — env-averaged truth-PR + structured-null region counts over `tau_C × l_min`.
- `fig3_manhattan.png` — C-score Manhattan with true-QTN crosses at the operating point.
- `results/pr_auc_sim.csv`, `heatmap_sim.rds`, `operating_point_sim.csv`, cached null bundles.

## Key `LDscnR` functions used

`compute_LD_decay` / `compute_ld_w` / `ld_complexity_reduction` / `ld_prune_and_eMLG`
/ `emmax` (regen) → `emmax_setup` / `emmax_fast` / `ld_cscore` (vector `alpha`) /
`ld_edges` / `ld_regions` / `null_fdr` / `calibrate_tauc` / `calibrate_lmin` /
`flag_true_qtns` / `qtn_ld_table` / `evaluate_ors` / `pr_auc` / `score_thresholds` /
`ld_manhattan`.

## GRM comparison (`grm_comparison.R`)

Which EMMAX GRM to use on the sims — the sim side of the 3sp GRM question.
Compares, replicate-averaged over env1–5, the **complexity-chain** GRM
(`ld_prune_and_eMLG`, ~45% of markers) against the 3sp-style **`ld_w<b`** GRM
(LD-independent markers below background LD). Reports per GRM: gif, pooled maxC,
and **PR-AUC vs truth on the adaptive tau grid** (the sweep = the distinct
observed C values, so it covers exactly where regions appear — a fixed `[0.05,1]`
grid is meaningless when a cell's C-score maxes out low). Also draws the
`ld_w_095` Manhattan (`figures/grm_ldw095_manhattan.png`) that explains the
mechanism.

**Findings (V2_c1, env1–5):**

| GRM | gif | PR-AUC l_min=1 | l_min=2 | l_min=4 | l_min=8 |
|-----|-----|-----|-----|-----|-----|
| chain | ~1.04 (calibrated) | 0.440 | 0.451 | 0.454 | 0.398 |
| ld_w<b | ~0.93 (deflated) | 0.342 | 0.390 | 0.403 | 0.395 |

- On the sims **~96% of markers fall below `b`** (the QTN included — sim QTN are
  *low*-LD), so `ld_w<b` keeps them in the GRM and mildly **deflates** it
  (gif ≈ 0.93). The chain, which physically thins every region to one
  representative, stays calibrated (gif ≈ 1).
- That deflation costs a little detection at **low l_min** (chain +0.06–0.10
  PR-AUC), but the gap **closes by l_min ≥ 4–8** — the region-size filter
  compensates. So the complexity chain is the (slightly) better sim GRM.
- This is the **opposite of 3sp**, where the signal is *high*-LD (inversions),
  so `ld_w<b` correctly *excludes* it and the chain (keeping 96%, gif 1.061)
  *over-corrects by regional absorption* and collapses the C-score.
- **gif ≈ 1** is a good check against *global* deflation (the sim failure mode)
  but blind to *regional* signal absorption (the 3sp failure mode). The
  transferable principle is to keep the GRM free of the signal — by physical
  thinning where the signal is low-LD (sims) and by an `ld_w` cutoff where it is
  high-LD (3sp).

Default GRM for the sims stays **complexity_chain** (`regen_sim_data.R`).

## Performance (I/O-bound, not RAM-bound)

The runs spend almost all their wall-clock **reading the regen bundles off the
external drive** (`/Volumes/Nemo`), serially, and we re-read the same files many
times. A single run is fine on RAM (~1–2 GB); RAM only bites when several heavy
jobs run at once. Two levers:

1. **Stage a local cache once** (`stage_cache.R`) — copies the V2_c1 files to
   `module_sim_LDscnR/data/cache/` as slim bundles (drops the bulky `LD_decay`
   edge lists, keeps everything the analysis reads). Then point the scripts at it:
   ```bash
   Rscript module_sim_LDscnR/stage_cache.R 2 1               # one-time, ~2 GB
   SIM_DATA=module_sim_LDscnR/data/cache SIM_CORES=5 \
     Rscript module_sim_LDscnR/grm_comparison.R 2 1          # SSD-speed + parallel
   ```
   The cache is regenerable — do not commit it (git-ignore `data/cache/`).
2. **`SIM_CORES`** — env-level parallelism over the 5 (independent) envs via
   `mclapply` (default 1). Only worthwhile once data is local — parallel reads off
   one USB drive just thrash the seek head. With `SIM_CORES>1` the inner
   `qtn_ld_table` drops to 1 thread to avoid oversubscription.

## Status

Written for review; **not yet run** end-to-end (the structured-null pooling section
mirrors the validated `module_sim/R/18d_null_bundle.R` pattern but should be
smoke-tested on one env first). The sim GRM choice (`complexity_chain` vs
`ld_w_threshold`) is pending the PR-AUC comparison.
