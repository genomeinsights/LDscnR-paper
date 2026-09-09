# module_sim_3sp53 -- canonical, live, 2026-09-08

The new NEMO production simulations ("3spC", 53 cM maps, finer resolution
`MAP_RES 3.3e-6` vs `bgs5`'s `3.3e-5`), processed on Petri's mini
(`192.168.1.150:/Volumes/Large_storage/prod_out`). Full 7-cell x 2-tag
design (`V0.5_c1`, `V0.5_c1.5`, `V0.5_c2`, `V1_c1`, `V1_c1.5`, `V2_c1`,
`V2_c1.5`), 1400 runs verified complete by the NEMO session (0 failures).

**PK's decision: this is the canonical dataset for reported values.**
`module_sim_bgs5` (4/7 cells, coarser map) is kept as a sensitivity check
for if a reviewer asks specifically -- not superseded, but not the primary
source either. Basis for the decision: a direct reproducibility comparison
restricted to the 4 cells `bgs5` covers showed the two datasets agree
strongly on everything that matters for the paper's core claims (Precision
r=0.93, Recall r=0.94, Fst r=1.00, method ranking identical), differing
mainly in *how strong* the BGS-recombination signal looks -- expected, and
already explained by the coarser map's recombination-distance artefact.

Unlike `bgs5`/`bgs2`, this is NOT a frozen snapshot -- `module_sim/` itself
still points at `bgs5` (kept working, unmodified, as its own sensitivity
line); this directory is a **separate, live** module that gets rerun as
the pipeline evolves, same reasoning as the `module_sim_bgs5` archive's own
"why a separate module" note (mini-local paths, doesn't destabilise
`module_sim/`).

## Same analysis set as module_sim_bgs5, plus the structured-null work

`results/` carries the same four canonical files (`pooled_pr.rds`,
`popgen_summary.rds`, `bgs_recomb_summary.rds`, `bgs_windows_summary.rds`)
as `module_sim_bgs5`, all of today's audit-driven fixes already applied
(single-SNP true-cluster-size scoring, matched FP-by-size denominator,
neutral-chromosome FP analysis, cluster bootstrap CI), plus two things
`bgs5` doesn't have:

- **`env_icc`** (in `popgen_summary.rds`, alongside Fst/Va): a one-way
  ANOVA variance decomposition of the environment across the 5 spatial
  population groups (kmeans(k=5) on population (x,y) -- confirmed
  directly, populations were sampled in the 4 corners + centre of the
  48x48 grid). ~0.887 universally, identical at every cell/tag (it depends
  only on population positions and env values, neither of which differs
  by cell or tag) -- ~89% of environmental variance sits BETWEEN the 5
  groups, only ~11% within. This is why a within-group phenotype
  permutation null reads conservative, and it's the same root cause behind
  the real false-positive rate itself: env is nearly a deterministic
  function of which genetically-differentiated deme an individual belongs
  to, so imperfectly-corrected structure leaks through as spurious
  environmental association.

- **`structured_null_summary.rds`**: a calibration check matching
  `module_3sp`/`module_9sp`'s manuscript methodology (see
  `LDscnR_manuscript/materials_and_methods.tex`, "Association models and
  empirical calibration" / "Nine-spined stickleback permutation stress
  test") as closely as the simulated design allows
  (`module_sim/R/12_structured_null.R`). Three null schemes, all built on
  `ld_outlier_test()`'s own `.ld_outlier_tested_units()` machinery:
  - `group`: population-level env-value permutation within the 5 spatial
    groups (the direct analog of module_3sp's regional-locality / module_
    9sp's lineage stratification) -- conservative, per the ICC above.
  - `mvn`: `s ~ MVN(0, GRM)`, orthogonalised against observed y -- a
    negative control.
  - `spatial`: `s ~ MVN(0, Gaussian kernel over individual (x,y))`,
    orthogonalised the same way -- isolates spatial confounding
    independent of the 5-group binning.

  Run 2026-09-08 22:53-23:50 (mini, concurrency 8) across the **full
  1400-combo grid** (7 cells x 2 tags x 10 reps x 10 envs), 0 failures --
  supersedes the earlier 14-combo representative-sample run (rep=1/env=1
  per cell x tag) that this section originally described. Pooled to 8400
  summary rows (1400 combos x 2 arms x 3 schemes) and 50400 sweep rows
  (x7 size floors) in `results/structured_null_summary.rds`.

  **[!] CORRECTED 2026-09-09** -- PK: "are you accounting for the fact
  that the number of false positives also drop by randomly removing
  clusters from the outlier list?" No, not adequately, and it materially
  changed two conclusions below. `realised_fdr = mean_surrogate /
  max(n_obs, 1)` uses `max(n_obs,1)` to avoid dividing by zero when the
  OBSERVED test finds no significant unit at all at some size floor. But
  the eligible pool of Stage-1 units genome-wide collapses fast with floor
  (mean ~2990 units at floor=2 down to ~2 at floor=50, some cells at 0 by
  floor=50, checked directly against `.ld_outlier_units()`'s own candidate
  table), so more and more combos hit `n_obs=0` as floor rises -- 36.5% at
  floor=2, **95.2% at floor=50**. For those combos `realised_fdr` isn't
  measuring an FDR at all; it collapses to the raw rate at which pure
  permutation noise clears an ever-emptying pool, which mechanically
  trends toward 0 regardless of whether real discoveries are more
  trustworthy at large sizes. All numbers below are now reported
  **conditional on `n_obs>0`** (the only combos where the ratio is a
  genuine FDR estimate) -- the pre-correction pooled numbers are kept
  alongside, in parentheses, to show how much they were inflated.

  **`group`**: mean `realised_fdr` 0.79 (was 0.70 pooled) at the base
  floor. By cell: 0.66-0.95, `V0.5_c1` highest, `V2_c1` lowest -- same
  ICC-conservatism conclusion as before, just a higher and more honest
  number.

  **`mvn` negative control**: unaffected in substance -- 0.025 corrected
  vs 0.035 pooled, still low and stable. The method itself is not a
  source of miscalibration.

  **`spatial`**: the dispersal-dependent gradient **survives** the
  correction, at somewhat smaller magnitude -- low in every high-
  dispersal (c1) cell (`V0.5_c1` 0.30, `V1_c1` 0.22, `V2_c1` 0.19, all
  comparable to `mvn`), 5.2-8.8 in every medium/low-dispersal (c1.5/c2)
  cell (was 6.8-13.5 pooled). Isolation-by-distance alone still outbids
  the real observed discoveries at low dispersal even correcting for pool
  shrinkage.

  **Minimum-cluster-size sweep**, corrected -- restricted at each floor to
  combos with a genuine discovery there (n=1777 at floor=2, falling to
  n=134 at floor=50, since that's also how fast real discoveries become
  rare, independent of the ratio issue): `group`-null `realised_fdr`
  declines from 0.79 (floor=2) to 0.52 (floor=10) and then **plateaus**
  around 0.52-0.57 through floor=50 -- a real, modest effect, NOT the
  "falls to ~0" collapse originally reported (that was floor 20-50 being
  95%+ dominated by zero-discovery combos). Monotone in 41.7% of
  individual trajectories (676/1622 with >=2 floors to compare), not 65%.

  **[!] SUPERSEDED 2026-09-09 (a few hours after the correction above) --
  PK: "if we randomly remove the same number of outlier clusters (but
  independently of size), would we also see a decline in FPs, and is the
  decline in fig_structured_null_sizesweep relative to that?"** This is
  the right control and the answer changes the conclusion again, this
  time in direction, not just magnitude. `R/14_random_removal_control.R`:
  instead of restricting the candidate pool to units with `n_markers >=
  floor`, draw `N_RANDOM=100` random subsets of the exact SAME
  cardinality from the full pool, ignoring size, and compute
  `realised_fdr` the same way (same `n_obs>0` conditioning -- this had to
  be re-verified per draw, not per-combo mean-then-floor, which silently
  reintroduces the identical `max(.,1)` artefact just fixed: checked
  directly, 88%/100% of combos had a mean random `n_obs` below 1 at
  floor=20/50, which would have manufactured a spurious "random looks
  better" result on its own before the per-draw fix caught it).

  Run on `bgs`, the full 10x10 rep x env grid (100 combos), for 4 of the 7
  cells (`V0.5_c1`, `V0.5_c2`, `V1_c1.5`, `V2_c1` -- a spread across the
  dispersal/variance grid), `emmax_consensus` only. **Identical at
  floor=2 in every cell, by construction** (floor=2 IS the full pool, so
  "random, same cardinality" trivially equals "size >= 2"). Beyond that,
  in **every one of the 4 cells checked**, random is equal or BETTER
  (lower `realised_fdr`) than size-based restriction, and the gap widens
  with floor:

  | floor | size-based (pooled) | random-matched (pooled) |
  |---|---|---|
  | 3  | 0.669 | 0.652 |
  | 5  | 0.630 | 0.526 |
  | 10 | 0.519 | 0.415 |
  | 20 | 0.489 | 0.297 |
  | 50 | 0.309 | 0.179 |

  (`fig_random_removal_control.pdf`; per-cell panels show the same
  direction in all 4, not just the pooled average.) **Size-based
  restriction shows no advantage over randomly restricting to the same
  number of clusters -- if anything it is consistently worse.** The
  plateau reported just above is therefore NOT evidence that big clusters
  are specifically more trustworthy; it is mostly (perhaps entirely) an
  artefact of restricting the candidate pool at all, something ANY
  same-sized restriction produces, size-based or not.

  A plausible mechanism, found while building the single-combo pilot for
  this check: the *surrogate* (pure-noise) significance rate per unit is
  NOT flat across sizes -- it is measurably higher among large units than
  a random-matched baseline (roughly 0.00014 at floor=2 rising to ~0.0012
  at floor=20 for the size-restricted pool, vs. a flat ~0.00014 for
  random). Large Stage-1 units are plausibly confounded with low-
  recombination regions, which this project has already shown
  independently have worse false-positive control (`fig_fp_neutral_chr`
  in `bgs5`) -- big clusters may concentrate both true AND false signal,
  not purify for true signal.

  This does **not** call `fig_fp_by_size` into question -- that is a
  separate, ground-truth (known-QTN) TP/FP analysis, not a permutation-
  null comparison, and is unaffected by anything here. What it DOES
  retract is the earlier claim that the permutation-based size-floor
  sweep independently corroborates `fig_fp_by_size` via "a truth-free
  route" -- it doesn't; once properly controlled, it points the other way.
  **Not yet run:** `nobgs`, the remaining 3 cells (`V0.5_c1.5`, `V1_c1`,
  `V2_c1.5`), `mvn`/`spatial`, `emmax_simes` -- the 4-cell/bgs/consensus
  result is consistent and the effect size is large, but treat it as
  strong-not-final until/unless the remaining scope is checked.

  **The `nobgs`/`V1_c1` anomaly was NOT resolved -- it was CONFIRMED, and
  the earlier "resolved" claim was itself an artefact of the same
  `max(n_obs,1)` problem.** Of its 200 replicate combos, only 16 ever
  reach floor=50 with a genuine discovery; restricted to those, mean
  `realised_fdr` = **0.937** -- the single worst cell x tag combination in
  the entire 1400-combo grid at that floor, not "in line with everything
  else." At every floor its corrected `realised_fdr` stays elevated
  (0.75-0.94), never declining the way most other cells do (contrast
  `bgs`/`V0.5_c2`, a comparably-powered cell, which runs 0.39-0.59 over
  the same range). Genuinely anomalous, still unexplained, still not
  investigated further -- carried forward again, this time with the right
  number.

  No meaningful `tag` effect on corrected `group`-null calibration (`bgs`
  0.78 vs `nobgs` 0.80) or `arm` effect (`emmax_consensus` 0.75 vs
  `emmax_simes` 0.83) -- both conclusions hold up under the correction.

## Not yet done

- `env_icc`'s implications haven't been folded back into `module_sim`'s
  own live pipeline or `module_sim_bgs5`'s archive.
- The `spatial`-null finding above (dispersal-dependent blowup) hasn't
  been cross-checked against an independent ground-truth signal the way
  the size-floor sweep was -- worth a look if it ends up load-bearing for
  the manuscript's low-dispersal discussion.
- `nobgs`/`V1_c1`'s calibration anomaly (confirmed, not resolved -- see
  the 2026-09-09 correction above) is still unexplained. Worth checking
  whether it's specific to that one cell's demography or a more general
  pattern once other cells get similarly scrutinised.
- Any future size-floor (or similarly conditioned) analysis on this data
  MUST report `realised_fdr` conditional on `n_obs>0` (or at minimum
  alongside the `n_obs>0` fraction) -- the `max(n_obs,1)` convention in
  `R/12_structured_null.R` silently makes the pooled ratio meaningless
  once most combos stop producing genuine discoveries, which happens
  faster than it looks (36.5% -> 95.2% zero-discovery across floor=2-50
  here). This was the mistake corrected above.
- `R/14_random_removal_control.R`'s matched-cardinality-random check
  (see above) has only been run on `bgs`/4-of-7-cells/`emmax_consensus`.
  Extending it to `nobgs`, the remaining 3 cells, `mvn`/`spatial`, and
  `emmax_simes` would either confirm this is general or localise it --
  worth doing before this result goes anywhere near the manuscript.
  `results/random_removal_control_summary.rds` holds the pooled 4-cell
  output; rerun `R/14_random_removal_control.R <tag> <cell>` per missing
  combo and `rbindlist()` the new `out/14_random_removal_control/rrc_*.rds`
  files in with it.
