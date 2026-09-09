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
  the right control, framed as a hypothesis test:

  **H0:** restricting outlier calls to larger Stage-1 clusters improves
  permutation calibration specifically *because* size signals real
  biology (extended LD around a true adaptive locus), not just because
  it shrinks the candidate list.

  **Test** (`R/14_random_removal_control.R`): instead of restricting the
  candidate pool to units with `n_markers >= floor`, draw `N_RANDOM=100`
  random subsets of the exact SAME cardinality from the full pool,
  ignoring size, and compute `realised_fdr` the same way (per-draw
  `n_obs>0` conditioning -- re-verifying this mattered: the first version
  computed it per-combo as `mean(n_obs_random)` then `max(.,1)`, which
  silently reintroduces the identical `max(.,1)` artefact just fixed
  above; checked directly, 88%/100% of combos had a mean random `n_obs`
  below 1 at floor=20/50, which would have manufactured a spurious
  "random looks better" result on its own -- caught before trusting it).

  **Result: H0 REJECTED**, and completely -- run on the **full grid**
  (all 7 cells x both tags, 10x10 rep x env each = 1400 combos,
  `emmax_consensus`/`group` only). Identical at floor=2 everywhere, by
  construction (floor=2 IS the full pool). Beyond that, random is equal
  or BETTER than size-based restriction at every floor, in every cell,
  in BOTH tags, gap widening with floor:

  | floor | size (bgs) | random (bgs) | size (nobgs) | random (nobgs) |
  |---|---|---|---|---|
  | 3  | 0.669 | 0.637 | 0.699 | 0.692 |
  | 5  | 0.651 | 0.499 | 0.618 | 0.538 |
  | 10 | 0.525 | 0.372 | 0.562 | 0.409 |
  | 20 | 0.470 | 0.269 | 0.511 | 0.346 |
  | 50 | 0.344 | 0.197 | 0.656 | 0.396 |

  (`fig_random_removal_control.pdf`.) **Size-based restriction shows no
  advantage over randomly restricting to the same number of clusters --
  if anything it is consistently worse.** The plateau reported just
  above is therefore NOT evidence that big clusters are specifically
  more trustworthy; it is mostly (perhaps entirely) an artefact of
  restricting the candidate pool at all, something ANY same-sized
  restriction produces, size-based or not.

  **Reconciling with `fig_fp_by_size` (ground truth), and a correction to
  the mechanism first proposed for this.** Checked directly rather than
  assumed: `fig_fp_by_size` genuinely does show FP proportion falling
  with cluster size here -- `bgs` 94.7% -> 80.0% (size-bin 2 -> 50+,
  modest, noisy), `nobgs` 91.9% -> 55.3% (size-bin 2 -> 50+, much
  stronger). So the ground-truth size-truth relationship is real, and
  real*ly* stronger in `nobgs`.

  The first mechanism proposed for the gap above (large units confounded
  with low-recombination regions, which independently have worse
  structure correction per `bgs5`'s `fig_fp_neutral_chr`) makes a
  falsifiable prediction: if that's right, the null-vs-size gap should be
  SMALLER in `nobgs`, where the true signal is much stronger relative to
  any structure-driven noise. **It isn't -- the table above shows `bgs`
  and `nobgs` behaving essentially identically.** That prediction failed,
  so recombination/BGS-specific structure is not the primary driver.

  A cleaner explanation is in the package's own documentation for
  `ld_unit_matrix()`'s `consensus_dosage` representation (what
  `emmax_consensus` tests): *"No size penalty -- a larger cluster is a
  BETTER-ESTIMATED consensus, so size helps rather than hurts."* A bigger
  unit's row-mean dosage is a single, lower-noise variable tested ONCE
  per unit (not a multiple-comparisons effect across markers). Lower
  noise means more statistical power, full stop -- power to detect a real
  QTN association (why `fig_fp_by_size` improves with size) AND power to
  detect any residual, imperfectly-corrected population structure that a
  permuted-but-still-structured phenotype still carries (why the
  permutation null's own hit rate rises with size). Both effects come
  from the same noise reduction, independent of BGS or recombination
  specifics -- which is exactly why `bgs` and `nobgs` look the same.

  This does **not** call `fig_fp_by_size` itself into question -- it is
  unaffected, a separate ground-truth analysis. What it retracts is the
  claim that the permutation-based size-floor sweep independently
  corroborates it via "a truth-free route": it doesn't. Size predicts
  truth (only checkable with ground truth) but does not validate itself
  against a permutation null, because the null gets the same power
  boost. **Still not run:** `mvn`/`spatial` schemes, `emmax_simes` arm --
  the `group`/`emmax_consensus` result is now complete (full grid) and
  the effect is large and consistent, but the mechanism explanation above
  (a property of `consensus_dosage` specifically) predicts `emmax_simes`
  -- a different combining rule, genuinely multi-comparison across a
  unit's markers -- could behave differently, and hasn't been checked.

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
- `R/14_random_removal_control.R`'s matched-cardinality-random check is
  now complete for `group`/`emmax_consensus` (full 7-cell x 2-tag grid),
  but still not run for `mvn`/`spatial` or `emmax_simes` --
  `emmax_simes` is the more important of the two to check, since the
  proposed mechanism (a `consensus_dosage`-specific noise/power effect)
  predicts it could behave differently (Simes IS a genuine multiple-
  comparisons combination across a unit's markers, unlike consensus's
  single averaged variable). `results/random_removal_control_summary.rds`
  holds the full pooled output; extend by rerunning
  `R/14_random_removal_control.R <tag> <cell>` with the arm changed to
  Simes (not currently a script argument -- would need a small edit) and
  `rbindlist()`-ing the new `out/14_random_removal_control/rrc_*.rds`
  files in with it.
