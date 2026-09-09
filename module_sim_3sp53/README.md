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

  **`group`**: mean `realised_fdr` 0.70 (median 0.49) at the base floor,
  cell means ranging 0.50 (`V2_c1`) - 0.80 (`V0.5_c1`) -- confirms the
  conservative-per-ICC prediction, now on 400 replicate combos per cell
  instead of 1.

  **`mvn` negative control**: holds cleanly at full scale -- mean 0.035
  (median 0.03), range 0-0.17 across all 2800 (tag x arm x combo) draws.
  The method itself is not a source of miscalibration anywhere in the grid.

  **`spatial`**: a genuinely new finding only visible at full-grid scale
  (the 14-combo sample had n=1 per cell, too noisy to trust) -- a sharp,
  dispersal-dependent split. `realised_fdr` stays low, same order as the
  `mvn` negative control, in every high-dispersal (c1) cell: `V0.5_c1`
  0.40, `V1_c1` 0.26, `V2_c1` 0.24. It explodes in every medium/low-
  dispersal (c1.5/c2) cell: 6.8-13.5. At low dispersal, pure spatial
  autocorrelation over individual coordinates -- no genotype, no genetic
  signal at all -- alone generates *more* "significant" surrogate hits
  than the real observed data contains. Independent of the group-ICC
  argument (which only speaks to the 5 coarse spatial bins), this shows
  isolation-by-distance itself, at fine spatial scale, is the dominant
  false-positive driver at low dispersal, not just coarse population
  membership.

  **Minimum-cluster-size sweep** (floors 2/3/5/10/20/50), now on 2800
  (tag x cell x rep x env x arm) trajectories instead of 13:
  `group`-null `realised_fdr` falls monotonically with floor in the
  pooled mean (0.70 at floor=2 -> 0.03 at floor=50) and strictly
  monotonically within 65% of individual trajectories (1826/2800) --
  same conclusion as the sample run, now on 200x the data: larger
  clusters are more trustworthy. Converges with the ground-truth-based
  `fig_fp_by_size` finding on the same conclusion via an independent,
  truth-free route, and the largest surviving clusters remain reasonable
  candidates pending independent evidence (GO enrichment, known genes)
  even where they don't clear a strict permutation FDR at the base floor.

  The `nobgs`/`V1_c1` elevation flagged as an unexplained exception in the
  14-combo sample (0.625-1.05 even at floor=50) is **resolved**, not
  confirmed, by the full grid: averaged over its 200 replicate combos (10
  reps x 10 envs x 2 arms), `realised_fdr` at floor=50 is 0.082, in line
  with every other cell/tag. The earlier number was noise from a single
  rep1/env1 draw, not a persistent property of that cell.

  No meaningful `tag` effect on `group`-null calibration (`bgs` 0.69 vs
  `nobgs` 0.71) or `arm` effect (`emmax_consensus` 0.68 vs `emmax_simes`
  0.72) -- expected, since these nulls calibrate the association-test
  procedure itself, not something BGS status or arm choice should move.

## Not yet done

- `env_icc`'s implications haven't been folded back into `module_sim`'s
  own live pipeline or `module_sim_bgs5`'s archive.
- The `spatial`-null finding above (dispersal-dependent blowup) hasn't
  been cross-checked against an independent ground-truth signal the way
  the size-floor sweep was -- worth a look if it ends up load-bearing for
  the manuscript's low-dispersal discussion.
