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
    negative control. Stays low and stable (~0.01-0.09) everywhere: the
    method itself is not the source of miscalibration.
  - `spatial`: `s ~ MVN(0, Gaussian kernel over individual (x,y))`,
    orthogonalised the same way -- isolates spatial confounding
    independent of the 5-group binning; often the most anti-conservative
    of the three, consistent with real spatial autocorrelation (not just
    coarse group membership) driving false positives.

  Run on a REPRESENTATIVE 14-combo sample (all 7 cells x 2 tags, rep=1
  env=1 fixed), matching module_3sp/9sp's own scale (one calibration check
  per dataset) rather than an exhaustive 1400-combo resweep. Also includes
  a truth-free minimum-stage-1-unit-size sweep (floors 2/3/5/10/20/50):
  `realised_fdr` under the `group` null falls monotonically with cluster
  size in 12 of 13 combos with real observed discoveries -- often from
  ~0.2-0.5 at floor=2 to 0 (the null never produces one) by floor=20-50.
  Converges with the ground-truth-based FP-by-size finding (`fig_fp_by_
  size`) on the same conclusion via an independent, truth-free route:
  larger clusters are more trustworthy, and the largest surviving ones are
  reasonable candidates pending independent evidence (GO enrichment,
  known genes) even where they don't clear a strict permutation FDR.
  One exception worth carrying forward, not smoothed over:
  `nobgs`/`V1_c1` stayed elevated (0.625-1.05) even at floor=50.

## Not yet done

- The 3 remaining structured-null combos' worth of coverage is thin (14
  combos, not the full grid) -- treat cell-level patterns as indicative,
  not final, until/unless a wider sweep is run.
- `env_icc`'s implications haven't been folded back into `module_sim`'s
  own live pipeline or `module_sim_bgs5`'s archive.
