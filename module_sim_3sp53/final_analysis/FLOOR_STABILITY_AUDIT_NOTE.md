# Floor-stability analysis: audit note

Written 2026-09-27, updated the same day after the real run completed on
mini1 (see "Execution status"). **Findings below are real, not placeholders.**

## What this analysis asks

In the 600 primary c=1 simulations (`V0.5_c1`, `V1_c1`, `V2_c1`; both BGS
treatments; 10 maps x 10 environmental continuations), does requiring a
reported Stage-2 region to persist across three Stage-1 size floors
(99.5%/99.7%/99.9% test-count reduction relative to all assayed markers)
increase precision, and how much recall is lost? PK's stated prediction --
lower recall, potentially substantially higher precision -- is a hypothesis
this analysis tests, not a result assumed going in.

This is entirely additive: `R/12_floor_decomposition.R`, `results/
simulation_stage2_region_details.rds`/`simulation_performance.tsv`, and every
manuscript value are untouched. New files only, under `floor_stability_*`
names.

## Execution status -- run on mini1, complete, 0 failures

This laptop cannot reach the per-combo bundles
(`stage_dir("02_build_ld_units", combo_id)/ld_units.rds`,
`stage_dir("03_emmax", combo_id)/emmax.rds`, optionally `04_lfmm.rds`) --
`00_config.R`'s `PATHS$parsed` points at an external volume
(`/Volumes/Large_storage/...`) not mounted here, and `out_final_v1/` on this
host contains only `18_bgs_marker_noise/`. Confirmed reachable instead on
`mini1` (SSH alias; `/Volumes/Large_storage/module_sim_3sp53_parsed_final_v1`
and the full `out_final_v1/` tree are both present there). `mini2` was also
checked and does NOT have this data -- its only relevant drive (`T9`) holds
the older, pre-correction-parser `module_sim` pipeline's cache/out/parsed
directories, not `module_sim_3sp53`/`final_analysis`.

The six new files (this note, `R/helpers_floor_stability.R`,
`R/21`-`R/23`, the figure script) were copied to mini1's existing
`LDscnR-paper` checkout by `scp` -- deliberately NOT via `git pull`, since
mini1's clone was several commits behind and had its own uncommitted local
edits to five files (including `ADDITIONAL_ANALYSES_AUDIT.md`); copying only
these six new, non-colliding filenames avoided that risk entirely. Confirmed
via `git status` on mini1 immediately after copying that nothing else
changed.

Run sequence `R/21` -> `R/22` -> `R/23` -> figure script, each gated on the
previous succeeding, executed 2026-09-27 14:16-14:20 (~4 minutes total, far
faster than `R/12`'s own ~10-hour full run -- this analysis evaluates far
fewer floor x arm combinations per combo). All four stages exited 0. `R/21`:
600/600 combos, 0 failures. `R/22`: 600/600 combos, 0 failures, 3,011 region
rows produced. `R/23`: 1,062 canonical (99.7%-floor) regions tiered across
600 combos x up to 3 methods.

**What WAS done to de-risk delivery before this run, on top of the real
run's own success (0 failures):**
- `select_stability_floors()` and `match_region_sets()` (the two genuinely
  new pieces of logic this analysis needed -- everything else reuses
  `R/12_floor_decomposition.R`'s and `R/helpers_stage2_truth.R`'s existing,
  already-reviewed machinery verbatim) were each directly unit-tested against
  hand-verified synthetic cases: floor selection against a bimodal cluster-
  size distribution with a known exact answer (confirmed matched); region
  matching against single-match, ambiguous one-to-many, no-overlap, and
  physical-overlap-only-vs-core_snp-required cases (all confirmed matched by
  hand).
- The full pooling/bootstrap/tier/call-set logic in `R/23` (the part with the
  most surface area for a subtle bug) was separately smoke-tested end to end
  against a synthetic 20-combo, 1-method dataset with a known built-in
  precision-by-tier gradient, confirming: recall values are finite, the
  incremental tier groups (`incr_eq1`/`incr_eq2`/`incr_eq3`) exactly partition
  the "all canonical" set's region count, and the recovered precision
  gradient across tiers matched what was built into the synthetic data.
- `R/22`'s association/assembly logic (`.region_set()`) is `R/12`'s own arms
  C/D machinery with the floor set changed and the return shape widened from
  summary counts to full region detail -- not independently re-verifiable
  without the bundles, but not new logic either.

## Headline finding

**PK's hypothesis is supported, and survives paired uncertainty, for all
three methods tested (EMMAX consensus, EMMAX Simes, LFMM Simes).** Precision
rises monotonically with tier (incr_eq1 -> incr_eq2 -> incr_eq3), e.g. EMMAX
consensus: 0.317 -> 0.524 -> 0.692; EMMAX Simes: 0.308 -> 0.574 -> 0.737; LFMM
Simes: 0.375 -> 0.514 -> 0.749 -- the same direction and a comparable
magnitude every time, not a fluke of one method. Recall falls in the same
direction and, at 3/3, substantially: pooled recall drops from ~0.22-0.24
(all canonical regions) to ~0.13-0.15 (3/3-stable regions only).

| method | tier | n_regions | precision | recall | FP reduction vs. all canonical |
|---|---|---:|---:|---:|---:|
| emmax_consensus | all (>=1/3) | 379 | 0.578 | 0.226 | -- |
| emmax_consensus | >=2/3 | 338 | 0.609 | 0.219 | 17.5% |
| emmax_consensus | 3/3 | 172 | 0.692 | 0.140 | 66.9% |
| emmax_simes | all (>=1/3) | 319 | 0.630 | 0.218 | -- |
| emmax_simes | >=2/3 | 293 | 0.659 | 0.212 | 15.3% |
| emmax_simes | 3/3 | 152 | 0.737 | 0.135 | 66.1% |
| lfmm_simes | all (>=1/3) | 364 | 0.613 | 0.239 | -- |
| lfmm_simes | >=2/3 | 340 | 0.629 | 0.233 | 10.6% |
| lfmm_simes | 3/3 | 167 | 0.749 | 0.149 | 70.2% |

### Does the precision gain survive paired uncertainty? Yes.

Map-cluster-bootstrapped (clustered by `rep`, 2,000 replicates) paired
contrasts, `results/floor_stability_contrasts.tsv`:

| method | contrast | diff precision (95% CI) | diff recall (95% CI) |
|---|---|---|---|
| emmax_consensus | ge2 - all | +0.032 [0.015, 0.052] | -0.008 [-0.014, -0.002] |
| emmax_consensus | eq3 - all | +0.114 [0.071, 0.156] | -0.086 [-0.103, -0.067] |
| emmax_simes | ge2 - all | +0.029 [0.012, 0.047] | -0.005 [-0.010, -0.001] |
| emmax_simes | eq3 - all | +0.107 [0.055, 0.154] | -0.083 [-0.095, -0.070] |
| lfmm_simes | ge2 - all | +0.017 [0.005, 0.034] | -0.007 [-0.011, -0.002] |
| lfmm_simes | eq3 - all | +0.136 [0.095, 0.179] | -0.090 [-0.107, -0.074] |

Every precision-difference CI excludes zero and is positive; every
recall-difference CI excludes zero and is negative. The `>=2/3` gain is real
but modest; the `3/3` gain is real and substantial (roughly +0.11 to +0.14
precision) at a real, substantial recall cost (roughly -0.08 to -0.09). This
pattern held up in every V level and both BGS treatments individually
(`results/floor_stability_summary_by_strata.tsv`; not separately
bootstrapped per stratum, but every one of the 18 strata x tier point
estimates moves the same direction as the pooled headline).

### Floor selection and collapse

`results/floor_stability_floor_selection_wide.tsv`, 600/600 combos, 0
failures. Achieved reductions sit very close to nominal targets (mean
0.9951/0.9970/0.9990 vs. 0.995/0.997/0.999 targets), and floor collapse is
rare: 99.5%==99.7% for only 2/600 combos (0.3%); 99.7%==99.9% and
99.5%==99.9% never collapse (0/600). 598/600 combos have three genuinely
distinct floors. Mean eligible Stage-1 units per combo: ~75 at the 99.5%
floor, ~45 at 99.7%, ~15 at 99.9% -- a stringent candidate set even at the
loosest of the three targets.

### Cross-floor matching ambiguities

Of 1,062 canonical regions, 28 (2.6%) have an ambiguous (>1) side-floor
match. This concentrates heavily at the STRICTER side: 24 ambiguous matches
at the 99.9% floor vs. only 4 at the 99.5% floor -- consistent with a
mechanism where dropping units at a stricter floor can fragment one
canonical region's territory into multiple, only-partially-overlapping
regions at that side, rather than the reverse. The physical-overlap-only
sensitivity variant (`tier_overlap_only`) agrees with the primary
core_snp+overlap criterion (`tier_effective`) for 1,030/1,062 regions
(97.0%); every disagreement moves the SAME direction (overlap-only assigns a
higher or equal tier, never lower -- expected, since overlap is a necessary
but not sufficient condition for the core_snp criterion). Neither the
ambiguity rate nor the sensitivity-variant divergence is large enough to
change the headline direction, but both are real, not negligible.

### Completeness

All 600 combos processed without a worker failure (`R/22`'s hard `stop()` on
any error never triggered). Zero-call rates (fraction of (combo, method)
cells with NO canonical region at the 99.7% floor) are HIGH and worth
flagging prominently: 43-76% depending on method/V/tag
(`results/floor_stability_zero_call_summary.tsv`) -- e.g. EMMAX consensus
ranges from 50% (V0.5, bgs) to 74% (V2, bgs). Zero-call rate rises with V
(weaker selection, presumably less signal to detect) for every method. This
means the tier table above describes the MINORITY of (combo, method) cells
that produce any canonical region at all at this floor -- the 99.7%-derived
floor is a genuinely aggressive filter, not a mild one. These zero-call
cells are correctly pooled as real zeros throughout (never dropped, never
imputed).

### Reference comparison (existing floor=2 canonical) -- separate, not blended

`results/floor_stability_reference_floor2.tsv`: at the manuscript's own
canonical floor=2 (same 600 c=1 combos, EMMAX methods only), precision is
0.564-0.580 and recall is 0.417-0.418 -- notably HIGHER recall than even this
analysis's "all canonical" (>=1/3, i.e. no stability filter) baseline at the
99.7%-derived floor (recall 0.218-0.226). This is expected and NOT part of
the stability-filter finding above: floor=2 admits far more eligible units
per combo (a much looser floor) than the 99.7%-derived floor does by design
(mean ~45 eligible units), so simply choosing this analysis's canonical
floor already trades substantial recall for modest precision, BEFORE any
cross-floor stability requirement is applied at all. Reported here as
context only, per PK's explicit instruction to keep "changing the canonical
floor" and "imposing a stability filter" as two separate, non-conflated
interventions.
