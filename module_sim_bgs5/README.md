# module_sim_bgs5 -- frozen snapshot, 2026-09-05

The first real BGS-comparison results, built against the correctly-
parameterized simulation set: `/Volumes/Nemo/Nemo_sim/bgs5/` ("scenario B",
1000 deleterious loci per chromosome, confirmed 1720 in a direct spot-check --
see `module_sim_bgs2/README.md` for the full story of why the earlier
`Nemo_out_bgs` archive was wrong, essentially zero deleterious load). Both
tags (`nobgs`, `bgs`) here are real; unlike `module_sim_bgs2`, nothing in this
snapshot needs to be treated as invalid.

**Coverage: 4 of the 7 target cells** -- `V0.5_c1`, `V0.5_c2`, `V1_c1.5`,
`V2_c1`. The remaining 3 (`V0.5_c1.5`, `V1_c1`, `V2_c1.5`) are still being
simulated; PK judged 4 cells enough for a first thorough comparison rather
than waiting. **Consequence for the high-dispersal (c=1) figures**
(`fig_fbeta.pdf`, `fig_precision_recall_scatter.pdf`, both restricted to
`V0.5_c1`/`V1_c1`/`V2_c1` per PK's "no point comparing" scoping): `V1_c1`
isn't one of the 4 available cells, so these are built from **2** high-
dispersal cells here (`V0.5_c1`, `V2_c1`), not 3. Re-run once `V1_c1` lands.

## What's new here vs. module_sim_bgs2

This snapshot was built AFTER a round of audit-driven fixes (external review
of the manuscript's simulation section), so it is not a like-for-like rerun
of bgs2's old code on new data:

- **Single-SNP arms score each significant marker as its own size-1 region**
  (not the enclosing multi-marker Stage-1 unit) -- R/04_score.R, PK: "Single
  SNPs go all the way to 1 not 2." Visibly changes single-SNP precision/FP
  counts (much lower precision, e.g. V2_c1 nobgs emmax_snp: 5067 FP vs 36 TP)
  compared to what bgs2's archived numbers would have shown under the old
  logic.
- FP-by-cluster-size uses a matched TP/(TP+FP) denominator at both the point
  estimate and its SE, plus a cell-stratified breakdown.
- Pooled Precision/Recall/PR carry both the across-env SE and a cluster
  bootstrap CI (resample envs, then reps within env).
- "EMMAX representative" -> "EMMAX consensus"; the PR-scatter's method-
  joining line (implied a trajectory 5 discrete methods don't have) is gone.

## Two independent BGS-effect estimators, now both showing a real signal

Unlike bgs2 (built on ~3 deleterious loci genome-wide -- functionally no BGS
at all), bgs5's real deleterious load produces a clean, mechanistically
sensible result on BOTH estimators built this session:

- **Fst-based** (R/08-09, `fig_bgs_recomb.pdf`): log2(Fst, bgs/nobgs)
  decreases monotonically from low recombination (0.197 +/- 0.016) to high
  recombination (0.126 +/- 0.011) -- the expected Hill-Robertson pattern
  (BGS interference stronger where recombination can't break up linkage to
  the deleterious background).
- **n_snp-window ratio** (R/10-11, `fig_bgs_windows.pdf`; the estimator
  established as the correct one -- a ratio like Fst can stay flat even as
  absolute diversity drops): B_obs ~0.77 genome-wide (BGS removes ~23% of
  segregating sites on average), falling to ~0.49-0.54 in the lowest-
  recombination quintile (Q1) vs ~0.87-0.88 in the highest (Q5) across all 4
  cells -- a much stronger, cleaner version of the same gradient.

`run_bgsrecomb_grid.sh` was missing the bundle-existence skip guard the
other grid drivers have when this first ran against the 4-cell grid (fixed
in the live `module_sim/`, ported here too) -- failed loudly instead of
skipping on the 3 not-yet-simulated cells, but did not affect the 800 valid
combos' results.

## Kept here as a snapshot

As with `module_sim_bgs2`: a permanent record of this comparison point, not
a directory that regenerates. The live pipeline (`module_sim/`) keeps
running against `/Volumes/Nemo/Nemo_sim/bgs5/` and will pick up the
remaining 3 cells automatically once their archives land (grid drivers skip
missing cells rather than fail).
