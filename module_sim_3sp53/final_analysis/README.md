# final_analysis -- clean manuscript-simulation pipeline (in progress)

Rebuild of `module_sim_3sp53` per `../CLAUDE_REANALYSIS_INSTRUCTIONS.md`, kept in
its own subtree so the existing exploratory work in `module_sim_3sp53/` is left
untouched (see that document's "Preserve existing work" section -- in
particular `R/14_random_removal_control.R` there currently has an uncommitted
change that must not be overwritten, moved, restored, or reformatted).

## Status: Phase 1 complete; Phase 2 primary EMMAX arm + LFMM portability
## arm both complete; Phase 3 gates 1-4 ALL PASS -- the full 1,400-
## combination primary EMMAX grid and the full 1,400-combination LFMM
## grid are both complete, `results/` has pooled precision/recall with
## bootstrap CIs for both arms, and 4 figures are built (2 named in the
## instructions, 2 illustrative supplementary Manhattan figures).

Implemented so far:

- `R/00_config.R` -- paths (all new roots, none overlapping the old
  `module_sim_3sp53` config), the full design grid (7 cells x 2 tags x 10
  reps x 10 envs), Stage-1/GRM/EMMAX/truth-scoring constants (rho=0.50,
  size floor=2, GCTA kinship on stage-1-pruned markers, alpha=0.05,
  rho_r2=0.75/rho_d=0.95 truth-matching, Va-share detectability >=5%),
  receipt/caching machinery (ported from `module_sim_3sp53`'s own, minus its
  hard LDscnR-version pin -- this session's own package fix moved that
  source, so a copied pin would fail spuriously), host-path validation.
- `R/01_parse_nemo.R` -- the corrected parser. Nemo locus identifiers are
  one-based; the inherited `idx := as.numeric(idx) + 1` (present in both
  `module_sim/R_parsing/01_parse_nemo.R` and `module_sim_3sp53/R_parsing/
  01_parse_nemo.R`, both still live and both still wrong) is removed. Runs a
  full assertion suite per combination (see `qc/`).
- `R/01_parse_nemo_qa.R` -- the three-run old-vs-corrected validation
  (Phase 3 gate 1, PASS): reparses the same three native runs under both the
  old (`+1`) and corrected (`+0`) offset, prints the QTN index mapping for
  human verification, and writes `qc/parser_qa_report.tsv` +
  `qc/three_run_offset_comparison.tsv`.
- `R/02_build_ld_units.R` -- Stage 1 (GDS, LD decay, `ld_complexity_
  reduction()` at rho=0.50, GCTA kinship on the stage-1-pruned basis). No
  association test (kept a separate stage on purpose, same reasoning
  `module_sim_3sp53/R/02_bundle.R` already gives for that split).
- `R/03_emmax.R` -- the three primary methods (`emmax_snp`, `emmax_simes`,
  `emmax_consensus`) plus the `emmax_snp_nonsingleton` diagnostic (same
  unrestricted marker BH decision as `emmax_snp`, singletons excluded from
  the scored set only -- NOT the old `*_snp_clustered` definition, which
  granted a marker truth credit from its whole enclosing unit). Stage 2
  (region assembly) is not run -- instructions: it cannot affect a p-value
  or the BH decision and "is not needed for the primary simulation score."
- `R/05_score_truth.R` -- Va/detectability (`flag_true_qtns()`, MAF>0.10 +
  >=5% per-chromosome Va share), decay-relative match thresholds
  (`score_thresholds()`, rho_r2=0.75/rho_d=0.95), and a from-scratch
  hypothesis-level precision/recall scorer -- NOT the package's own
  `classify_ors()`/`evaluate_ors()`, which implement a *different*,
  dedup-neutral region-level metric (each region claims at most one QTN;
  duplicate claims dropped from both TP and FP) used elsewhere in this
  project for exploratory C-score work. The instructions explicitly reject
  that shape here ("do not silently discard duplicate links... from this
  denominator"), so precision/recall are computed directly per spec:
  precision = truth-linked significant hypotheses / all significant
  hypotheses (no dedup); recall = unique detectable QTN recovered / all
  detectable QTN; conditional recall = same numerator / detectable QTN
  covered by >=1 eligible Stage-1 unit (a coverage diagnostic).
- `R/gate2_verify_marker_order.R` -- Phase 3 gate 2 (PASS): confirms
  identical marker order/record counts through parsing -> Stage 1/GRM ->
  association -> truth scoring for one combination
  (`bgs/V0.5_c2/rep1/env1`). Report: `qc/gate2_validation_report.txt`,
  scores: `qc/gate2_truth_scores_bgs_V0.5_c2_rep1_env1.tsv`.
- `qc/offset_propagation_check.R` -- PK: "Don't we only have to fix what
  has to do with TP/FPs LD/distance to QTN... would only require one to
  update the map data, no need to rerun the EMMAX/LFMM analyses?" Checked
  directly rather than assumed: **no** -- confirmed the parser fix requires
  a full Stage-1/GRM/association rerun, not just a truth-column remap.
  Genotypes are unaffected by the bug, but a QTN's ~100kb+ position shift
  moves its sort-order position relative to many neighbouring markers,
  which `ld_complexity_reduction()`'s LD-decay windowing uses -- Stage-1
  clustering (and hence the GRM basis, `stage1_pruned`) genuinely differs
  between old and corrected parsing. Result on `bgs/V0.5_c2/rep1/env1`:
  GRM differs (max abs diff 1.15e-02, 6534 vs 6569 markers selected for
  kinship); `emmax_snp`'s p-value differs for ALL 17,713 markers common to
  both parses (max abs diff 3.15e-02) -- not just markers near a QTN,
  since EMMAX's kinship correction is a global adjustment. Log:
  `qc/offset_propagation_check.txt`.
- `R/run_combo.R` -- chains parse -> Stage 1/GRM -> EMMAX -> truth-scoring
  for one combination; found and fixed a real bug immediately on first use
  (see below).
- `run_gate3_grid.sh` -- Phase 3 gate 3 (PASS): all 7 cells x both tags,
  rep=1/env=1 (14 combinations). 14/14 complete, all 4 methods present, all
  expected summary fields present. Report: `qc/gate3_validation_report.tsv`.
  Two real bugs found and fixed running this at 14-combo scale for the
  first time (neither visible at gate 2's single-combination scale):
  1. **STAGE-clobbering.** `02_build_ld_units.R`/`03_emmax.R`/`05_score_
     truth.R` each set a top-level `STAGE <- "..."` used inside their own
     functions as a free variable resolved at CALL time. Invisible when
     each script runs in its own fresh session (Phase 1/gate 2), but
     `run_combo.R` sources all four into ONE process -- by the time any
     function is actually called, every script's STAGE assignment has
     already run, leaving the LAST-sourced value for all of them. Found
     immediately: `build_ld_units()` wrote its output under
     `05_score_truth/`'s directory instead of its own. Fixed by making
     STAGE a local binding inside each function.
  2. **Zero-QTN runs.** A raw run can save literally zero QTN loci at all
     (each run saves only a small, variable subset of the reference map's
     101 QTN -- see `qc/offset_propagation_check.R`'s header). When that
     happens, `qtn_ld_table()` has nothing to compute against and returns
     a genuinely columnless empty table, not just an empty-but-correctly-
     shaped one -- `qtn_lut[r2 > ...]` then errors instead of matching
     zero rows. Guarded in `05_score_truth.R`.
- `run_full_grid.sh` -- Phase 3 gate 4 (PASS): the full manifest, 7 cells x
  2 tags x 10 reps x 10 envs = 1,400 combinations, through the complete
  primary-EMMAX pipeline. Missing raw inputs are treated as failures, not
  silent skips (checked explicitly before each combo runs). Launched as a
  monitored background job (2026-09-09 19:07-~19:50, concurrency 7) --
  **1,400/1,400 combinations complete, 0 failures** (14 already done from
  gate 3, 1,386 newly run). Independently cross-checked against the actual
  per-combo output files (not just the manifest's own shell-exit status):
  `qc/full_grid_validation_report.tsv` confirms all 1,400 have a valid
  `truth_scores.rds` with all 4 methods and all expected fields present.
  `qc/validate_full_grid.R` is the check script.
  `qc/full_grid_manifest.tsv` is a tracked copy of the run manifest
  (`out_final_v1/full_grid_manifest.tsv` itself is gitignored, per
  `out_final_v1/`'s regenerable-output convention) -- satisfies the
  instructions' "Final outputs" item 11 (`results/run_manifest.tsv`) for
  now; will move under `results/` once that directory exists.

  **First pooled look** (raw hypothesis-level pooled counts across the
  whole grid, `qc/full_grid_pooled_summary_preview.tsv` -- NOT yet the
  real primary result: no cluster-bootstrap uncertainty, no map/burn-in
  pairing, this is a sanity check, not a reportable number):

  | method | n_significant | TP | FP | pooled precision |
  |---|---|---|---|---|
  | `emmax_snp` | 36,027 | 6,274 | 29,753 | 0.174 |
  | `emmax_snp_nonsingleton` | 35,007 | 6,190 | 28,817 | 0.177 |
  | `emmax_simes` | 5,016 | 974 | 4,042 | 0.194 |
  | `emmax_consensus` | 4,483 | 961 | 3,522 | 0.214 |

  Precision rises monotonically size-conscious-method -> Stage-1 Simes ->
  Stage-1 consensus, exactly the qualitative pattern the reanalysis's
  central question asks about (does phenotype-blind Stage-1 LD complexity
  reduction improve precision relative to unrestricted marker-wise
  testing) -- encouraging, but this raw pooled-count preview is not the
  primary result: it doesn't yet use the map-cluster bootstrap the
  instructions require for uncertainty, and recall isn't shown here since
  it needs unique-QTN deduplication across the whole grid, not a per-
  combo sum.

- `R/06_summarise.R` -- the real pooled result (supersedes the raw preview
  above). Point estimates are ratios of pooled counts, never means of
  per-combo precision/recall (`precision = sum(TP)/(sum(TP)+sum(FP))`,
  `recall = sum(n_recovered)/sum(n_detectable_qtn)`, "unique QTN" scoped
  WITHIN each combo -- QTN identity isn't comparable across reps, each is
  an independent simulated genome with its own reference map). Does NOT
  use `R/05_pool.R`'s across-environment SE -- the ten environmental
  continuations for a map share its burn-in, not ten independent
  replicates (see `R/14_random_removal_control.R`'s own confirmation of
  this against the `.ini` files). Primary CI: cluster-bootstrap the ten
  map/burn-in IDs (B=2000), retaining all ten paired environmental
  continuations per resampled map, all 4 methods paired within each
  replicate so contrasts vs `emmax_snp` are genuine paired comparisons.
  Sensitivity: a crossed two-way bootstrap (reps and envs resampled
  INDEPENDENTLY, Cartesian product), run at the grand-pooled level, per
  instructions' "if time permits" framing for this specific check. Both
  bootstraps are vectorised as matrix algebra over pre-aggregated sums
  (never re-subsetting 1,400 rows per replicate) -- full run (56 strata x
  2000 replicates x 4 methods, plus the crossed sensitivity) takes well
  under a second.

  **Grand-pooled result** (all cells/tags together, `results/simulation_
  bootstrap_sensitivity.tsv`): the primary (map-cluster) CIs for
  `emmax_snp` (precision 0.174, CI [0.156, 0.194]) and `emmax_consensus`
  (0.214, CI [0.204, 0.225]) do NOT overlap -- a real, defensible
  precision improvement from phenotype-blind Stage-1 clustering at the
  grand-pooled level. The crossed sensitivity bootstrap widens every CI
  substantially (`emmax_consensus`: [0.177, 0.281]) without reversing the
  ordering -- conclusions are not an artefact of conditioning on the ten
  observed environmental surfaces, though the crossed CIs do overlap more,
  as expected once the environmental axis is also treated as sampled.
  Recall shows the expected trade-off in the other direction (`emmax_snp`
  0.216 > `emmax_consensus` 0.176).

  **Per-cell result is genuinely heterogeneous** (`results/simulation_
  method_contrasts.tsv`, 56 rows) -- NOT a uniform "Stage-1 always helps."
  Some cells show a large, CI-excludes-zero precision GAIN for
  `emmax_consensus` over `emmax_snp` (`nobgs/V0.5_c1`: +0.226, CI [0.143,
  0.309]; `nobgs/V1_c1`: +0.185 [0.102, 0.240]); others show a genuine,
  CI-excludes-zero LOSS (`nobgs/V1_c1.5`: -0.076 [-0.138, -0.024];
  `nobgs/V2_c1.5`: -0.107 [-0.193, -0.038]); several show no detectable
  difference. Kept as-is, not smoothed into the grand-pooled number --
  per instructions, "if some regimes fail, that is part of the method's
  operating range and belongs in the result."

  `results/simulation_performance.tsv` (56 rows, one per cell x tag x
  method): the primary output table -- pooled TP/FP/FN, precision/recall
  with CIs, coverage, test counts. `results/simulation_summary_full.rds`
  caches the full pooled/point/contrast/sensitivity objects plus the raw
  5,600-row (1,400 combos x 4 methods) per-combo table, for rescoring
  without rereading 1,400 individual files.

- `R/04_lfmm.R` -- the LFMM portability check (`lfmm_snp`, `lfmm_simes`
  only, no consensus, per instructions). Full 1,400-combination grid run
  (2026-09-09, `run_lfmm_grid.sh`): 1,400/1,400 complete, 0 failures.
  `R/06_summarise.R` gained a second `bootstrap_arm()` call for
  `LFMM_METHODS`, with its own within-engine reference (`lfmm_simes` vs
  `lfmm_snp`, NOT vs `emmax_snp` -- LFMM is "retain only as a portability
  analysis", not a candidate for the primary EMMAX contrast). Writes
  `results/simulation_performance_lfmm.tsv` and `results/simulation_
  lfmm_portability_contrast.tsv`. Grand-pooled: the same qualitative
  Stage-1-helps-precision pattern survives the engine swap (`lfmm_snp`
  0.134 -> `lfmm_simes` 0.150), though LFMM runs systematically
  higher-recall/lower-precision than EMMAX at matched restriction (e.g.
  unrestricted: 0.134/0.325 vs EMMAX's 0.174/0.216) -- expected from
  `lfmm2`'s genomic-control correction behaving less conservatively than
  EMMAX's kinship correction here, not a new finding.

- `R_figures/figure_simulation_performance.R` / `figureS_simulation_
  absolute_performance.R` -- the two figures named in the instructions'
  "Final outputs" list (paired change vs `emmax_snp`, all 7 cells; and
  absolute precision/recall, all 4 EMMAX methods).

- `R_figures/manhattan_example_data.R` (shared helper, not a standalone
  figure) + `R_figures/figureS_simulation_manhattan_example.R` /
  `figureS_simulation_manhattan_tpfp.R` -- two illustrative supplementary
  Manhattan figures, NOT named in the instructions' output list, built on
  request. Both concatenate all 10 reps of ONE representative combo
  (`nobgs/V0.5_c1`/env=3 -- the only one of the 10 environments where
  every rep has >=1 significant EMMAX marker, found by a one-time scan)
  into a 20-"chromosome" illustrative genome (rep r's Chr1/Chr2 become
  chromosome 2r-1/2r, odd=real/even=near-neutral), EMMAX top / LFMM
  bottom row. `..._example.R` colours every marker by which Stage-2
  assembled outlier region it physically falls in (illustrative only --
  Stage 2 is not what precision/recall are scored from); `..._tpfp.R`
  instead colours by TP/FP status of the significant Stage-1
  emmax_simes/lfmm_simes unit a marker belongs to, reproducing
  `05_score_truth.R`'s hypothesis-level truth-linkage logic exactly (same
  primitives/PARAMS) -- i.e. precisely what IS counted towards the pooled
  precision/recall above, visualised directly.

Not implemented: Stage 2 / reported-region output as a scored quantity
(Stage 2 remains illustrative-only, per instructions), the truth-threshold
sensitivity grid, the BGS validation figure/table, or any final
table/manuscript macro.

`LFMM_K=5`'s justification stop point IS resolved (2026-09-09, PK): the 80
sampled populations fall into 5 discrete spatial groups (4 grid corners +
1 centre) BY DESIGN -- "deliberate to avoid the discussion of what K to use
if the samples were... randomly sampled across the landscape." This is
verifiable from population (x,y) coordinates alone, independent of
genotypes or phenotype (exactly the instructions' "structure diagnostic
independent of association truth"), and was already confirmed
independently during this project's structured-null work (`kmeans(k=5)`
on coordinates recovers the identical grouping). See `00_config.R`'s
`LFMM_K` comment. LFMM was subsequently built and run in full (see
`R/04_lfmm.R` above).

## Output roots (new, none overlapping old `module_sim_3sp53` paths)

- Parsed bundles: `/Volumes/Large_storage/module_sim_3sp53_parsed_final_v1/`
  (mini-local, mirrors the old `..._parsed` root's location on the same
  volume, new name).
- `final_analysis/qc/` -- QA/validation reports, git-tracked, meant to be
  human-reviewed directly.
- `final_analysis/out_final_v1/` -- per-stage, per-combination outputs
  (`02_build_ld_units/`, `03_emmax/`, `05_score_truth/`, each keyed by
  `<tag>_<cell>_rep<rep>_env<env>/`, with a `_receipt.rds` alongside every
  stage's `.rds` output for the caching/rescoring machinery). Gitignored --
  regenerable from raw NEMO output + this pipeline, same convention as
  `module_sim_3sp53/out/` (never committed; only a later `results/`
  equivalent will be). Reviewable copies of gate-2's own single-combination
  output are saved under `qc/` instead (small, git-tracked).

Raw NEMO output (`/Volumes/Large_storage/prod_out`, `~/LDscnR-NEMO/
params_3spC_53cM/`) is read-only throughout; nothing here writes to it.
