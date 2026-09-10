# final_analysis -- clean manuscript-simulation pipeline (in progress)

Rebuild of `module_sim_3sp53` per `../CLAUDE_REANALYSIS_INSTRUCTIONS.md`, kept in
its own subtree so the existing exploratory work in `module_sim_3sp53/` is left
untouched (see that document's "Preserve existing work" section -- in
particular `R/14_random_removal_control.R` there currently has an uncommitted
change that must not be overwritten, moved, restored, or reformatted).

## Status: Phase 1 complete; Phase 2 complete under the 2026-09-10 revised
## primary estimand -- STAGE-2 ASSEMBLED REGIONS are now the primary
## scoring/reporting unit for every method, including the unrestricted
## marker-wise comparator (emmax_snp_region/lfmm_snp_region, made
## comparable via post hoc Stage-2 assembly of every discovered
## phenotype-blind cluster). Marker-/Stage-1-unit-level scores are
## retained as diagnostics only. Phase 3 gates 1-4 ALL PASS -- the full
## 1,400-combination grid has pooled region precision/recall with
## bootstrap CIs, and 9 figures are built (2 named in the instructions, 7
## illustrative/supplementary). The BGS validation, truth-threshold
## sensitivity grid, and both remaining "Final outputs" tables
## (simulation_qc.tsv, simulation_truth_sensitivity.tsv) are also done
## (2026-09-10, per ~/gitlab/LDscnR_manuscript/AUDIT.md's outstanding-
## items list). Remaining: any final table/manuscript macro.

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

  **PRIMARY ESTIMAND REVISED 2026-09-10 (PK)**: Stage-2 ASSEMBLED REGIONS,
  not Stage-1 units or markers, are the primary scoring/reporting unit --
  see `05_score_truth.R`'s entry below for the full rationale and the real
  bug this fixed. The grand-pooled numbers below are the CURRENT
  (region-level) results; the marker-/unit-level numbers this section
  reported before 2026-09-10 are superseded and now live only in
  `results/simulation_performance_diagnostic.tsv`.

  **Grand-pooled result** (all cells/tags together): `emmax_snp_region`
  precision 0.124 -> `emmax_simes_region` 0.163 -> `emmax_consensus_region`
  0.175 -- monotonic, matching the central question, with `emmax_snp_region`
  now a genuinely comparable reported-call baseline (post hoc Stage-2
  assembly of every discovered phenotype-blind cluster, not a raw
  significant-marker count). Recall shows the expected trade-off in the
  other direction (0.230 -> 0.190 -> 0.176). LFMM region arm: `lfmm_snp_
  region` 0.110 -> `lfmm_simes_region` 0.129 precision, 0.341 -> 0.296
  recall -- same qualitative pattern survives the engine swap.

  **Per-cell result** (`results/simulation_method_contrasts.tsv`, 42
  rows): precision CI excludes zero (a real gain) in 7/14 strata for
  `emmax_simes_region` vs `emmax_snp_region`, with no CI-excludes-zero
  LOSSES -- e.g. `nobgs/V1_c1`: +0.343 [0.262, 0.418]; `nobgs/V0.5_c1`:
  +0.264 [0.194, 0.327]. Recall CI excludes zero (a real cost) in 11/14
  strata, consistently negative. This per-cell picture is now MORE
  consistent/one-directional than the old unit-level story was -- region
  assembly does not just fix the unrestricted comparator's scoring
  artefact, it changes the substantive per-cell conclusion. Kept as-is,
  not smoothed into the grand-pooled number -- per instructions, "if some
  regimes fail, that is part of the method's operating range and belongs
  in the result."

  `results/simulation_performance.tsv` (42 rows, one per cell x tag x
  region method): the PRIMARY output table -- pooled TP/FP/FN region
  counts, region precision/recall with CIs, coverage, test counts.
  `results/simulation_performance_diagnostic.tsv` / `simulation_
  diagnostic_contrasts.tsv`: the marker-/Stage-1-unit-level DIAGNOSTIC
  arm (`emmax_snp`, `emmax_snp_nonsingleton`, `emmax_simes`,
  `emmax_consensus`, plus `lfmm_snp`/`lfmm_simes`), retained only to show
  why region assembly is necessary, never the primary estimand.
  `results/simulation_performance_region.tsv` / `simulation_region_
  granularity_contrast.tsv`: paired region-vs-unit contrasts for the SAME
  method (same bootstrap draw) -- the direct evidence for that "why."
  `results/simulation_summary_full.rds` caches everything, including the
  raw per-combo table, for rescoring without rereading 1,400 files.

- `R/04_lfmm.R` -- the LFMM portability check (`lfmm_snp`, `lfmm_simes`
  only, no consensus, per instructions). Full 1,400-combination grid run
  (2026-09-09, `run_lfmm_grid.sh`): 1,400/1,400 complete, 0 failures.
  `R/06_summarise.R` gained a second `bootstrap_arm()` call for
  `LFMM_METHODS` (now `lfmm_snp_region`/`lfmm_simes_region`, region-level
  since the 2026-09-10 update), with its own within-engine reference
  (`lfmm_simes_region` vs `lfmm_snp_region`, NOT vs `emmax_snp_region` --
  LFMM is "retain only as a portability analysis", not a candidate for
  the primary EMMAX contrast). Writes `results/simulation_performance_
  lfmm.tsv` and `results/simulation_lfmm_portability_contrast.tsv`. See
  the grand-pooled numbers above.

- `R_figures/figure_simulation_performance.R` / `figureS_simulation_
  absolute_performance.R` -- the two figures named in the instructions'
  "Final outputs" list (paired STAGE-2 REGION change vs `emmax_snp_region`,
  all 7 cells; and absolute Stage-2 region precision/recall, all 3 primary
  region methods -- updated 2026-09-10 for the revised primary estimand).

- `R/05_score_truth.R` -- Stage-2 assembled regions became the PRIMARY
  scoring unit here (2026-09-10, PK, after inspecting the illustrative
  Manhattan figures below). `emmax_snp_region`/`emmax_simes_region`/
  `emmax_consensus_region`/`lfmm_snp_region`/`lfmm_simes_region` each
  score one WHOLE Stage-2 region as one hypothesis. Two seeding routes
  (instructions' "Testing unit versus reported-call unit"): Simes/
  consensus seed Stage 2 with their BH-significant Stage-1 UNITS, as
  before; the unrestricted marker-wise comparators instead mark EVERY
  phenotype-blind Stage-1 cluster -- INCLUDING SINGLETONS -- as
  discovered when it contains >=1 BH-significant marker, then run the
  SAME assembly, so `emmax_snp_region`/`lfmm_snp_region` are genuinely
  comparable reported-call units, not "every significant SNP" scored
  against "one assembled region." Both routes call `.run_stage2()`,
  which reproduces `ld_outlier_test()`'s own `"stage2_discovered"` branch
  directly (same `cl_sig`/`mk_sig`/`ms_sig`/`sub` construction, same
  `ld_prune_and_eMLG()` call) to get `pr$groups$members` -- the ACTUAL
  constituent discovered-cluster membership -- rather than approximating
  a region's markers from its `[Chr,from,to]` bounds.

  **Real bug fixed by this same update**: the first region-scoring pass
  (still in git history) assigned every marker PHYSICALLY BETWEEN a
  region's bounds to that region, including untested/non-significant
  intervening markers that were never part of any discovered cluster --
  "an intervening, untested marker could otherwise lend truth credit to
  the region." Fixed by using `pr$groups$members` throughout. Validated
  before the full rescore on the illustrative `nobgs/V0.5_c1/env3` combo,
  all 10 reps x 5 methods (`qc/validate_stage2_v2_PK.R`, not committed):
  every region's members are an exact subset of its seed clusters' union,
  every seed cluster's members land in exactly one output region (never
  split across two) -- `subset_ok`/`partition_ok` TRUE and
  `clusters_split == 0` in all 50 checked rows.

  Marker-/Stage-1-unit-level scores (`emmax_snp`, `emmax_snp_
  nonsingleton`, `emmax_simes`, `emmax_consensus`, `lfmm_snp`,
  `lfmm_simes`) are retained as DIAGNOSTICS ONLY, to show why region
  assembly is necessary -- never the primary estimand.

- `R_figures/manhattan_example_data.R` (shared helper, not a standalone
  figure) + five illustrative supplementary Manhattan figures, NOT named
  in the instructions' output list, built on request. All five concatenate
  all 10 reps of ONE representative combo (`nobgs/V0.5_c1`/env=3 -- the
  only one of the 10 environments where every rep has >=1 significant
  EMMAX marker, found by a one-time scan) into a 20-"chromosome"
  illustrative genome (rep r's Chr1/Chr2 become chromosome 2r-1/2r,
  odd=real/even=near-neutral), EMMAX top / LFMM bottom row:
  - `figureS_simulation_manhattan_example.R` -- every marker coloured by
    which Stage-2 assembled outlier region it physically falls in
    (purely visual grouping, not a truth-scoring claim -- unlike the
    scoring below, span-based colouring is fine here).
  - `figureS_simulation_manhattan_tpfp.R` -- coloured by TP/FP of the
    significant Stage-1 emmax_simes/lfmm_simes UNIT a marker belongs to
    (diagnostic-level scoring, for comparison against the region-level
    figure below).
  - `figureS_simulation_manhattan_tpfp_unrestricted.R` -- same TP/FP
    colouring for the UNRESTRICTED marker-wise engines (emmax_snp/
    lfmm_snp) instead -- shows the FP freckling Stage-1 restriction
    removes (this combo, diagnostic level: TP 1331/FP 394 restricted vs
    TP 667/FP 1052 unrestricted).
  - `figureS_simulation_manhattan_tpfp_stage2.R` -- coloured by TP/FP of
    the whole STAGE-2 REGION a marker's constituent cluster belongs to
    (the PRIMARY scoring, `score_stage2_tpfp()`, using the same
    corrected `.run_stage2()` logic as `05_score_truth.R`). Built
    because the unit-level figure can show mixed TP/FP inside a single
    visual peak when Stage 2 merges a truth-linked unit with an adjacent
    non-linked one (PK) -- this fixes that. The subtitle reports TP/FP
    REGION counts (this combo: 45 TP / 24 FP regions, EMMAX+LFMM
    combined), not coloured-member-marker totals, per instructions
    ("print TP/FP region counts prominently").
  `manhattan_example_data.R`'s `.truth_linkage_for_rep()` reproduces
  `05_score_truth.R`'s truth definition exactly (same primitives/PARAMS)
  throughout.

- `R/08_bgs_validation.R` (2026-09-10, per `~/gitlab/LDscnR_manuscript/
  AUDIT.md`'s outstanding-items list) -- the one compact BGS validation.
  Computed BOTH pre- and post-MAF-filter (new `parse_nemo_run(...,
  keep_prefilter=TRUE)` option, purely additive to that function's return
  shape). CONFIRMED, not just disclosed: post-filter He alone is
  BACKWARDS (bgs looks more diverse, +0.020 overall, CI excludes 0). The
  pre-filter marker-RETENTION ratio is the real, validated signature --
  bgs keeps ~17% fewer segregating markers overall, steepest at low
  recombination (ratio 0.738 vs 0.931 at high), the classic BGS
  signature -- reproducing a MAF-trim measurement trap this project had
  already flagged elsewhere.`figureS_simulation_bgs_qc.R` is the compact
  figure (retention ratio | post-filter He, side by side, so the trap is
  shown, not hidden).
- `R/09_truth_sensitivity.R` -- rescores truth (never association
  statistics) across the prespecified Va-share threshold grid
  (0.01-0.20). Cheap by construction: Stage-2 region assembly doesn't
  depend on the truth threshold, so it is computed once per combo and
  reused across all 5 thresholds. The primary 5% row reproduces the
  grand-pooled numbers exactly; method ranking has zero crossovers across
  the whole grid -- the primary threshold choice is not doing special
  work. `figureS_simulation_truth_sensitivity.R` is the supplementary
  figure.
- `R/10_export_qc.R` -- assembles `results/simulation_qc.tsv` (parser
  checks, input counts, compact BGS validation headline numbers) entirely
  from already-committed small artifacts, no new computation.
- Minor fix (AUDIT.md's code-audit note): `01_parse_nemo.R`'s `counts` QA
  table had `n_monomorphic` mislabeled under its `n_QTN` column;
  `n_monomorphic` is now its own column. `counts` turns out to be fully
  dead code (computed, never returned or printed) -- left as-is per the
  instructions' "correct the label," not asked to remove it.
- `R/00_config.R` gained host-aware path resolution (`.first_existing()`)
  for `raw_nemo`/`raw_recmap_dir`/`raw_env_dir`: the same raw NEMO output,
  recmaps and env files also live on Petri's laptop (verified byte-
  identical to the mini's copy) at different mount points -- picks
  whichever exists on the current host. Most of this pipeline now runs on
  either machine; `out_final_v1`/parsed bundles remain mini-only for now.

Not implemented: any final table/manuscript macro.

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
