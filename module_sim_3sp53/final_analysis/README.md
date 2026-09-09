# final_analysis -- clean manuscript-simulation pipeline (in progress)

Rebuild of `module_sim_3sp53` per `../CLAUDE_REANALYSIS_INSTRUCTIONS.md`, kept in
its own subtree so the existing exploratory work in `module_sim_3sp53/` is left
untouched (see that document's "Preserve existing work" section -- in
particular `R/14_random_removal_control.R` there currently has an uncommitted
change that must not be overwritten, moved, restored, or reformatted).

## Status: Phase 1 complete; Phase 2 partial (primary EMMAX arm only);
## Phase 3 gates 1-3 pass

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

Not implemented: Phase 3 gate 4 (full 1,400-grid),
LFMM itself (deliberately deferred, split from the primary EMMAX arm so
that can finish and be audited independently), Stage 2 / reported-region
output, pooling across combinations, the cluster-bootstrap uncertainty
machinery, or any final table/figure.

`LFMM_K=5`'s justification stop point IS resolved (2026-09-09, PK): the 80
sampled populations fall into 5 discrete spatial groups (4 grid corners +
1 centre) BY DESIGN -- "deliberate to avoid the discussion of what K to use
if the samples were... randomly sampled across the landscape." This is
verifiable from population (x,y) coordinates alone, independent of
genotypes or phenotype (exactly the instructions' "structure diagnostic
independent of association truth"), and was already confirmed
independently during this project's structured-null work (`kmeans(k=5)`
on coordinates recovers the identical grouping). See `00_config.R`'s
`LFMM_K` comment. This unblocks LFMM's own stop point but does not by
itself put LFMM in scope -- that's still a separate build/run decision.

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
