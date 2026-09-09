# final_analysis -- clean manuscript-simulation pipeline (in progress)

Rebuild of `module_sim_3sp53` per `../CLAUDE_REANALYSIS_INSTRUCTIONS.md`, kept in
its own subtree so the existing exploratory work in `module_sim_3sp53/` is left
untouched (see that document's "Preserve existing work" section -- in
particular `R/14_random_removal_control.R` there currently has an uncommitted
change that must not be overwritten, moved, restored, or reformatted).

## Status: Phase 1 complete; Phase 2 partial (primary EMMAX arm only);
## Phase 3 gates 1-2 pass

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

Not implemented: Phase 3 gates 3-4 (14-combination grid; full 1,400-grid),
LFMM (deliberately deferred -- instructions: verify K=5 with a structure
diagnostic before running it, split from the primary EMMAX arm so that can
finish and be audited independently), Stage 2 / reported-region output,
pooling across combinations, the cluster-bootstrap uncertainty machinery, or
any final table/figure. `LFMM_K <- 5L` in `00_config.R` is flagged there as
the still-unjustified legacy value the instructions ask to label as such.

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
