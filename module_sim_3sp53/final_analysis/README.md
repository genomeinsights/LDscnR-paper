# final_analysis -- clean manuscript-simulation pipeline (in progress)

Rebuild of `module_sim_3sp53` per `../CLAUDE_REANALYSIS_INSTRUCTIONS.md`, kept in
its own subtree so the existing exploratory work in `module_sim_3sp53/` is left
untouched (see that document's "Preserve existing work" section -- in
particular `R/14_random_removal_control.R` there currently has an uncommitted
change that must not be overwritten, moved, restored, or reformatted).

## Status: Phase 1 only (parser correction + validation)

Implemented so far:

- `R/00_config.R` -- paths (all new roots, none overlapping the old
  `module_sim_3sp53` config), design constants (7 cells x 2 tags x 10 reps x
  10 envs), host-path validation that fails loudly on a missing mount.
- `R/01_parse_nemo.R` -- the corrected parser. Nemo locus identifiers are
  one-based; the inherited `idx := as.numeric(idx) + 1` (present in both
  `module_sim/R_parsing/01_parse_nemo.R` and `module_sim_3sp53/R_parsing/
  01_parse_nemo.R`, both still live and both still wrong) is removed. Runs a
  full assertion suite per combination (see `qc/`).
- `R/01_parse_nemo_qa.R` -- the three-run old-vs-corrected validation
  (Phase 3 gate 1): reparses the same three native runs under both the old
  (`+1`) and corrected (`+0`) offset, prints the QTN index mapping for
  human verification, and writes `qc/parser_qa_report.tsv` +
  `qc/three_run_offset_comparison.tsv`.

Not implemented: Phase 2 (`02_build_ld_units.R` onward), Phase 3 gates 2-4,
LFMM, the full 1,400-combination grid, pooling/bootstrap, or any final
figure/table. Stopped deliberately after the three-run QA report, per
instruction.

## Output roots (new, none overlapping old `module_sim_3sp53` paths)

- Parsed bundles: `/Volumes/Large_storage/module_sim_3sp53_parsed_final_v1/`
  (mini-local, mirrors the old `..._parsed` root's location on the same
  volume, new name).
- Everything else this subtree produces: `final_analysis/qc/` and (once
  later phases exist) `final_analysis/out_final_v1/`, both git-tracked
  here rather than on the external volume, since Phase 1's outputs are
  small and meant to be human-reviewed.

Raw NEMO output (`/Volumes/Large_storage/prod_out`, `~/LDscnR-NEMO/
params_3spC_53cM/`) is read-only throughout; nothing here writes to it.
