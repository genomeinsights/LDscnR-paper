# Simulation half — claim-to-code inventory

Built **backwards**, from the 21 `\claim{}` entries in `synthesis/sec_sims.tex`
to the script that produced each number and the cached object it read. Anything
not reached by this trace is dead by construction, so nothing has to be argued
about. A forward catalogue of 90 scripts would have produced a directory listing
rather than the map a referee wants.

Method taken from the panel session's `inventory/PANEL_INVENTORY.md`.

---

## Load-bearing scripts: 21 of 90

| # | Claim (sec_sims.tex) | Script | Reads |
|---|---|---|---|
| 1 | Mechanism: false positives are isolated | `isolation_mechanism.R` | bundles, `structure_null/` |
| 2 | …and the mechanism is not a filter | `floor_sweep.R` | bundles |
| 3 | Testing stage-2 clusters instead of SNPs | `snp_vs_cluster_dedup.R` | bundles |
| 4 | Five summarisation rules compared | `cluster_summary_test.R` | bundles |
| 5 | $F_\beta$ crossover | `fbeta_tradeoff.R` | `filter_then_test/` |
| 6 | Engine replication EMMAX/LFMM | `sensitivity_grid_analyse.R` | `sensitivity_grid/panels/` |
| 7 | Precision by span quintile | `span_matched_precision.R` | bundles |
| 8 | Region-level precision not comparable | `occupancy_and_merge_decomposition.R` | bundles, `chaining_vs_dcap.csv` |
| 9 | Mechanism survives the degenerate regime | `mechanism_robustness.R` | bundles |
| 10 | Four geometry filters rejected | `filter_then_test.R`, `filter_gain_audit.R` | `filter_then_test/` |
| 11 | Minimum cluster size swept | `floor_sweep.R` | bundles |
| 12 | Neutral chromosomes as FP control | `neutral_chr_control.R` | `structure_null/` |
| 13 | Locating the miscalibration | `neutral_chr_control.R` | `structure_null/` |
| 14 | Tail calibration explains the cells | `neutral_chr_control.R` | `structure_null/` |
| 15 | Structure-aware nulls and what breaks them | `structure_null.R`, `null_tail_calibration.R`, `region_locking_check.R`, `mixed_basis_null.R`, `env_noise_null.R` | bundles |
| 16 | Kinship effective rank orders the cells | `kinship_rank_diagnostic.R` | bundles |
| 17 | Pre-flight REML diagnostic | `vg_share_calibration.R` | bundles |
| 18 | What an anti-conservative null is good for | (interpretation of 15, 17) | — |
| 19 | Full factorial sensitivity grid | `sensitivity_grid.R` + `_analyse.R` | bundles |
| 20 | Stability score over the grid | `sensitivity_grid_analyse.R` | `qtn_bins_true.rds` |
| 21 | Joint admissibility of the two thresholds | `admissible_region.R` | bundles |

Supporting, reached but not cited directly: `engine_x_statistic.R`,
`test_then_cluster.R`, `stage1_vs_stage2_units.R`, `chaining_vs_dcap.R`,
`region_definition_comparison.R`, `decay_window_test.R`, `proxy_T1_local_decay.R`,
`proxy_T1_by_window.R`, `cap_cost.R`, `operating_points.R`, `signal_fragility.R`,
`env_structure_alignment.R`, `ldw_dcap_coupling.R`.

**So ~34 of 90 scripts are reached; 56 are not load-bearing.**

---

## Data the trace actually needs

- `regen_sim_data_bgs5/` — 800 bundles, 9.9 GB. Every script reads these.
  Each carries `GTs`, `map` (with `emx_p`, `lfmm_p`, `ld_w_095`, `chr_type`,
  `true_QTN`, `true_pos_QTN`), `LD_decay`, `GRM`, `complexity_reduction$stage1`.
- `results/sensitivity_grid/panels/` — 80 panel objects (1.6 MB)
- `results/structure_null/` — per-panel null objects + `cl_chrtype.rds`
- `results/chaining_vs_dcap.csv`, `results/engine_x_statistic*.csv`,
  `results/test_then_cluster*.csv`
- Nothing else. The other 365 result files are from superseded analyses.

---

## First finding: what the inventory caught

**Two load-bearing numbers had no committed script until this trace was run.**
The occupancy-by-span-quartile table and the true-positive decomposition of the
merge sweep were both computed in scratch during conversation, and the
manuscript's *withdrawal* of the consolidation claim — the strongest correction
in the section — rested entirely on them. Now `occupancy_and_merge_decomposition.R`
(commit `00f58ca`).

The panel session predicted this exact gap before I looked, on the grounds that
the numbers most likely to be missing are the ones computed quickly in answer to
a question. Both of mine were from this morning.

---

## Second finding: three claims cite numbers from superseded runs

| Claim | Cites | Now available at |
|---|---|---|
| 6, engine replication | 10-panel run | 80 panels (`exs_full_*.csv`) |
| 3, precision advantage | `snp_vs_cluster_dedup` | consistent, but predates the stage-1 comparison |
| 7, span quintiles | 6 panels | not re-run at 80 |

None of the directions changed when replicated at 80 panels, but the **numbers
in the text are from the smaller runs** and should be refreshed before
submission rather than at proof stage.

---

## Records are histories, not statements of belief

`module_sim_LDscnR/doc/development_history.pdf` (23 pp) contains every claim made
in this work including those later withdrawn — the C-score, four filters, the
shrinkage transfer function, the feature-extent cap, the effective-rank
generalisation, the consolidation reading. Withdrawals are marked in the commit
log and in script headers, not systematically in the document.

**The synthesis is the statement of belief. The history is evidence to consult
and never to summarise from.** A read-only session cannot run the check that
would settle which version is current.
