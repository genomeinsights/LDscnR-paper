# Simulation half — six analysis blocks

Granularity fixed by PK to match the panel half: **six blocks, not per-claim or
per-script mappings**. A block is the smallest set of analyses a reader can check
**without opening a second file** — so where a claim needed two files to verify,
those two files belong together.

The two halves do **not** contain the same six, and should not. This side has no
decay/recombination counterpart, and the region work is inseparable from the
merge-distance sweep because the same measurement both characterises regions and
overturned the sweep's interpretation. What matches is the granularity.

Built backwards, from `synthesis/sec_sims.tex` to code.

---

## Claim-to-block map

Every claim lands in exactly one block; no block is empty; no claim needs two.

| Block | Claims (sec_sims.tex) |
|---|---|
| `01_units_and_scan` | mechanism is isolation; mechanism is not a filter; stage-2 clusters vs SNPs; five summarisation rules; $F_\beta$ crossover; engine replication; mechanism survives the degenerate regime |
| `02_nulls` | structure-aware nulls and what breaks them; what an anti-conservative null is good for; pre-flight REML diagnostic |
| `03_regions_and_merging` | precision by span quintile; region-level precision not comparable across partitions |
| `04_parameters` | full factorial grid; stability score; joint admissibility; four geometry filters rejected; minimum cluster size |
| `05_calibration` | neutral chromosomes as FP control; locating the miscalibration; tail calibration explains the cells; kinship effective rank |
| `06_validity` | **dependency block** — no claims; holds the truth definitions and scoring conventions the other five depend on |

### Block types

The rule is: **a block either holds claims or is an explicit shared dependency,
and must be labelled as one or the other.** What is forbidden is a block that is
neither. `01`–`05` here are claim blocks; `06` is a dependency block, holding the
dedup convention, the QTN-tagging rule and the neutral-chromosome definition.
Folding it into `05` would not remove it — it would make `05` the quiet owner of
conventions five blocks rest on, which is worse than owning them openly.

**The two halves' `06`s are different kinds.** The panel half's `06` is a claim
block (region locking). This one is a dependency block. Same number, same
granularity, different type — stated here rather than left for a curator to
notice.

---

## Blocks and their contents

| Block | Scripts | Reads |
|---|---|---|
| `01_units_and_scan` | `isolation_mechanism.R`, `floor_sweep.R`, `snp_vs_cluster_dedup.R`, `cluster_summary_test.R`, `fbeta_tradeoff.R`, `mechanism_robustness.R`, `engine_x_statistic.R`, `test_then_cluster.R`, `stage1_vs_stage2_units.R` | bundles |
| `02_nulls` | `structure_null.R`, `null_tail_calibration.R`, `region_locking_check.R`, `mixed_basis_null.R`, `env_noise_null.R`, `vg_share_calibration.R` | bundles |
| `03_regions_and_merging` | `occupancy_and_merge_decomposition.R`, `span_matched_precision.R`, `chaining_vs_dcap.R`, `region_definition_comparison.R` | bundles, `chaining_vs_dcap.csv` |
| `04_parameters` | `sensitivity_grid.R`, `sensitivity_grid_analyse.R`, `admissible_region.R`, `operating_points.R`, `ldw_dcap_coupling.R`, `filter_then_test.R`, `filter_gain_audit.R` | bundles, `sensitivity_grid/panels/` |
| `05_calibration` | `neutral_chr_control.R`, `kinship_rank_diagnostic.R`, `env_structure_alignment.R`, `signal_fragility.R` | bundles, `structure_null/` |
| `06_validity` | `decay_window_test.R`, `proxy_T1_local_decay.R`, `proxy_T1_by_window.R`, `cap_cost.R` | bundles |

**34 of 90 scripts reached; 56 not load-bearing.**

## Data the trace needs

`regen_sim_data_bgs5/` (800 bundles, 9.9 GB) — every script reads these.
Then `results/sensitivity_grid/panels/`, `results/structure_null/` (incl.
`cl_chrtype.rds`), and four CSVs: `chaining_vs_dcap`, `engine_x_statistic*`,
`ttc_full_*`, `span_full`. The other ~360 result files are superseded.

---

## Findings from building it

**1. Two load-bearing numbers had no committed script.** The occupancy-by-span
table and the true-positive decomposition of the merge sweep were scratch files,
and the manuscript's *withdrawal* of the consolidation claim rested on them. Now
`occupancy_and_merge_decomposition.R` (`00f58ca`). The panel session predicted
this gap before I looked: the numbers most likely to be missing are the ones
computed quickly in answer to a question.

**2. One of my own three "superseded" findings was wrong.** I recorded the
engine-replication claim as citing a 10-panel run. It cites the sensitivity grid,
which was 80 panels × 27 cells all along, and the figures in the text match the
committed output exactly. **An inventory can be wrong in the direction of alarm
as well as omission.**

**3. Where the text and the code named different quantities.** Applying the panel
session's rewrite test — write the script from scratch and the mismatch has
nowhere to hide — the claim comparing single-SNP against cluster-level flags used
a "51 of 67 / 7 of 10" construction that cannot be reproduced from any committed
80-panel output. Restated from `ttc_full_*` in quantities the code computes.

---

## Records are histories, not statements of belief

`doc/development_history.pdf` (23 pp) contains every claim made here including
the withdrawn ones — the C-score, four filters, two transfer functions, the
effective-rank generalisation, the consolidation reading. Withdrawals are marked
in commits and script headers, not systematically in the document.

**The synthesis is the statement of belief. The history is evidence to consult
and never to summarise from.**
