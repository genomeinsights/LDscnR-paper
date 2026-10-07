# Additional simulation analyses: audit report

Audit date: 2026-09-17. Revised 2026-09-18 after PK's first-pass review ("Second pass" section). **Revised again 2026-09-18 (later the same day)** after PK's second review, which found two silent-failure bugs in the second-pass code itself ("Third pass" section) -- **all numbers in the "Analysis 1" and "Analysis 2" RESULT sections below are from the fully corrected, complete final rerun and supersede every earlier number in this document**, including the "Second pass" section's own point-2 and point-3 tables (kept below for the historical record of what changed and why, but explicitly marked superseded where relevant).

## Output provenance investigation (2026-09-18, third round)

PK flagged that locally-inspected `simulation_summary_full.rds` and `simulation_stage2_region_details.rds` contained older marker-wise totals than the canonical TSVs -- correct, and traced to ground before building the primary-`c=1` summaries requested in the same message.

**Independent re-aggregation, ground truth.** Every `05_score_truth/*/truth_scores.rds` file (all 1,400) was read directly and pooled from scratch, bypassing every cached summary object: `emmax_snp_region` 6,623 regions / 813 TP / 5,810 FP, `emmax_simes_region` 3,936 / 640 / 3,296, `emmax_consensus_region` 3,450 / 604 / 2,846, `lfmm_snp_region` 11,581 / 1,258 / 10,323, `lfmm_simes_region` 7,810 / 1,011 / 6,799 -- an exact match to `results/simulation_performance.tsv` on the machine where the pipeline actually runs (the mini). All 1,400 input files have mtimes within an 87-second window on 2026-09-17 (15:11:06-15:12:33) -- static, unchanged since the ordering-unification rescore; there is no possibility they were touched again afterward.

**A confusing false alarm along the way, run down before concluding anything.** `R/06_summarise.R`'s own aggregation logic (`load_all_scores()`, a plain sequential loop + `rbindlist()`, no parallelism) was confirmed fully deterministic by inspection. A fresh rerun of it nonetheless appeared to disagree with per-`(tag,cell)`-stratum numbers when diffed via `git diff` **run on the mini**. Chased down rather than assumed benign: the mini has its own independent git clone of this repository that was never `git pull`ed at any point in this multi-day session (only `rsync` was ever used to move files there) -- its `HEAD` was still sitting at `7452266`, the *original pre-session* commit, entirely unrelated to any fix made since. Comparing the mini's (correct, current) working files against the mini's own ancient git history produced a spurious "these disagree" signal that had nothing to do with data correctness. Confirmed by checking `git diff --stat HEAD -- results/` **on the local Mac clone** (the one actually pushed to origin) after pulling fresh copies from the mini: every TSV was already byte-for-byte identical to what's committed -- the ordering-unification fix's numbers were correctly carried through the "Second pass" commit for every TSV output.

**The one genuine gap, narrower than first feared.** Two files *are* tracked in git and *were* stale: `results/simulation_stage2_region_details.rds` and `results/simulation_summary_full.rds`, both last committed in `7452266` -- the very first commit of this module, before the ordering-unification fix, before any bug fix in this document. They were never re-synced and re-committed after later reruns because my own results-sync commands repeatedly used `rsync --exclude='*.rds'` (to avoid pulling large binary files unnecessarily for TSV-only checks) and I never went back to specifically re-pull these two tracked exceptions. R/11's own on-disk RDS and TSV agree with each other exactly on the mini (verified directly: pooled `n`/TP/FP identical for all 5 region methods) and both match the independent truth-scores aggregation -- so **R/11 itself never needed rerunning**; this was purely a local-sync/commit gap on my end, now fixed by pulling the current RDS files fresh and re-committing them (no upstream stage rerun involved).

**Conclusion for the remaining steps:** no rerun of `05_score_truth`, `06_summarise`, `R/11`, or `R/13` is needed -- all are confirmed mutually consistent and correct on the execution host, and now correctly reflected in the committed repository too. The `c=1` primary summary below is built entirely from this now-verified-correct state.

## What was implemented

Three connected analyses, all scored at the **Stage-2 reported-region** scale (never Stage-1 units), sharing one Stage-2-assembly/truth-scoring implementation (`R/helpers_stage2_truth.R`) with the canonical pipeline so none of them can silently diverge from it:

1. **Null-region calibration** (`R/13_region_null_calibration.R`) -- does a structure-aware null's Stage-2 discovery burden predict the known false-positive burden?
2. **Stage-2 region size vs. truth** (`R/11_stage2_region_details.R` + `R/14_summarise_additional_analyses.R`) -- are smaller reported regions more likely false positives, and does that survive an opportunity-effect (genomic-coverage) control?
3. **Floor decomposition** (`R/12_floor_decomposition.R`) -- how much of the precision/recall change comes from filtering small Stage-1 units versus from LD-aggregation itself?

New files: `R/helpers_stage2_truth.R`, `R/11_stage2_region_details.R`, `R/12_floor_decomposition.R`, `R/13_region_null_calibration.R`, `R/14_summarise_additional_analyses.R`, `R_figures/figure_simulation_null_truth_calibration.R`, `R_figures/figure_simulation_floor_decomposition.R`, `R_figures/figure_simulation_floor_decomposition_by_cell.R`, `R_figures/figureS_simulation_stage2_size_truth.R`, plus the one-off full-grid rerun driver `R/rerun_score_truth_grid.R` and orchestration script `run_rerun_sequence.sh`. `05_score_truth.R` was modified only to call the new shared helper instead of its own private copy of the same logic (`region_scoring_version` bumped 2->3 on first pass, then 3->4 on the second pass).

## Which old exploratory results were NOT reused, and why

`module_sim_3sp53/R/12_structured_null.R` and `R/14_random_removal_control.R` (plus their pooled `results/*.rds`) were **not** used as a numeric source for Analysis 1 or 3. Per `module_sim_3sp53/CLAUDE_REANALYSIS_INSTRUCTIONS.md`, those scripts read bundles built by the pre-correction parser, so their counts are not comparable to anything in `final_analysis/`. What *was* reused, as design/method choices rather than numbers:

- The three null-construction recipes (`group`/`mvn`/`spatial`) -- ported verbatim from `R/12_structured_null.R`'s validated draw code, applied to `final_analysis/`'s own corrected-parser inputs.
- The `n_obs > 0`-conditioning discipline for `realised_fdr`/`R_null/obs` -- that module's own hard-won 2026-09-09 correction, carried forward here as a hard requirement (never `max(n_observed, 1)`).
- `R/14_random_removal_control.R`'s matched-cardinality-random-subset design was **read but not rerun**, and the file itself, including its uncommitted local WIP, was **not modified** at any point (explicit standing instruction).
- `results/simulation_stage2_region_details.rds`/`.tsv` (R/11) and everything downstream deliberately does **not** reuse the earlier `module_sim_3sp53` "cluster size vs FP" figures -- those scored Stage-1 units; this scores Stage-2 reported regions, a different estimand.

## Shared Stage-2 helper: bugs found and fixed while building and validating it (first pass)

`R/helpers_stage2_truth.R` refactors the previously-duplicated Stage-2 assembly/scoring logic into `assemble_stage2()` / `score_stage2_regions()` / `stage2_seed_from_units()` / `stage2_seed_from_markers()` / `bootstrap_rep_matrix()`. Building and validating this, and the floor sweep and null calibration built on top of it, surfaced four genuine, previously-latent bugs:

1. **Order-sensitivity of `ld_prune_and_eMLG()`'s input.** Its distance-restricted dynamic cut merges *adjacent* input clusters, so its output partition depends on input order, not just the set. First-pass fix kept the marker-seeded and unit-seeded routes on different orderings to match historical output. **Superseded by the second-pass fix below.**
2. **`unit_id` is not a stable identifier across `size_floor` values.** `LDscnR:::.ld_outlier_units()` assigns `unit_id` before its own internal `setorder(Chr, from)`, so it is a permutation of `1:N`, reassigned from scratch at every call. Fixed by keying on `core_snp` instead of `unit_id` in `R/12_floor_decomposition.R`.
3. **`floor == floor` self-comparison in the grand-pooled contrast loop.** The bare loop variable `floor` collided with the `floor` *column*, silently pooling all floors together. Fixed by renaming to `this_floor`.
4. **`..env` silently failed to resolve** inside `R/14`'s opportunity-effect combo loop, sidestepped by copying loop arguments to distinctly-prefixed locals (`q_tag`/`q_cell`/`q_rep`/`q_env`).

Also fixed: `.genomic_sort_clusters()` originally called `match()` once *per cluster*, re-hashing the full marker universe every time -- fixed to hash once via `LDscnR:::.marker_positions()`.

None of these first-pass bugs affected any existing, already-published result.

## Second pass (2026-09-18): response to PK's first review

PK's first review raised four issues. All four were addressed. **Points 2 and 3's numeric results below turned out to still be wrong** -- silently affected by two further bugs PK caught in the *second* review (see "Third pass"); the numbers here are kept only as a record of what the redesign changed structurally.

### 1. Stage-2 assembly ordering, unified

**Issue:** the first-pass fix kept the marker-seeded and unit-seeded assembly routes on different orderings. PK: reproducing historical output is not sufficient reason for two logically-equivalent paths to behave inconsistently.

**Decision, with PK's explicit sign-off:** a full 1,400-combo side-by-side comparison (native vs. genomic order) was run before committing: aggregate change small (3rd decimal of precision) but non-trivial per-combo (8.5%/16.3% of individual EMMAX/LFMM combos get a different region count). PK confirmed: **apply everywhere, update manuscript macros.**

**Implementation:** `stage2_seed_from_markers()` now calls `.genomic_sort_clusters()`, matching `stage2_seed_from_units()`. `region_scoring_version` bumped 3->4, full 1,400-combo rescore (1.1 min), `06_summarise.R` regenerated canonical results. **This numeric result is unaffected by the Third-pass bugs (which are confined to R/13 and R/14) and stands as final.**

**Measured impact (matches the pre-decision prediction almost exactly, and is the current, final canonical result):**

| method | regions (before) | regions (after) | FP (after) | precision (before) | precision (after) |
|---|---:|---:|---:|---:|---:|
| `emmax_snp_region` | 6,658 | **6,623** | 5,810 | 0.1236 | **0.1228** |
| `lfmm_snp_region` | 11,569 | **11,581** | 10,323 | 0.1104 | **0.1086** |
| `emmax_simes_region` | 3,936 | 3,936 (unchanged) | 3,296 | 0.1626 | 0.1626 |
| `emmax_consensus_region` | 3,450 | 3,450 (unchanged) | 2,846 | 0.1751 | 0.1751 |
| `lfmm_simes_region` | 7,810 | 7,810 (unchanged) | 6,799 | 0.1294 | 0.1294 |

Only the two marker-seeded methods are affected. Relative to the marker-wise baseline, LD-aggregation reduces false-positive *regions* by **43.3%** (Simes: 5,810 -> 3,296) and **51.0%** (consensus: 5,810 -> 2,846). **`values_simulation.tex`'s affected macros have NOT yet been updated** -- no manuscript file has been touched; this table is the input for that.

### 2. Opportunity-effect control replaced with exact per-region relocation -- [superseded, see Third pass]

**Issue:** the first-pass opportunity control drew K=30 generic contiguous windows per (combo, chromosome, size-bin) at the bin's median size, and compared exactly two point-estimate slopes with no uncertainty on their difference.

**Redesign (`R/14_summarise_additional_analyses.R`):** every observed Stage-2 region is relocated by an exact circular shift of its own constituent-marker index-spacing pattern to a uniformly random start point on the same chromosome, R=200 independent replicates, giving a full null distribution of the opportunity-only TP-rate-vs-size slope.

**[!] The numbers first reported here (obs_slope 0.253 vs. null_slope 0.528 [0.4655,0.5960]) were computed on a silently incomplete null -- see "Third pass" point 2 below for the bug and the corrected result (0.253 vs. 0.512 [0.470, 0.556], essentially the same conclusion but now on the complete grid).** The qualitative finding -- that the observed slope is significantly *shallower* than the opportunity-null slope, not steeper -- is unchanged and, if anything, more firmly established on the complete data. See Analysis 2's RESULT section below for the final, authoritative numbers.

### 3. Null-calibration redesigned with paired map identifiers and cell/tag-adjusted correlation -- [superseded, see Third pass]

**Issue:** the first-pass "balanced" design sampled reps independently per tag, so BGS/no-BGS treatments rarely shared a map. PK also asked for the association to be checked after adjusting for demographic cell.

**Redesign:** `R/13`'s balanced-mode combo selection now samples reps **once per cell**, shared across both `tag` values, `N_REPS_PER_CELL <- 5L` -- 7 cells x 5 reps x 2 tags x 10 envs = 700 combos. A cell(+tag)-adjusted Spearman correlation is now reported alongside the raw one.

**[!] The design and code changes described above are correct and final.** But the FIRST run under this redesign only realised 52 of the intended 70 map/burn-in groups, and the audit at the time wrongly attributed this to "gaps in the raw simulation grid." **That diagnosis was wrong** -- PK's second review found a cache-read bug (see "Third pass" point 1) that silently dropped 1,134/4,200 jobs, and confirmed by direct check that the raw grid has no such gaps. The corrected rerun realises all 70/70 groups. See Analysis 1's RESULT section below for the final, authoritative numbers -- the table originally printed here is superseded in full, not just refined.

### 4. Bootstrap pairing corrected in both R/13 and R/14

**Issue:** two bootstraps resampled `(tag, cell, rep)` as the clustering unit, resampling BGS and no-BGS rows independently despite them being a paired map/burn-in draw.

**Fixed** in `R/14`'s `.fit_size_model()` and `R/13`'s own summary bootstrap, both now clustering by `(cell, rep)`. This fix was correct as far as it went (BGS/no-BGS pairing), but incomplete.

**[!] FURTHER CORRECTED (PK, 2026-09-26 third review): `(cell, rep)` is still too fine.** Per `~/gitlab/LDscnR-NEMO/make_prod.sh`'s STEP 1, the genetic maps/QTN positions/environmental values/dispersal template are built ONCE per `rep` and reused identically across ALL 7 cells and both tags -- so rows sharing a `rep` are not independent even across different cells, and any bootstrap pooling multiple cells (R/13's grand-pooled summary; R/14's `.fit_size_model()`, which includes `cell` as a covariate across all 7 cells) needs to cluster by `rep` alone, not `(cell, rep)`. This is also what materials_and_methods.tex:60 already describes for the manuscript's PRIMARY performance numbers ("2,000 bootstrap samples of the ten map/burn-in histories... the three selection strengths were also retained together") -- R/13/R/14's own bootstraps had simply drifted from that already-correct design description. Point estimates in both scripts are unaffected (refit on the same full data regardless of resampling scheme); only the CI columns changed. Rerun summary-only (not the expensive upstream stages) via `R/13b_refresh_null_summary_after_cluster_fix.R` and `R/14b_refresh_size_models_after_cluster_fix.R`, both reading already-saved per-combo/per-region tables rather than recomputing them. R/15's alignment-model bootstrap had the identical `(cell,rep)` bug (inherited by a since-corrected `R/20_fst_alignment_joint_null_model.R` too) -- refreshed the same way (`R/15b_refresh_alignment_models_after_cluster_fix.R`); see the updated numbers throughout the Analysis-4 section below. `results/simulation_stage2_size_models.tsv`'s refreshed CIs are also now quoted correctly in `LDscnR_manuscript/Supplementary.tex:614` (point estimates $-0.628$/$-1.090$ unchanged; CIs updated from $[-0.706,-0.566]$/$[-1.230,-0.967]$ to $[-0.723,-0.559]$/$[-1.204,-0.984]$).

## Third pass (2026-09-18, later the same day): two silent-failure bugs in the second-pass code

PK's second review found the second-pass fixes were logically correct but two **silent-failure bugs elsewhere in the same scripts** had corrupted their outputs -- both confirmed with direct log evidence before fixing, not taken on faith:

### 1. R/13's cache-read returned the wrong object shape

`run_one()`'s cache-hit branch (`stage_stale()` == FALSE) returned the *entire saved list* (`tag`, `cell`, ..., `counts`), while the fresh-compute branch returned `counts` directly. The caller's `mean(counts)` on a list just warns and returns `NA`; `max(counts)` on a list throws `"invalid 'type' (list) of argument"`, caught by the job's own `tryCatch` and silently dropped. **Confirmed in the overnight log: exactly 1,134/4,200 jobs failed this way**, explaining the "52/70 groups" the previous audit wrongly attributed to missing raw data (confirmed by direct check: raw `ld_units.rds` exists for every rep of every affected cell/tag).

**Fix:** cache-hit branch now returns `readRDS(...)$counts`; `NULL_CALIB_VERSION <- 2L` added to `PARAMS` to force every job to recompute cleanly (the on-disk values themselves were never corrupted -- only the reader was -- but a clean version bump removes any doubt); a hard `stop()` if `nrow(dt) != nrow(jobs)`; an explicit design-completeness assertion on the balanced-mode combo table (exactly `N_REPS_PER_CELL` reps/cell, both tags, all 10 envs); and a payload-shape validation (`is.numeric`, `length(counts)==B`, all finite non-negative integers) so a malformed cache entry can't silently produce one bad row that still passes the row-count check.

### 2. R/14's per-combo seed silently dropped ~half the opportunity-relocation grid

`strtoi()` on a full 8-hex-digit CRC32 hash overflows R's signed 32-bit integer for any hash >= `0x80000000` (roughly half of all values) and returns `NA`. `NA` propagated to `set.seed(NA)`, which throws `"supplied seed is not a valid integer,"` caught by the combo's `tryCatch` and silently dropped. **Confirmed: 619/1,223 combos (50.6%) failed this way**, while `observed_by_bin` (used for `obs_slope`) was still built from the *full* grid -- comparing a full-grid observed slope against a roughly-half-grid null slope.

A second, independent issue in the same function: combos with no detectable QTN, or no marker passing the LD/distance criteria, `return(NULL)` -- removing their observed regions (necessarily all false positives) from the null denominator entirely, rather than correctly contributing zero-hit relocation draws.

**Fix:** the CRC32 hash is now split into two 4-hex-digit halves (each safely `< 2^16`) and combined with double-precision arithmetic before ever coercing to a 32-bit integer for `set.seed()` (verified: 0 collisions across 1,223 real combo IDs). `is_qtn_linked` now starts empty and only gets populated when QTNs/linked markers exist, so every relocation draw for a QTN-less combo correctly scores as a miss rather than the combo being dropped. A hard `stop()` was added if the number of combos represented in the pooled relocation table doesn't exactly match `nrow(combos_present)`.

**A related, cosmetic bug** was also found (by self-review, not by PK) while inspecting the validation run: the printed interpretation text only branched on "slope-difference CI entirely positive" vs. an else covering both "includes zero" and "entirely negative" -- so the CI-entirely-negative case (exactly this analysis's actual result) printed the wrong narrative ("does NOT exclude zero"). Only the console/log text was affected, never the numeric TSV outputs. Fixed to a proper three-way branch.

### 3. Additional safeguards PK requested before the corrected rerun

- **R/11**: worker results now wrapped as `{ok, data, error}` so a genuine worker failure (crash, OOM, corrupted read) is distinguishable from `region_details_one_combo()`'s own legitimate `NULL` (a combo with zero reported regions across every method) -- the plain `tryCatch(..., error=function(e) NULL)` pattern made these indistinguishable and a real failure would have silently produced apparently-complete output. Hard `stop()` if any `ok == FALSE`.
- **R/12**: `errs` was already counted and printed but never gated on. Added `if (errs > 0) stop(...)` before `saveRDS()` (`floor_decomposition_one_combo()` has no legitimate `NULL`-return case, so any `NULL` here really is a worker error).
- **`run_rerun_sequence.sh`**: added the by-cell floor figure as its own stage (previously regenerated manually, not part of the automated sequence).
- `null_calib_version` added to R/13's final summary receipt; R/14's slope-difference bounds are now explicitly documented (in code comments and printed text) as a **95% relocation interval**, not a general confidence interval -- they describe variation among relocations of the fixed, observed region set, not sampling uncertainty across simulated maps.
- `mc.cores` raised from 7 to 12 across R/11-R/14 (the mini has 14 physical cores; 2 left for the system, per PK).

**Corrected rerun, launched only after all of the above were independently verified** (unit tests of the seed function, the payload-shape predicate, and the `{ok,data,error}` logic; a real end-to-end run of the fixed R/14 on the complete grid, confirming 1,223/1,223 combos and the same qualitative result as before the fix) and after PK's own independent code check found no remaining logic error. Full sequence (06->11->12->13->14->figures) completed cleanly overnight: **R/13: 4,200/4,200 jobs (design check: 7 cells x 5 reps/cell x 2 tags x 10 envs = 700 combos, all realised), R/14: 1,223/1,223 combos, no failures anywhere.** Total R/13 wall time 594.9 min (~9.9h) at `mc.cores=12` -- slower in wall-clock terms than the (invalid, partial) first attempt because every one of the 4,200 jobs now genuinely computes B=200 null draws from scratch, versus a mix of real computation and near-instant silent failures before.

## Analysis 3: floor decomposition -- RESULT

Full 1,400-combo grid, floor grid `{1, 2, 3, 5, 10, 20}`, on the final ordering-unified pipeline (unaffected by the Third-pass bugs, which are confined to R/13/R/14). **Validation check "floor=1 arm B == arm A" passed for all 1,400 combos.**

**Finding (grand-pooled across all 1,400 combos).** At floor=2, filtering out small Stage-1 units (arm B, floor-filtered marker) changes precision by only **+0.0015** relative to the unrestricted marker-wise baseline (95% CI [-0.0023, 0.0054], spans zero) -- essentially no effect. LD-aggregation (Simes/consensus, arms C/D) relative to floor-filtered marker at floor=2 gives **+0.0384** (Simes) and **+0.0508** (consensus) precision, at a recall cost of **-0.0475** and **-0.0607** respectively (all CIs exclude zero). **This directly answers the analysis question: the precision/recall trade-off this manuscript reports is attributable to combining correlated markers into fewer, better-powered tests (LD-aggregation), not to discarding small Stage-1 units as a filtering step.**

**Correction: filtering-alone is not negligible at every floor.** At floor=10, filtering alone (B vs. A) increases EMMAX precision by **+0.0223** (95% CI [0.0114, 0.0344], excludes zero) -- real but modest, and it comes with a substantial recall cost (-0.0656, CI [-0.0830, -0.0473]) and, per the phenotype-blind floor-choice table, a severe genomic-coverage cost. An earlier draft of this section described "almost no precision gain through floor 10" -- that overstated the null result at floor=2 into a claim about the whole floor range; corrected here.

**By-cell result (`figure_simulation_floor_decomposition_by_cell.pdf`) -- corrected.** An earlier draft of this section claimed the LD-aggregation advantage was "largest in the low-gene-flow/c2 cells." **That was wrong, checked directly against the by-cell data and reversed.** At floor=2, EMMAX Precision x Recall in the `V0.5, c=2` cell (low gene flow) is **0.00322 for floor-filtered marker, 0.00174 for Simes, 0.00073 for consensus** -- aggregation is *worse* than floor-filtering alone in this cell, not better. The consistent aggregation advantage (Simes and consensus both clearly exceeding floor-filtered marker) occurs **primarily in the `c=1` (high-gene-flow) cells** -- all three (`V0.5_c1`, `V1_c1`, `V2_c1`) show Simes/consensus precision x recall clearly above floor-filtered marker; the `c=1.5`/`c=2` cells show weak, inconsistent, or reversed performance.

**Manuscript-hierarchy recommendation (new, per PK):** given this cell-dependent pattern, the manuscript should **not** present a single pooled 1,400-run mean as though it describes one operating regime. Recommended: use the **600 high-gene-flow (`c=1`) simulations** (3 cells x 2 tags x 10 reps x 10 envs) as the primary method-performance summary, and retain the remaining **800 simulations (`c=1.5`/`c=2`)** as structured stress tests demonstrating behaviour under stronger genetic structure -- not pooled into the primary estimate.

Precision for the aggregation arms is U-shaped across the floor grid; recall declines monotonically; Precision x Recall (its own panel, replacing "markers retained," per PK's request) also dips through the middle of the grid, i.e. neither floor extreme is jointly optimal.

See `results/simulation_floor_sweep.tsv`, `results/simulation_floor_contrasts.tsv`, `results/simulation_floor_choice_blind.tsv`, `results/figure_simulation_floor_decomposition_by_cell_data.tsv`, `figure_simulation_floor_decomposition.pdf` (4 panels: tests %, precision, recall, precision x recall), and `figure_simulation_floor_decomposition_by_cell.pdf` (Precision x Recall only, `scales="free_y"`, faceted by cell, labelled by selection intensity/gene flow, ordered gene-flow-high-to-low then selection-strong-to-weak).

## Analysis 1: null-region calibration -- RESULT

**Design (final, corrected).** 7 cells x 5 reps/cell x both tags x 10 envs x 2 methods x 3 schemes, B=200. **All 70/70 nominal map/burn-in groups realised** (the earlier "52/70" was the R/13 cache-read bug, not a data gap -- see "Third pass"). 4,200/4,200 jobs succeeded.

**Two different questions, deliberately kept separate (per PK's explicit framing -- this replaces the first-pass audit's treatment of the between-cell pattern as a "confound" to explain away):**

1. *Within-cell, map-level*: does the null-to-observed ratio predict the exact FDP among maps sharing the same demographic regime? **Answer: no, or only weakly.**
2. *Between-cell*: does the null warn that an entire demographic regime is untrustworthy? **This is not a nuisance confound -- it is a central, and apparently real, finding.**

**Question 1 (within-cell): raw vs. cell+tag-adjusted Spearman rho (first value = EMMAX consensus, second = EMMAX Simes):**

| scheme | role | raw rho | cell+tag-adjusted rho |
|---|---|---:|---:|
| group | candidate null | -0.409 / -0.351 | -0.053 / -0.190 |
| mvn (kinship-matched) | **negative control** | -0.625 / -0.599 | +0.058 / -0.170 |
| spatial | candidate null | **+0.680 / +0.639** | **-0.058 / -0.108** |

`spatial`'s raw positive correlation -- the only scheme with the theoretically expected sign -- is **fully explained by demographic cell**: the cell+tag-adjusted correlation is slightly *negative* for both methods, not merely reduced. `group` and `mvn` show weak, inconsistent, sign-flipping adjusted correlations. **None of the three schemes' null-to-observed ratio is a usable within-cell FDP estimator.**

**Question 2 (between-cell): `spatial`-null ratio and known FDP, pooled by cell (this is the finding that matters):**

| gene flow | cells | spatial null/obs ratio | known FDP |
|---|---|---:|---:|
| high (`c=1`) | `V0.5_c1`, `V1_c1`, `V2_c1` | **0.25 - 0.38** | 0.37 - 0.45 |
| medium/low (`c=1.5`, `c=2`) | `V0.5_c1.5`, `V1_c1.5`, `V2_c1.5`, `V0.5_c2` | **2.77 - 5.29** | 0.80 - 0.98 |

In the `c=1` cells, the spatial null under-predicts the observed region count (ratio < 1) and the true FDP is moderate. In the `c=1.5`/`c=2` cells, the spatial null predicts **2.8x to 5.3x more regions than were even observed** -- structure alone could reproduce or exceed the actual discovery burden -- and the true FDP is very high (80-98%). This is exactly the separation PK's own review predicted as the most likely outcome, and it survives on the complete, bug-corrected data.

**Conclusion for the manuscript, stated the way PK asked it to be stated:**
- The null-to-observed ratio is **not** an FDP estimator.
- Within-cell, map-level calibration is weak.
- The spatial null **may nevertheless identify demographic regimes in which structure alone can reproduce or exceed the observed discovery burden** -- a qualitative trust diagnostic at the regime level, distinct from (and not reducible to) a quantitative within-regime FDP estimate. Do not describe the between-cell separation as a confound that the cell-adjustment "removes" -- removing it answers question 1, not question 2, and question 2 is the more directly useful result for flagging untrustworthy regimes.

`mvn` is labelled `null_role = "negative_control"` (no spatial/environmental signal by construction; its role is to confirm the pipeline doesn't manufacture spurious calibration from kinship structure alone, not to be judged as an FDP proxy). `group`/`spatial` are labelled `"candidate_null"`.

See `results/simulation_null_truth_calibration.tsv` (per-map/burn-in-level rows, 70 unique `(tag,cell,rep)` groups), `results/simulation_null_truth_summary.tsv` (raw and cell+tag-adjusted rho, `null_role`, pooled ratios with `rep`-clustered bootstrap CIs -- see item 4's further correction above; point estimates unchanged), `figure_simulation_null_truth_calibration.pdf` / `_labelled.pdf`.

## Analysis 2: Stage-2 region size and false-positive status -- RESULT

**Descriptive.** FP proportion declines monotonically with both size measures, consistently across all 5 region methods and both BGS treatments: e.g. for `emmax_consensus_region`, 2-marker regions are 90.1% FP (bgs) / 86.5% FP (nobgs), falling to 64.8% FP (bgs) / 22.4% FP (nobgs) for 51+-marker regions -- the decline is real in both treatments, but BGS retains a substantially higher FP rate at large sizes than no-BGS.

**Model.** `FP ~ log2(size) + method + cell + tag`, map/burn-in-cluster bootstrap (B=2,000, clustered by `rep` -- corrected 2026-09-26 from `(cell,rep)`, see item 4): `log2(n_markers)` coefficient **-0.628** (95% CI -0.723 to -0.559), `log2(n_units)` coefficient **-1.090** (95% CI -1.204 to -0.984) -- both stable after controlling for method, cell and BGS treatment; both still clearly exclude zero, same direction as before the fix. Manual VIF confirms span/density are severely collinear with marker count and correctly excluded from the primary model.

**Opportunity-effect sensitivity (final, corrected numbers -- see "Third pass" point 2 for the bug that affected the first version of this result).** Every observed region relocated by an exact circular shift of its own marker-spacing pattern on its own chromosome, R=200 replicates, **all 1,223/1,223 combos with >=1 detectable QTN represented** (the previous run silently dropped 619 of them):

| quantity | value |
|---|---:|
| observed log-odds slope (TP-rate vs. log2 size) | **0.2528** |
| opportunity-null slope, mean [95% relocation interval] | **0.5115** [0.4696, 0.5555] |
| slope difference (observed − null) [95% relocation interval] | **−0.2587** [−0.3027, −0.2168] |
| P(null slope ≥ observed slope), 200 replicates | **1.0** (200/200) |

The interval is a **95% relocation interval**, not a sampling-uncertainty confidence interval: `obs_slope` is a single fixed value from the real, observed regions; the interval describes variation among 200 independent relocations of that same fixed set, conditional on it -- it does not describe uncertainty in the observed slope across simulated maps.

The observed slope is significantly *shallower* than the opportunity-null slope. The excess ratio (observed TP rate / opportunity hit rate) is large everywhere (3.2x-16.1x) but *decreases* monotonically with size -- the relative excess over chance is concentrated in small regions, not large ones, even though large regions still have a higher absolute observed TP rate.

**What this result can and cannot support.** It is a **conditional relocation test**: its interval reflects relocation-to-relocation variability of the observed regions, not map-to-map sampling uncertainty, and it evaluates **Stage-2 reported regions**, not the abstract's own claim about **Stage-1 unit size and QTN enrichment among tested units** -- a related but different estimand. It therefore cannot *directly* invalidate the abstract's Stage-1-unit sentence. What it does support: removing the abstract's mechanistic claim that selection-generated local LD explains size enrichment, while retaining an appropriately scaled, purely descriptive statement (regions with more markers/units have a higher absolute rate of being linked to a causal variant) is the safest course until a Stage-1-unit-scale version of this same check exists.

See `results/simulation_stage2_size_truth.tsv`, `results/simulation_stage2_size_models.tsv`, `results/simulation_stage2_opportunity_null.tsv`, `results/simulation_stage2_opportunity_slope.tsv`, `figureS_simulation_stage2_size_truth.pdf`.

## Manuscript consistency

The canonical `values_simulation.tex` macros are stale relative to the current, ordering-unified full-grid rerun (see "Second pass" point 1's table above for the complete set). Current canonical EMMAX values:

- marker-wise (`emmax_snp_region`): 6,623 regions, 5,810 false positives, precision 0.123.
- Simes (`emmax_simes_region`): 3,936 regions, 3,296 false positives -- a **43.3%** reduction in false-positive regions relative to marker-wise.
- consensus (`emmax_consensus_region`): 3,450 regions, 2,846 false positives -- a **51.0%** reduction.

**The provisional paragraph already placed in the manuscript's conclusion should remain marked as provisional.** In particular, any high- vs. low-gene-flow statement there should not be finalised on the strength of this document alone until PK has reviewed the between-cell null-calibration result and the corrected by-cell floor-decomposition finding above together -- both now point the same direction (high gene flow / `c=1` is the well-behaved regime; low gene flow / `c=1.5`,`c=2` is where structure dominates and LD-aggregation's advantage weakens or reverses), which is a stronger, more specific claim than what the provisional paragraph currently says, but it is PK's call whether and how to fold it in.

## Computational cost (actual)

- 05_score_truth full-grid rescore (ordering fix only): 1,400/1,400, 1.1 min.
- 06_summarise: 2s. R/11: 37s (full grid, `mc.cores=12`, `{ok,data,error}` wrapper). R/12: 651s (~11 min).
- **R/13** (final, corrected, `mc.cores=12`): 4,200/4,200 jobs, **594.9 min (~9.9h)** -- every job now a genuine fresh compute (the version bump invalidated all prior caches), versus the earlier invalid run's mix of real computation and near-instant silent failures.
- **R/14** (final, corrected): 296s (~5 min), 1,223/1,223 combos.
- Figures: ~1-3s each.

## Validation results (mandatory checklist)

1. Floor-1 filtered-SNP == unrestricted-SNP -- **PASS**, all 1,400 combos.
2. Canonical floor-2 Simes/consensus reproduce final results -- **PASS** (unaffected by the ordering fix); marker-wise canonical numbers **intentionally changed**, with PK's sign-off.
3. BH recomputed fresh after marker filtering, never reused -- **PASS**.
4. GRM identical across floors -- **PASS by construction**.
5. Truth scoring uses actual constituent markers only -- **PASS**.
6. Observed Stage-2 counts in the null analysis equal R/11's truth-scoring counts -- **PASS by construction**.
7. No null ratio uses `max(n_observed, 1)` -- **PASS**.
8. Zero-discovery cases retained and counted -- **PASS**, 0 occurred in the final run.
9. Pooled precision/recall from pooled counts, never means of ratios -- **PASS**.
10. Bootstrap resampling respects shared map/burn-in structure -- **PASS, corrected on the second pass, then further corrected 2026-09-26** (clustered by `rep` alone -- the actual shared map/burn-in unit, reused across all 7 cells per `make_prod.sh`'s STEP 1 -- not `(cell,rep)`, which was itself an improvement on the original `(tag,cell,rep)` but still too fine; see item 4).
11. Methods/BGS treatments paired within bootstrap replicates -- **PASS, same fix as item 10.**
12. All figures recreated from saved tables, no model rerun -- **PASS**.
13. Every output has a receipt with inputs/params/version/seed -- **PASS**, including `null_calib_version` now in R/13's final receipt.
14. Smoke-tested before the full run -- **PASS**, both on the second pass and again on the third pass (unit tests of the seed function, payload-shape predicate, and `{ok,data,error}` logic; a full-grid real run of the fixed R/14 before committing to the multi-hour R/13 run).
15. **(New, third pass) Worker failures cannot silently produce apparently-complete output** -- **PASS**: R/11, R/12, R/13, R/14 all now hard-stop on any job/combo failure rather than pooling a partial result set. This is the checklist item the two Third-pass bugs would have been caught by immediately, had it existed on the second pass.

## Abstract's cluster-size statement: verdict

**Remove the mechanistic claim, retain a scaled descriptive statement -- confirmed on the corrected, complete data.** The properly-quantified, complete-grid opportunity control shows the size gradient's steepness is significantly *below* what genomic-coverage opportunity alone predicts (final numbers above). This evaluates Stage-2 reported regions, not the abstract's own Stage-1-unit estimand, so it cannot *directly* overturn that sentence -- but pending a Stage-1-scale version of the same check, removing the "selection-generated local LD" mechanistic language and keeping only the descriptive claim (regions/units with more markers have a higher absolute rate of true-positive linkage) is the safer course.

## Recommended figures/tables for the manuscript

**Main text:** none of these three analyses is proposed as a main-text figure.

**Supplementary Material:**
- `figure_simulation_floor_decomposition.pdf` and `figure_simulation_floor_decomposition_by_cell.pdf` -- together support the LD-aggregation finding *and* its cell-dependence; the by-cell figure is now important, not merely a robustness check, given the manuscript-hierarchy recommendation above.
- `figureS_simulation_stage2_size_truth.pdf` -- with the opportunity-null slope-difference result as a co-equal finding, not a caveat.
- `figure_simulation_null_truth_calibration.pdf` -- recommended with BOTH the within-cell (weak) and between-cell (strong, regime-level) results shown or cited, per Analysis 1's two-questions framing above; showing only the cell-adjusted numbers would discard the between-cell finding, which PK specifically flagged as the more useful one.
- `results/simulation_floor_choice_blind.tsv`, `results/simulation_null_truth_summary.tsv`, `results/simulation_stage2_opportunity_slope.tsv` as supplementary tables.

Not recommended: `results/simulation_stage2_region_details.rds` (too large/granular), the `_labelled` calibration figure (exploratory only).

## Follow-up: environment-genetic-structure alignment (new, 2026-09-18)

PK's hypothesis: false positives should become common when environmental variation follows the same spatial pattern as genetic relatedness -- **alignment**, since relatedness is a matrix, not a single variable (see `env-structure-alignment-hypothesis.md` memory for the full design rationale). Implemented in `R/15_env_structure_alignment.R`, no new permutations or association reruns (GRM/env already in each combo's bundle; joined against R/11's observed FP and a new per-env, unpooled export from R/13).

**Design:** per `(tag, cell, rep, env)` -- individual env continuations, not pooled across the 10 per map/burn-in (pooling would average away the variation of interest). Primary measure: population-level R² of env ~ top 5 eigenvectors of the population-level (block-mean) relationship matrix. Sensitivity measure: a Mantel-style correlation between the flattened relationship matrix and a Gaussian-kernel env-similarity matrix. Verified against synthetic aligned/unaligned data before running on the real grid (aligned case: R²=0.999, Mantel r=0.998; unaligned: R²=0.30, Mantel r=-0.03). All 700/700 combos succeeded (hard `{ok,data,error}`-gated, same convention as R/11/R/12/R/13).

**Result, revised after PK's second review of this figure/analysis: more nuanced than "alignment adds nothing beyond cell identity."** That was the first draft's headline; it is only correct for the two *discovery-conditional* outcomes below, and does not hold once the analysis is corrected to include zero-discovery combos.

**A. Discovery-conditional outcomes (FP proportion and null/obs ratio among reported regions) -- adjusted effect vanishes, matching the original headline:**

| model | outcome | r2_axes slope (log-odds / log2) | 95% cluster-bootstrap CI |
|---|---|---:|---:|
| unconditional | FP proportion | **+12.92** | [11.02, 14.63] -- excludes zero |
| adjusted for cell+tag+method | FP proportion | -1.37 | [-4.96, 0.71] -- includes zero |
| unconditional | null/obs ratio (log2) | **+19.76** | [14.53, 24.17] -- excludes zero |
| adjusted for cell+tag+method | null/obs ratio (log2) | -2.84 | [-5.82, 0.31] -- includes zero |

*(CIs above refreshed 2026-09-26 for the `rep`-alone clustering fix -- see item 4. Point estimates unchanged; both rows' zero-crossing conclusions are unaffected.)*

**B. Complementary full-grid outcomes (PK's point 3, added on second review) -- adjusted effect does NOT vanish:** `gp` above is conditioned on >=1 reported region existing (884/1,400 (combo,method) pairs; the other 516 had zero reported regions and are silently excluded, since FP proportion is undefined there). `results/simulation_env_structure_alignment_full_grid.tsv` keeps all 1,400 pairs (zero-filled, a real zero not a missing value) and models whether a false positive occurs at all, how many occur, and the expected spatial-null burden -- without conditioning on a discovery having happened:

| model | outcome | r2_axes slope | 95% cluster-bootstrap CI |
|---|---|---:|---:|
| unconditional | any FP occurred (logit) | +14.86 | [11.34, 18.53] -- excludes zero |
| **adjusted** for cell+tag+method | any FP occurred (logit) | **+9.35** | **[5.38, 13.60] -- excludes zero** |
| unconditional | # FP regions (log, Poisson) | +15.93 | [14.42, 17.47] -- excludes zero |
| **adjusted** for cell+tag+method | # FP regions (log, Poisson) | **+9.38** | **[7.66, 10.46] -- excludes zero** |
| unconditional | expected null region count (log2) | +17.52 | [15.44, 18.80] -- excludes zero |
| **adjusted** for cell+tag+method | expected null region count (log2) | **+0.99** | **[0.69, 1.31] -- excludes zero** |

*(CIs above refreshed 2026-09-26 for the `rep`-alone clustering fix -- see item 4. Point estimates unchanged; all six rows' zero-crossing conclusions are unaffected -- this section's "adjusted effect does NOT vanish" finding stands.)*

**Once zero-discovery combos are correctly retained, alignment DOES predict both whether a false positive occurs and how many occur, net of demographic cell.** This contradicts the "collinear with cell, no independent effect" reading of A alone.

**C. Sensitivity check (PK's point 4): the Mantel-style measure, run through the identical models, does not agree with r2_axes in the adjusted models -- several flip sign:**

| outcome | adjusted r2_axes slope | adjusted mantel_r slope |
|---|---:|---:|
| FP proportion | -1.37 [-4.96, 0.71] (incl. zero) | **-2.89 [-4.62, -1.35]** (excl. zero, negative) |
| any FP occurred | **+9.35 [5.38, 13.60]** (excl. zero, positive) | **-2.26 [-3.68, -0.52]** (excl. zero, **negative**) |
| # FP regions | **+9.38 [7.66, 10.46]** (excl. zero, positive) | **-1.35 [-2.30, -0.80]** (excl. zero, **negative**) |
| expected null count | **+0.99 [0.69, 1.31]** (excl. zero, positive) | **-0.68 [-0.96, -0.36]** (excl. zero, **negative**) |

*(CIs above refreshed 2026-09-26 for the `rep`-alone clustering fix -- see item 4. Point estimates unchanged; all four rows' zero-crossing conclusions, and the "measures disagree on direction" finding, are unaffected. Note: the unconditional -- not adjusted -- `fp_prop x mantel_r` row, not shown in this table, is the one row anywhere in this section whose zero-crossing status changed: [-2.87, 0.28] (incl. zero) -> [-2.20, -0.20] (excl. zero); it does not appear in any of this document's stated conclusions.)*

The two alignment measures **disagree on direction** for 3 of 4 adjusted comparisons where both are significant. This means "does alignment predict FP risk beyond cell" does not have one clean answer -- it depends on which of the two measures is used, and R² on 5 arbitrarily-chosen axes vs. a full-matrix Mantel correlation are not simply two views of the same thing; they can rank the same combos differently. **This divergence, not either measure's result taken alone, is the honest summary of the sensitivity check** -- reported per PK's request, not silently resolved in favour of one measure.

**D. Within-cell (r2_axes only, FP proportion and null-ratio only, per the original spec's scope): 12 of 14 CIs include zero, but two do not, both in the negative direction (opposite the hypothesis) -- corrected from the first draft, which wrongly stated all 14 include zero:**
- Spatial-null ratio, `V0.5_c2`: **-11.70 [-16.99, -7.38]**.
- FP proportion, `V2_c1.5`: **-15.00 [-51.76, -4.32]**.

No cell shows a positive within-cell relationship; the correct summary is "no consistent *positive* within-cell relationship," not "all intervals include zero." Both exceptions are single-cell, wide-CI, and in the direction that argues against the hypothesis, not for it -- not treated as a second finding requiring its own explanation, but not hidden either.

*(The `rep`-alone clustering fix -- see item 4 -- does not change any number in this section: restricting to one cell already makes `(cell,rep)` and `rep`-alone identical clustering, verified byte-identical against the pre-fix values.)*

**Overall:** the between-cell separation (Analysis 1's spatial-null regime-level result) remains the dominant, most robust pattern. Alignment adds real, adjusted-significant information about false-positive *occurrence and count* beyond cell identity (panel B), but not about false-positive *proportion among discoveries* or the *null/obs ratio* (panel A) -- and even where it is significant, the two ways of measuring alignment disagree on direction (panel C). This is a genuinely mixed result, not a clean confirmation or refutation of PK's hypothesis, and should be reported as such.

See `results/simulation_env_structure_alignment.tsv` (discovery-conditional, 884 rows), `results/simulation_env_structure_alignment_full_grid.tsv` (all 1,400 combo-method pairs, zero-filled), `results/simulation_env_structure_alignment_models.tsv` (all 34 slope models: 2 alignment measures x 5 outcomes x 2 variants, plus 14 within-cell), `results/simulation_null_truth_calibration_by_env.tsv` (the per-env, unpooled export from R/13 this analysis depends on), `figureS_simulation_env_structure_alignment.pdf` / `_labelled.pdf`.

**Bugs found and fixed during this analysis (both self-caught during PK's second review of it, confirmed against my own code, not assumed from PK's report alone):**
1. **Cluster-bootstrap index mismatch.** `cl_rows`/`n_cl` were built once from the full `gp` table and reused for the ratio model, which is fit on `ratio_d` (`gp` filtered to `ratio_null_obs>0`, a smaller table with its own 1..nrow(ratio_d) row numbering). A bootstrap index built from `gp`'s numbering routinely exceeded `ratio_d`'s row count; `data.table` silently returns an all-`NA` row for an out-of-range index rather than erroring, corrupting that replicate instead of failing loudly. Confirmed with a small reproduction (a 9-row/7-row table pair) before fixing. Fixed by making `.boot_slope()` always build its own cluster info fresh from whatever `data` it is actually given -- eliminates this whole class of mismatch and replaces the previous separate, duplicated `.boot_slope_c()` with one helper used everywhere. Corrected numbers above match PK's own independent rerun exactly (19.76 [14.45,24.04] and -2.84 [-6.26,1.05]).
2. **`intersect()` in an earlier draft of the same helper** silently deduplicated a cluster drawn more than once in the same bootstrap resample, breaking resampling-with-replacement. Fixed to direct list-indexing (`cl_rows[as.character(draw)]`), verified with a small synthetic index test.

**Figure fixes (PK's third review of the figure specifically):** confidence ribbons removed from both panels (they used ordinary model/OLS standard errors, which treat every environmental-continuation/BGS-pair row as independent -- directly contradicting the figure's own "only 5 independent histories per cell" caption; a correctly map-cluster-bootstrapped ribbon was judged not worth the added complexity for a supplementary figure, so the fitted line is now shown without one, per PK's own suggested resolution). Both panels' fitted lines now come from one model per cell with a shared alignment slope across methods (`r2_axes + method_label`), matching the saved within-cell models exactly -- the first draft fit fully independent slopes per method in panel A and left panel B unfixed even after panel A was corrected (caught and fixed together, not incrementally). "not simulated" labels added to the two (V,c) combinations this design excludes (`V1_c2`, `V2_c2`). Selection levels spelled out in full. The x-axis superscript (`R²`) rendered as "R…" in the PNG device on the mini; replaced with a plain-ASCII two-line label. The manuscript version (`figureS_simulation_env_structure_alignment.pdf`, no title) now also has no in-plot caption -- that text belongs in the manuscript's own LaTeX `\caption{}` and is written out verbatim to `figures/figureS_simulation_env_structure_alignment_caption.txt` by the same script; the `_labelled` version keeps title+caption for internal review. Caption also now states that 3 observations with `ratio_null_obs == 0` cannot appear on panel B's log axis and are omitted there only (not from panel A or the models).

## c=1 primary manuscript summary (2026-09-18, fourth round)

Per PK's request: the manuscript's primary method-performance analysis uses only the 600 high-gene-flow simulations at `c=1` (`V0.5_c1`, `V1_c1`, `V2_c1`); the remaining 800 (`c=1.5`/`c=2`) stay in the full-grid outputs as demographic stress tests. **Downstream summarisation only** -- no clustering, association testing, truth scoring, or null permutation was rerun; see the provenance section above for why none was needed. Implemented as new functions appended to `R/06_summarise.R` (`summarise_c1_primary()`, plus a small backward-compatible generalisation of `crossed_bootstrap_stratum()` to take a `methods` argument instead of a hardcoded global), reusing the existing `map_cluster_bootstrap_stratum()`/`crossed_bootstrap_stratum()` primitives unchanged -- fed pre-pooled-over-`V` data so a resampled rep pulls in all three `c=1` selection settings, both BGS treatments, and all ten envs together, exactly as specified, with no changes to the bootstrap machinery itself.

**Point estimates -- exact match to PK's own independently-derived checkpoint table (not hard-coded, computed fresh from pooled counts):**

| Method | Regions | TP | FP | Precision | Recall | Precision x Recall |
|---|---:|---:|---:|---:|---:|---:|
| Marker-wise EMMAX | 1,447 | 530 | 917 | 0.36628 | 0.47225 | 0.17297 |
| EMMAX Simes | 753 | 437 | 316 | 0.58035 | 0.41785 | 0.24249 |
| EMMAX consensus | 776 | 438 | 338 | 0.56443 | 0.41676 | 0.23523 |
| Marker-wise LFMM | 1,792 | 697 | 1,095 | 0.38895 | 0.56366 | 0.21923 |
| LFMM Simes | 1,032 | 563 | 469 | 0.54554 | 0.51251 | 0.27960 |

Reductions relative to each engine's own marker-wise baseline: tests -80.30% (both Simes/consensus, identical unit-set); regions Simes -47.96%, consensus -46.37%; FP regions Simes -65.54%, consensus -63.14% -- all matching PK's checkpoint figures (80.3%, 48.0%/46.4%, 65.5%/63.1%) to the reported precision.

**Required paired contrasts (2,000-replicate map-cluster bootstrap, percentile 95% CI, all four from the identical bootstrap draws -- verified via `stopifnot(identical(ba$b, bb$b))` in code, not assumed):**

| Contrast | diff precision | 95% CI | diff recall | 95% CI |
|---|---:|---:|---:|---:|
| EMMAX Simes - marker-wise EMMAX | +0.2141 | [0.1718, 0.2562] | -0.0544 | [-0.0714, -0.0397] |
| EMMAX consensus - marker-wise EMMAX | +0.1982 | [0.1724, 0.2218] | -0.0555 | [-0.0757, -0.0328] |
| **EMMAX Simes - EMMAX consensus** | **+0.0159** | **[-0.0050, 0.0369] -- includes zero** | **+0.0011** | **[-0.0170, 0.0194] -- includes zero** |
| LFMM Simes - marker-wise LFMM | +0.1566 | [0.1359, 0.1746] | -0.0511 | [-0.0682, -0.0362] |

**The direct Simes-vs-consensus contrast PK specifically flagged (point estimates 0.580 vs 0.564 precision, 0.418 vs 0.417 recall) is not statistically distinguishable**: both the precision and recall differences have bootstrap intervals spanning zero. Per PK's explicit instruction, neither representation is described as superior on the strength of this data.

**Sensitivity (labelled, not primary): crossed two-way bootstrap (reps and envs resampled independently), restricted to `c=1`.** Wider than the primary map-cluster intervals for every method (as expected -- it adds env-surface sampling uncertainty on top of map/burn-in uncertainty), but does not change which contrasts exclude/include zero. E.g. `emmax_snp_region` precision: primary [0.323, 0.418] vs crossed [0.297, 0.453].

**Mandatory checks, verified not assumed:**

1. Exactly 600 combos per `c=1` method -- **PASS**, hard-coded `stop()` gate in `summarise_c1_primary()`, printed and confirmed for all 5 methods including both LFMM methods (not just EMMAX).
2. All three intended cells present, no `c1.5`/`c2` leakage -- **PASS**, exact `cell %chin% c("V0.5_c1","V1_c1","V2_c1")` selection (never substring), `setequal()`-gated.
3. Point estimates equal direct sums from the current per-combination truth-score objects -- **PASS**, exact match to PK's independently-derived checkpoint table above.
4. Methods, BGS treatments and the three `V` settings stay paired in every bootstrap replicate -- **PASS by construction**: `rep_dt_c1` pools over `V`/tag/env *before* the bootstrap ever runs, so one resampled rep pulls all of them together; `map_cluster_bootstrap_stratum()` is called once for all 5 methods with one seed.
5. All ten environments remain together in the primary bootstrap -- **PASS**, same pre-pooling as check 4.
6. Simes and consensus compared using identical bootstrap draws -- **PASS**, `stopifnot(identical(ba$b, bb$b))` in the contrast loop.
7. Existing complete-grid outputs remain byte-identical -- **PASS**, verified directly: `results/simulation_performance.tsv`'s pooled totals unchanged (6,623/813/5,810 etc.) after adding the `c=1` code; the two RDS files in the "Output provenance investigation" section above were regenerated because they were demonstrably stale (last committed before the ordering-unification fix), not because of this `c=1` addition.
8. New results reproducible without rerunning association analyses -- **PASS**: `results/simulation_bootstrap_replicates_c1.rds` saves the actual bootstrap replicate draws (both primary and crossed), so the intervals can be recomputed from that file alone.
9. Log and this document report exact inputs/totals/validation -- this section.

**New outputs** (full-grid files listed elsewhere in this document are unchanged): `results/simulation_performance_c1.tsv`, `results/simulation_method_contrasts_c1.tsv`, `results/simulation_bootstrap_sensitivity_c1.tsv`, `results/simulation_bootstrap_replicates_c1.rds`. `R/06_summarise.R` now also writes a receipt (it previously had none) recording `n_bootstrap`, `seed_bootstrap`, and the primary-manuscript-scope note that `c=1` is primary while the complete 7-cell grid is retained.

**No manuscript, figure, or LaTeX value has been edited** -- per the standing instruction and this round's explicit request, this section is for review before any of that.

## Bootstrapped c=1 floor decomposition (2026-09-19, fifth round)

Per PK's request: quantify, within the 600 primary high-gene-flow simulations, how much of the precision/recall change across the floor grid comes from (1) filtering small Stage-1 units below the floor vs. (2) Simes/consensus aggregation. **Downstream summary only** -- implemented as a new standalone script, `R/16_floor_decomposition_c1.R`, that reads `out_final_v1/12_floor_decomposition_raw.rds` (already computed, unchanged) and never sources or re-executes `R/12_floor_decomposition.R` itself, which has no per-combo caching and would otherwise recompute the full 1,400-combo x 6-floor x 4-arm Stage-2 assembly grid from scratch. Chromosome pooling is satisfied by construction: each row of the raw table is already a per-combo total summed across both simulated chromosomes.

**Bootstrap design, exactly as specified:** 2,000 replicates, the ten map/burn-in identities as the sole resampling unit. All three `V` settings, both BGS treatments and all ten envs are pooled to `(method, arm, floor, rep)` *before* the bootstrap runs, so a resampled rep pulls all of them together automatically. Reuses `bootstrap_rep_matrix()` (already used throughout this pipeline) called with the **identical fixed seed for every one of the 42 (method, arm, floor) combinations** -- since that function reseeds internally before drawing, one shared seed guarantees byte-identical resampling multiplicities everywhere, so every contrast (including across floors, not just within one) is paired by construction.

### Two bugs found and fixed before trusting any output (both caught by the mandatory checks, not assumed correct from the script running without error)

1. **Wrong region-count column.** The raw table has two distinct counts: `n_significant` (significant Stage-1 units/markers *before* Stage-2 assembly) and `n_regions` (the actual assembled Stage-2 region count -- the one that belongs with TP/FP and matches `simulation_performance_c1.tsv`). The first draft summed `n_significant`, silently reporting a region count 7-15x too large (e.g. 10,861 instead of 1,447 at floor=2/arm A/EMMAX) while TP/FP were correct throughout. Caught immediately by mandatory check 5 (below), not forced past: fixed by summing `n_regions` instead.
2. **Missing count columns in the bootstrap helper.** `.derived()` (this script's per-(method,arm,floor) point-estimate/bootstrap-ratio helper) originally returned only precision/recall/fdp/frac_markers_retained/coverage, omitting the raw `n_tests`/`n_regions`/`FP` counts the count-based contrasts (`diff_n_tests`, `diff_n_regions`, `diff_n_fp`) need. This produced an all-`NA` (auto-typed `logical`) CI column rather than an error -- the first pass of mandatory check 8 (finiteness) did not catch it because it only inspected the per-arm bootstrap tables, which never had those columns to begin with, not the downstream contrasts table. Fixed by adding the three counts to `.derived()`'s output, and check 8 was strengthened with a second pass that explicitly inspects every `_ci_lo`/`_ci_hi` column of the final contrasts table for non-numeric/non-finite values -- this second pass is what would have caught the original bug, and is now permanent.

### Mandatory checks -- all PASS, verified not assumed

1. Exactly 600 combos per applicable (method, arm, floor) -- **PASS**, printed and gated for all 42 combinations.
2. Cells exactly `V0.5_c1`/`V1_c1`/`V2_c1` -- **PASS**.
3. No `c1.5`/`c2` leakage -- **PASS**.
4. Floor=1 arm B == arm A per combo, before pooling -- **PASS**, `TRUE` for both EMMAX (600/600) and LFMM (600/600), checked on TP/FP/`n_significant`/`n_regions`/`n_tests` together.
5. Floor=2 arms A/C/D reproduce `simulation_performance_c1.tsv` -- **PASS after the bug-1 fix above**: `n_regions`/TP/FP identical for EMMAX A/C/D and LFMM A/C (1,447/530/917; 753/437/316; 776/438/338; 1,792/697/1,095; 1,032/563/469) -- confirms both paths (the floor-sweep's own Stage-2 assembly and `05_score_truth.R`'s canonical scoring) are the same code on the same ordering-unified pipeline, not two implementations that happen to agree by luck.
6. Methods/BGS/`V`/envs/floors/arms paired within every bootstrap draw -- **PASS by construction**: pooling happens before the bootstrap; one seed for all 42 strata.
7. Simes/consensus-vs-marker contrasts use identical draws -- **PASS**, asserted via `stopifnot(nrow(ba)==nrow(bb))` plus the shared-seed construction.
8. All bootstrap replicates finite -- **PASS after the bug-2 fix above**, both the per-arm tables and (now explicitly checked) the final contrasts table.
9. Existing full-grid floor outputs byte-for-byte unchanged -- **PASS**, verified via a SHA-256 checksum of `results/simulation_floor_sweep.tsv` taken before and after this script ran, not merely "the script didn't write to that path."

### Point estimates at the canonical floor (2), EMMAX -- ratios of pooled counts, 95% percentile bootstrap CI

| Arm | Regions | TP | FP | Precision | Recall | FDP | Precision x Recall |
|---|---:|---:|---:|---:|---:|---:|---:|
| A: Unrestricted marker | 1,447 | 530 | 917 | 0.3663 [0.324, 0.419] | 0.4723 [0.435, 0.508] | 0.6337 [0.581, 0.676] | 0.1730 [0.150, 0.203] |
| B: Floor-filtered marker | 1,372 | 535 | 837 | 0.3899 [0.347, 0.445] | 0.4625 [0.425, 0.498] | 0.6101 [0.555, 0.653] | 0.1803 [0.157, 0.209] |
| C: Stage-1 Simes | 753 | 437 | 316 | 0.5803 [0.544, 0.628] | 0.4178 [0.378, 0.454] | 0.4197 [0.372, 0.456] | 0.2425 [0.216, 0.267] |
| D: Stage-1 consensus | 776 | 438 | 338 | 0.5644 [0.526, 0.611] | 0.4168 [0.376, 0.451] | 0.4356 [0.389, 0.474] | 0.2352 [0.214, 0.259] |

### Required paired contrasts at floor=2, EMMAX

| Contrast | diff precision | 95% CI | diff recall | 95% CI | diff FP regions | 95% CI |
|---|---:|---:|---:|---:|---:|---:|
| B - A (filtering alone) | +0.0237 | [0.0088, 0.0395] excl. 0 | -0.0098 | [-0.0250, 0.0022] incl. 0 | -80 | [-121, -45] excl. 0 |
| C - B (Simes after matching floor) | +0.1904 | [0.1394, 0.2377] excl. 0 | -0.0446 | [-0.0509, -0.0378] excl. 0 | -521 | [-696, -364] excl. 0 |
| C - A (filtering + Simes, total) | +0.2141 | [0.1712, 0.2583] excl. 0 | -0.0544 | [-0.0716, -0.0397] excl. 0 | -601 | [-764, -442] excl. 0 |
| D - B (consensus after matching floor) | +0.1745 | [0.1375, 0.2064] excl. 0 | -0.0457 | [-0.0636, -0.0270] excl. 0 | -499 | [-659, -358] excl. 0 |
| D - A (filtering + consensus, total) | +0.1982 | [0.1729, 0.2228] excl. 0 | -0.0555 | [-0.0754, -0.0334] excl. 0 | -579 | [-727, -442] excl. 0 |

### Does filtering alone explain a meaningful part of the improvement? -- Yes, within c=1, unlike the full 7-cell grid

At floor=2, B-A's precision CI **excludes zero** (+0.0237 [0.0088, 0.0395]) -- small but real. This is a genuine refinement of the full-grid finding (`ADDITIONAL_ANALYSES_AUDIT.md`'s Analysis 3 section, where the pooled-across-all-7-cells B-A contrast at floor=2 was reported as "close to zero, CI spans zero"): restricted to the high-gene-flow regime specifically, filtering alone is doing real, non-trivial, statistically supported work even at the canonical floor, not none. The effect **grows monotonically with floor** -- B-A precision diff is +0.0424 [0.0220, 0.0640] at floor=3, +0.0574 [0.0309, 0.0901] at floor=5, +0.1175 [0.0635, 0.1815] at floor=10, +0.1718 [0.0714, 0.3047] at floor=20 -- but so does the recall cost (diff recall -0.0098 at floor=2 down to -0.351 at floor=20) and the FP-region reduction (-80 at floor=2 down to -826 at floor=20, all CIs excluding zero from floor=2 onward). **None of this should be read as evidence that filtering explains most of the LD-aggregation-vs-filtering gap** -- C-B and D-B (the aggregation-after-matching-floor contrasts) remain roughly 4-8x larger than B-A in absolute precision terms at every floor -- but "filtering alone contributes nothing" is no longer an accurate summary for the c=1 regime specifically.

### Are Simes and consensus distinguishable after floor matching? -- Not at floor=2, but yes from floor=3 upward

A direct C-D (Simes minus consensus) paired contrast was added (not one of the 5 required contrasts, but needed to answer this question at every floor, not only floor=2 where the separate `R/06` c=1 primary summary already checked it once): reuses the same paired bootstrap draws as everything else here.

| Floor | diff precision (Simes - consensus) | 95% CI |
|---|---:|---:|
| 1 | +0.0109 | [-0.0090, 0.0342] includes 0 |
| 2 | +0.0159 | [-0.0047, 0.0379] includes 0 |
| 3 | +0.0420 | [0.0108, 0.0724] **excludes 0** |
| 5 | +0.0509 | [0.0342, 0.0690] **excludes 0** |
| 10 | +0.0405 | [0.0195, 0.0604] **excludes 0** |
| 20 | +0.0316 | [0.0034, 0.0676] **excludes 0** |

At floor=2 this exactly reproduces the earlier, independently-built `R/06` c=1 primary summary's result (+0.0159 [-0.0050, 0.0369] there vs. +0.0159 [-0.0047, 0.0379] here -- same point estimate, matching CIs within bootstrap noise, cross-validating both independently-written implementations). **At floor=2, the manuscript's canonical value, Simes and consensus remain statistically indistinguishable on precision, consistent with the earlier finding.** But at every higher floor tested (3, 5, 10, 20), Simes becomes **significantly more precise** than consensus, with recall differences staying small and mostly non-significant. This is a genuinely new finding from this analysis, not previously reported: the choice between Simes and consensus is closer to immaterial only near the floor actually used in the manuscript; at more aggressive floors it is not.

### Precision x Recall across the floor grid -- U-shaped for the aggregation arms, peaking around floor 2-3

Unrestricted marker's Precision x Recall is constant across floors by construction (0.1730). Floor-filtered marker peaks at floor=2-3 (~0.180) before declining to 0.065 at floor=20. Stage-1 Simes and consensus both peak at floor=2 (0.2425 and 0.2352 respectively) and decline similarly at high floors, with Simes consistently at or above consensus from floor=2 onward -- matching the C-D contrast above. Full table in `results/simulation_floor_sweep_c1.tsv`.

### Phenotype-blind floor-choice table, c=1 (no association testing, structural only)

Pooled as **ratios of pooled counts** (a deliberate methodological difference from the full-grid `results/simulation_floor_choice_blind.tsv`, which averages per-combo fractions -- noted here explicitly, not an oversight): `frac_markers_retained = sum(markers retained across 600 combos) / sum(total markers across 600 combos)`.

| Floor | Mean eligible units | Mean tests | Fraction markers retained | Mean chromosomes represented |
|---|---:|---:|---:|---:|
| 1 | 8,971.4 | 8,971.4 | 1.0000 | 2.000 |
| 2 | 2,994.9 | 2,994.9 | 0.6068 | 2.000 |
| 3 | 1,261.7 | 1,261.7 | 0.3788 | 2.000 |
| 5 | 357.4 | 357.4 | 0.1822 | 2.000 |
| 10 | 56.8 | 56.8 | 0.0626 | 2.000 |
| 20 | 11.5 | 11.5 | 0.0251 | 1.952 |

At floor=20, mean chromosome representation drops below 2 (1.952) -- a minority of c=1 combos lose an entire chromosome's worth of eligible units at the most aggressive floor tested, consistent with the severe coverage cost already flagged for floor=20 elsewhere in this document.

### Outputs, provenance, seeds

New files (existing full-grid files listed elsewhere in this document are unchanged, confirmed by checksum): `results/simulation_floor_sweep_c1.tsv`, `results/simulation_floor_contrasts_c1.tsv` (the 5 required contrasts plus the supplementary C-D), `results/simulation_floor_choice_blind_c1.tsv`, `results/simulation_floor_bootstrap_replicates_c1.rds` (the actual bootstrap draws, so intervals can be recomputed without rerunning anything upstream). `figures/figureS_simulation_floor_decomposition_c1.pdf`/`.png` (no title, no in-plot caption; the caption text is written verbatim to `figures/figureS_simulation_floor_decomposition_c1_caption.txt` for the manuscript's own LaTeX `\caption{}`) -- four panels (tests as % of unrestricted, precision, recall, precision x recall), arms A-D, EMMAX solid/LFMM dashed (no LFMM consensus arm), 95% bootstrap ribbons on the three bootstrapped panels (the tests panel is a deterministic percentage, not bootstrapped, and has no ribbon). `R/16_floor_decomposition_c1.R` writes its own receipt recording `c1_cells`, `floor_grid`, `n_bootstrap`, and the fixed `boot_seed` (`SEEDS[["bootstrap"]] + 30000`). The existing `figure_simulation_floor_decomposition_by_cell.pdf` (all 7 cells, showing where methods fail at `c=1.5`/`c=2`) is untouched and remains the complementary view, per PK's explicit instruction to retain it.

**No manuscript file has been touched.** This section is for review before any of it is incorporated.

## Status (2026-10-07 audit)

The earlier "Still outstanding" list is resolved: `values_simulation.tex` is
generated by `LDscnR_manuscript/generate_simulation_manuscript.R`, the abstract's
cluster-size sentence is gone, and the manuscript has since been edited against
these results. All generators and outputs the manuscript reads (R/16-R/18, the
floor-stability tables, `module_9sp/R/02b`) are now committed.

Caveats and open items:

- **Rep-only bootstraps have only 9 clusters.** The null-calibration and
  alignment inputs use 5 randomly chosen reps per cell (`R/13`), and rep 7 is
  never drawn, so the rep-clustered bootstraps in R/13-R/15 and R/20 resample 9
  units. Percentile intervals from 9 clusters are likely anti-conservative; 17
  of 20 alignment-model intervals narrowed after the fix. Treat these CIs as
  descriptive. The one row that newly excludes zero (fp_prop x mantel_r,
  unconditional, [-2.20, -0.20]) is not quoted in the manuscript.
- **Fst adds no out-of-sample prediction in c=1.** `R/20` leave-rep-out
  delta-R^2 (joint vs alignment-only) is 0.028 [-0.008, 0.063] in
  `c1_primary` (`results/exploratory_fst_alignment_joint_delta_r2_oos.tsv`); it
  is small but positive in the stress cells (0.0075 [0.0031, 0.0125]) and pooled
  (0.0026 [0.0013, 0.0037]). This qualifies commit 7409c69's reading that Fst
  carries independent information: in the primary regime it does not generalise
  across reps.
- **R/19 LD-decay provenance mismatch is unresolved.** R/19 still reads
  `module_sim/out/02_bundle` (pre-correction parser). It is exploratory and not
  cited, but `RESTRUCTURE.sh` would move `module_sim/` and break it. Recompute
  decay from the final_analysis bundles or retire R/19 before restructuring.
- **R/13-R/15 should be rerun end to end** once the bundles are mounted, so the
  13b/14b/15b summary-only refresh scripts (which duplicate their parents' code)
  can be retired.
