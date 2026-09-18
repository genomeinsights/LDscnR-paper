# Additional simulation analyses: audit report

Audit date: 2026-09-17. Revised 2026-09-18 after PK's first-pass review ("Second pass" section). **Revised again 2026-09-18 (later the same day)** after PK's second review, which found two silent-failure bugs in the second-pass code itself ("Third pass" section) -- **all numbers in the "Analysis 1" and "Analysis 2" RESULT sections below are from the fully corrected, complete final rerun and supersede every earlier number in this document**, including the "Second pass" section's own point-2 and point-3 tables (kept below for the historical record of what changed and why, but explicitly marked superseded where relevant).

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

**Fixed** in `R/14`'s `.fit_size_model()` and `R/13`'s own summary bootstrap, both now clustering by `(cell, rep)` only. This fix is correct and unaffected by the Third-pass bugs.

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

See `results/simulation_null_truth_calibration.tsv` (per-map/burn-in-level rows, 70 unique `(tag,cell,rep)` groups), `results/simulation_null_truth_summary.tsv` (raw and cell+tag-adjusted rho, `null_role`, pooled ratios with `(cell,rep)`-clustered bootstrap CIs), `figure_simulation_null_truth_calibration.pdf` / `_labelled.pdf`.

## Analysis 2: Stage-2 region size and false-positive status -- RESULT

**Descriptive.** FP proportion declines monotonically with both size measures, consistently across all 5 region methods and both BGS treatments: e.g. for `emmax_consensus_region`, 2-marker regions are 90.1% FP (bgs) / 86.5% FP (nobgs), falling to 64.8% FP (bgs) / 22.4% FP (nobgs) for 51+-marker regions -- the decline is real in both treatments, but BGS retains a substantially higher FP rate at large sizes than no-BGS.

**Model.** `FP ~ log2(size) + method + cell + tag`, map/burn-in-cluster bootstrap (B=2,000, clustered by `(cell,rep)`): `log2(n_markers)` coefficient **-0.628** (95% CI -0.706 to -0.566), `log2(n_units)` coefficient **-1.090** (95% CI -1.230 to -0.967) -- both stable after controlling for method, cell and BGS treatment. Manual VIF confirms span/density are severely collinear with marker count and correctly excluded from the primary model.

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
10. Bootstrap resampling respects shared map/burn-in structure -- **PASS, corrected on the second pass** (clustered by `(cell,rep)`, not `(tag,cell,rep)`).
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
| unconditional | FP proportion | **+12.92** | [8.33, 16.33] -- excludes zero |
| adjusted for cell+tag+method | FP proportion | -1.37 | [-5.52, 1.88] -- includes zero |
| unconditional | null/obs ratio (log2) | **+19.76** | [14.45, 24.04] -- excludes zero |
| adjusted for cell+tag+method | null/obs ratio (log2) | -2.84 | [-6.26, 1.05] -- includes zero |

**B. Complementary full-grid outcomes (PK's point 3, added on second review) -- adjusted effect does NOT vanish:** `gp` above is conditioned on >=1 reported region existing (884/1,400 (combo,method) pairs; the other 516 had zero reported regions and are silently excluded, since FP proportion is undefined there). `results/simulation_env_structure_alignment_full_grid.tsv` keeps all 1,400 pairs (zero-filled, a real zero not a missing value) and models whether a false positive occurs at all, how many occur, and the expected spatial-null burden -- without conditioning on a discovery having happened:

| model | outcome | r2_axes slope | 95% cluster-bootstrap CI |
|---|---|---:|---:|
| unconditional | any FP occurred (logit) | +14.86 | [11.66, 18.38] -- excludes zero |
| **adjusted** for cell+tag+method | any FP occurred (logit) | **+9.35** | **[5.75, 13.45] -- excludes zero** |
| unconditional | # FP regions (log, Poisson) | +15.93 | [12.34, 19.46] -- excludes zero |
| **adjusted** for cell+tag+method | # FP regions (log, Poisson) | **+9.38** | **[7.13, 11.60] -- excludes zero** |
| unconditional | expected null region count (log2) | +17.52 | [14.28, 20.25] -- excludes zero |
| **adjusted** for cell+tag+method | expected null region count (log2) | **+0.99** | **[0.30, 1.75] -- excludes zero** |

**Once zero-discovery combos are correctly retained, alignment DOES predict both whether a false positive occurs and how many occur, net of demographic cell.** This contradicts the "collinear with cell, no independent effect" reading of A alone.

**C. Sensitivity check (PK's point 4): the Mantel-style measure, run through the identical models, does not agree with r2_axes in the adjusted models -- several flip sign:**

| outcome | adjusted r2_axes slope | adjusted mantel_r slope |
|---|---:|---:|
| FP proportion | -1.37 [-5.52, 1.88] (incl. zero) | **-2.89 [-4.85, -1.28]** (excl. zero, negative) |
| any FP occurred | **+9.35 [5.75, 13.45]** (excl. zero, positive) | **-2.26 [-4.08, -0.41]** (excl. zero, **negative**) |
| # FP regions | **+9.38 [7.13, 11.60]** (excl. zero, positive) | **-1.35 [-2.63, -0.33]** (excl. zero, **negative**) |
| expected null count | **+0.99 [0.30, 1.75]** (excl. zero, positive) | **-0.68 [-1.03, -0.29]** (excl. zero, **negative**) |

The two alignment measures **disagree on direction** for 3 of 4 adjusted comparisons where both are significant. This means "does alignment predict FP risk beyond cell" does not have one clean answer -- it depends on which of the two measures is used, and R² on 5 arbitrarily-chosen axes vs. a full-matrix Mantel correlation are not simply two views of the same thing; they can rank the same combos differently. **This divergence, not either measure's result taken alone, is the honest summary of the sensitivity check** -- reported per PK's request, not silently resolved in favour of one measure.

**D. Within-cell (r2_axes only, FP proportion and null-ratio only, per the original spec's scope): 12 of 14 CIs include zero, but two do not, both in the negative direction (opposite the hypothesis) -- corrected from the first draft, which wrongly stated all 14 include zero:**
- Spatial-null ratio, `V0.5_c2`: **-11.70 [-16.99, -7.38]**.
- FP proportion, `V2_c1.5`: **-15.00 [-51.76, -4.32]**.

No cell shows a positive within-cell relationship; the correct summary is "no consistent *positive* within-cell relationship," not "all intervals include zero." Both exceptions are single-cell, wide-CI, and in the direction that argues against the hypothesis, not for it -- not treated as a second finding requiring its own explanation, but not hidden either.

**Overall:** the between-cell separation (Analysis 1's spatial-null regime-level result) remains the dominant, most robust pattern. Alignment adds real, adjusted-significant information about false-positive *occurrence and count* beyond cell identity (panel B), but not about false-positive *proportion among discoveries* or the *null/obs ratio* (panel A) -- and even where it is significant, the two ways of measuring alignment disagree on direction (panel C). This is a genuinely mixed result, not a clean confirmation or refutation of PK's hypothesis, and should be reported as such.

See `results/simulation_env_structure_alignment.tsv` (discovery-conditional, 884 rows), `results/simulation_env_structure_alignment_full_grid.tsv` (all 1,400 combo-method pairs, zero-filled), `results/simulation_env_structure_alignment_models.tsv` (all 34 slope models: 2 alignment measures x 5 outcomes x 2 variants, plus 14 within-cell), `results/simulation_null_truth_calibration_by_env.tsv` (the per-env, unpooled export from R/13 this analysis depends on), `figureS_simulation_env_structure_alignment.pdf` / `_labelled.pdf`.

**Bugs found and fixed during this analysis (both self-caught during PK's second review of it, confirmed against my own code, not assumed from PK's report alone):**
1. **Cluster-bootstrap index mismatch.** `cl_rows`/`n_cl` were built once from the full `gp` table and reused for the ratio model, which is fit on `ratio_d` (`gp` filtered to `ratio_null_obs>0`, a smaller table with its own 1..nrow(ratio_d) row numbering). A bootstrap index built from `gp`'s numbering routinely exceeded `ratio_d`'s row count; `data.table` silently returns an all-`NA` row for an out-of-range index rather than erroring, corrupting that replicate instead of failing loudly. Confirmed with a small reproduction (a 9-row/7-row table pair) before fixing. Fixed by making `.boot_slope()` always build its own cluster info fresh from whatever `data` it is actually given -- eliminates this whole class of mismatch and replaces the previous separate, duplicated `.boot_slope_c()` with one helper used everywhere. Corrected numbers above match PK's own independent rerun exactly (19.76 [14.45,24.04] and -2.84 [-6.26,1.05]).
2. **`intersect()` in an earlier draft of the same helper** silently deduplicated a cluster drawn more than once in the same bootstrap resample, breaking resampling-with-replacement. Fixed to direct list-indexing (`cl_rows[as.character(draw)]`), verified with a small synthetic index test.

**Figure fixes (PK's third review of the figure specifically):** confidence ribbons removed from both panels (they used ordinary model/OLS standard errors, which treat every environmental-continuation/BGS-pair row as independent -- directly contradicting the figure's own "only 5 independent histories per cell" caption; a correctly map-cluster-bootstrapped ribbon was judged not worth the added complexity for a supplementary figure, so the fitted line is now shown without one, per PK's own suggested resolution). Both panels' fitted lines now come from one model per cell with a shared alignment slope across methods (`r2_axes + method_label`), matching the saved within-cell models exactly -- the first draft fit fully independent slopes per method in panel A and left panel B unfixed even after panel A was corrected (caught and fixed together, not incrementally). "not simulated" labels added to the two (V,c) combinations this design excludes (`V1_c2`, `V2_c2`). Selection levels spelled out in full. The x-axis superscript (`R²`) rendered as "R…" in the PNG device on the mini; replaced with a plain-ASCII two-line label. The manuscript version (`figureS_simulation_env_structure_alignment.pdf`, no title) now also has no in-plot caption -- that text belongs in the manuscript's own LaTeX `\caption{}` and is written out verbatim to `figures/figureS_simulation_env_structure_alignment_caption.txt` by the same script; the `_labelled` version keeps title+caption for internal review. Caption also now states that 3 observations with `ratio_null_obs == 0` cannot appear on panel B's log axis and are omitted there only (not from panel A or the models).

## Still outstanding

- `values_simulation.tex`'s macros for `emmax_snp_region`/`lfmm_snp_region` have **not** been updated.
- The abstract's cluster-size sentence has **not** been edited.
- The manuscript's provisional conclusion paragraph has **not** been edited or finalised.
- No manuscript file has been touched at any point in this work.
