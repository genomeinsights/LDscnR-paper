## =============================================================================
## final_analysis/R/16_floor_decomposition_c1.R
##
## Bootstrapped, c=1-pooled version of the floor-decomposition analysis
## (PK, 2026-09-19): how much of the precision/recall change across the
## floor grid, WITHIN the primary high-gene-flow simulations, comes from
## (1) filtering small Stage-1 units below the floor vs. (2) aggregating
## the remaining LD-correlated markers (Simes/consensus)?
##
## DOWNSTREAM SUMMARY ONLY -- reads out_final_v1/12_floor_decomposition_raw.rds
## (already computed, unchanged by this script) and pools/bootstraps from it.
## No clustering, association testing, truth scoring or null permutation is
## rerun; R/12_floor_decomposition.R itself is never sourced or re-executed
## (it has no per-combo stage_stale() caching -- running it again would
## recompute the full 1,400-combo x 6-floor x 4-arm Stage-2 assembly grid
## from scratch, which is exactly what this script must NOT do).
##
## Scope: the 600 primary c=1 simulations (V0.5_c1, V1_c1, V2_c1) only,
## exact cell selection (never a substring match). Pooled over all three V
## settings, both BGS treatments, all ten reps, all ten envs. Chromosome
## pooling is already satisfied by construction -- each row of the raw
## table is a (combo, arm, floor, method) total already summed across both
## simulated chromosomes within that combo; there is no finer-grained
## per-chromosome table to additionally pool.
##
## Bootstrap: 2,000 map-cluster replicates, REP ALONE as the resampling
## unit (matching R/06's summarise_c1_primary() exactly, for the same
## reason: V/tag/env are pooled here, not held fixed per stratum, so a
## resampled rep identity must pull in that rep's rows from all three V
## settings, both tags, and all ten envs together). Uses the SAME
## bootstrap_rep_matrix() helper already used throughout this pipeline,
## called with the IDENTICAL seed for every (method, arm, floor)
## combination -- since bootstrap_rep_matrix() calls set.seed(seed) itself
## before drawing, reusing one seed everywhere guarantees byte-identical
## resampling multiplicities across methods, arms AND floors, so every
## contrast (including across floors) is paired by construction, not by
## convention.
## =============================================================================
suppressMessages({library(data.table)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "R", "00_config.R"))
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "R", "helpers_stage2_truth.R"))
say("=== 16_floor_decomposition_c1 ===\n\n")

C1_CELLS <- c("V0.5_c1", "V1_c1", "V2_c1")
FLOOR_GRID_C1 <- c(1, 2, 3, 5, 10, 20)
ARM_LABELS_C1 <- c(A_unrestricted_marker = "Unrestricted marker", B_floor_filtered_marker = "Floor-filtered marker",
                   C_stage1_simes = "Stage-1 Simes", D_stage1_consensus = "Stage-1 consensus")
C1_BOOT_SEED <- SEEDS[["bootstrap"]] + 30000   ## one fixed seed for every (method,arm,floor) bootstrap_rep_matrix() call

RAW_FILE <- "out_final_v1/12_floor_decomposition_raw.rds"
if (!file.exists(RAW_FILE)) stop("R/12_floor_decomposition.R has not produced ", RAW_FILE, " -- run it first (this script never reruns it).")
raw <- readRDS(RAW_FILE)
raw_before_checksum <- digest::digest("results/simulation_floor_sweep.tsv", algo = "sha256", file = TRUE)   ## mandatory check 9 (paired with the check after writing outputs, below)

d1 <- raw[cell %chin% C1_CELLS]

## ---- mandatory checks 2-3: exact cell scope, no c1.5/c2 leakage ----
if (!setequal(unique(d1$cell), C1_CELLS))
  stop("floor c=1 summary: unexpected cell set -- got ", paste(sort(unique(d1$cell)), collapse = ", "),
      "; expected exactly ", paste(C1_CELLS, collapse = ", "))
if (any(grepl("c1\\.5|c2$", unique(d1$cell))))
  stop("floor c=1 summary: a c1.5/c2 cell leaked into the c=1 subset -- ", paste(unique(d1$cell), collapse = ", "))

## ---- mandatory check 1: exactly 600 combos per applicable (method, arm, floor) ----
combo_check <- d1[, .(n_combo = uniqueN(paste(tag, cell, rep, env))), by = .(method, arm, floor)]
say("[1] combo count per (method, arm, floor) -- expect 600 everywhere\n")
print(combo_check[order(method, arm, floor)])
if (any(combo_check$n_combo != 600))
  stop("floor c=1 summary: expected exactly 600 combos for every (method,arm,floor), got: ",
      paste(sprintf("%s/%s/floor%d=%d", combo_check$method, combo_check$arm, combo_check$floor, combo_check$n_combo)[combo_check$n_combo != 600], collapse = "; "))

## ---- mandatory check 4: floor=1 arm B == arm A per combo, BEFORE pooling ----
say("\n[2] validation: floor=1 arm B == arm A, per combo, before any pooling\n")
for (m in unique(d1$method)) {
  a1 <- d1[method == m & arm == "A_unrestricted_marker" & floor == 1]
  b1 <- d1[method == m & arm == "B_floor_filtered_marker" & floor == 1]
  key <- c("tag", "cell", "rep", "env")
  chk <- merge(a1[, .(tag, cell, rep, env, TP, FP, n_significant, n_regions, n_tests)],
              b1[, .(tag, cell, rep, env, TP, FP, n_significant, n_regions, n_tests)], by = key, suffixes = c("_A", "_B"))
  ok <- nrow(chk) == 600 && all(chk$TP_A == chk$TP_B & chk$FP_A == chk$FP_B & chk$n_significant_A == chk$n_significant_B &
                                chk$n_regions_A == chk$n_regions_B & chk$n_tests_A == chk$n_tests_B)
  say("    %s: %s (%d combos compared)\n", m, ok, nrow(chk))
  if (!ok) stop("floor c=1 summary: floor=1 arm A != arm B for method '", m, "' -- per-combo mismatch found, not a pooling artefact")
}

## ---- pooled-ratio denominator for "fraction of markers retained": ----
## pool the NUMERATOR (n_eligible_markers at this floor) and the
## DENOMINATOR (n_eligible_markers at floor=1, same arm+method+combo)
## separately across combos, then divide -- a ratio of pooled counts, never
## a mean of per-combo ratios (frac_markers_retained itself is not used
## directly for exactly this reason).
floor1_markers <- d1[floor == 1, .(tag, cell, rep, env, arm, method, n_eligible_markers_floor1 = n_eligible_markers)]
d1 <- merge(d1, floor1_markers, by = c("tag", "cell", "rep", "env", "arm", "method"))

## ---- pool to (method, arm, floor, rep) -- the unit the bootstrap resamples ----
## [!] BUG FIXED before first real run (PK's mandatory check 5 caught this
## immediately): the raw table has TWO distinct counts -- n_significant
## (significant Stage-1 units/markers BEFORE Stage-2 assembly) and
## n_regions (the actual assembled Stage-2 region count, the one that
## belongs with TP/FP and matches simulation_performance_c1.tsv). An
## earlier draft of this line summed n_significant here, silently reporting
## a ~7-15x-too-large "region count" that still happened to leave TP/FP
## correct (they were never wrong) -- confirmed via direct comparison:
## sum(n_regions) at floor=2/arm=A/emmax = 1447 = TP+FP = canonical exactly;
## sum(n_significant) for the same slice = 10861, a different quantity
## entirely. Never silently forced past; see check 5 below.
rep_dt <- d1[, .(n_tests = sum(n_tests), n_regions = sum(n_regions), TP = sum(TP), FP = sum(FP),
                 n_detectable_qtn = sum(n_detectable_qtn), n_recovered = sum(n_recovered),
                 n_eligible_markers = sum(n_eligible_markers), n_eligible_markers_floor1 = sum(n_eligible_markers_floor1),
                 n_qtn_covered_by_eligible = sum(n_qtn_covered_by_eligible)),
             by = .(method, arm, floor, rep)]

.derived <- function(mat) {
  ## mat: B x k (or 1 x k for the point estimate) matrix of pooled counts.
  ## Same never-divide-by-zero convention as pooled_point()/06_summarise.R
  ## throughout this pipeline: max(denominator, 1), never a bare 0.
  ## [!] BUG FIXED before this was caught downstream: the raw counts
  ## (n_tests, n_regions, FP) were not passed through here, so the count-
  ## based paired contrasts (diff_n_tests, diff_n_regions, diff_n_fp) tried
  ## to read a column that did not exist, silently producing an all-NA
  ## (auto-typed logical) CI column rather than erroring -- caught by
  ## inspecting the actual output table, not assumed correct from the
  ## script running without error. Now included explicitly so every
  ## required contrast has a real bootstrap CI.
  data.table(n_tests = mat[, "n_tests"], n_regions = mat[, "n_regions"], FP = mat[, "FP"],
            precision = mat[, "TP"] / pmax(mat[, "TP"] + mat[, "FP"], 1),
            recall = mat[, "n_recovered"] / pmax(mat[, "n_detectable_qtn"], 1),
            fdp = mat[, "FP"] / pmax(mat[, "TP"] + mat[, "FP"], 1),
            frac_markers_retained = mat[, "n_eligible_markers"] / pmax(mat[, "n_eligible_markers_floor1"], 1),
            coverage = mat[, "n_qtn_covered_by_eligible"] / pmax(mat[, "n_detectable_qtn"], 1))
}

say("\n[3] c=1 map-cluster bootstrap (B=%d, one fixed seed for every method x arm x floor)\n", N_BOOTSTRAP)
strata <- unique(rep_dt[, .(method, arm, floor)])
point_rows <- list(); boot_store <- list()
for (i in seq_len(nrow(strata))) {
  m <- strata$method[i]; a <- strata$arm[i]; fl <- strata$floor[i]
  rd <- rep_dt[method == m & arm == a & floor == fl][order(rep)]
  stopifnot("rep_dt must have exactly REPS_ALL rows per (method,arm,floor)" = nrow(rd) == length(REPS_ALL) && identical(rd$rep, REPS_ALL))
  cnt_cols <- c("n_tests", "n_regions", "TP", "FP", "n_detectable_qtn", "n_recovered", "n_eligible_markers", "n_eligible_markers_floor1", "n_qtn_covered_by_eligible")
  point_mat <- matrix(colSums(rd[, ..cnt_cols]), nrow = 1, dimnames = list(NULL, cnt_cols))
  boot_mat <- bootstrap_rep_matrix(rd[, ..cnt_cols], N_BOOTSTRAP, seed = C1_BOOT_SEED)

  pt_d <- .derived(point_mat); bt_d <- .derived(boot_mat)
  pt_d[, prec_x_rec := precision * recall]; bt_d[, prec_x_rec := precision * recall]

  key <- sprintf("%s__%s__%d", m, a, fl)
  boot_store[[key]] <- bt_d   ## kept for the paired contrasts below

  point_rows[[length(point_rows) + 1]] <- data.table(
    method = m, arm = a, arm_label = ARM_LABELS_C1[[a]], floor = fl,
    n_tests = point_mat[, "n_tests"], n_regions = point_mat[, "n_regions"], TP = point_mat[, "TP"], FP = point_mat[, "FP"],
    n_detectable_qtn = point_mat[, "n_detectable_qtn"], n_recovered = point_mat[, "n_recovered"],
    precision = pt_d$precision, precision_ci_lo = ci_quantile(bt_d$precision)[1], precision_ci_hi = ci_quantile(bt_d$precision)[2],
    recall = pt_d$recall, recall_ci_lo = ci_quantile(bt_d$recall)[1], recall_ci_hi = ci_quantile(bt_d$recall)[2],
    fdp = pt_d$fdp, fdp_ci_lo = ci_quantile(bt_d$fdp)[1], fdp_ci_hi = ci_quantile(bt_d$fdp)[2],
    prec_x_rec = pt_d$prec_x_rec, prec_x_rec_ci_lo = ci_quantile(bt_d$prec_x_rec)[1], prec_x_rec_ci_hi = ci_quantile(bt_d$prec_x_rec)[2],
    frac_markers_retained = pt_d$frac_markers_retained, coverage = pt_d$coverage)
}
performance_c1 <- rbindlist(point_rows)

## ---- mandatory check 8: every bootstrap replicate produced finite values ----
all_finite <- all(vapply(boot_store, function(x) all(vapply(x, function(col) all(is.finite(col)), logical(1))), logical(1)))
say("\n[4] all bootstrap replicates finite: %s\n", all_finite)
if (!all_finite) stop("floor c=1 summary: some bootstrap replicate produced a non-finite value despite the max(.,1) safeguard -- investigate before trusting any interval")

## ---- required paired contrasts, per floor, from the SAME boot_store draws ----
say("\n[5] required paired contrasts per floor\n")
CONTRAST_SPECS <- list(
  emmax = list(list(a = "B_floor_filtered_marker", b = "A_unrestricted_marker", label = "B - A (filtering alone)"),
              list(a = "C_stage1_simes", b = "B_floor_filtered_marker", label = "C - B (Simes aggregation after matching floor)"),
              list(a = "C_stage1_simes", b = "A_unrestricted_marker", label = "C - A (filtering + Simes, total)"),
              list(a = "D_stage1_consensus", b = "B_floor_filtered_marker", label = "D - B (consensus aggregation after matching floor)"),
              list(a = "D_stage1_consensus", b = "A_unrestricted_marker", label = "D - A (filtering + consensus, total)"),
              ## Not one of the 5 REQUIRED contrasts, but needed to directly
              ## answer the explicit reporting question "are Simes and
              ## consensus distinguishable after floor matching" at every
              ## floor, not just floor=2 (already answered once, at floor=2
              ## only, by the separate c=1 primary summary in R/06). Uses the
              ## same paired bootstrap draws as everything else here.
              list(a = "C_stage1_simes", b = "D_stage1_consensus", label = "C - D (Simes vs consensus, SUPPLEMENTARY, not one of the 5 required contrasts)")),
  lfmm = list(list(a = "B_floor_filtered_marker", b = "A_unrestricted_marker", label = "B - A (filtering alone)"),
             list(a = "C_stage1_simes", b = "B_floor_filtered_marker", label = "C - B (Simes aggregation after matching floor)"),
             list(a = "C_stage1_simes", b = "A_unrestricted_marker", label = "C - A (filtering + Simes, total)"))
)
contrast_rows <- list()
for (m in names(CONTRAST_SPECS)) for (fl in FLOOR_GRID_C1) {
  for (spec in CONTRAST_SPECS[[m]]) {
    key_a <- sprintf("%s__%s__%d", m, spec$a, fl); key_b <- sprintf("%s__%s__%d", m, spec$b, fl)
    if (is.null(boot_store[[key_a]]) || is.null(boot_store[[key_b]])) next   ## e.g. lfmm has no D arm
    ba <- boot_store[[key_a]]; bb <- boot_store[[key_b]]
    pa <- performance_c1[method == m & arm == spec$a & floor == fl]
    pb <- performance_c1[method == m & arm == spec$b & floor == fl]
    ## mandatory check 7: identical bootstrap draws for Simes/consensus-vs-marker
    ## contrasts -- guaranteed by construction (one fixed seed everywhere), but
    ## asserted explicitly rather than only assumed: same number of rows, and
    ## both derived from bootstrap_rep_matrix() calls seeded identically.
    stopifnot(nrow(ba) == nrow(bb))
    contrast_rows[[length(contrast_rows) + 1]] <- data.table(
      method = m, floor = fl, arm_a = spec$a, arm_b = spec$b, label = spec$label,
      diff_precision = pa$precision - pb$precision, diff_precision_ci_lo = ci_quantile(ba$precision - bb$precision)[1], diff_precision_ci_hi = ci_quantile(ba$precision - bb$precision)[2],
      diff_recall = pa$recall - pb$recall, diff_recall_ci_lo = ci_quantile(ba$recall - bb$recall)[1], diff_recall_ci_hi = ci_quantile(ba$recall - bb$recall)[2],
      diff_fdp = pa$fdp - pb$fdp, diff_fdp_ci_lo = ci_quantile(ba$fdp - bb$fdp)[1], diff_fdp_ci_hi = ci_quantile(ba$fdp - bb$fdp)[2],
      diff_prec_x_rec = pa$prec_x_rec - pb$prec_x_rec, diff_prec_x_rec_ci_lo = ci_quantile(ba$prec_x_rec - bb$prec_x_rec)[1], diff_prec_x_rec_ci_hi = ci_quantile(ba$prec_x_rec - bb$prec_x_rec)[2],
      diff_n_tests = pa$n_tests - pb$n_tests, diff_n_tests_ci_lo = ci_quantile(ba$n_tests - bb$n_tests)[1], diff_n_tests_ci_hi = ci_quantile(ba$n_tests - bb$n_tests)[2],
      diff_n_regions = pa$n_regions - pb$n_regions, diff_n_regions_ci_lo = ci_quantile(ba$n_regions - bb$n_regions)[1], diff_n_regions_ci_hi = ci_quantile(ba$n_regions - bb$n_regions)[2],
      diff_n_fp = pa$FP - pb$FP, diff_n_fp_ci_lo = ci_quantile(ba$FP - bb$FP)[1], diff_n_fp_ci_hi = ci_quantile(ba$FP - bb$FP)[2])
  }
}
contrasts_c1 <- rbindlist(contrast_rows)

## ---- mandatory check 8, part 2: the CONTRASTS table itself, not just the ----
## underlying per-arm boot tables (the n_tests/n_regions/FP omission from
## .derived() was only caught here, downstream of the first check-8 pass --
## an all-NA column silently types as logical rather than erroring, so this
## check must inspect the actual numeric CI columns explicitly, not just
## run without error).
ci_cols <- grep("_ci_lo$|_ci_hi$", names(contrasts_c1), value = TRUE)
non_finite <- vapply(ci_cols, function(cn) !is.numeric(contrasts_c1[[cn]]) || anyNA(contrasts_c1[[cn]]) || any(!is.finite(contrasts_c1[[cn]])), logical(1))
say("    contrasts table CI columns all finite: %s\n", !any(non_finite))
if (any(non_finite)) stop("floor c=1 summary: non-finite/non-numeric contrast CI column(s): ", paste(ci_cols[non_finite], collapse = ", "))

## ---- mandatory check 5: floor=2 arms A/C/D reproduce simulation_performance_c1.tsv ----
say("\n[6] validation: floor=2 arms A/C/D vs. results/simulation_performance_c1.tsv\n")
perf_c1_file <- "results/simulation_performance_c1.tsv"
if (!file.exists(perf_c1_file)) stop("results/simulation_performance_c1.tsv not found -- run R/06_summarise.R's c=1 addition first.")
canon <- fread(perf_c1_file)
check5_map <- list(list(arm = "A_unrestricted_marker", method = "emmax", canon_method = "emmax_snp_region"),
                   list(arm = "C_stage1_simes", method = "emmax", canon_method = "emmax_simes_region"),
                   list(arm = "D_stage1_consensus", method = "emmax", canon_method = "emmax_consensus_region"),
                   list(arm = "A_unrestricted_marker", method = "lfmm", canon_method = "lfmm_snp_region"),
                   list(arm = "C_stage1_simes", method = "lfmm", canon_method = "lfmm_simes_region"))
check5_ok <- TRUE
for (spec in check5_map) {
  fl2 <- performance_c1[method == spec$method & arm == spec$arm & floor == 2]
  cn <- canon[method == spec$canon_method]
  match <- fl2$n_regions == cn$n_significant && fl2$TP == cn$TP && fl2$FP == cn$FP
  say("    %s/%s @ floor=2 vs %s: n_regions %d vs %d, TP %d vs %d, FP %d vs %d -- %s\n",
      spec$method, spec$arm, spec$canon_method, fl2$n_regions, cn$n_significant, fl2$TP, cn$TP, fl2$FP, cn$FP, match)
  if (!match) check5_ok <- FALSE
}
if (!check5_ok)
  stop("floor c=1 summary: floor=2 arms A/C/D do not reproduce results/simulation_performance_c1.tsv exactly -- ",
      "this is NOT expected to differ (both paths call the same assemble_stage2()/score_stage2_regions() on the ",
      "same ordering-unified pipeline), so this is treated as a real discrepancy to resolve, not forced past.")

## =============================================================================
## ---- phenotype-blind floor-choice table, c=1 only (no association testing) ----
## Mirrors R/12's own blind table exactly in spirit (structural only: eligible
## units/tests/marker retention from the ALREADY-COMPUTED Stage-1 clustering
## in each combo's cached ld_units.rds bundle, re-filtered by floor via
## LDscnR:::.ld_outlier_units() -- a cheap deterministic size filter over an
## existing clustering object, not a re-clustering). Pooled as ratios of
## pooled counts (not per-combo means, unlike R/12's own mean()-based
## version) for consistency with the rest of this c=1 analysis's discipline;
## noted explicitly in the audit as a deliberate methodological difference
## from the full-grid table, not an oversight.
## =============================================================================
say("\n[7] phenotype-blind floor-choice table, c=1 only (structural, no test statistics)\n")
combos_c1 <- CJ(tag = TAGS_ALL, cell = C1_CELLS, rep = REPS_ALL, env = ENVS_ALL)
stopifnot(nrow(combos_c1) == 600)
blind_rows <- list()
for (i in seq_len(nrow(combos_c1))) {
  combo_id <- sprintf("%s_%s_rep%d_env%d", combos_c1$tag[i], combos_c1$cell[i], combos_c1$rep[i], combos_c1$env[i])
  bf <- file.path(stage_dir("02_build_ld_units", combo_id), "ld_units.rds")
  if (!file.exists(bf)) stop("floor c=1 blind table: missing bundle for ", combo_id, " -- expected all 600 c=1 bundles to exist")
  bb <- readRDS(bf)
  u1 <- LDscnR:::.ld_outlier_units(bb$stage1, bb$map, 1L)
  n_markers_total <- nrow(bb$map)
  for (fl in FLOOR_GRID_C1) {
    el <- u1[n_markers >= fl]
    mk <- unique(unlist(el$members, use.names = FALSE))
    cov_by_chr <- el[, .(n_units = .N, n_markers = sum(n_markers)), by = Chr]
    blind_rows[[length(blind_rows) + 1]] <- data.table(
      tag = combos_c1$tag[i], cell = combos_c1$cell[i], rep = combos_c1$rep[i], env = combos_c1$env[i], floor = fl,
      n_eligible_units = nrow(el), n_tests = nrow(el), n_markers_retained = length(mk), n_markers_total = n_markers_total,
      n_chr_represented = nrow(cov_by_chr))
  }
}
blind_c1 <- rbindlist(blind_rows)
blind_pooled_c1 <- blind_c1[, .(n_combo = .N, sum_n_eligible_units = sum(n_eligible_units), sum_n_tests = sum(n_tests),
                                sum_n_markers_retained = sum(n_markers_retained), sum_n_markers_total = sum(n_markers_total),
                                mean_n_chr_represented = mean(n_chr_represented)),
                            by = floor]
blind_pooled_c1[, `:=`(mean_n_eligible_units = sum_n_eligible_units / n_combo, mean_n_tests = sum_n_tests / n_combo,
                       frac_markers_retained = sum_n_markers_retained / sum_n_markers_total)]
setorder(blind_pooled_c1, floor)
say("    n_combo per floor (expect 600 everywhere): %s\n", paste(unique(blind_pooled_c1$n_combo), collapse = ", "))
if (!all(blind_pooled_c1$n_combo == 600)) stop("floor c=1 blind table: expected 600 combos at every floor")
print(blind_pooled_c1)

## ---- write outputs (NEW, separately-named files only) ----
fwrite(performance_c1, "results/simulation_floor_sweep_c1.tsv", sep = "\t")
fwrite(contrasts_c1, "results/simulation_floor_contrasts_c1.tsv", sep = "\t")
fwrite(blind_pooled_c1, "results/simulation_floor_choice_blind_c1.tsv", sep = "\t")
saveRDS(list(boot = boot_store, rep_dt = rep_dt, performance = performance_c1, contrasts = contrasts_c1),
       "results/simulation_floor_bootstrap_replicates_c1.rds", compress = "xz")
say("\n[8] wrote results/simulation_floor_sweep_c1.tsv, results/simulation_floor_contrasts_c1.tsv, results/simulation_floor_choice_blind_c1.tsv, results/simulation_floor_bootstrap_replicates_c1.rds\n")

## ---- mandatory check 9: existing full-grid floor outputs untouched ----
raw_after_checksum <- digest::digest("results/simulation_floor_sweep.tsv", algo = "sha256", file = TRUE)
say("[9] full-grid results/simulation_floor_sweep.tsv unchanged: %s\n", identical(raw_before_checksum, raw_after_checksum))
if (!identical(raw_before_checksum, raw_after_checksum))
  stop("floor c=1 summary: results/simulation_floor_sweep.tsv changed during this run -- this script must never touch it")

write_receipt("16_floor_decomposition_c1", inputs = RAW_FILE,
              params = list(c1_cells = C1_CELLS, floor_grid = FLOOR_GRID_C1, n_bootstrap = N_BOOTSTRAP, boot_seed = C1_BOOT_SEED),
              outputs = c("results/simulation_floor_sweep_c1.tsv", "results/simulation_floor_contrasts_c1.tsv",
                         "results/simulation_floor_choice_blind_c1.tsv", "results/simulation_floor_bootstrap_replicates_c1.rds"))
cat("FLOOR_C1_DONE\n")
