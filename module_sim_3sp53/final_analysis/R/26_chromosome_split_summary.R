## =============================================================================
## final_analysis/R/26_chromosome_split_summary.R
##
## PK (2026-10-09): results split into the QTN chromosome (Chr1) and the
## neutral chromosome (Chr2; no QTN of non-zero effect -- see
## neutral_chromosomes() in helpers_stage2_truth.R). Every region or
## significant marker on the neutral chromosome is a false positive by
## construction, independent of the truth-matching distance rule; on the QTN
## chromosome some false regions are in LD with a QTN but scored false
## because they are too distant from it. The main results keep both
## chromosomes exactly as before; this script writes NEW files only.
##
## Primary c=1 runs only (600). Inputs:
##   results/simulation_stage2_region_details.rds      (R/11; canonical floor 2)
##   results/simulation_summary_full.rds               (R/06; per-run recall inputs)
##   results/floor_stability_density_half_double_region_audit.tsv (R/23)
##   out_final_v1/27_floor_decomposition_raw_c1_chrsplit.rds      (R/27, optional)
##   results/simulation_marker_qtn_redundancy_raw_c1_chrsplit.tsv (R/27, optional)
## Uncertainty: 2,000 bootstrap resamples of the 10 replicate maps, keeping all
## cells, BGS treatments and environmental continuations of a sampled map
## together (as in the primary analysis); paired contrasts use the same draws.
## =============================================================================
suppressMessages(library(data.table))
MODULE_ROOT <- path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis")
setwd(MODULE_ROOT)
source(file.path(MODULE_ROOT, "R", "00_config.R"))
say("=== 26_chromosome_split_summary ===\n\n")

C1_CELLS <- c("V0.5_c1", "V1_c1", "V2_c1")
NTRL_CHR <- "Chr2"; QTN_CHR <- "Chr1"   ## verified for every local c=1 run: all detectable QTN on Chr1
B <- 2000L; SEED <- SEEDS[["bootstrap"]] + 26000L
R <- function(x) file.path(MODULE_ROOT, "results", x)
q <- function(x, p) unname(stats::quantile(x, p, na.rm = TRUE))

## ---- bootstrap helper: rep-level sums -> B resampled pooled statistics -------
## `rep_dt`: one row per rep (all 10 present) with numeric count columns;
## `stat_fun(sums)` maps a named list of pooled counts to a named vector.
set.seed(SEED)
DRAWS <- matrix(sample.int(10L, 10L * B, replace = TRUE), nrow = 10L)
MULT <- apply(DRAWS, 2, tabulate, nbins = 10L)   ## 10 x B
boot_counts <- function(rep_dt, cols) {
  rep_dt <- rep_dt[order(rep)]; stopifnot(identical(rep_dt$rep, 1:10))
  stats::setNames(lapply(cols, function(cc) as.numeric(crossprod(MULT, rep_dt[[cc]]))), cols)
}

## ---- 1. primary performance by chromosome ------------------------------------
rd <- as.data.table(readRDS(R("simulation_stage2_region_details.rds")))[cell %chin% C1_CELLS]
stopifnot(all(rd[Chr == NTRL_CHR, TP] == FALSE))
raw <- as.data.table(readRDS(R("simulation_summary_full.rds"))$raw)[cell %chin% C1_CELLS & grepl("_region$", method)]
per_run <- rd[, .(TP_qtnchr = sum(TP & Chr == QTN_CHR), FP_qtnchr = sum(!TP & Chr == QTN_CHR),
                  FP_ntrl = sum(Chr == NTRL_CHR)), by = .(tag, cell, rep, env, method)]
per_run <- merge(raw[, .(tag, cell, rep, env, method, TP, FP, n_detectable_qtn, n_recovered)], per_run,
                 by = c("tag", "cell", "rep", "env", "method"), all.x = TRUE)
for (cc in c("TP_qtnchr", "FP_qtnchr", "FP_ntrl")) per_run[is.na(get(cc)), (cc) := 0L]
stopifnot(nrow(per_run) == 600L * 5L, all(per_run$TP == per_run$TP_qtnchr),
          all(per_run$FP == per_run$FP_qtnchr + per_run$FP_ntrl))
say("[1] per-run counts for %d runs x %d methods; region-level split reproduces the saved TP/FP exactly\n",
    uniqueN(per_run[, .(tag, cell, rep, env)]), uniqueN(per_run$method))

perf_rows <- list(); boot_store <- list()
for (m in unique(per_run$method)) {
  rp <- per_run[method == m, lapply(.SD, sum), by = rep,
                .SDcols = c("TP", "FP", "TP_qtnchr", "FP_qtnchr", "FP_ntrl", "n_detectable_qtn", "n_recovered")]
  bc <- boot_counts(rp, setdiff(names(rp), "rep"))
  n_runs <- per_run[method == m, .N]
  b <- data.table(precision_all = bc$TP / pmax(bc$TP + bc$FP, 1),
                  precision_qtnchr = bc$TP_qtnchr / pmax(bc$TP_qtnchr + bc$FP_qtnchr, 1),
                  fp_ntrl_per_run = bc$FP_ntrl / n_runs, fp_qtnchr_per_run = bc$FP_qtnchr / n_runs,
                  share_fp_ntrl = bc$FP_ntrl / pmax(bc$FP, 1),
                  recall = bc$n_recovered / pmax(bc$n_detectable_qtn, 1))
  boot_store[[m]] <- b
  s <- rp[, lapply(.SD, sum), .SDcols = -"rep"]
  est <- c(precision_all = s$TP / (s$TP + s$FP), precision_qtnchr = s$TP_qtnchr / (s$TP_qtnchr + s$FP_qtnchr),
           fp_ntrl_per_run = s$FP_ntrl / n_runs, fp_qtnchr_per_run = s$FP_qtnchr / n_runs,
           share_fp_ntrl = s$FP_ntrl / s$FP, recall = s$n_recovered / s$n_detectable_qtn)
  perf_rows[[m]] <- rbindlist(lapply(names(est), function(k) data.table(method = m, statistic = k,
    estimate = est[[k]], ci_lo = q(b[[k]], 0.025), ci_hi = q(b[[k]], 0.975))))[,
    `:=`(TP = s$TP, FP = s$FP, TP_qtnchr = s$TP_qtnchr, FP_qtnchr = s$FP_qtnchr, FP_ntrl = s$FP_ntrl, n_runs = n_runs)]
}
perf <- rbindlist(perf_rows)
saved <- fread(R("simulation_performance_c1.tsv"))
chk <- merge(perf[statistic == "precision_all", .(method, estimate)], saved[, .(method, precision)], by = "method")
stopifnot(isTRUE(all.equal(chk$estimate, chk$precision)))
say("[2] pooled precision over both chromosomes reproduces simulation_performance_c1.tsv\n")
fwrite(perf, R("simulation_chrsplit_performance_c1.tsv"), sep = "\t")

## paired contrasts: aggregation vs marker-wise, same bootstrap draws
pairs <- list(c("emmax_simes_region", "emmax_snp_region"), c("emmax_consensus_region", "emmax_snp_region"),
              c("lfmm_simes_region", "lfmm_snp_region"))
contr <- rbindlist(lapply(pairs, function(pr) {
  a <- boot_store[[pr[1]]]; z <- boot_store[[pr[2]]]
  ea <- perf[method == pr[1]]; ez <- perf[method == pr[2]]
  g <- function(d, k) d[statistic == k, estimate]
  rbindlist(list(
    data.table(method = pr[1], versus = pr[2], statistic = "ratio_fp_ntrl_per_run",
               estimate = g(ea, "fp_ntrl_per_run") / g(ez, "fp_ntrl_per_run"),
               ci_lo = q(a$fp_ntrl_per_run / z$fp_ntrl_per_run, 0.025), ci_hi = q(a$fp_ntrl_per_run / z$fp_ntrl_per_run, 0.975)),
    data.table(method = pr[1], versus = pr[2], statistic = "ratio_fp_qtnchr_per_run",
               estimate = g(ea, "fp_qtnchr_per_run") / g(ez, "fp_qtnchr_per_run"),
               ci_lo = q(a$fp_qtnchr_per_run / z$fp_qtnchr_per_run, 0.025), ci_hi = q(a$fp_qtnchr_per_run / z$fp_qtnchr_per_run, 0.975)),
    data.table(method = pr[1], versus = pr[2], statistic = "diff_precision_qtnchr",
               estimate = g(ea, "precision_qtnchr") - g(ez, "precision_qtnchr"),
               ci_lo = q(a$precision_qtnchr - z$precision_qtnchr, 0.025), ci_hi = q(a$precision_qtnchr - z$precision_qtnchr, 0.975))))
}))
fwrite(contr, R("simulation_chrsplit_contrasts_c1.tsv"), sep = "\t")
say("\n[3] performance by chromosome (c=1, floor 2):\n")
print(dcast(perf, method ~ statistic, value.var = "estimate"), digits = 3)
print(contr, digits = 3)

## ---- 2. floor-stability tiers by chromosome ------------------------------------
aud <- fread(R("floor_stability_density_half_double_region_audit.tsv"))
stopifnot(all(aud[Chr == NTRL_CHR, TP] == FALSE))
tiers <- aud[, .(n = .N, n_ntrl = sum(Chr == NTRL_CHR), share_ntrl = mean(Chr == NTRL_CHR),
                 precision_all = mean(TP), precision_qtnchr = mean(TP[Chr == QTN_CHR])),
             by = .(method, tier = tier_effective)][order(method, tier)]
fwrite(tiers, R("simulation_chrsplit_tiers_c1.tsv"), sep = "\t")
say("\n[4] floor-stability tiers by chromosome:\n"); print(tiers, digits = 3)

## ---- 3. floor decomposition by chromosome (R/27) --------------------------------
f27 <- file.path(MODULE_ROOT, "out_final_v1", "27_floor_decomposition_raw_c1_chrsplit.rds")
if (file.exists(f27)) {
  fl <- as.data.table(readRDS(f27))
  fd_rows <- list()
  for (g in split(fl, by = c("method", "arm", "floor"))) {
    rp <- g[, lapply(.SD, sum), by = rep, .SDcols = c("TP_qtnchr", "FP_qtnchr", "FP_ntrl", "TP", "FP")]
    bc <- boot_counts(rp, setdiff(names(rp), "rep")); s <- rp[, lapply(.SD, sum), .SDcols = -"rep"]; n_runs <- nrow(g)
    fd_rows[[length(fd_rows) + 1]] <- data.table(method = g$method[1], arm = g$arm[1], floor = g$floor[1], n_runs = n_runs,
      TP = s$TP, FP = s$FP, FP_ntrl = s$FP_ntrl, FP_qtnchr = s$FP_qtnchr,
      precision_qtnchr = s$TP_qtnchr / max(s$TP_qtnchr + s$FP_qtnchr, 1),
      precision_qtnchr_ci_lo = q(bc$TP_qtnchr / pmax(bc$TP_qtnchr + bc$FP_qtnchr, 1), 0.025),
      precision_qtnchr_ci_hi = q(bc$TP_qtnchr / pmax(bc$TP_qtnchr + bc$FP_qtnchr, 1), 0.975),
      fp_ntrl_per_run = s$FP_ntrl / n_runs,
      fp_ntrl_per_run_ci_lo = q(bc$FP_ntrl / n_runs, 0.025), fp_ntrl_per_run_ci_hi = q(bc$FP_ntrl / n_runs, 0.975))
  }
  fd <- rbindlist(fd_rows)[order(method, floor, arm)]
  fwrite(fd, R("simulation_chrsplit_floor_decomposition_c1.tsv"), sep = "\t")
  say("\n[5] floor decomposition by chromosome, floor 2:\n"); print(fd[floor == 2], digits = 3)
} else say("\n[5] R/27 floor decomposition output not found -- skipped\n")

## ---- 4. marker-wise significant markers by chromosome (R/27) --------------------
r27 <- R("simulation_marker_qtn_redundancy_raw_c1_chrsplit.tsv")
if (file.exists(r27)) {
  red <- fread(r27)
  stopifnot(all(red$n_qtn_linked_ntrl == 0L))
  rs <- red[, .(n_runs = .N, n_significant = sum(n_significant), n_significant_ntrl = sum(n_significant_ntrl),
                n_false_positive = sum(n_false_positive), n_false_positive_ntrl = sum(n_false_positive_ntrl)), by = method]
  rs[, `:=`(share_sig_on_ntrl = n_significant_ntrl / pmax(n_significant, 1),
            share_fp_on_ntrl = n_false_positive_ntrl / pmax(n_false_positive, 1),
            ntrl_sig_per_run = n_significant_ntrl / n_runs)]
  fwrite(rs, R("simulation_chrsplit_marker_redundancy_c1.tsv"), sep = "\t")
  say("\n[6] significant markers by chromosome:\n"); print(rs, digits = 3)
} else say("\n[6] R/27 redundancy output not found -- skipped\n")

say("\nwrote results/simulation_chrsplit_*_c1.tsv\n")
