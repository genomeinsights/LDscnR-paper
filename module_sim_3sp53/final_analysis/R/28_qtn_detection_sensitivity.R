## =============================================================================
## final_analysis/R/28_qtn_detection_sensitivity.R
##
## Review request (Frederic, 2026-09-25): instead of only the "detectable QTN"
## used for recall (MAF > 0.10 and >= 5% of the chromosome's additive
## variance), keep every analysed QTN and show how detection depends on its
## variance share and MAF, plus the distribution of variance shares.
##
## Primary c=1 runs (600). For every QTN in the analysed marker set (MAF >
## 0.10; QTN below that were removed by the marker filter and cannot be
## scored), the share of additive variance is computed exactly as
## flag_true_qtns() does. A QTN counts as detected by a method when any
## constituent marker of that method's reported regions (canonical floor 2,
## results/simulation_stage2_region_details.rds) matches it under the same
## LD and distance rule used for truth scoring (score_thresholds()), now
## applied to all QTN rather than only the detectable ones.
## Output: results/simulation_qtn_detection_c1.tsv (one row per QTN per run).
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(parallel)})
MODULE_ROOT <- file.path(Sys.getenv("LDSCNR_PAPER_ROOT", "."), "module_sim_3sp53/final_analysis")
setwd(MODULE_ROOT)
source(file.path(MODULE_ROOT, "R", "00_config.R"))
say("=== 28_qtn_detection_sensitivity ===\n\n")

C1_CELLS <- c("V0.5_c1", "V1_c1", "V2_c1")
METHODS <- c("emmax_snp_region", "emmax_simes_region", "emmax_consensus_region", "lfmm_snp_region", "lfmm_simes_region")
rd_all <- as.data.table(readRDS("results/simulation_stage2_region_details.rds"))[cell %chin% C1_CELLS]
combos <- CJ(tag = TAGS_ALL, cell = C1_CELLS, rep = REPS_ALL, env = ENVS_ALL)

one_combo <- function(tag, cell, rep, env) {
  combo_id <- sprintf("%s_%s_rep%d_env%d", tag, cell, rep, env)
  b <- readRDS(file.path(stage_dir("02_build_ld_units", combo_id), "ld_units.rds"))
  GTs <- b$GTs; map <- copy(as.data.table(b$map))
  p <- colSums(GTs) / nrow(GTs) / 2
  map[, p_freq := p[match(marker, colnames(GTs))]]
  map[, Va := 2 * p_freq * (1 - p_freq) * allelic_values^2]
  map <- flag_true_qtns(map, va_col = "Va", maf_col = "MAF", maf_min = MAF_KEEP, p_va_min = VA_SHARE_DETECTABLE)
  q <- map[type == "QTN" & MAF > MAF_KEEP & !is.na(Va) & Va > 0,
           .(qtn_marker = marker, Chr, MAF, Va, p_Va = p_Va_ok, detectable = true_pos_QTN)]
  if (!nrow(q)) return(NULL)
  thr <- score_thresholds(b$LD_decay$decay_sum, rho_r2 = TRUTH_RHO_R2, rho_d = TRUTH_RHO_D, dmax_cap = TRUTH_DMAX_CAP)
  lut <- qtn_ld_table(GTs, map, candidate_markers = map$marker, cores = 1)
  lut <- lut[r2 > thr$r2min & dist_bp < thr$dmax & qtn_marker %chin% q$qtn_marker]
  tg <- tag; ce <- cell; rp <- rep; ev <- env
  rd <- rd_all[tag == tg & cell == ce & rep == rp & env == ev]
  for (m in METHODS) {
    mk <- unique(unlist(rd[method == m, member_markers], use.names = FALSE))
    hit <- unique(lut[marker %chin% mk, qtn_marker])
    q[, (m) := qtn_marker %chin% hit]
  }
  q[, `:=`(tag = tag, cell = cell, rep = rep, env = env)][]
}

t0 <- Sys.time()
res <- mclapply(seq_len(nrow(combos)), function(i) {
  cb <- combos[i]
  tryCatch(one_combo(cb$tag, cb$cell, cb$rep, cb$env),
           error = function(e) stop(sprintf("%s_%s_rep%d_env%d: %s", cb$tag, cb$cell, cb$rep, cb$env, conditionMessage(e))))
}, mc.cores = as.integer(Sys.getenv("R28_CORES", "10")))
bad <- vapply(res, inherits, logical(1), "try-error")
if (any(bad)) stop(paste(unlist(res[bad]), collapse = "\n"))
dt <- rbindlist(res, fill = TRUE)
say("[1] %d QTN across %d runs in %.1f min\n", nrow(dt), uniqueN(dt[, .(tag, cell, rep, env)]),
    as.numeric(difftime(Sys.time(), t0, units = "mins")))

## validation: among detectable QTN, recall must reproduce the primary analysis
perf <- fread("results/simulation_performance_c1.tsv")
chk <- rbindlist(lapply(METHODS, function(m) data.table(method = m, recall_here = mean(dt[detectable == TRUE][[m]]))))
chk <- merge(chk, perf[, .(method, recall)], by = "method")
print(chk)
stopifnot(isTRUE(all.equal(chk$recall_here, chk$recall, tolerance = 1e-9)))
say("[2] recall among detectable QTN reproduces simulation_performance_c1.tsv exactly\n")

fwrite(dt, "results/simulation_qtn_detection_c1.tsv", sep = "\t")
say("[3] wrote results/simulation_qtn_detection_c1.tsv\n")
