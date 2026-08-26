## =====================================================================
## module_sim_LDscnR / operating_point_vs_truth.R
##
## THE SIM-TO-EMPIRICAL BRIDGE. On simulations you can see the truth-PR surface
## over (tau_C, l_min) and read off the best operating point. On real data you
## cannot -- you must choose without truth, from the structure-aware null. This
## asks whether the null-chosen point lands where truth would have put it.
##
## If it does, the sims license the empirical operating point and the method is
## usable on real data. If it does not, the empirical operating point is
## unjustified, which matters more than any PR-AUC margin.
##
## Three candidate rules, all truth-free:
##   null_calib   calibrate_lmin() then calibrate_tauc() -- the pooled count-FDR
##                route. NOTE: framework section 4 demoted this to a diagnostic
##                because it returned tau_C = NA on the stickleback data. The open
##                question is whether it works on SIMS, where the null is quiet.
##   gate_qR      the operating point the empirical pipeline actually uses:
##                tau = 0.05, l_min = 3, regions with q_R < 0.05.
##   fixed        tau = 0.05, l_min = 2, no null at all -- the naive default.
##
## Each is compared against the per-genome truth optimum. Compared PER GENOME,
## because the truth-optimal l_min varies a lot across genomes and an average
## optimum would hide that.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/operating_point_vs_truth.R <scan_dir> [outdir]
## =====================================================================
suppressMessages({ library(data.table); library(LDscnR) })
`%||%` <- function(a, b) if (is.null(a)) b else a

a <- commandArgs(trailingOnly = TRUE)
if (!length(a)) stop("usage: operating_point_vs_truth.R <scan_dir> [outdir]")
SCAN_DIR <- a[1]; OUT <- if (length(a) >= 2) a[2] else SCAN_DIR
PANEL_DIR <- Sys.getenv("PANEL_DIR", "/Volumes/Nemo/Nemo_sim/analysis_inputs")
ENGINE <- Sys.getenv("ENGINE", "emmax"); BASIS <- Sys.getenv("BASIS", "env_orth")
RHO_LD <- 0.75; RHO_D <- 0.95; DCAP <- 5e5; FDR <- 0.05
TAUS  <- seq(0.02, 0.50, by = 0.02)
LMINS <- c(1L, 2L, 3L, 5L, 10L, 20L)
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)

scans <- list.files(SCAN_DIR, pattern = sprintf("^scan_.*_%s_%s_.*[.]rds$", ENGINE, BASIS), full.names = TRUE)
cat(sprintf("[1] %d genome(s), engine %s, basis %s\n", length(scans), ENGINE, BASIS)); flush.console()

res <- list(); k <- 0L
for (f in scans) {
  cell <- sub(sprintf("^scan_(V[0-9.]+_c[0-9.]+_env[0-9]+)_%s_.*$", ENGINE), "\\1", basename(f))
  pf <- file.path(PANEL_DIR, sprintf("panel_%s.rds", cell)); if (!file.exists(pf)) next
  panel <- readRDS(pf); map <- as.data.table(panel$map)
  if (!"true_pos_QTN" %in% names(map)) map <- flag_true_qtns(map)
  if (!sum(map$true_pos_QTN %in% TRUE)) next
  th <- score_thresholds(as.data.table(panel$decay_sum), rho_r2 = RHO_LD, rho_d = RHO_D, dmax_cap = DCAP)
  x <- readRDS(f); C <- x$null$C_obs
  if (!any(C > 0)) { cat(sprintf("  [skip] %s: no C > 0\n", cell)); next }
  qtab <- qtn_ld_table(panel$GTs, map, names(C)[C > 0], 2e6, cores = 1)

  ## --- the truth surface -------------------------------------------------
  surf <- rbindlist(lapply(TAUS, function(tau) {
    mk <- names(C)[which(C >= tau)]
    r_all <- if (length(mk)) ld_regions(mk, x$edges) else list()
    rbindlist(lapply(LMINS, function(L) {
      r <- r_all[lengths(r_all) >= L]
      ev <- if (length(r)) evaluate_ors(r, map, qtab, th$r2min, th$dmax)
            else list(TP=0L, FP=0L, Precision=NA_real_, Recall=0, PR=0)
      data.table(tau, l_min = L, TP = ev$TP, FP = ev$FP,
                 precision = ev$Precision, recall = ev$Recall,
                 PR = ifelse(is.na(ev$PR %||% NA), 0, ev$PR %||% 0))
    }))
  }))
  best <- surf[which.max(PR)]

  ## --- truth-free rules ---------------------------------------------------
  lm_auto <- tryCatch(calibrate_lmin(x$null, x$edges, tau = 0.05, q = 0.99), error = function(e) NA_integer_)
  tau_auto <- tryCatch(calibrate_tauc(x$null, x$edges, l_min = if (is.na(lm_auto)) 2L else lm_auto,
                                      fdr = FDR, tau_grid = TAUS), error = function(e) NA_real_)
  ## look the rule's (tau, l_min) up on the surface. The argument is named
  ## tau_v, not tau: data.table resolves a bare `tau` inside [ ] to the COLUMN,
  ## and the `..` prefix does not disambiguate it in `i`.
  ## calibrate_lmin() can return an l_min off the grid, so it is snapped to the
  ## nearest evaluated value and the snapped value is recorded.
  pick <- function(tau_v, L) {
    na <- data.table(tau = NA_real_, l_min = NA_integer_, PR = NA_real_,
                     precision = NA_real_, recall = NA_real_)
    if (is.na(tau_v) || is.na(L)) return(na)
    L_use <- LMINS[which.min(abs(LMINS - L))]
    s <- surf[l_min == L_use]
    if (!nrow(s)) return(na)
    s[which.min(abs(s$tau - tau_v))][, .(tau, l_min, PR, precision, recall)]
  }
  r_null  <- pick(tau_auto, lm_auto)
  r_gate  <- pick(0.05, 3L)
  r_fixed <- pick(0.05, 2L)

  k <- k + 1L
  res[[k]] <- data.table(cell,
    best_tau = best$tau, best_lmin = best$l_min, best_PR = round(best$PR, 4),
    null_tau = tau_auto, null_lmin = lm_auto,
    null_tau_used = r_null$tau, null_lmin_used = r_null$l_min, null_PR = round(r_null$PR, 4),
    gate_PR = round(r_gate$PR, 4), fixed_PR = round(r_fixed$PR, 4))
  cat(sprintf("  %-14s best (tau %.2f, l_min %2d) PR %.3f | null (tau %s, l_min %s) PR %s | gate PR %.3f\n",
              cell, best$tau, best$l_min, best$PR,
              ifelse(is.na(tau_auto), "NA", sprintf("%.2f", tau_auto)),
              ifelse(is.na(lm_auto), "NA", as.character(lm_auto)),
              ifelse(is.na(r_null$PR), "NA", sprintf("%.3f", r_null$PR)), r_gate$PR)); flush.console()
}
out <- rbindlist(res)
fwrite(out, file.path(OUT, "operating_point_vs_truth.csv"))
cat(sprintf("\n=== how much PR does each truth-free rule give up vs the per-genome optimum? (%d genomes) ===\n", nrow(out)))
print(out[, .(rule = c("null_calib","gate_qR(0.05,3)","fixed(0.05,2)"),
              mean_PR = round(c(mean(null_PR, na.rm=TRUE), mean(gate_PR, na.rm=TRUE), mean(fixed_PR, na.rm=TRUE)), 3),
              mean_gap = round(c(mean(best_PR - null_PR, na.rm=TRUE), mean(best_PR - gate_PR, na.rm=TRUE),
                                 mean(best_PR - fixed_PR, na.rm=TRUE)), 3),
              n_usable = c(sum(!is.na(null_PR)), sum(!is.na(gate_PR)), sum(!is.na(fixed_PR))))])
cat(sprintf("\n  per-genome optimum mean PR: %.3f\n", mean(out$best_PR, na.rm=TRUE)))
cat(sprintf("  truth-optimal l_min varies: %s\n", paste(sort(unique(out$best_lmin)), collapse=", ")))
