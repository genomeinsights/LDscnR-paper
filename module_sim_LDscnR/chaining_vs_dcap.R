## =============================================================================
## module_sim_LDscnR / chaining_vs_dcap.R
##
## THE DEFECT THIS TARGETS was found today and has not been characterised:
## median discovered-region span is 2.2 Mb, and the widest two span quintiles
## (8.5 and 14.6 Mb) run at 0.10-0.13 precision against 0.28-0.36 below. A
## substantial share of discoveries are megabase-scale chained objects that tag
## nothing. The panel session hit the same thing and fixed it by changing the
## test unit; the untested alternative is that the CHAINING ITSELF is
## controllable, because distance_threshold is the parameter that builds runs.
##
## Sweeping it asks a question no other run here has: does a tighter merge
## distance remove the chained objects WITHOUT costing recall? If it does, the
## defect is a parameter choice. If precision and recall move together, chaining
## is intrinsic to the partition and only the reordering can address it.
##
## Reported per dcap: region spans, the precision of the widest quintile
## specifically, overall precision and recall, and the occupancy of discovered
## regions -- so the mechanism is visible rather than inferred from the outcome.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/chaining_vs_dcap.R
## Env: SIM_DATA, CELLS, TAGS, ENVS, FILES, DCAPS, ALPHA, CORES, OUT
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(parallel)})
SIM   <- Sys.getenv("SIM_DATA", "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5")
CELLS <- strsplit(Sys.getenv("CELLS", "V0.5_c1,V0.5_c2,V1_c1.5,V2_c1"), ",")[[1]]
TAGS  <- strsplit(Sys.getenv("TAGS", "nobgs,bgs"), ",")[[1]]
ENVS  <- as.integer(strsplit(Sys.getenv("ENVS", "1,2,3,4,5"), ",")[[1]])
FILES <- as.integer(strsplit(Sys.getenv("FILES", "1,2,3,4,5,6,7,8,9,10"), ",")[[1]])
DCAPS <- as.numeric(strsplit(Sys.getenv("DCAPS", "1e4,3e4,1e5,3e5"), ",")[[1]])
ALPHA <- as.numeric(Sys.getenv("ALPHA", "0.05")); CORES <- as.integer(Sys.getenv("CORES", "8"))
OUTF  <- Sys.getenv("OUT", "module_sim_LDscnR/results/chaining_vs_dcap.csv")
dir.create(dirname(OUTF), recursive = TRUE, showWarnings = FALSE)
simes <- function(p) { p <- sort(p[is.finite(p)]); n <- length(p)
  if (!n) return(NA_real_); min(n * p / seq_len(n)) }

one <- function(CELL, TAG, ENV) {
  per <- list()
  for (dc in DCAPS) {
    D <- list(); LK <- list()
    for (i in FILES) {
      f <- file.path(SIM, sprintf("adapt_%s_chr%d_%s_env%d.rds", TAG, i, CELL, ENV))
      if (!file.exists(f)) next
      x <- readRDS(f); m <- flag_true_qtns(as.data.table(x$map))
      pr <- ld_prune_and_eMLG(GTs = x$GTs, stage1 = x$complexity_reduction$stage1,
            LD_decay = x$LD_decay, ld_w_col = "ld_w_095", ld_w_threshold = 0.025,
            score_threshold = 0.80, min_r2_rho = 0.5, distance_threshold = dc,
            compute_unflagged_eMLG = TRUE, cores = 1)
      s2 <- ld_group_map(pr, prefix = i)[, .(marker, CL2 = group_id)]
      mm <- merge(m, s2, by = "marker", all.x = TRUE)[!is.na(CL2)]
      th <- score_thresholds(as.data.table(x$LD_decay$decay_sum),
                             rho_r2 = 0.75, rho_d = 0.95, dmax_cap = 1e5)
      drv <- mm[true_pos_QTN %in% TRUE]
      if (nrow(drv)) LK[[length(LK)+1]] <- rbindlist(lapply(seq_len(nrow(drv)), function(j) {
        near <- mm[as.character(Chr) == as.character(drv$Chr[j]) & abs(Pos - drv$Pos[j]) < th$dmax]
        if (!nrow(near)) return(NULL)
        r2 <- suppressWarnings(cor(x$GTs[, drv$marker[j]], x$GTs[, near$marker],
                                   use = "pairwise.complete.obs")^2)
        d <- data.table(CL2 = near$CL2, r2 = as.numeric(r2))[is.finite(r2) & r2 >= th$r2min]
        if (!nrow(d)) return(NULL)
        d[, .(r2 = max(r2)), by = CL2][, qtn := paste0(i, "_", drv$marker[j])][] }))
      D[[length(D)+1]] <- mm[, .(marker, CL2, Pos, Chr = as.character(Chr), p = emx_p)]
    }
    DD <- rbindlist(D); lk <- rbindlist(LK, fill = TRUE); nq <- uniqueN(lk$qtn)
    if (!nq) next
    U <- DD[, .(p = simes(p), n = .N, span = max(Pos) - min(Pos)), by = CL2]
    ok <- is.finite(U$p); q <- rep(NA_real_, nrow(U)); q[ok] <- p.adjust(U$p[ok], "BH")
    sig <- U[which(!is.na(q) & q < ALPHA)]
    if (!nrow(sig)) next
    sub <- lk[CL2 %in% sig$CL2]
    best <- if (nrow(sub)) sub[order(-r2)][, .SD[1], by = qtn] else sub
    kept <- setdiff(sig$CL2, setdiff(sub$CL2, best$CL2))
    sig[, tp := CL2 %in% best$CL2]
    wq <- if (nrow(sig) >= 5) sig[span >= quantile(span, 0.8)] else sig
    per[[length(per)+1]] <- data.table(cell = CELL, tag = TAG, env = ENV, dcap = dc,
      n_units = nrow(U), regions = length(kept),
      precision = uniqueN(best$CL2) / max(length(kept), 1),
      recall = uniqueN(best$qtn) / nq,
      med_span_kb = median(sig$span) / 1e3, max_span_Mb = max(sig$span) / 1e6,
      widest_q_precision = mean(wq$tp), pct_over_1Mb = 100 * mean(sig$span > 1e6),
      n_qtn = nq)
  }
  rbindlist(per)
}
G <- CJ(cell = CELLS, tag = TAGS, env = ENVS, sorted = FALSE)
cat(sprintf("%d panels x %d dcap levels | CORES %d\n", nrow(G), length(DCAPS), CORES))
R <- rbindlist(Filter(Negate(is.null), mclapply(seq_len(nrow(G)),
  function(k) tryCatch(one(G$cell[k], G$tag[k], G$env[k]),
    error = function(e) { cat("FAIL", G$cell[k], G$tag[k], G$env[k], conditionMessage(e), "\n"); NULL }),
  mc.cores = CORES)))
fwrite(R, OUTF)
cat("\n== chaining against merge distance\n")
print(R[, .(panels = .N, regions = round(mean(regions), 1),
            precision = round(mean(precision), 3), recall = round(mean(recall), 3),
            med_span_kb = round(mean(med_span_kb)), max_span_Mb = round(mean(max_span_Mb), 1),
            widest_q_prec = round(mean(widest_q_precision), 3),
            pct_over_1Mb = round(mean(pct_over_1Mb), 1)), by = dcap][order(dcap)])
