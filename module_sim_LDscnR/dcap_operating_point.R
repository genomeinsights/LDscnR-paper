## =============================================================================
## module_sim_LDscnR / dcap_operating_point.R
##
## CORRECTED gap analysis. The earlier version (operating_points.R commentary and
## the paper outline) computed gaps between ADJACENT MARKERS OVER THE WHOLE MAP.
## That is the wrong quantity twice over, as 2c found by reading the source:
##
##   split_by_distance(pos_min, pos_max, dt):
##     gap <- pos_min[ord] - c(NA, pos_max[ord][-length(ord)])
##
## it is the gap between consecutive STAGE-1 CLUSTERS, and it is applied to the
## FLAGGED subset only -- clusters with any member ld_w > ld_w_threshold
## (ld_prune_and_eMLG.R:389, strict >). Flagged clusters are sparse, so gaps
## between them are far larger than gaps between adjacent markers, and the
## earlier numbers understate them by a wide margin.
##
## THE CLAIM AT RISK IS MINE, NOT 2c's. Theirs was "the derived threshold is
## active on the panel" -- direction unchanged by the correction. Mine was "both
## settings sit past the top of the gap distribution, so switching is INERT",
## which rested entirely on the derived value exceeding the largest observed gap
## on 160 of 160 chromosomes. On the correct quantity that may no longer hold.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/dcap_operating_point.R
## Env: SIM_DATA, CELLS, TAGS, ENV, FILES, LDW, CORES, OUT
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(parallel)})
SIM   <- Sys.getenv("SIM_DATA", "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5")
OUT   <- Sys.getenv("OUT", "module_sim_LDscnR/results/operating_points")
CELLS <- strsplit(Sys.getenv("CELLS", "V0.5_c1,V0.5_c2,V1_c1.5,V2_c1"), ",")[[1]]
TAGS  <- strsplit(Sys.getenv("TAGS", "nobgs,bgs"), ",")[[1]]
ENV   <- as.integer(Sys.getenv("ENV", "1"))
FILES <- as.integer(strsplit(Sys.getenv("FILES", "1,2,3,4,5,6,7,8,9,10"), ",")[[1]])
LDW   <- as.numeric(Sys.getenv("LDW", "0.025"))
CORES <- as.integer(Sys.getenv("CORES", "8"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

jobs <- CJ(CELL = CELLS, TAG = TAGS, i = FILES, sorted = FALSE)
one <- function(k) {
  r <- jobs[k]
  f <- sprintf("%s/adapt_%s_chr%d_%s_env%d.rds", SIM, r$TAG, r$i, r$CELL, ENV)
  if (!file.exists(f)) return(NULL)
  x  <- readRDS(f); s1 <- x$complexity_reduction$stage1
  ms <- as.data.table(s1$map_snp); ds <- as.data.table(x$LD_decay$decay_sum)
  ## flagged = any member above the threshold, matching line 389 exactly (strict >)
  flag_ids <- ms[ld_w_095 > LDW, unique(CL_id)]
  ## one flat aggregation, not a per-cluster lookup -- 2c's first attempt hung
  ## on a named-vector match() re-run over the whole map for every cluster.
  ext <- ms[CL_id %in% flag_ids, .(Chr = Chr[1], pmin = min(Pos), pmax = max(Pos)),
            by = CL_id]
  rbindlist(lapply(ds$Chr, function(ch) {
    e <- ext[Chr == ch][order(pmin)]
    if (nrow(e) < 2) return(NULL)
    g <- e$pmin[-1] - e$pmax[-nrow(e)]
    g <- g[is.finite(g)]
    data.table(CELL = r$CELL, TAG = r$TAG, i = r$i, Chr = ch,
               n_flagged = nrow(e), n_gap = length(g),
               med = median(g), q99 = quantile(g, .99), mx = max(g),
               over_1e5 = sum(g > 1e5), over_5e5 = sum(g > 5e5),
               d_der = d_from_rho(ds[Chr == ch, a_pred], 0.95),
               over_der = sum(g > d_from_rho(ds[Chr == ch, a_pred], 0.95))) }))
}
G <- rbindlist(Filter(Negate(is.null), mclapply(seq_len(nrow(jobs)), one, mc.cores = CORES)))
fwrite(G, file.path(OUT, "dcap_operating_point_by_chr.csv"))

S <- G[, .(chr = .N, flagged = sum(n_flagged), gaps = sum(n_gap),
           med_gap = median(med), max_gap = max(mx),
           pct_over_1e5 = 100 * sum(over_1e5) / sum(n_gap),
           pct_over_5e5 = 100 * sum(over_5e5) / sum(n_gap),
           d_derived_kb = median(d_der) / 1e3,
           pct_over_derived = 100 * sum(over_der) / sum(n_gap),
           chr_derived_over_max = sum(d_der > mx)), by = .(CELL, TAG)]
fwrite(S, file.path(OUT, "dcap_operating_point.csv"))
print(S[order(CELL, TAG)])
cat(sprintf("\nPOOLED: %d flagged clusters, %d gaps | median %.0f bp | max %.0f bp\n",
            sum(G$n_flagged), sum(G$n_gap), median(G$med), max(G$mx)))
cat(sprintf("  gaps over 1e5      : %.3f%%\n", 100*sum(G$over_1e5)/sum(G$n_gap)))
cat(sprintf("  gaps over derived  : %.3f%%  (median derived %.0f kb)\n",
            100*sum(G$over_der)/sum(G$n_gap), median(G$d_der)/1e3))
cat(sprintf("  derived exceeds the largest observed gap on %d of %d chromosomes\n",
            sum(G$d_der > G$mx), nrow(G)))
