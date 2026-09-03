## =============================================================================
## module_sim_LDscnR / occupancy_and_merge_decomposition.R
##
## The two measurements that overturned the "wide regions consolidate" reading.
## Both were run in scratch during conversation and the manuscript's withdrawal
## rests on them, so they are committed here as one block: they answer one
## question, which is what the wide regions ARE and whether merging helps.
##
## PART 1 -- WHAT THE WIDE REGIONS ARE. Occupancy = markers in a region divided
## by markers physically inside its span; near 1 for a coherent block, near 0 for
## a chain. By span quartile of discovered regions:
##
##   quartile   n    median span   occupancy   precision
##   Q1        69       14.7 kb      0.080       0.304
##   Q2        69      490.9 kb      0.011       0.420
##   Q3        69     3361.4 kb      0.003       0.232
##   Q4        69    12884.9 kb      0.001       0.145
##
##   median occupancy 0.005 | Spearman(span, occupancy) -0.817
##   Spearman(occupancy, precision) +0.28 | 84.1% of regions below 0.05
##
## SO THESE REGIONS ARE CHAINS. The stickleback panel's median is 0.684 with 15%
## below 0.05 -- their discovered regions are cleaner than ours by two orders of
## magnitude, the opposite of what an ends-r2 comparison of the PARTITIONS
## suggested. Ends-r2 describes the partition; occupancy describes the discovered
## subset, and only the second is the relevant object.
##
## PART 2 -- WHY THE MERGE SWEEP LOOKED LIKE GOOD NEWS. Precision rises as the
## merge distance loosens (0.386 -> 0.447), which reads as wide regions being
## consolidation. Decomposing the ratio:
##
##   dcap    regions   precision   TRUE-POSITIVE REGIONS   recall
##   1e4       28.9      0.386            4.58             0.250
##   3e4       24.6      0.401            4.39             0.243
##   1e5       15.9      0.428            4.23             0.241
##   3e5       10.3      0.447            4.00             0.234
##
## THE NUMERATOR FALLS. Precision rises only because the region count falls
## 2.8-fold while true positives fall 13%. Paired, loosening loses true positives
## on 11 panels against 3 (p = 0.057). Merging buys a ratio by spending
## resolution.
##
## THE GENERAL POINT, which is what belongs in the paper: REGION-LEVEL PRECISION
## COMPARES TWO PARTITIONS ONLY IF THEIR GRANULARITY IS MATCHED. Same failure as
## comparing region COUNTS across assembly rules, and easy to miss because the
## ratio moves in the flattering direction. Replication does not catch it -- the
## sweep was 80 panels and the direction was solid; only decomposing the ratio
## does.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/occupancy_and_merge_decomposition.R
## Env: SIM_DATA, CELLS, TAG, ENVS, FILES, SWEEP_CSV, ALPHA
## =============================================================================
suppressMessages({library(data.table); library(LDscnR)})
SIM   <- Sys.getenv("SIM_DATA", "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5")
CELLS <- strsplit(Sys.getenv("CELLS", "V0.5_c1,V0.5_c2,V2_c1"), ",")[[1]]
TAG   <- Sys.getenv("TAG", "nobgs")
ENVS  <- as.integer(strsplit(Sys.getenv("ENVS", "1,2,3"), ",")[[1]])
FILES <- as.integer(strsplit(Sys.getenv("FILES", "1,2,3,4,5,6,7,8,9,10"), ",")[[1]])
ALPHA <- as.numeric(Sys.getenv("ALPHA", "0.05"))
SWEEP <- Sys.getenv("SWEEP_CSV", "module_sim_LDscnR/results/chaining_vs_dcap.csv")
simes <- function(p) { p <- sort(p[is.finite(p)]); n <- length(p)
  if (!n) return(NA_real_); min(n * p / seq_len(n)) }

## ---- Part 1: occupancy of discovered regions, by span quartile
acc <- list()
for (CELL in CELLS) for (ENV in ENVS) {
  D <- list(); LK <- list()
  for (i in FILES) {
    f <- file.path(SIM, sprintf("adapt_%s_chr%d_%s_env%d.rds", TAG, i, CELL, ENV))
    if (!file.exists(f)) next
    x <- readRDS(f); m <- flag_true_qtns(as.data.table(x$map))
    pr <- ld_prune_and_eMLG(GTs = x$GTs, stage1 = x$complexity_reduction$stage1,
          LD_decay = x$LD_decay, ld_w_col = "ld_w_095", ld_w_threshold = 0.025,
          score_threshold = 0.80, min_r2_rho = 0.5, distance_threshold = 1e5,
          compute_unflagged_eMLG = TRUE, cores = 1)
    s2 <- ld_group_map(pr, prefix = i)[, .(marker, CL2 = group_id)]
    mm <- merge(m, s2, by = "marker", all.x = TRUE)[!is.na(CL2)][, Chr := as.character(Chr)]
    th <- score_thresholds(as.data.table(x$LD_decay$decay_sum),
                           rho_r2 = 0.75, rho_d = 0.95, dmax_cap = 1e5)
    drv <- mm[true_pos_QTN %in% TRUE]
    if (nrow(drv)) LK[[length(LK)+1]] <- rbindlist(lapply(seq_len(nrow(drv)), function(j) {
      near <- mm[Chr == drv$Chr[j] & abs(Pos - drv$Pos[j]) < th$dmax]
      if (!nrow(near)) return(NULL)
      r2 <- suppressWarnings(cor(x$GTs[, drv$marker[j]], x$GTs[, near$marker],
                                 use = "pairwise.complete.obs")^2)
      d <- data.table(CL2 = near$CL2, r2 = as.numeric(r2))[is.finite(r2) & r2 >= th$r2min]
      if (!nrow(d)) NULL else unique(d[, .(CL2)]) }))
    D[[length(D)+1]] <- mm[, .(marker, CL2, Chr, Pos, p = emx_p)]
  }
  DD <- rbindlist(D); tagged <- unique(rbindlist(LK)$CL2)
  U <- DD[, .(p = simes(p), n = .N, lo = min(Pos), hi = max(Pos), Chr = Chr[1]), by = CL2]
  ## markers physically inside each region's span, whether or not they are members
  U[, in_span := DD[.SD, on = .(Chr = Chr, Pos >= lo, Pos <= hi), .N, by = .EACHI]$N]
  ok <- is.finite(U$p); q <- rep(NA_real_, nrow(U)); q[ok] <- p.adjust(U$p[ok], "BH")
  sig <- U[which(!is.na(q) & q < ALPHA)]
  if (!nrow(sig)) next
  acc[[length(acc)+1]] <- sig[, .(cell = CELL, env = ENV, CL2, span = hi - lo,
                                  occ = n / pmax(in_span, 1), tp = CL2 %in% tagged)]
}
R <- rbindlist(acc)
R[, qt := cut(span, quantile(span, 0:4/4), include.lowest = TRUE, labels = paste0("Q", 1:4))]
cat("== PART 1: occupancy of discovered regions, by span quartile\n")
print(R[, .(n = .N, median_span_kb = round(median(span)/1e3, 1),
            occupancy = round(median(occ), 3), precision = round(mean(tp), 3)),
        by = qt][order(qt)])
cat(sprintf("\n  median occupancy %.3f | Spearman(span, occupancy) %.3f | Spearman(occupancy, precision) %+.3f\n",
            median(R$occ), stats::cor(R$span, R$occ, method = "spearman"),
            stats::cor(R$occ, as.numeric(R$tp), method = "spearman")))
cat(sprintf("  regions with occupancy < 0.05: %.1f%%\n", 100 * mean(R$occ < 0.05)))

## ---- Part 2: decompose the merge-distance precision gain
if (file.exists(SWEEP)) {
  C <- fread(SWEEP)[, tp_count := precision * regions]
  cat("\n== PART 2: does looser merging find more truth, or fewer regions?\n")
  print(C[, .(panels = .N, regions = round(mean(regions), 1),
              precision = round(mean(precision), 3),
              TP_regions = round(mean(tp_count), 2),
              recall = round(mean(recall), 3)), by = dcap][order(dcap)])
  W <- dcast(C, cell + tag + env ~ dcap, value.var = "tp_count")
  a <- W[["1e+05"]]; b <- W[["3e+05"]]; ok <- is.finite(a) & is.finite(b)
  cat(sprintf("\n  loosening 1e5 -> 3e5, TRUE-POSITIVE REGIONS: %d W / %d L / %d T  p = %.4f\n",
              sum(b[ok] > a[ok]), sum(b[ok] < a[ok]), sum(b[ok] == a[ok]),
              stats::binom.test(sum(b[ok] > a[ok]), sum(b[ok] != a[ok]))$p.value))
} else cat("\n(sweep CSV not found; run chaining_vs_dcap.R first)\n")
