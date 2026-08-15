## module_sim/06a_run_caller.R
## EXPENSIVE step (run once per condition). Everything that needs the genotype
## matrix happens here; the result is a small cache that 06b (scoring) and 07
## (plots) read with NO genotypes, so design iteration is seconds.
##   1. pool the 10 replicates -> 20-chromosome genome
##   2. outlier SNP sets: single-SNP EMMAX/LFMM (BH-FDR, no null) and the ld_w
##      method (ld_w-filter + cluster null; its significant-cluster driver SNPs)
##   3. shared region frame: cluster the UNION of all outlier SNPs once
##   4. QTN LD table + TP/FP labels for the regions
## Cache (git-ignored): module_sim/cache_V{V}_c{c}_env{env}.rds
## Run from LDscnR-paper/:  Rscript module_sim/06a_run_caller.R [V c env]

source("module_sim/R/_config.R")
a   <- commandArgs(trailingOnly = TRUE)
V   <- if (length(a) >= 1) a[1] else "2"
cc  <- if (length(a) >= 2) a[2] else "1"
env <- if (length(a) >= 3) a[3] else "1"
FDR <- 0.05
RMSC_GRID <- seq(0, 0.98, by = 0.02)          # coarser than 0.005: ~4x fewer FDR passes

P   <- pool_group(sprintf("^adapt_bgs_chr[0-9]+_V%s_c%s_env%s\\.rds$", V, cc, env))
map <- flag_true_positive_QTNs(as.data.table(P$map))
th  <- score_thresholds(P$decay$decay_sum)
totQTN <- map[true_pos_QTN == TRUE, .N]
cat(sprintf("V%s_c%s_env%s: %d SNPs, %d chr, true_pos_QTN=%d | TP-match r2>%.2f dist<%.0fkb\n",
    V, cc, env, nrow(map), uniqueN(map$Chr), totQTN, th$r2min, th$dmax / 1e3))

## ld_w method: outlier SNPs = FDR-significant SNPs that land in a SIGNIFICANT
## cluster (its actual calls; NOT the LD-recruited neighbours).
ldw_outliers <- function(pcol) {
  r <- ld_outlier_clusters(setNames(map[[pcol]], map$marker), P$ldw[map$marker],
       map[, .(marker, Chr, Pos)], P$GTs, null = "background",
       rho_ld = 0.95, rho_d = 0.95, LD_decay = P$decay,
       rmsc_grid = RMSC_GRID, B = 1000, cores = 4, verbose = FALSE)
  sig_cl <- as.data.table(r$clusters)[significant == TRUE, cluster]
  as.data.table(r$candidates)[sig == TRUE & cluster %in% sig_cl, marker]
}
single_outliers <- function(pcol) map$marker[p.adjust(map[[pcol]], "BH") < FDR]

sets <- list(
  "EMMAX single" = single_outliers("emx_p"),
  "LFMM single"  = single_outliers("lfmm_p"),
  "EMMAX ld_w"   = ldw_outliers("emx_p"),
  "LFMM ld_w"    = ldw_outliers("lfmm_p"))
cat("outlier SNPs per method:",
    paste(sprintf("%s=%d", names(sets), lengths(sets)), collapse = "  "), "\n")

## shared region frame (cluster the union once) + QTN LD table + TP/FP labels
union_mk <- unique(unlist(sets))
dcap     <- d_from_rho(median(P$decay$decay_sum$a_pred), 0.95)
regions  <- cluster_regions(union_mk, map, P$GTs, r2_link = 0.5, dcap = dcap)
qtab     <- precompute_QTN_LD(P$GTs, map, union_mk, 2e6, cores = 4)
lab      <- classify_ORs(regions, map, qtab, th$r2min, th$dmax)
cat(sprintf("shared frame: %d outlier SNPs -> %d regions (%d TP matching %d distinct QTN)\n",
    length(union_mk), length(regions), sum(lab$is_TP), uniqueN(lab[is_TP == TRUE, qtn])))

## compact per-SNP map (positions + p-values + truth) -- lets 06b/07 skip
## genotypes. `type` + `true_pos_QTN` are needed by classify_ORs()/evaluate_ORs().
mapc <- map[, .(marker, Chr, Pos, emx_p, lfmm_p, type, true_pos_QTN)]
saveRDS(list(sets = sets, regions = regions, lab = lab, qtab = qtab, map = mapc,
             th = th, dcap = dcap, total_true_QTN = totQTN,
             params = list(V = V, c = cc, env = env, FDR = FDR, rmsc_grid = RMSC_GRID)),
        file.path(dir_data, sprintf("cache_V%s_c%s_env%s.rds", V, cc, env)))
cat("wrote cache_V", V, "_c", cc, "_env", env, ".rds\n", sep = "")
