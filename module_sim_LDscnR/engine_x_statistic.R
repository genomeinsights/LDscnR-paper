## =============================================================================
## module_sim_LDscnR / engine_x_statistic.R
##
## THE CELL THE PANEL SESSION CANNOT FILL. On the stickleback panel LFMM keeps
## Simes because its p-values are precomputed and consensus refits the
## association on a new genotype column, so their comparison is EMMAX/Simes vs
## EMMAX/consensus vs LFMM/Simes -- engine and statistic confounded in the third
## cell. Here both engines can be run on both statistics, because LFMM can be
## refitted on the eMLG consensus matrix (lfmm2 on a written .lfmm file, 4 s per
## chromosome).
##
## Six arms, all scored at the stage-2 region level with the module's usual
## dedup convention so they are directly comparable:
##   stage2 Simes     x {EMMAX, LFMM}     the current default
##   stage2 consensus x {EMMAX, LFMM}     no n penalty
##   stage1 Simes     x {EMMAX, LFMM}     PK's reordering
##
## ALSO CHECKS THE PANEL'S OCCUPANCY DIAGNOSTIC. They found their stage-2
## discoveries were anchored by distance-run objects spanning whole chromosomes:
## occupancy = markers in a unit / markers physically inside its span, which
## approaches 1 for a coherent block and 0 for a chain. Whether bgs5 has the same
## objects decides whether their failure mode is reproducible here at all.
##
## RESULT, 10 panels, paired within panel:
##
##   STAGE 1 vs STAGE 2 (Simes)      recall            precision
##     EMMAX                       5 W / 0 L  p=0.063   2 W / 6 L  ns
##     LFMM                        4 W / 0 L  p=0.125   0 W / 9 L  p=0.004
##   Combined 9 W / 0 L on recall (p = 0.004): THE REORDERING NEVER LOSES RECALL
##   UNDER EITHER ENGINE, at a precision cost that is real under LFMM.
##
##   CONSENSUS vs SIMES (stage 2)    recall            precision
##     EMMAX                       4 W / 3 L  ns        1 W / 9 L  p=0.021
##     LFMM                        9 W / 1 L  p=0.021   0 W /10 L  p=0.002
##   CONSENSUS IS NOT FREE HERE. It buys power only under LFMM and costs
##   precision under BOTH engines. This contradicts the panel's finding of more
##   power at equal precision, and the likely reason is measurement: their
##   precision proxy is EcoPeak overlap, which a larger region satisfies more
##   easily, while QTN-tagging is not size-sensitive in that way.
##
##   ENGINE at fixed statistic      recall 9 W / 1 L p=0.021, precision ns.
##   SO THE FREE LUNCH IS THE ENGINE, NOT THE STATISTIC: LFMM gains power at
##   equal precision, consensus does not.
##
## AND THE PANEL'S DISTANCE-RUN ARTEFACT EXISTS HERE, contrary to what I told
## them earlier. The widest discovered region spans 29.22 Mb on a ~30 Mb
## chromosome, with median occupancy 0.13-0.16 among discoveries. bgs5 does
## produce chromosome-scale low-occupancy objects, so it can speak to their
## failure mode after all.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/engine_x_statistic.R
## Env: SIM_DATA, CELLS, TAG, ENVS, FILES, ALPHA, LFMM_K, CORES
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(LEA); library(parallel)})
SIM   <- Sys.getenv("SIM_DATA", "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5")
CELLS <- strsplit(Sys.getenv("CELLS", "V0.5_c1,V0.5_c2"), ",")[[1]]
TAG   <- Sys.getenv("TAG", "nobgs")
ENVS  <- as.integer(strsplit(Sys.getenv("ENVS", "1,2,3,4,5"), ",")[[1]])
FILES <- as.integer(strsplit(Sys.getenv("FILES", "1,2,3,4,5,6,7,8,9,10"), ",")[[1]])
ALPHA <- as.numeric(Sys.getenv("ALPHA", "0.05")); LFMM_K <- as.integer(Sys.getenv("LFMM_K", "5"))
CORES <- as.integer(Sys.getenv("CORES", "3"))
simes <- function(p) { p <- sort(p[is.finite(p)]); n <- length(p)
  if (!n) return(NA_real_); min(n * p / seq_len(n)) }

one <- function(CELL, ENV) {
  D <- list(); LK <- list(); CE <- list(); CL <- list(); OC <- list()
  for (i in FILES) {
    f <- file.path(SIM, sprintf("adapt_%s_chr%d_%s_env%d.rds", TAG, i, CELL, ENV))
    if (!file.exists(f)) next
    x <- readRDS(f); m <- flag_true_qtns(as.data.table(x$map))
    s1 <- as.data.table(x$complexity_reduction$stage1$map_snp)[, .(marker, CL1 = CL_id)]
    pr <- ld_prune_and_eMLG(GTs = x$GTs, stage1 = x$complexity_reduction$stage1,
          LD_decay = x$LD_decay, ld_w_col = "ld_w_095", ld_w_threshold = 0.025,
          score_threshold = 0.80, min_r2_rho = 0.5, distance_threshold = 1e5,
          compute_unflagged_eMLG = TRUE, cores = 1)
    g  <- as.data.table(pr$groups)
    s2 <- ld_group_map(pr, prefix = i)[, .(marker, CL2 = group_id)]
    mm <- merge(merge(m, s1, by = "marker", all.x = TRUE), s2, by = "marker", all.x = TRUE)
    mm <- mm[!is.na(CL1) & !is.na(CL2)][, CL1 := paste0(i, "_", CL1)]
    ## occupancy: unit markers / markers physically inside the unit's span
    oc <- mm[, .(n = .N, lo = min(Pos), hi = max(Pos), Chr = as.character(Chr)[1]), by = CL2]
    oc[, in_span := mm[.SD, on = .(Chr = Chr, Pos >= lo, Pos <= hi), .N, by = .EACHI]$N]
    OC[[length(OC)+1]] <- oc[, .(CL2, n, span = hi - lo, occ = n / pmax(in_span, 1))]
    th  <- score_thresholds(as.data.table(x$LD_decay$decay_sum),
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
    D[[length(D)+1]] <- mm[, .(marker, CL1, CL2, chr_type,
                               p_emx = emx_p, p_lfm = lfmm_p)]
    ## consensus under BOTH engines on the same eMLG matrix
    E <- pr$eMLG
    if (!is.null(E) && ncol(E)) {
      pe <- tryCatch(as.numeric(emmax_fast(emmax_setup(E, x$GRM), as.numeric(x$env$env))),
                     error = function(e) rep(NA_real_, ncol(E)))
      tmp <- tempfile(); dir.create(tmp)
      pl <- tryCatch({
        write.lfmm(E, file.path(tmp, "g.lfmm")); write.env(as.numeric(x$env$env), file.path(tmp, "e.env"))
        pj <- lfmm2(file.path(tmp, "g.lfmm"), file.path(tmp, "e.env"), K = LFMM_K)
        as.numeric(suppressWarnings(lfmm2.test(pj, file.path(tmp, "g.lfmm"),
          file.path(tmp, "e.env"), full = FALSE, genomic.control = TRUE))$pvalues)
      }, error = function(e) rep(NA_real_, ncol(E)))
      unlink(tmp, recursive = TRUE)
      CE[[length(CE)+1]] <- data.table(CL2 = paste0(i, "_", colnames(E)),
                                       p_emx = pe, p_lfm = pl)
    }
  }
  DD <- rbindlist(D); lk <- rbindlist(LK, fill = TRUE); nq <- uniqueN(lk$qtn)
  if (!nq) return(NULL)
  CC <- rbindlist(CE); OCC <- rbindlist(OC); map12 <- unique(DD[, .(CL1, CL2)])
  bh <- function(p) { ok <- is.finite(p); q <- rep(NA_real_, length(p))
                      q[ok] <- p.adjust(p[ok], "BH"); which(!is.na(q) & q < ALPHA) }
  score <- function(regions, n_test, part, stat, eng) {
    sub  <- lk[CL2 %in% regions]
    best <- if (nrow(sub)) sub[order(-r2)][, .SD[1], by = qtn] else sub
    kept <- setdiff(regions, setdiff(sub$CL2, best$CL2))
    tp   <- if (nrow(best)) uniqueN(best$CL2) else 0L
    o    <- OCC[CL2 %in% regions]
    data.table(cell = CELL, env = ENV, partition = part, statistic = stat, engine = eng,
               n_test = n_test, regions = length(kept), tp = tp,
               precision = tp / max(length(kept), 1),
               recall = (if (nrow(best)) uniqueN(best$qtn) else 0L) / nq,
               widest_Mb = if (nrow(o)) max(o$span)/1e6 else NA_real_,
               med_occupancy = if (nrow(o)) median(o$occ) else NA_real_, n_qtn = nq)
  }
  R <- list()
  for (eng in c("EMMAX","LFMM")) {
    pc <- if (eng == "EMMAX") "p_emx" else "p_lfm"
    U2 <- DD[, .(p = simes(get(pc))), by = CL2]
    R[[length(R)+1]] <- score(U2$CL2[bh(U2$p)], nrow(U2), "stage 2", "Simes", eng)
    U1 <- DD[, .(p = simes(get(pc))), by = CL1]
    s1sig <- U1$CL1[bh(U1$p)]
    R[[length(R)+1]] <- score(unique(map12[CL1 %in% s1sig]$CL2), nrow(U1), "stage 1", "Simes", eng)
    if (nrow(CC)) R[[length(R)+1]] <-
      score(CC$CL2[bh(CC[[pc]])], nrow(CC), "stage 2", "consensus", eng)
  }
  rbindlist(R)
}
G <- CJ(cell = CELLS, env = ENVS, sorted = FALSE)
R <- rbindlist(Filter(Negate(is.null), mclapply(seq_len(nrow(G)),
       function(k) tryCatch(one(G$cell[k], G$env[k]), error = function(e) {
         cat("FAIL", G$cell[k], G$env[k], conditionMessage(e), "\n"); NULL }),
       mc.cores = CORES)))
fwrite(R, "module_sim_LDscnR/results/engine_x_statistic.csv")
cat(sprintf("\n%d panels\n\n", uniqueN(R[, .(cell, env)])))
print(R[, .(panels = .N, regions = round(mean(regions), 1), tp = round(mean(tp), 1),
            precision = round(mean(precision), 3), recall = round(mean(recall), 3),
            widest_Mb = round(mean(widest_Mb, na.rm = TRUE), 2),
            occupancy = round(mean(med_occupancy, na.rm = TRUE), 3)),
        by = .(partition, statistic, engine)][order(partition, statistic, engine)])
