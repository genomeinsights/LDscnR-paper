## =============================================================================
## module_sim_LDscnR / test_then_cluster.R
##
## PK's proposed reordering. Current pipeline clusters to stage 2 and THEN tests;
## the proposal is to TEST STAGE-1 CLUSTERS FIRST and apply stage-2 clustering
## only to the survivors, using it to assemble regions for reporting rather than
## to define the units that get tested.
##
## THE MOTIVATION IS A FAILURE MODE SEEN ON THE STICKLEBACK PANEL: sparse
## structure-driven clusters spanning whole chromosomes came back significant
## BECAUSE they were sparse, while the Chr1 inversion was missed because it was
## penalised for being large and dense.
##
## THAT MECHANISM IS CHECKABLE AND IT IS ARITHMETIC. Simes is
## min_i (n * p_(i) / i), so a cluster's smallest p-value is multiplied by its
## member count: a 500-marker inversion carrying one strong signal is penalised
## 100-fold harder than a 5-marker sparse artefact carrying the same one. Testing
## at stage 1 evens that out, because stage-1 clusters are far more uniform in
## size than stage-2 groups.
##
## ARMS, both scored at the REGION level so they are comparable:
##   cluster-then-test   stage-2 units -> Simes -> BH.  Regions are the
##                       significant stage-2 clusters. (current)
##   test-then-cluster   stage-1 units -> Simes -> BH.  Regions are the distinct
##                       stage-2 groups containing at least one significant
##                       stage-1 cluster. (proposed)
##   consensus           stage-2 eMLG consensus genotypes -> one test per group
##                       -> BH. PK's point: consensus CARRIES NO n PENALTY, so it
##                       sidesteps the size effect that motivates the reordering.
##   single SNP          markers -> BH, regions likewise assembled by stage-2
##                       membership, as a floor.
##
## RESULT, 15 panels (V0.5_c1, V0.5_c2, V2_c1 x 5 environments), paired:
##
##   test-then-cluster vs cluster-then-test    recall     6 W / 0 L / 9 T  p=0.031
##                                             true pos   6 W / 0 L / 9 T  p=0.031
##                                             precision  2 W / 7 L / 6 T  p=0.18
##
## THE REORDERING NEVER LOSES RECALL AND SOMETIMES GAINS, at a precision cost
## that does not reach significance. And THE GAIN IS CONCENTRATED WHERE THE
## CLUSTERING IS WORST:
##
##   cell       arm                 regions   TP   precision   recall
##   V0.5_c1    cluster-then-test     14.8    9.4    0.654      0.508
##   V0.5_c1    test-then-cluster     17.0    9.8    0.592      0.529
##   V0.5_c2    cluster-then-test     36.6    3.4    0.118      0.101
##   V0.5_c2    test-then-cluster     63.8    6.4    0.114      0.194
##   V2_c1      cluster-then-test      7.8    3.2    0.594      0.264
##   V2_c1      test-then-cluster      8.2    3.4    0.570      0.278
##
## In V0.5_c2 it nearly DOUBLES true positives and recall at unchanged precision.
## That is the cell with the lowest kinship effective rank, the heaviest tail and
## the worst calibration -- so the reordering pays exactly where stage-2 grouping
## is most distorted, which is the shape of PK's stickleback motivation.
##
## THE CONSENSUS ARM ISOLATES THE n PENALTY AND DOES NOT RESCUE ANYTHING.
## Consensus carries no n multiplier at all, so if the Simes size penalty were
## the binding constraint it should win. It does not: recall 5 W / 4 L (p = 1.0),
## true positives 5 W / 4 L (p = 1.0), and precision 1 W / 12 L (p = 0.0034). It
## finds more regions with the same true positives. SO THE PENALTY IS REAL --
## Simes/min-p reaches 13.7x on 21-100-marker clusters -- BUT IT IS NOT WHAT
## LIMITS DETECTION HERE, and removing it costs precision instead.
##
## CAUTION ON TRANSFER. These simulations have no 430 kb inversion. The failure
## mode motivating the reordering appears only faintly: detection of QTN-bearing
## clusters rises with size up to 21-100 markers and falls only above 100, where
## median min-p is also weaker, so the drop is confounded. A weak result on data
## that barely exhibits the problem is weak evidence against the proposal.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/test_then_cluster.R
## Env: SIM_DATA, CELL, TAG, ENVS, FILES, ALPHA
## =============================================================================
suppressMessages({library(data.table); library(LDscnR)})
SIM   <- Sys.getenv("SIM_DATA", "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5")
CELLS <- strsplit(Sys.getenv("CELLS", "V0.5_c1"), ",")[[1]]
TAG   <- Sys.getenv("TAG", "nobgs")
ENVS  <- as.integer(strsplit(Sys.getenv("ENVS", "1,2,3,4,5"), ",")[[1]])
FILES <- as.integer(strsplit(Sys.getenv("FILES", "1,2,3,4,5,6,7,8,9,10"), ",")[[1]])
ALPHA <- as.numeric(Sys.getenv("ALPHA", "0.05"))
simes <- function(p) { p <- sort(p[is.finite(p)]); n <- length(p)
  if (!n) return(NA_real_); min(n * p / seq_len(n)) }

one_env <- function(CELL, ENV) {
  MM <- list(); LK <- list(); CN <- list()
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
    s2 <- ld_group_map(g, prefix = i)[, .(marker, CL2 = group_id)]
    mm <- merge(merge(m, s1, by = "marker", all.x = TRUE), s2, by = "marker", all.x = TRUE)
    mm <- mm[!is.na(CL1) & !is.na(CL2)][, CL1 := paste0(i, "_", CL1)]
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
    MM[[length(MM)+1]] <- mm[, .(marker, CL1, CL2, p = emx_p, chr_type,
                                 qtn = true_pos_QTN %in% TRUE)]
    E <- pr$eMLG
    if (!is.null(E) && ncol(E)) {
      pe <- tryCatch(as.numeric(emmax_fast(emmax_setup(E, x$GRM), as.numeric(x$env$env))),
                     error = function(e) NULL)
      if (!is.null(pe)) CN[[length(CN)+1]] <-
        data.table(CL2 = paste0(i, "_", colnames(E)), p = pe)
    }
  }
  D <- rbindlist(MM); lk <- rbindlist(LK, fill = TRUE); nq <- uniqueN(lk$qtn)
  if (!nq) return(NULL)
  map12 <- unique(D[, .(CL1, CL2)])

  ## region-level scoring, identical for every arm
  score <- function(regions, n_test, lab) {
    sub  <- lk[CL2 %in% regions]
    best <- if (nrow(sub)) sub[order(-r2)][, .SD[1], by = qtn] else sub
    kept <- setdiff(regions, setdiff(sub$CL2, best$CL2))
    tp   <- if (nrow(best)) uniqueN(best$CL2) else 0L
    ntr  <- unique(D[CL2 %in% regions, .(CL2, chr_type)])[chr_type == "ntrl", .N]
    fn   <- unique(D[, .(CL2, chr_type)])[, mean(chr_type == "ntrl")]
    data.table(arm = lab, cell = CELL, tag = TAG, env = ENV, n_test = n_test, regions = length(kept), tp = tp,
               precision = tp / max(length(kept), 1),
               recall = (if (nrow(best)) uniqueN(best$qtn) else 0L) / nq,
               on_ntrl = ntr,
               FDP_ntrl = (ntr / fn) / max(length(regions), 1), n_qtn = nq)
  }
  bh <- function(p) { ok <- is.finite(p); q <- rep(NA_real_, length(p))
                      q[ok] <- p.adjust(p[ok], "BH"); which(!is.na(q) & q < ALPHA) }

  U2 <- D[, .(p = simes(p)), by = CL2]
  a  <- score(U2$CL2[bh(U2$p)], nrow(U2), "cluster-then-test")
  U1 <- D[, .(p = simes(p)), by = CL1]
  sig1 <- U1$CL1[bh(U1$p)]
  b  <- score(unique(map12[CL1 %in% sig1]$CL2), nrow(U1), "test-then-cluster")
  sigS <- D$marker[bh(D$p)]
  c_ <- score(unique(D[marker %in% sigS]$CL2), nrow(D), "single SNP")
  res <- list(a, b, c_)
  if (length(CN)) { CC <- rbindlist(CN)
    res[[length(res)+1]] <- score(CC$CL2[bh(CC$p)], nrow(CC), "consensus (no n penalty)") }
  rbindlist(res)
}
G <- CJ(cell = CELLS, env = ENVS, sorted = FALSE)
R <- rbindlist(Filter(Negate(is.null),
       lapply(seq_len(nrow(G)), function(k) one_env(G$cell[k], G$env[k]))))
OUT <- Sys.getenv("OUT", "module_sim_LDscnR/results/test_then_cluster.csv")
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
source("module_sim_LDscnR/prov.R")
fwrite(R, OUT)
write_prov(OUT, list(SIM_DATA = SIM, CELLS = CELLS, TAG = TAG, ENVS = ENVS,
                    FILES = FILES, ALPHA = ALPHA))
cat(sprintf("%s %s | %d envs\n\n", paste(CELLS, collapse=","), TAG, uniqueN(R$env)))
print(R[, .(envs = .N, n_test = round(mean(n_test)), regions = round(mean(regions), 1),
            tp = round(mean(tp), 1), precision = round(mean(precision), 3),
            recall = round(mean(recall), 3), FDP_ntrl = round(mean(FDP_ntrl), 3)), by = .(cell, arm)][order(cell, arm)])
W <- dcast(R, cell + env ~ arm, value.var = c("precision", "recall", "regions"))
cat("\npaired, test-then-cluster vs cluster-then-test:\n")
for (mm in c("precision", "recall", "regions")) {
  a <- W[[paste0(mm, "_cluster-then-test")]]; b <- W[[paste0(mm, "_test-then-cluster")]]
  cat(sprintf("  %-10s %d win / %d loss / %d tie\n", mm, sum(b > a), sum(b < a), sum(b == a)))
}
