## =============================================================================
## module_sim_LDscnR / stage1_vs_stage2_units.R
##
## PK: run the method on STAGE-1 clusters instead of stage-2. They are clustered
## into regions at the end anyway, so what does the choice of test unit do to the
## multiple testing?
##
## Stage 1 (ld_complexity_reduction) groups markers into LD clusters. Stage 2
## (ld_prune_and_eMLG) merges stage-1 clusters that are both linked and
## physically contiguous, so it is strictly coarser: fewer, larger units.
##
## Testing at stage 1 therefore raises the multiplicity but tests purer units;
## testing at stage 2 lowers it but tests units that may mix signal with
## neighbouring background. Scored identically at both levels -- Simes over
## members, BH across units, then the module's usual dedup convention (one region
## per QTN, best tagger, satellites removed) -- plus the convention-free
## neutral-chromosome false-positive rate.
##
## RESULT (V0.5_c1 nobgs env2, 20 detectable QTN):
##
##   units              n_test   discoveries   TP   precision   recall   FDP_ntrl
##   single SNP        298,741       564       13     0.023      0.65      0.11
##   stage-1 clusters  101,983        17        8     0.471      0.45      0.34
##   stage-2 clusters  101,569        15        8     0.533      0.45      0.38
##
## STAGE 1 DOES ESSENTIALLY ALL THE MULTIPLICITY REDUCTION, AND STAGE 2 DOES NONE.
## 298,741 markers -> 101,983 stage-1 clusters is a 2.93x cut; stage 2 removes a
## further 414 units, 0.41%. The multiple-testing benefit is stage 1's alone.
##
## WHAT STAGE 2 DOES INSTEAD IS RESHAPE THE FLAGGED SUBSET. It merges 1,407
## flagged stage-1 clusters into 993 groups -- a 1.42-fold merge -- while leaving
## the ~100,000 unflagged clusters untouched as their own units. So it changes
## which units exist where the signal is, not how many are tested, and it buys
## precision 0.471 -> 0.533 at identical recall and identical TP: two of the 17
## discoveries were separate pieces of one locus or a false unit absorbed into a
## true one.
##
## TWO METRICS DISAGREE ABOUT THE SINGLE-SNP SCAN AND BOTH ARE REPORTED. By
## tagging, its precision is 0.023 against 0.53. By the convention-free
## neutral-chromosome rate it looks BETTER (0.11 against 0.38). They measure
## different failures: only 45 of its 564 discoveries are on a chromosome with no
## QTN, so the other ~500 non-tagging hits sit on QTN chromosomes -- they are
## REDUNDANT hits around real signal rather than structure artefacts. That is the
## redundancy argument restated: 564 significant SNPs are not 564 loci, and a
## metric that only asks "is this near truth" cannot see it.
##
## Single panel, as requested. Given how many single-panel numbers have failed to
## replicate in this project, treat the magnitudes as provisional.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/stage1_vs_stage2_units.R
## Env: SIM_DATA, CELL, TAG, ENV, FILES, ALPHA
## =============================================================================
suppressMessages({library(data.table); library(LDscnR)})
SIM   <- Sys.getenv("SIM_DATA", "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5")
CELL  <- Sys.getenv("CELL", "V0.5_c1"); TAG <- Sys.getenv("TAG", "nobgs")
ENV   <- as.integer(Sys.getenv("ENV", "2"))
FILES <- as.integer(strsplit(Sys.getenv("FILES", "1,2,3,4,5,6,7,8,9,10"), ",")[[1]])
ALPHA <- as.numeric(Sys.getenv("ALPHA", "0.05"))
simes <- function(p) { p <- sort(p[is.finite(p)]); n <- length(p)
  if (!n) return(NA_real_); min(n * p / seq_len(n)) }

u1 <- list(); u2 <- list(); sn <- list(); lk1 <- list(); lk2 <- list(); lks <- list()
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
  mm <- mm[!is.na(CL1) & !is.na(CL2)]
  mm[, CL1 := paste0(i, "_", CL1)]
  th  <- score_thresholds(as.data.table(x$LD_decay$decay_sum),
                          rho_r2 = 0.75, rho_d = 0.95, dmax_cap = 1e5)
  drv <- mm[true_pos_QTN %in% TRUE]
  mklk <- function(col) if (!nrow(drv)) NULL else rbindlist(lapply(seq_len(nrow(drv)), function(j) {
    near <- mm[as.character(Chr) == as.character(drv$Chr[j]) & abs(Pos - drv$Pos[j]) < th$dmax]
    if (!nrow(near)) return(NULL)
    r2 <- suppressWarnings(cor(x$GTs[, drv$marker[j]], x$GTs[, near$marker],
                               use = "pairwise.complete.obs")^2)
    d <- data.table(CL = near[[col]], r2 = as.numeric(r2))[is.finite(r2) & r2 >= th$r2min]
    if (!nrow(d)) return(NULL)
    d[, .(r2 = max(r2)), by = CL][, qtn := paste0(i, "_", drv$marker[j])][] }))
  lk1[[length(lk1)+1]] <- mklk("CL1"); lk2[[length(lk2)+1]] <- mklk("CL2")
  ## the single-SNP arm needs its OWN link table keyed on marker -- passing it a
  ## cluster-keyed one silently scores every discovery as false.
  lks[[length(lks)+1]] <- mklk("marker")
  u1[[length(u1)+1]] <- mm[, .(p = simes(emx_p), n_loci = .N,
                               chr_type = chr_type[1]), by = .(CL = CL1)]
  u2[[length(u2)+1]] <- mm[, .(p = simes(emx_p), n_loci = .N,
                               chr_type = chr_type[1]), by = .(CL = CL2)]
  sn[[length(sn)+1]] <- mm[, .(CL = marker, p = emx_p, n_loci = 1L, chr_type)]
}
L1 <- rbindlist(lk1, fill = TRUE); L2 <- rbindlist(lk2, fill = TRUE)
LS <- rbindlist(lks, fill = TRUE)
nq <- uniqueN(c(L1$qtn, L2$qtn, LS$qtn))

score <- function(U, lk, lab) {
  ok  <- is.finite(U$p); q <- rep(NA_real_, nrow(U)); q[ok] <- p.adjust(U$p[ok], "BH")
  fl  <- U$CL[which(!is.na(q) & q < ALPHA)]
  sub <- lk[CL %in% fl]
  best <- if (nrow(sub)) sub[order(-r2)][, .SD[1], by = qtn] else sub
  kept <- setdiff(fl, setdiff(sub$CL, best$CL))
  tp   <- if (nrow(best)) uniqueN(best$CL) else 0L
  ntr  <- U[CL %in% fl & chr_type == "ntrl", .N]
  fntr <- U[ok & chr_type == "ntrl", .N] / sum(ok)
  data.table(units = lab, n_test = sum(ok), median_size = as.double(median(U$n_loci)),
             discoveries = length(kept), tp = tp,
             precision = round(tp / max(length(kept), 1), 3),
             recall = round((if (nrow(best)) uniqueN(best$qtn) else 0L) / nq, 3),
             on_neutral_chr = ntr,
             FDP_neutral = round((ntr / fntr) / max(length(fl), 1), 3))
}
U1 <- rbindlist(u1); U2 <- rbindlist(u2); SN <- rbindlist(sn)
R <- rbind(score(SN, LS, "single SNP"), score(U1, L1, "stage-1 clusters"),
           score(U2, L2, "stage-2 clusters"))
cat(sprintf("%s %s env%d | %d detectable QTN\n\n", CELL, TAG, ENV, nq))
print(R)
cat(sprintf("\nmultiplicity: stage 1 is %.1fx stage 2, and %.3fx the SNP count\n",
            R[units=="stage-1 clusters"]$n_test / R[units=="stage-2 clusters"]$n_test,
            R[units=="stage-1 clusters"]$n_test / R[units=="single SNP"]$n_test))
