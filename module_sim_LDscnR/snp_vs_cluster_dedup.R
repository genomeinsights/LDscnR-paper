## =============================================================================
## snp_vs_cluster_dedup.R -- single-SNP and cluster-level analyses on ONE footing.
##
## The two have been scored differently and are not comparable as reported. Here
## both are reduced to the same units (stage-2 clusters) and the same satellite
## treatment, so the ONLY difference is how significance is decided:
##
##   snp      BH over all SNPs; a cluster is flagged if any member is significant
##   snp_n2   the same, but only SNPs in clusters of >= 2 members. Isolated
##            markers are the most obvious false positives, and dropping them is
##            a structural form of ld_w filtering that keeps SNP-level testing.
##   snp_n5   the same at >= 5 members
##   rep      BH over the LD-central representative p-values
##   simes    BH over Simes-combined member p-values
##   emlg     BH over EMMAX on the consensus genotype
##   simes_n2 Simes restricted to clusters of >= 2 members
##   simes_ns Simes excluding LONG-AND-SPARSE clusters: span > 100 kb with fewer
##            than 5 markers per 100 kb. That class is 3.4% of units and holds
##            2.5% of QTN, while long-and-DENSE units are 2.3% of units and hold
##            70% of QTN. Span alone would be destructive (dropping everything
##            over 100 kb removes 29 of 40 QTN clusters); span conditioned on
##            density is nearly free. The criterion is phenotype-blind -- span and
##            marker count come from the map and the clustering.
##
## snp_n2 vs simes_n2 is the decisive contrast: SAME marker set, so the only
## difference is whether member p-values are aggregated. It separates what
## clustering gains by REMOVING isolated markers from what it gains by
## COMBINING evidence within a block.
##
## SATELLITE REMOVAL. A flagged cluster that tags a QTN (r2 >= r2min within dmax)
## but is not the BEST tagger of it is a satellite of a locus already found. Two
## accountings are reported:
##   raw     every flagged cluster counts; satellites are false positives
##   dedup   satellites are removed from the flagged set entirely, which is the
##           evaluate_ors convention -- one region per QTN
##
## Env: SIM_DATA, CELL, TAG, ENV, FILES, ALPHA
## =============================================================================
suppressMessages({library(data.table); library(LDscnR)})
SIM  <- Sys.getenv("SIM_DATA", "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5")
OUT  <- Sys.getenv("OUT", "module_sim_LDscnR/results/filter_then_test")
CELLS <- strsplit(Sys.getenv("CELLS", "V0.5_c1,V0.5_c2,V1_c1.5,V2_c1"), ",")[[1]]
TAGS  <- strsplit(Sys.getenv("TAGS", "nobgs,bgs"), ",")[[1]]
ENVS  <- as.integer(strsplit(Sys.getenv("ENVS", "1,2,3,4,5,6,7,8,9,10"), ",")[[1]])
FILES <- as.integer(strsplit(Sys.getenv("FILES", "1,2,3,4,5,6,7,8,9,10"), ",")[[1]])
ALPHA <- as.numeric(Sys.getenv("ALPHA", "0.05"))
## The cluster partition is built from GENOTYPES and never sees a p-value, so
## swapping the association method changes exactly one column. The eMLG arm is
## the exception -- it refits the association on consensus genotypes -- and is
## skipped unless the engine is emmax.
ENG   <- Sys.getenv("ENGINE", "emmax")
PCOL  <- if (ENG == "emmax") "emx_p" else "lfmm_p"
CORES <- as.integer(Sys.getenv("CORES", "1"))
suppressMessages(library(parallel))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
simes <- function(p) { p <- sort(p[is.finite(p)]); n <- length(p)
  if (!n) return(NA_real_); min(n * p / seq_len(n)) }

one_panel <- function(CELL, TAG, ENV) {
clu <- list(); snp <- list(); lnk <- list()
for (i in FILES) {
  f <- sprintf("%s/adapt_%s_chr%d_%s_env%d.rds", SIM, TAG, i, CELL, ENV)
  if (!file.exists(f)) next
  x <- readRDS(f); m <- flag_true_qtns(as.data.table(x$map))
  pr <- ld_prune_and_eMLG(GTs = x$GTs, stage1 = x$complexity_reduction$stage1,
          LD_decay = x$LD_decay, ld_w_col = "ld_w_095", ld_w_threshold = 0.025,
          score_threshold = 0.80, min_r2_rho = 0.5, distance_threshold = 1e5,
          compute_unflagged_eMLG = TRUE, cores = 1)
  g  <- as.data.table(pr$groups)
  ms <- rbindlist(lapply(seq_len(nrow(g)), function(k)
          data.table(marker = g$members[[k]], CL = paste0(i, "_", g$group_id[k]))))
  mm <- merge(m, ms, by = "marker", all.x = TRUE)[!is.na(CL)]
  th  <- score_thresholds(as.data.table(x$LD_decay$decay_sum),
                          rho_r2 = 0.75, rho_d = 0.95, dmax_cap = 1e5)
  drv <- mm[true_pos_QTN %in% TRUE]
  ## per (cluster, QTN) best r2, for satellite identification
  if (nrow(drv)) lnk[[length(lnk)+1]] <- rbindlist(lapply(seq_len(nrow(drv)), function(j) {
    ch <- as.character(drv$Chr[j])
    near <- mm[as.character(Chr) == ch & abs(Pos - drv$Pos[j]) < th$dmax]
    if (!nrow(near)) return(NULL)
    r2 <- suppressWarnings(cor(x$GTs[, drv$marker[j]], x$GTs[, near$marker],
                               use = "pairwise.complete.obs")^2)
    d <- data.table(CL = near$CL, r2 = as.numeric(r2))[is.finite(r2) & r2 >= th$r2min]
    if (!nrow(d)) return(NULL)
    d[, .(r2 = max(r2)), by = CL][, qtn := paste0(i, "_", drv$marker[j])][]
  }))
  ## cluster-level summaries
  su <- mm[, .(p_rep = get(PCOL)[marker %in% pr$pruned][1], p_simes = simes(get(PCOL)),
               has_qtn = any(true_pos_QTN %in% TRUE)), by = CL]
  E <- if (ENG == "emmax") pr$eMLG else NULL
  if (!is.null(E) && ncol(E)) {
    pp <- tryCatch(emmax(Y = as.numeric(x$env$env), X = E, K = x$GRM, cores = 1)$pval,
                   error = function(e) NULL)
    if (!is.null(pp) && length(pp) == ncol(E)) {
      em <- data.table(CL = paste0(i, "_", colnames(E)), p_emlg = as.numeric(pp))
      su <- merge(su, em, by = "CL", all.x = TRUE)
    }
  }
  if (!"p_emlg" %in% names(su)) su[, p_emlg := NA_real_]
  geo <- mm[, .(n_loci = .N, span = max(Pos) - min(Pos)), by = CL]
  su <- merge(su, geo, by = "CL", all.x = TRUE)
  su[, dens := n_loci / pmax(span/1e5, 1e-9)]
  su[n_loci == 1, dens := NA_real_]
  su[, long_sparse := span > 1e5 & (is.na(dens) | dens < 5)]
  clu[[length(clu)+1]] <- su
  snp[[length(snp)+1]] <- mm[, .(CL, p = get(PCOL), n_loci = .N), by = CL][, .(CL, p, n_loci)]
}
su  <- rbindlist(clu, fill = TRUE)
sn  <- rbindlist(snp, fill = TRUE)
lk  <- rbindlist(lnk, fill = TRUE)
nq  <- sum(su$has_qtn)
cat(sprintf("  %d clusters, %d contain a QTN, %d tag one\n", nrow(su), nq, uniqueN(lk$CL)))

## flagged sets
flag <- list()
q_snp <- p.adjust(sn$p, "BH"); flag$snp <- unique(sn$CL[which(q_snp < ALPHA)])
for (mn in c(2, 5)) {
  ok <- sn$n_loci >= mn
  q <- rep(NA_real_, nrow(sn)); q[ok] <- p.adjust(sn$p[ok], "BH")
  flag[[paste0("snp_n", mn)]] <- unique(sn$CL[which(!is.na(q) & q < ALPHA)])
}
## INVERSE ARM: keep ONLY long-and-sparse, the class that holds 1 QTN of 40.
## This is a falsification test. The claim drawn from the other filters is that
## enrichment for CONTAINING a QTN does not predict detection performance. If
## that is right, restricting to a QTN-poor class should not be catastrophic.
## If it is catastrophic, enrichment does predict detection and the claim fails.
okls <- su$long_sparse & is.finite(su$p_simes)
qls <- rep(NA_real_, nrow(su)); qls[okls] <- p.adjust(su$p_simes[okls], "BH")
flag$simes_only_sparse <- su$CL[which(!is.na(qls) & qls < ALPHA)]
## and the dense counterpart, for symmetry
okld <- (su$span > 1e5) & !su$long_sparse & is.finite(su$p_simes)
qld <- rep(NA_real_, nrow(su)); qld[okld] <- p.adjust(su$p_simes[okld], "BH")
flag$simes_only_longdense <- su$CL[which(!is.na(qld) & qld < ALPHA)]
oknp <- !su$long_sparse & is.finite(su$p_simes)
qnp <- rep(NA_real_, nrow(su)); qnp[oknp] <- p.adjust(su$p_simes[oknp], "BH")
flag$simes_nosparse <- su$CL[which(!is.na(qnp) & qnp < ALPHA)]
oknb <- !su$long_sparse & su$n_loci >= 2 & is.finite(su$p_simes)
qnb <- rep(NA_real_, nrow(su)); qnb[oknb] <- p.adjust(su$p_simes[oknb], "BH")
flag$simes_n2_nosparse <- su$CL[which(!is.na(qnb) & qnb < ALPHA)]
ok2 <- su$n_loci >= 2
qs2 <- rep(NA_real_, nrow(su)); qs2[ok2] <- p.adjust(su$p_simes[ok2], "BH")
flag$simes_n2 <- su$CL[which(!is.na(qs2) & qs2 < ALPHA)]
for (cl in c("p_rep","p_simes","p_emlg")) {
  v <- su[[cl]]; ok <- is.finite(v)
  q <- rep(NA_real_, length(v)); q[ok] <- p.adjust(v[ok], "BH")
  flag[[sub("^p_", "", cl)]] <- su$CL[which(!is.na(q) & q < ALPHA)]
}

score <- function(fl) {
  contain <- sum(su$has_qtn[su$CL %in% fl])
  sub <- lk[CL %in% fl]
  ## dedup: per QTN keep only the best-tagging flagged cluster
  best <- if (nrow(sub)) sub[order(-r2)][, .SD[1], by = qtn] else sub
  sats <- setdiff(sub$CL, best$CL)
  kept <- setdiff(fl, sats)
  ## A region can be the best tagger of MORE THAN ONE QTN. Precision must count
  ## REGIONS that are a best tagger; recall must count QTN recovered. An earlier
  ## version used the QTN count for both and produced precision above 1.
  tp_reg <- if (nrow(best)) uniqueN(best$CL)  else 0L
  tp_qtn <- if (nrow(best)) uniqueN(best$qtn) else 0L
  data.table(flagged = length(fl), contain = contain,
             tag_not_contain = uniqueN(setdiff(sub$CL, su[has_qtn == TRUE]$CL)),
             neither = length(setdiff(fl, unique(c(sub$CL, su[has_qtn == TRUE]$CL)))),
             satellites = length(sats),
             raw_prec = contain/length(fl), raw_rec = contain/nq,
             dedup_regions = length(kept), dedup_tp = tp_reg, dedup_qtn = tp_qtn,
             dedup_prec = tp_reg/max(length(kept),1), dedup_rec = tp_qtn/nq)
}
res <- rbindlist(lapply(names(flag), function(n) cbind(analysis = n, score(flag[[n]]))))
res[, `:=`(cell = CELL, tag = TAG, env = ENV, n_qtn = nq, n_clusters = nrow(su))]
return(res)
}
grid <- CJ(cell = CELLS, tag = TAGS, env = ENVS, sorted = FALSE)
cat(sprintf("  %d panels, CORES=%d\n", nrow(grid), CORES))
runp <- function(z) tryCatch(one_panel(grid$cell[z], grid$tag[z], grid$env[z]),
                             error = function(e) { message("  FAIL ", grid$cell[z], " ", grid$tag[z],
                             " env", grid$env[z], ": ", conditionMessage(e)); NULL })
out <- if (CORES > 1) {
  mclapply(seq_len(nrow(grid)), runp, mc.cores = CORES, mc.preschedule = FALSE)
} else {
  lapply(seq_len(nrow(grid)), runp)
}
allr <- rbindlist(Filter(Negate(is.null), out), fill = TRUE)
stopifnot(nrow(allr) > 0)
fwrite(allr, file.path(OUT, sprintf("snp_vs_cluster_dedup_allpanels_%s.csv", ENG)))
cat(sprintf("\n  written: %d rows from %d panels\n", nrow(allr), uniqueN(allr[, .(cell,tag,env)])))
