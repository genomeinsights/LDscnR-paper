## How many stage-2 clusters does the ld_w filter ADD that the conventional
## single-SNP scan did not already find, and how many does it LOSE?
## Filtering relaxes the BH threshold, so it can gain clusters; it also removes
## markers from eligibility, so it can lose them. Both are counted here, and the
## gained/lost sets are checked for whether they contain a driving QTN.
suppressMessages({library(data.table); library(LDscnR)})
SIM  <- Sys.getenv("SIM_DATA", "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5")
CELL <- Sys.getenv("CELL", "V0.5_c1"); TAG <- Sys.getenv("TAG", "nobgs")
ENV  <- as.integer(Sys.getenv("ENV", "3")); FILES <- 1:10
KS   <- as.integer(strsplit(Sys.getenv("KS", "1000,5000,20000,50000"), ",")[[1]])
ALPHA <- 0.05

per_file <- function(i) {
  f <- sprintf("%s/adapt_%s_chr%d_%s_env%d.rds", SIM, TAG, i, CELL, ENV)
  if (!file.exists(f)) return(NULL)
  x <- readRDS(f); m <- as.data.table(x$map)
  pr <- ld_prune_and_eMLG(GTs = x$GTs, stage1 = x$complexity_reduction$stage1,
          LD_decay = x$LD_decay, ld_w_col = "ld_w_095", ld_w_threshold = 0.025,
          score_threshold = 0.80, min_r2_rho = 0.5, distance_threshold = 1e5,
          compute_unflagged_eMLG = FALSE, cores = 1)
  stopifnot(identical(sort(pr$pruned), sort(x$grm_markers)))
  g <- as.data.table(pr$groups)
  ms <- rbindlist(lapply(seq_len(nrow(g)), function(k)
          data.table(marker = g$members[[k]], CL = paste0(i, "_", g$group_id[k]))))
  m <- merge(m, ms, by = "marker", all.x = TRUE)
  m[, driving := true_QTN %in% TRUE & MAF > 0.1 & p_Va > 0.05]
  m[, .(marker, p = emx_p, ld_w = ld_w_095, CL, driving)]
}
m <- rbindlist(lapply(FILES, per_file), fill = TRUE)[is.finite(p) & !is.na(CL)]
qtn_cl <- unique(m[driving == TRUE]$CL)

flagged <- function(idx) {
  q <- p.adjust(m$p[idx], "BH")
  unique(m$CL[idx][which(q < ALPHA)])
}
## TWO ranking units. ld_w is a LOCAL-LD statistic and so autocorrelated along
## the genome: ranking SNPs by it returns a few large blocks many times over
## (the top 1,000 SNPs here are 13 clusters), which spends the budget on
## redundancy. Ranking CLUSTERS by their median ld_w gives one value per block,
## which is the level at which ld_w carries independent information.
cl <- m[, .(ld_w_med = median(ld_w, na.rm = TRUE)), by = CL]
sel_snp <- function(k) head(order(-m$ld_w), k)
sel_clu <- function(k) which(m$CL %in% head(cl[order(-ld_w_med)]$CL, k))
conv <- flagged(seq_len(nrow(m)))
cat(sprintf("  conventional scan: %d clusters flagged, %d contain a driving QTN\n",
            length(conv), sum(conv %in% qtn_cl)))
cat(sprintf("  (%d driving QTN in %d distinct clusters)\n\n", sum(m$driving), length(qtn_cl)))
audit <- function(unit) rbindlist(lapply(KS, function(kk) {
  idx <- if (unit == "cluster") sel_clu(kk) else sel_snp(kk)
  fl <- flagged(idx)
  gained <- setdiff(fl, conv); lost <- setdiff(conv, fl)
  data.table(unit = unit, k = kk, n_sel_snp = length(idx),
             n_sel_clu = uniqueN(m$CL[idx]), n_filtered = length(fl),
             gained = length(gained), gained_with_qtn = sum(gained %in% qtn_cl),
             lost = length(lost), lost_with_qtn = sum(lost %in% qtn_cl),
             kept = length(intersect(fl, conv)),
             qtn_filt = sum(fl %in% qtn_cl))
}))
res <- rbind(audit("cluster"), audit("snp"))
cat("  qtn_filt is out of", length(qtn_cl), "QTN-bearing clusters;",
    "the conventional scan finds", sum(conv %in% qtn_cl), "\n\n")
print(res[unit=="cluster"]); cat("\n"); print(res[unit=="snp"])
cat("\n  gained = flagged after filtering but NOT by the conventional scan\n")
cat("  lost   = flagged by the conventional scan but not after filtering\n")
