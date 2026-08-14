## module_sim/04_structure_null.R
## Structure-aware null experiment, run on the 3-chromosome stickleback DEMO as a
## sanity check. Idea: draw each cluster's null from the genome-wide HIGH-ld_w
## (candidate) pool instead of the low-ld_w tail, so ordinary structure inflation
## is the expectation and only exceptional peaks survive.
## RESULT: it OVER-CORRECTS. On the demo the candidate pool is ~70% significant
## (there candidates ARE the signal), so a cluster must beat 70% to survive and
## Eda cannot -- only the inversion is recovered. Same failure mode as the
## rotation null. Conclusion: no single null handles both the localized-signal
## (c1/demo) and structure-dominated (c2) regimes; keep the low-ld_w background
## null as THE method and use ld_outlier_clusters()'s structure diagnostics
## (ldw_tracks_structure, whole_chr_clusters) to flag when calls are unreliable.
## Needs the demo cache built in the LDscnR package repo (dev/outlier/cache.rds)
## and the stickleback genotypes shipped with the package.
## Run from LDscnR-paper/:  Rscript module_sim/04_structure_null.R

suppressMessages({ library(LDscnR); library(data.table); library(parallel) })
CACHE <- Sys.getenv("LDSCNR_DEMO_CACHE",
                    unset = "/Users/petrikem/gitlab/LDscnR/dev/outlier/cache.rds")
o <- readRDS(CACHE); map <- copy(o$map); data(stickleback); GTs <- stickleback$GTs
setnames(map, c("ld_w_095", "p_loco"), c("ld_w", "p")); setorder(map, Chr, Pos)
DCLUST <- 2e5; B <- 1000                       # tight 200kb distance cap
pos_of <- setNames(map$Pos, map$marker); chr_of <- setNames(map$Chr, map$marker)

## single-linkage on 1 - r^2 within chromosome, split at DCLUST gaps
sl <- function(mk) {
  if (length(mk) <= 1) return(if (length(mk)) list(mk) else list())
  out <- list()
  for (ch in unique(chr_of[mk])) {
    m <- mk[chr_of[mk] == ch]; if (length(m) == 1) { out <- c(out, list(m)); next }
    pos <- pos_of[m]; R <- suppressWarnings(cor(GTs[, m])^2); R[!is.finite(R)] <- 0
    cl <- cutree(hclust(as.dist(1 - R), "single"), h = 0.5); ordr <- order(pos)
    run <- integer(length(m)); run[ordr] <- cumsum(c(TRUE, diff(pos[ordr]) > DCLUST))
    out <- c(out, split(m, paste(cl, run)))
  }
  out
}

## structure-aware null: draw n from the CANDIDATE pool (consecutive => captures LD)
struct_p <- function(clusters, cand_mk, sig_of) {
  cm <- map[marker %in% cand_mk][order(Chr, Pos), marker]
  cand_sig <- sig_of[cm]; nb <- length(cand_sig)
  cat(sprintf("  candidate-pool sig-rate=%.3f (n=%d)\n", mean(cand_sig), nb))
  unlist(mclapply(clusters, function(mk) {
    n <- length(mk); k <- sum(sig_of[mk])
    nk <- vapply(seq_len(B), function(b) {
      st <- sample.int(nb, 1); sum(cand_sig[((st - 1 + seq_len(n) - 1) %% nb) + 1]) }, numeric(1))
    (1 + sum(nk >= k)) / (B + 1)
  }, mc.cores = 4))
}

q <- 0.87; cand <- map$ld_w >= quantile(map$ld_w, q)
qsnp <- rep(NA_real_, nrow(map)); qsnp[cand] <- p.adjust(map$p[cand], "fdr")
pcut <- max(map$p[!is.na(qsnp) & qsnp < 0.05]); sig_of <- setNames(map$p <= pcut & cand, map$marker)
cl <- sl(map$marker[cand])
emp_q <- p.adjust(struct_p(cl, map$marker[cand], sig_of), "fdr"); keep <- cl[emp_q < 0.05]
tab <- if (length(keep)) rbindlist(lapply(keep, function(mk)
  data.table(Chr = chr_of[mk[1]], n = length(mk),
             Mb = sprintf("%.2f-%.2f", min(pos_of[mk]) / 1e6, max(pos_of[mk]) / 1e6)))) else data.table()
cat(sprintf("=== DEMO structure-aware null + tight clustering (200kb): %d clusters -> %d sig ===\n",
    length(cl), length(keep)))
if (nrow(tab)) print(tab[order(Chr, -n)]) else cat("  NONE\n")
cat(">> only the inversion survives; Eda is lost => the structure-aware null over-corrects.\n")
