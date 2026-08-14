## module_sticklebacks/06_subcluster.R
## Prototype: post-process the loose (single-linkage r0.5) candidate clusters by
##   (1) splitting on PHYSICAL DISTANCE  (consecutive gap > distance_threshold)
##   (2) optionally a 1-r^2 hclust cut WITHIN each physically-contiguous run
## Goal: split the chained clusters (Eda 12.8->19.9 Mb, Chr7 2->18 Mb) into tight
## units WITHOUT shattering a coherent block (Chr1 inversion 21.49-21.94) into
## singletons the way global r2_link=0.8 did.
## Run from LDscnR-paper/:  Rscript module_sticklebacks/06_subcluster.R

suppressMessages({ library(data.table); library(ggplot2) })
mod <- "/Users/petrikem/gitlab/LDscnR-paper/module_sticklebacks"
R5 <- readRDS(file.path(mod, "results_three_nulls.rds"))
cand <- as.data.table(R5$emx_perm$candidates)   # marker, Chr, Pos, ld_w, pval, q_snp, sig, cluster
setkey(cand, cluster)

## candidate genotypes only (small)
e <- new.env(); load("/Users/petrikem/gitlab/LDscnR-paper/LFMM_3sp/data/3sp_data.RData", envir = e)
G <- e$GTs_3sp; colnames(G) <- e$map_3sp$marker; G <- G[, cand$marker]; rm(e); gc()

## ---- sub-clustering of one cluster ---------------------------------------
sub_one <- function(mk, pos, distance_threshold, r2_sub, linkage = "complete") {
  n <- length(mk); if (n == 1L) return(rep(1L, 1L))
  ord <- order(pos)
  run <- integer(n); run[ord] <- cumsum(c(TRUE, diff(pos[ord]) > distance_threshold))
  out <- integer(n); nxt <- 0L
  for (r in unique(run)) {
    idx <- which(run == r)
    if (length(idx) == 1L || is.null(r2_sub)) { out[idx] <- nxt + 1L; nxt <- nxt + 1L; next }
    R2 <- suppressWarnings(stats::cor(G[, mk[idx], drop = FALSE], use = "pairwise.complete.obs")^2)
    R2[!is.finite(R2)] <- 0
    cl <- stats::cutree(stats::hclust(stats::as.dist(1 - R2), linkage), h = 1 - r2_sub)
    out[idx] <- nxt + cl; nxt <- nxt + max(cl)
  }
  out
}
refine <- function(distance_threshold, r2_sub, linkage = "complete") {
  cand[, sub := sub_one(marker, Pos, distance_threshold, r2_sub, linkage), by = cluster]
  cand[, rcl := paste0(cluster, ".", sub)]
  cand[, .(n = .N, n_sig = sum(sig), lo = min(Pos), hi = max(Pos), span_kb = (max(Pos)-min(Pos))/1e3),
       by = .(Chr, rcl)]
}

show_region <- function(cl, chr, lab) {
  x <- cl[Chr == chr][order(-n)][1:min(6,.N)]
  x[, Mb := sprintf("%.2f-%.2f", lo/1e6, hi/1e6)]
  cat(sprintf("  [%s] %s: %d clusters (n>=3: %d)\n", lab, chr, nrow(cl[Chr==chr]), nrow(cl[Chr==chr & n>=3])))
  print(x[, .(rcl = sub("Chr[0-9]+_", "", rcl), Mb, n, n_sig, span_kb = round(span_kb))])
}

cat("baseline r0.5 (no sub-clustering):\n")
base <- cand[, .(n=.N, n_sig=sum(sig), lo=min(Pos), hi=max(Pos)), by=.(Chr,cluster)][, rcl := cluster]
for (ch in c("Chr4","Chr7","Chr1")) { x<-base[Chr==ch][order(-n)][1:min(4,.N)]
  x[,Mb:=sprintf("%.2f-%.2f",lo/1e6,hi/1e6)]; cat(sprintf("  %s:\n",ch)); print(x[,.(cluster,Mb,n,n_sig)]) }

for (cfg in list(list(d=5e5, r=NULL, l="-",        tag="dist 500kb only"),
                 list(d=2e5, r=NULL, l="-",        tag="dist 200kb only"),
                 list(d=5e5, r=0.5,  l="complete", tag="dist 500kb + r2>=0.5 complete"),
                 list(d=5e5, r=0.5,  l="average",  tag="dist 500kb + r2>=0.5 average"))) {
  cl <- refine(cfg$d, cfg$r, cfg$l)
  cat(sprintf("\n===== %s : total clusters=%d (median sig-cluster size=%s) =====\n",
              cfg$tag, nrow(cl),
              ifelse(any(cl$n_sig>0), round(median(cl[n_sig>0,n])), "-")))
  show_region(cl, "Chr4", cfg$tag); show_region(cl, "Chr7", cfg$tag); show_region(cl, "Chr1", cfg$tag)
}
