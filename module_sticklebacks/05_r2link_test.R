## module_sticklebacks/05_r2link_test.R
## Re-run the three analyses with tighter single-linkage (r2_link = 0.8) and
## compare cluster structure to r2_link = 0.5, focusing on the chained clusters
## (Chr4 Eda 12.8-19.9 Mb, Chr7 2-18 Mb) and whether tightening rescues Eda for
## LFMM. Candidates (=RMSC q*) are unchanged; only clustering + null differ.
## Run from LDscnR-paper/:  Rscript module_sticklebacks/05_r2link_test.R

suppressMessages({ library(LDscnR); library(data.table) })
set.seed(1)
mod <- "/Users/petrikem/gitlab/LDscnR-paper/module_sticklebacks"
B <- 1000L; CORES <- 4L; GRID <- seq(0, 0.99, by = 0.005); R2LINK <- 0.8

sr <- readRDS(file.path(mod, "snp_stats_aligned.rds")); setDT(sr)
e <- new.env(); load("/Users/petrikem/gitlab/LDscnR-paper/LFMM_3sp/data/3sp_data.RData", envir = e)
GTs_all <- e$GTs_3sp; colnames(GTs_all) <- e$map_3sp$marker
GTs <- GTs_all[, sr$marker]; rm(GTs_all); gc()
eco <- as.integer(e$pheno_3sp$ecotype == "Marine")
K <- readRDS("/Users/petrikem/gitlab/LDscnR-paper/3sp_data/grm_null.rds")$K
map <- sr[, .(marker, Chr, Pos)]
pv_emx <- setNames(sr$emx_p, sr$marker); pv_lfmm <- setNames(sr$lfmm_p, sr$marker)
ldw <- setNames(sr$ld_w, sr$marker)

run <- function(pv, null, Y = NULL, Kk = NULL)
  ld_outlier_clusters(pv, ldw, map, GTs, null = null, Y = Y, K = Kk, B = B,
                      cores = CORES, r2_link = R2LINK, rmsc_grid = GRID, verbose = FALSE)

R8 <- list(emx_perm = run(pv_emx, "permutation", Y = eco, Kk = K),
           emx_bg   = run(pv_emx, "background"),
           lfmm_bg  = run(pv_lfmm, "background"))
saveRDS(R8, file.path(mod, "results_three_nulls_r08.rds"))

R5 <- readRDS(file.path(mod, "results_three_nulls.rds"))
for (RR in list(R5, R8)) for (k in names(RR)) RR[[k]]$clusters <- as.data.table(RR[[k]]$clusters)

hd <- function(R, tag) for (k in names(R)) {
  cl <- as.data.table(R[[k]]$clusters)
  cat(sprintf("  [%s %-8s] clusters=%d significant=%d\n", tag, k, nrow(cl), sum(cl$significant)))
}
cat("=== cluster/significant counts ===\n"); hd(R5, "r0.5"); hd(R8, "r0.8")

## focus: Chr4 (Eda ~12.8 Mb) and Chr7 clusters, both linkages
focus <- function(R, chr, tag) {
  cat(sprintf("\n-- %s  %s clusters (n>=3), sig flags [emx_perm/emx_bg/lfmm_bg] --\n", tag, chr))
  cl <- as.data.table(R$emx_perm$clusters)[Chr == chr & n >= 3]
  ep <- as.data.table(R$emx_perm$clusters); eb <- as.data.table(R$emx_bg$clusters); lb <- as.data.table(R$lfmm_bg$clusters)
  cl[, `:=`(Mb = sprintf("%.2f-%.2f", start/1e6, end/1e6),
            sig_ep = significant,
            sig_eb = eb$significant[match(cluster, eb$cluster)],
            sig_lb = lb$significant[match(cluster, lb$cluster)])]
  print(cl[order(-n)][1:min(8,.N), .(cluster, Mb, n, n_sig, sig_ep, sig_eb, sig_lb)])
}
focus(R5, "Chr4", "r0.5"); focus(R8, "Chr4", "r0.8")
focus(R5, "Chr7", "r0.5"); focus(R8, "Chr7", "r0.8")
cat("\nsaved results_three_nulls_r08.rds\n")
