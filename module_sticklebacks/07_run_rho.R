## module_sticklebacks/07_run_rho.R
## Re-run the three analyses with distance-restricted clustering, thresholds
## derived per chromosome from the LD-decay fit: rho_ld=0.95 (r2 linkage) and
## rho_d=0.99 (distance cap). Compare to the un-capped r0.5 baseline.
## Run from LDscnR-paper/:  Rscript module_sticklebacks/07_run_rho.R

suppressMessages({ library(LDscnR); library(data.table) })
set.seed(1)
mod <- "/Users/petrikem/gitlab/LDscnR-paper/module_sticklebacks"
B <- 1000L; CORES <- 4L; GRID <- seq(0, 0.99, by = 0.005)
RHO_LD <- 0.95; RHO_D <- 0.99

sr <- readRDS(file.path(mod, "snp_stats_aligned.rds")); setDT(sr)
ds <- readRDS(file.path(mod, "decay_sum_3sp.rds"))
LDd <- list(decay_sum = ds)
e <- new.env(); load("/Users/petrikem/gitlab/LDscnR-paper/LFMM_3sp/data/3sp_data.RData", envir = e)
GTs_all <- e$GTs_3sp; colnames(GTs_all) <- e$map_3sp$marker
GTs <- GTs_all[, sr$marker]; rm(GTs_all); gc()
eco <- as.integer(e$pheno_3sp$ecotype == "Marine")
K <- readRDS("/Users/petrikem/gitlab/LDscnR-paper/3sp_data/grm_null.rds")$K
map <- sr[, .(marker, Chr, Pos)]
pv_emx <- setNames(sr$emx_p, sr$marker); pv_lfmm <- setNames(sr$lfmm_p, sr$marker)
ldw <- setNames(sr$ld_w, sr$marker)

run <- function(pv, null, Y = NULL, Kk = NULL)
  ld_outlier_clusters(pv, ldw, map, GTs, null = null, Y = Y, K = Kk, B = B, cores = CORES,
                      rho_ld = RHO_LD, rho_d = RHO_D, LD_decay = LDd,
                      rmsc_grid = GRID, verbose = FALSE)

R <- list(emx_perm = run(pv_emx, "permutation", Y = eco, Kk = K),
          emx_bg   = run(pv_emx, "background"),
          lfmm_bg  = run(pv_lfmm, "background"))
saveRDS(R, file.path(mod, "results_three_nulls_rho.rds"))

for (k in names(R)) {
  cl <- as.data.table(R[[k]]$clusters)
  cat(sprintf("[%s] clusters=%d significant=%d median_sig_size=%s\n", k, nrow(cl),
              sum(cl$significant), ifelse(any(cl$significant), round(median(cl[significant==TRUE,n])), "-")))
}
## focus regions: Eda (Chr4), Chr7, inversion (Chr1)
ep <- as.data.table(R$emx_perm$clusters); eb <- as.data.table(R$emx_bg$clusters); lb <- as.data.table(R$lfmm_bg$clusters)
for (ch in c("Chr4","Chr7","Chr1")) {
  x <- ep[Chr==ch][order(-n)][1:min(6,.N)]
  x[, `:=`(Mb=sprintf("%.2f-%.2f",start/1e6,end/1e6),
           sig_eb=eb$significant[match(cluster,eb$cluster)], sig_lb=lb$significant[match(cluster,lb$cluster)])]
  cat(sprintf("\n-- %s (rho: r2~0.1, dist~%dkb) --\n", ch, round(d_from_rho(ds[Chr==ch,a_pred],0.99)/1e3)))
  print(x[, .(cluster, Mb, n, n_sig, sig_ep=significant, sig_eb, sig_lb)])
}
cat("\nsaved results_three_nulls_rho.rds\n")
