## module_sticklebacks/02_run_three_nulls.R
## Candidate-first LD-aware outlier clustering on the 3sp (three-spine
## stickleback) MAF>=0.1 data, three ways:
##   (1) EMMAX  + permutation null   (fast EMMAX; structure-aware)
##   (2) EMMAX  + background null    (method-agnostic)
##   (3) LFMM   + background null    (LFMM is non-separable -> background only)
## Observed per-SNP statistics come from SNP_res_3sp.rds (the aligned source);
## ld_w is the rho=0.95 column of ld_ws_3sp_MAF01.rds.
## Run from LDscnR-paper/:  Rscript module_sticklebacks/02_run_three_nulls.R

suppressMessages({ library(LDscnR); library(data.table) })
set.seed(1)
mod <- "/Users/petrikem/gitlab/LDscnR-paper/module_sticklebacks"
B <- 1000L; CORES <- 4L; RHO <- "0.95"; R2LINK <- 0.5

## ---- aligned per-SNP stats (marker, Chr, Pos, maf, ld_w, emx_p, lfmm_p) ----
sr <- readRDS(file.path(mod, "snp_stats_aligned.rds"))
setDT(sr)
cat(sprintf("SNPs: %d | chromosomes: %d\n", nrow(sr), uniqueN(sr$Chr)))

## ---- genotypes (ML dosages, MAF>=0.1) + GRM + phenotype --------------------
e <- new.env()
load("/Users/petrikem/gitlab/LDscnR-paper/LFMM_3sp/data/3sp_data.RData", envir = e)
GTs_all <- e$GTs_3sp
colnames(GTs_all) <- e$map_3sp$marker
GTs <- GTs_all[, sr$marker]                 # 117 x 790k, aligned to sr
rm(GTs_all); gc()
eco <- as.integer(e$pheno_3sp$ecotype == "Marine")
K <- readRDS("/Users/petrikem/gitlab/LDscnR-paper/3sp_data/grm_null.rds")$K
stopifnot(nrow(K) == length(eco), ncol(GTs) == nrow(sr))
cat(sprintf("GTs: %d x %d | eco Marine/Fresh = %d/%d | K: %dx%d\n",
            nrow(GTs), ncol(GTs), sum(eco), sum(1 - eco), nrow(K), ncol(K)))

map <- sr[, .(marker, Chr, Pos)]
pv_emx  <- setNames(sr$emx_p,  sr$marker)
pv_lfmm <- setNames(sr$lfmm_p, sr$marker)
ldw     <- setNames(sr$ld_w,   sr$marker)
GRID <- seq(0, 0.99, by = 0.005)            # extend past 0.98 to locate the peak

run <- function(pv, null, Y = NULL, Kk = NULL, tag) {
  t0 <- proc.time()[3]
  r <- ld_outlier_clusters(pv, ldw, map, GTs, null = null, Y = Y, K = Kk,
                           B = B, cores = CORES, r2_link = R2LINK,
                           rmsc_grid = GRID, verbose = TRUE)
  cat(sprintf("\n### %s : %.1fs\n", tag, proc.time()[3] - t0))
  print(r)
  r
}

res_emx_perm <- run(pv_emx,  "permutation", Y = eco, Kk = K, tag = "EMMAX + permutation")
res_emx_bg   <- run(pv_emx,  "background",                    tag = "EMMAX + background")
res_lfmm_bg  <- run(pv_lfmm, "background",                    tag = "LFMM  + background")

saveRDS(list(emx_perm = res_emx_perm, emx_bg = res_emx_bg, lfmm_bg = res_lfmm_bg),
        file.path(mod, "results_three_nulls.rds"))

## ---- comparison: significant clusters per chromosome ----------------------
sig_by_chr <- function(r) r$clusters[significant == TRUE, .N, by = Chr][order(-N)]
cat("\n================ significant clusters per chromosome ================\n")
for (nm in c("emx_perm", "emx_bg", "lfmm_bg")) {
  r <- get(paste0("res_", nm))
  s <- sig_by_chr(r)
  cat(sprintf("\n[%s]  q*=%.3f  candidates=%d  clusters=%d  significant=%d\n",
              nm, r$ld_w_threshold, nrow(r$candidates), nrow(r$clusters),
              sum(r$clusters$significant)))
  print(s)
}
cat("\nsaved results_three_nulls.rds\n")
