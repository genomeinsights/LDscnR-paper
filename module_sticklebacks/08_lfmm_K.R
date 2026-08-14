## module_sticklebacks/08_lfmm_K.R
## Definitive test of the Eda / LFMM story: rerun lfmm2 across K (number of
## latent factors) with the ecotype as env, and watch whether Eda's signal is
## SPECIFICALLY suppressed as K grows (latent-factor absorption) while a non-Eda
## strong locus (Chr20) stays stable. Same genotypes/env as the paper's LFMM.
## Run from LDscnR-paper/:  Rscript module_sticklebacks/08_lfmm_K.R

suppressMessages({ library(data.table); library(LEA) })
mod <- "/Users/petrikem/gitlab/LDscnR-paper/module_sticklebacks"
scratch <- file.path(tempdir(), "lfmmK")
dir.create(scratch, showWarnings = FALSE, recursive = TRUE); setwd(scratch)

sr <- readRDS(file.path(mod, "snp_stats_aligned.rds")); setDT(sr); sr[, idx := .I]
e <- new.env(); load("/Users/petrikem/gitlab/LDscnR-paper/LFMM_3sp/data/3sp_data.RData", envir = e)
G <- e$GTs_3sp; colnames(G) <- e$map_3sp$marker; G <- G[, sr$marker]      # 117 x 790k, aligned
eco <- as.integer(e$pheno_3sp$ecotype == "Marine")
write.lfmm(G, "geno.lfmm"); write.env(eco, "eco.env")

## regions to track
eda_idx <- sr[Chr == "Chr4" & Pos >= 12.80e6 & Pos <= 12.82e6, idx]     # Eda core
c20_idx <- sr[Chr == "Chr20" & Pos >= 0.17e6 & Pos <= 0.27e6, idx]      # strong non-Eda control
cat(sprintf("Eda core SNPs: %d | Chr20 control SNPs: %d\n\n", length(eda_idx), length(c20_idx)))

res <- rbindlist(lapply(c(1, 2, 4, 7, 10), function(K) {
  proj <- lfmm2("geno.lfmm", "eco.env", K = K)
  pv <- suppressWarnings(lfmm2.test(proj, "geno.lfmm", "eco.env",
                                    genomic.control = TRUE, full = TRUE))$pvalues
  if (is.matrix(pv)) pv <- pv[, 1]
  cat(sprintf("K=%2d done: Eda minP=%.2g (nsig %d) | Chr20 minP=%.2g (nsig %d)\n",
              K, min(pv[eda_idx], na.rm = TRUE), sum(pv[eda_idx] < 0.01, na.rm = TRUE),
              min(pv[c20_idx], na.rm = TRUE), sum(pv[c20_idx] < 0.01, na.rm = TRUE)))
  data.table(K = K,
             eda_minP  = min(pv[eda_idx], na.rm = TRUE),
             eda_nsig  = sum(pv[eda_idx] < 0.01, na.rm = TRUE),
             c20_minP  = min(pv[c20_idx], na.rm = TRUE),
             c20_nsig  = sum(pv[c20_idx] < 0.01, na.rm = TRUE),
             gif_note  = "gc=TRUE")
}))
setwd(mod)
cat("Eda vs Chr20 signal across K (lfmm2, ecotype env, genomic control):\n")
res[, `:=`(eda_minP = signif(eda_minP, 2), c20_minP = signif(c20_minP, 2))]
print(res)
cat("\nInterpretation: if eda_minP/nsig weakens sharply as K rises while Chr20 holds,\n",
    "the latent factors specifically absorb Eda.\n")
saveRDS(res, file.path(mod, "lfmm_K_scan.rds"))
