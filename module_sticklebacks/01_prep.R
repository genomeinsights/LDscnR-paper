## module_sticklebacks/01_prep.R
## Build the aligned per-SNP input for the candidate-first outlier analysis of
## the 3sp (three-spine stickleback) MAF>=0.1 data, and extract the compact
## LD-decay summary used to derive per-chromosome rho thresholds.
## Outputs (regenerable; git-ignored):
##   module_sticklebacks/snp_stats_aligned.rds  (marker, Chr, Pos, maf, ld_w, emx_F/p, lfmm_F/p)
##   module_sticklebacks/decay_sum_3sp.rds       (per-chr a/a_pred/b/c decay params)
## Run from LDscnR-paper/:  Rscript module_sticklebacks/01_prep.R

suppressMessages({ library(data.table) })
mod  <- "/Users/petrikem/gitlab/LDscnR-paper/module_sticklebacks"
base <- "/Users/petrikem/gitlab/LD-scaling-genome-scans/empirical_data/3sp"

## aligned source: both methods' per-SNP stats, MAF>=0.1, keyed by marker="Chr:Pos"
sr <- as.data.table(readRDS(file.path(base, "SNP_res_3sp.rds")))    # 790,578 x 11
ldw_mat <- readRDS("/Users/petrikem/gitlab/LDscnR-paper/3sp_data/ld_ws_3sp_MAF01.rds")
stopifnot(all(sr$marker %in% rownames(ldw_mat)))
sr[, ld_w := ldw_mat[marker, "0.95"]]              # rho=0.95 local LD support
## raw F -> p, df = (1, N-2) = (1, 115); nulls self-calibrate, so no genomic control here
sr[, emx_p  := pf(emx_F,  1, 115, lower.tail = FALSE)]
sr[, lfmm_p := pf(lfmm_F, 1, 115, lower.tail = FALSE)]
sr <- sr[is.finite(ld_w) & is.finite(emx_p) & is.finite(lfmm_p)]
saveRDS(sr[, .(marker, Chr, Pos, maf, ld_w, emx_F, emx_p, lfmm_F, lfmm_p)],
        file.path(mod, "snp_stats_aligned.rds"))
cat(sprintf("snp_stats_aligned.rds: %d SNPs, %d chromosomes\n", nrow(sr), uniqueN(sr$Chr)))

## compact decay summary (avoid loading the 318 MB decay object downstream)
ds <- as.data.table(readRDS("/Users/petrikem/gitlab/LDscnR-paper/3sp_data/ld_decay_3sp_MAF01.rds")$decay_sum)
saveRDS(ds, file.path(mod, "decay_sum_3sp.rds"))
cat(sprintf("decay_sum_3sp.rds: %d chromosomes (a_pred, b, c)\n", nrow(ds)))
