## =====================================================================
## module_sticklebacks_LDscnR / kingman_cscore_enrichment.R
##
## SNP-level positive-control: is a SNP with a higher C-score (the count of
## rho x q* grid cells in which it was a valid local-LD outlier candidate, as a
## fraction) more likely to sit inside a Kingman EcoPeak? Bin SNPs by C-score and
## compute the fold enrichment of global-specific EcoPeak membership per bin,
## relative to the genome-wide (all-SNP) rate, for EMMAX and LFMM.
##
## Reads kingman2021's lifted peak BED READ-ONLY; writes only to this module.
## Run from the LDscnR-paper root:
##   Rscript module_sticklebacks_LDscnR/kingman_cscore_enrichment.R
## Writes figures/kingman_cscore_enrichment.png + results/kingman_cscore_enrichment.csv
## =====================================================================
suppressMessages({ library(data.table); library(LDscnR); library(ggplot2); library(patchwork) })

BND  <- "module_sticklebacks_LDscnR/data/3sp_LDscnR_data.rds"
PEAK <- "kingman2021/data/liftover/pv_c155.specific.bed"       # global-specific EcoPeaks (gasAcu1)
OUT  <- "module_sticklebacks_LDscnR"; ROMAN <- c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI",
        "XII","XIII","XIV","XV","XVI","XVII","XVIII","XIX","XX","XXI")
BREAKS <- c(-1e-9, 0, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0)
LABS   <- c("0", "(0,.05]", "(.05,.1]", "(.1,.2]", "(.2,.3]", "(.3,.5]", "(.5,1]")

d  <- readRDS(BND); m3 <- as.data.table(d$map); m3[, chr_num := as.integer(gsub("Chr","",Chr))]
C_emx  <- ld_cscore(emmax_fast(emmax_setup(d$GTs, d$GRM), d$eco), d$ld_ws, alpha = 0.05, qstar = seq(0,.95,.05))
C_lfmm <- ld_cscore(m3$lfmm_p, d$ld_ws, alpha = 0.05, qstar = seq(0,.95,.05))

## per-SNP membership in a global-specific EcoPeak (gasAcu1)
Pk <- fread(PEAK, header = FALSE, col.names = c("chr","start","end","pv"))
Pk[, chr_num := match(sub("^chr","",chr), ROMAN)]; Pk <- Pk[!is.na(chr_num), .(chr_num,start,end)]
setkey(Pk, chr_num, start, end)
m3[, in_peak := !is.na(foverlaps(m3[, .(chr_num, start = Pos, end = Pos)], Pk,
                                 by.x = c("chr_num","start","end"), type = "any", mult = "first", nomatch = NA)$start)]

enrich <- function(C, method) {
  dt <- data.table(C = as.numeric(C), in_peak = m3$in_peak)
  base <- mean(dt$in_peak)
  dt[, Cbin := cut(C, breaks = BREAKS, labels = LABS, right = TRUE)]
  a <- dt[, .(n = .N, n_peak = sum(in_peak), prop = mean(in_peak)), by = Cbin][order(Cbin)]
  a[, `:=`(fold = prop / base, method = method)]; a
}
res <- rbind(enrich(C_emx, "EMMAX"), enrich(C_lfmm, "LFMM"))
base_rate <- mean(m3$in_peak)
sp_e <- cor(as.numeric(C_emx),  as.integer(m3$in_peak), method = "spearman")
sp_l <- cor(as.numeric(C_lfmm), as.integer(m3$in_peak), method = "spearman")
cat(sprintf("genome-wide SNP-in-EcoPeak rate = %.4f\n", base_rate))
cat(sprintf("Spearman(C, in_peak): EMMAX=%.3f  LFMM=%.3f\n", sp_e, sp_l))
cat("\n=== fold enrichment of EcoPeak membership by C-score bin ===\n")
print(res[, .(method, Cbin, n, n_peak, prop = round(prop,4), fold = round(fold,1))])
fwrite(res, file.path(OUT, "results", "kingman_cscore_enrichment.csv"))

## ---- plot: fold vs C-score bin, per method -----------------------------------------
res[, Cbin := factor(Cbin, levels = LABS)]
g <- ggplot(res, aes(Cbin, fold, fill = method)) +
  geom_col(position = position_dodge(0.8), width = 0.72, color = "grey30", linewidth = 0.25) +
  geom_hline(yintercept = 1, linetype = 2, color = "grey40") +
  geom_text(aes(label = n_peak), position = position_dodge(0.8), vjust = -0.4, size = 2.6) +
  scale_fill_manual(values = c(EMMAX = "#E1AF00", LFMM = "#3B9AB2"), name = NULL) +
  labs(x = "per-SNP C-score bin", y = "fold enrichment of EcoPeak membership\n(vs genome-wide rate)",
       title = "Higher C-score -> higher chance of overlapping a Kingman global-specific EcoPeak",
       subtitle = sprintf("genome-wide rate %.3f%%;  Spearman(C, in-peak): EMMAX %.3f, LFMM %.3f;  labels = # SNPs in a peak",
                          100*base_rate, sp_e, sp_l)) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.x = element_blank(), legend.position = c(0.12, 0.88),
        plot.title = element_text(face = "bold", size = 11))
ggsave(file.path(OUT, "figures", "kingman_cscore_enrichment.png"), g, width = 10, height = 5.5, dpi = 180)
cat("wrote figures/kingman_cscore_enrichment.png\n")
