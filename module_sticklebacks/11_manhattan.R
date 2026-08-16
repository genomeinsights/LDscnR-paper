## module_sticklebacks/11_manhattan.R
## C-score Manhattan for the empirical 3sp data (EMMAX vs LFMM), from the observed
## alpha=0.05 C landscape (09_cscore.R -> cscore_obs.rds). Shows WHERE consistency
## concentrates: Eda (Chr4) highlighted; EMMAX panel carries its structured-null
## tau_C (FDR<=0.05) reference line. Only C>0 SNPs are drawn (C=0 is the baseline).
## Run from LDscnR-paper/:  Rscript module_sticklebacks/11_manhattan.R
## Output: module_sticklebacks/fig_cscore_manhattan_3sp.png

suppressMessages({ library(data.table); library(ggplot2) })
mod <- "/Users/petrikem/gitlab/LDscnR-paper/module_sticklebacks"
o  <- readRDS(file.path(mod, "cscore_obs.rds")); sr <- as.data.table(o$sr)
sn <- tryCatch(readRDS(file.path(mod, "structured_null_3sp.rds")), error = function(e) NULL)
tau_emx <- if (!is.null(sn)) sn$tau05 else NA_real_

## chromosome order + cumulative x
chr_lev <- paste0("Chr", c(1:18, 20, 21))
sr <- sr[Chr %in% chr_lev]; sr[, Chr := factor(Chr, levels = chr_lev)]
setorder(sr, Chr, Pos)
chr_len <- sr[, .(mx = max(Pos)), by = Chr]; chr_len[, off := cumsum(shift(mx, fill = 0))]
sr <- merge(sr, chr_len[, .(Chr, off)], by = "Chr"); sr[, gx := Pos + off]
ax <- sr[, .(center = (min(gx) + max(gx)) / 2), by = Chr]

## long format, C>0 only; flag Eda/Chr4
L <- rbind(sr[C_emx  > 0, .(Chr, gx, C = C_emx,  method = "EMMAX")],
           sr[C_lfmm > 0, .(Chr, gx, C = C_lfmm, method = "LFMM")])
L[, method := factor(method, levels = c("EMMAX", "LFMM"))]
L[, shade := ifelse(Chr == "Chr4", "Eda (Chr4)", ifelse(as.integer(Chr) %% 2 == 0, "even", "odd"))]

p <- ggplot(L, aes(gx, C, color = shade)) +
  geom_point(size = 0.5) +
  facet_wrap(~method, ncol = 1, scales = "free_y") +
  scale_color_manual(values = c("Eda (Chr4)" = "#D62828", "even" = "grey60", "odd" = "grey30"), name = NULL) +
  scale_x_continuous(breaks = ax$center, labels = sub("Chr", "", ax$Chr), expand = c(0.01, 0)) +
  labs(title = "3sp C-score Manhattan (observed, alpha=0.05): EMMAX vs LFMM",
       subtitle = "C = frac of (rho,q*) cells candidate & FDR<0.05 | EMMAX = precomputed emx_p | Eda = Chr4 (red). No null line: EMMAX engine not yet reconciled with the null.",
       x = "Chromosome", y = "C-score") +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(), legend.position = "top")
ggsave(file.path(mod, "fig_cscore_manhattan_3sp.png"), p, width = 12, height = 6.5, dpi = 150)
cat(sprintf("wrote fig_cscore_manhattan_3sp.png | EMMAX C>0: %d SNPs (max %.3f), LFMM C>0: %d (max %.3f)\n",
            sr[C_emx > 0, .N], max(sr$C_emx), sr[C_lfmm > 0, .N], max(sr$C_lfmm)))
## top loci per method
cat("\nEMMAX top chromosomes by max C:\n"); print(sr[, .(maxC = round(max(C_emx),3), nSNP = sum(C_emx > 0)), by = Chr][order(-maxC)][1:6])
cat("\nLFMM top chromosomes by max C:\n");  print(sr[, .(maxC = round(max(C_lfmm),3), nSNP = sum(C_lfmm > 0)), by = Chr][order(-maxC)][1:6])
