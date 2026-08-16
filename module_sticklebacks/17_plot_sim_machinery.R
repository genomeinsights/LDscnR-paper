## module_sticklebacks/17_plot_sim_machinery.R
## Regenerate the two-panel Manhattan from the saved full run (16) without re-running the
## null: rebuild the edge cache from the saved C vectors, recompute survivors, plot.
## Top: EMMAX at its null-calibrated tau_C (incl Eda). Bottom: LFMM at the
## genomic-control-mapped tau_C (EMMAX's operating point on LFMM's C-quantile scale).
## Run: Rscript module_sticklebacks/17_plot_sim_machinery.R

suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
source("module_sim/R/_config.R")
mod <- "/Users/petrikem/gitlab/LDscnR-paper/module_sticklebacks"
R2LINK <- 0.1; DC <- 5e5; LMIN <- 10L
x <- readRDS(file.path(mod, "sim_machinery_3sp.rds"))
C_emx <- x$C_emx; C_lfmm <- x$C_lfmm; tc <- x$tau_c; tc_lfmm <- x$tau_lfmm
sr <- readRDS(file.path(mod, "snp_stats_aligned.rds")); setDT(sr)
decs <- as.data.table(readRDS(file.path(mod, "decay_sum_3sp.rds")))
e <- new.env(); load("/Users/petrikem/gitlab/LDscnR-paper/LFMM_3sp/data/3sp_data.RData", envir = e)
GTs <- e$GTs_3sp; colnames(GTs) <- e$map_3sp$marker; GTs <- GTs[, sr$marker]
uni <- unique(c(names(C_emx)[C_emx > 0], names(C_lfmm)[C_lfmm > 0]))
edges <- build_edge_cache(uni, sr[, .(marker, Chr, Pos)], GTs, decs, r2_link = R2LINK, dcap = DC); rm(GTs); gc()
surv <- function(Cv, tau){ mk <- names(Cv)[Cv >= tau]; reg <- cluster_from_cache(mk, edges); unlist(reg[lengths(reg) >= LMIN]) }
s_emx <- surv(C_emx, tc); s_lf <- surv(C_lfmm, tc_lfmm)

chr_lev <- paste0("Chr", c(1:18,20,21)); sr[, Chr := factor(Chr, levels = chr_lev)]; setorder(sr, Chr, Pos)
clen <- sr[, .(mx = max(Pos)), by = Chr]; clen[, off := cumsum(shift(mx, fill = 0))]; xoff <- setNames(clen$off, as.character(clen$Chr))
sr[, gx := Pos + xoff[as.character(Chr)]]; ax <- clen[, .(center = off + mx/2, Chr)]
panel <- function(Cv, survmk, tau, ttl) {
  d <- data.table(marker = names(Cv), C = as.numeric(Cv))[sr[, .(marker, gx, Chr)], on = "marker"][C > 0]
  d[, is_surv := marker %in% survmk]
  ggplot() + geom_point(data = d[is_surv == FALSE], aes(gx, C), color = "grey72", size = 0.4) +
    geom_point(data = d[is_surv == TRUE], aes(gx, C, color = Chr), size = 1.2) +
    geom_hline(yintercept = tau, linetype = 2, color = "#D62828") +
    scale_x_continuous(breaks = ax$center, labels = sub("Chr", "", ax$Chr), expand = c(0.01, 0)) +
    scale_color_viridis_d(guide = "none") + labs(title = ttl, x = NULL, y = "C-score") +
    theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()) }
p <- panel(C_emx, s_emx, tc, sprintf("EMMAX (low-ld_w GRM), l_min=10, null-calibrated tau_C=%.2f -> %d loci (incl Eda/Chr4)", tc, nrow(x$emx))) /
     panel(C_lfmm, s_lf, tc_lfmm, sprintf("LFMM @ genomic-control-mapped tau_C=%.2f, l_min=10 -> %d loci across 14 chr (naive tau_C=%.2f -> %d)", tc_lfmm, nrow(x$lfmm_gc), tc, nrow(x$lfmm_naive)))
ggsave(file.path(mod, "fig_sim_machinery_3sp.png"), p, width = 13, height = 7, dpi = 150)
cat(sprintf("wrote figure: EMMAX %d loci, LFMM-GC %d loci\n", nrow(x$emx), nrow(x$lfmm_gc)))
