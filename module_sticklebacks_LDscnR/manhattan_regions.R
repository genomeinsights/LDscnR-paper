## =====================================================================
## module_sticklebacks_LDscnR / manhattan_regions.R
##
## Genome-wide -log10(q) Manhattans for 3sp, per method (EMMAX + LFMM), with the
## C-score outlier regions coloured -- in the original "sim-machinery" plotting
## style (single-row chromosome facets, grey background, no x-axis, bold strips,
## the col_vec palette). Uses the migrated pipeline end to end: the bundle's
## ld_w + GRM, the structure-aware null (EMMAX), ld_cscore for LFMM, and
## ld_edges/ld_regions for clustering.
##
## Parameterised on the three clustering knobs (defaults reproduce the reference
## figure): tau_C (C-score threshold), l_min (min cluster size), rho_ld (decay-
## relative r^2 link; rho_ld=0.60 ~ r2=0.4 on 3sp). Pass "auto" for tau_C / l_min
## to take the null-calibrated values instead.
##
## Run from the LDscnR-paper root:
##   Rscript module_sticklebacks_LDscnR/manhattan_regions.R [tau_C] [l_min] [rho_ld]
##   Rscript module_sticklebacks_LDscnR/manhattan_regions.R 0.05 10 0.60   # reference look
##   Rscript module_sticklebacks_LDscnR/manhattan_regions.R auto auto 0.75 # null-calibrated
## =====================================================================

suppressMessages({ library(data.table); library(LDscnR); library(ggplot2); library(patchwork) })

## ---- 0. config -------------------------------------------------------
BUNDLE     <- "module_sticklebacks_LDscnR/data/3sp_LDscnR_data.rds"
NULL_CACHE <- "module_sticklebacks_LDscnR/results/null_uncapped_3sp.rds"   # reused if present
OUTFIG     <- "module_sticklebacks_LDscnR/figures"
a <- commandArgs(trailingOnly = TRUE)
PAR <- list(
  tau_C  = if (length(a) >= 1) a[1] else "0.05",     # numeric, or "auto" (calibrate_tauc)
  l_min  = if (length(a) >= 2) a[2] else "10",       # integer, or "auto" (calibrate_lmin)
  rho_ld = if (length(a) >= 3) as.numeric(a[3]) else 0.60,
  alpha  = 0.05, qstar = seq(0, 0.95, by = 0.05), dcap = 5e5, B = 100L, seed = 1L,
  fdr = 0.05, lmin_q = 0.99, lmin_tau = 0.05)
if (!dir.exists(OUTFIG)) dir.create(OUTFIG, recursive = TRUE)

## reference palette + style (from module_sticklebacks/15_sim_vs_poster_manhattan.R)
COL_VEC <- rep(c("#B2DF8A","#FFD92F","firebrick","#33A02C","#7FC97F","#CAB2D6","#FB8072","#E6AB02","steelblue",
  "#FB9A99","#1B9E77","#BC80BD","#E31A1C","#7570B3","#A6761D","#A6CEE3","salmon","forestgreen","#BF5B17"), 30)

## ---- 1. data + structure-aware null ----------------------------------
d   <- readRDS(BUNDLE); map <- as.data.table(d$map)
if (file.exists(NULL_CACHE)) { null <- readRDS(NULL_CACHE); cat("[1] loaded cached null\n") } else {
  cat("[1] computing structured null (B=", PAR$B, ") ...\n", sep = ""); utils::flush.console()
  prep <- emmax_setup(d$GTs, d$GRM)
  null <- structured_null(d$eco, d$GTs, d$GRM, d$ld_ws, basis = "genetic", B = PAR$B,
                          alpha = PAR$alpha, qstar = PAR$qstar, prep = prep, seed = PAR$seed)
  saveRDS(null, NULL_CACHE)
}
C_emx  <- null$C_obs
C_lfmm <- ld_cscore(map$lfmm_p, d$ld_ws, alpha = PAR$alpha, qstar = PAR$qstar); names(C_lfmm) <- map$marker
uni    <- unique(c(null$universe, names(C_lfmm)[C_lfmm > 0]))
edges  <- ld_edges(uni, d$GTs, map[, .(marker, Chr, Pos)], as.data.table(d$LD_decay$decay_sum),
                   rho_ld = PAR$rho_ld, dcap = PAR$dcap)

## ---- 2. resolve tau_C / l_min (fixed, or null-calibrated) -------------
l_min <- if (identical(PAR$l_min, "auto")) calibrate_lmin(null, edges, tau = PAR$lmin_tau, q = PAR$lmin_q) else as.integer(PAR$l_min)
tau_C <- if (identical(PAR$tau_C, "auto")) calibrate_tauc(null, edges, l_min = l_min, fdr = PAR$fdr,
                                                          tau_grid = seq(0.05, 1, 0.05)) else as.numeric(PAR$tau_C)
r2 <- stats::median(ld_from_rho(as.data.table(d$LD_decay$decay_sum)$b, as.data.table(d$LD_decay$decay_sum)$c, PAR$rho_ld))
cat(sprintf("[2] tau_C=%.3f  l_min=%d  rho_ld=%.2f (~ r2=%.2f)\n", tau_C, l_min, PAR$rho_ld, r2))
regs_of <- function(C) { r <- ld_regions(names(C)[C >= tau_C], edges); r[lengths(r) >= l_min] }

## ---- 3. plot one method in the reference style -----------------------
chr_lev <- paste0("Chr", sort(as.integer(gsub("Chr", "", unique(map$Chr)))))
plot_one <- function(regs, qcol, tag) {
  mh <- copy(map[, .(marker, Chr, Pos, q = get(qcol))]); mh[, CL_id := "ns"]
  for (i in seq_along(regs)) mh[marker %in% regs[[i]], CL_id := as.character(i)]
  mh[, Chr := factor(Chr, levels = chr_lev)]
  ggplot(mh, aes(Pos, -log10(q))) +
    geom_point(data = mh[CL_id == "ns"], color = "grey35", size = 0.5, alpha = 0.3) +
    geom_point(data = mh[CL_id != "ns"], aes(color = CL_id), size = 1.8) +
    geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.6) +
    facet_wrap(~ Chr, nrow = 1, scales = "free_x") +
    scale_color_manual(values = rep(COL_VEC, 10)) +
    labs(x = "Genomic position", y = bquote(-log[10](italic(q)) ~ "(" * .(tag) * ")"),
         title = sprintf("3sp %s : C>=%.2f, >=%d-SNP clusters (rho_ld=%.2f ~ r2=%.2f) -> %d regions",
                         tag, tau_C, l_min, PAR$rho_ld, r2, length(regs))) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none", panel.grid = element_blank(),
          strip.text = element_text(face = "bold", size = 7),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.spacing.x = unit(0.05, "lines"))
}
regE <- regs_of(C_emx); regL <- regs_of(C_lfmm)
cat(sprintf("[3] EMMAX regions: %d ; LFMM regions: %d\n", length(regE), length(regL)))
gE <- plot_one(regE, "emx_q",  "EMMAX")
gL <- plot_one(regL, "lfmm_q", "LFMM")
tag <- sprintf("tau%.2f_lmin%d_rho%.2f", tau_C, l_min, PAR$rho_ld)
ggsave(file.path(OUTFIG, paste0("manhattan_EMMAX_", tag, ".png")), gE, width = 18, height = 4.5, dpi = 180)
ggsave(file.path(OUTFIG, paste0("manhattan_LFMM_",  tag, ".png")), gL, width = 18, height = 4.5, dpi = 180)
ggsave(file.path(OUTFIG, paste0("manhattan_both_",  tag, ".png")), gE / gL, width = 18, height = 9, dpi = 180)
cat(sprintf("[3] wrote manhattan_{EMMAX,LFMM,both}_%s.png to %s\n", tag, OUTFIG))
