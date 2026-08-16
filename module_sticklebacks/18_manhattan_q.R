## module_sticklebacks/18_manhattan_q.R
## Poster-style Manhattan for the full sim-machinery run (16): y = -log10(q), faceted by
## chromosome, the surviving outlier REGIONS coloured one distinct colour each (poster's
## discriminating palette, not viridis). Top = EMMAX at its null-calibrated tau_C (incl
## Eda); LFMM at the genomic-control-mapped tau_C. q = BH-FDR of the same engine's p
## (EMMAX = low-ld_w-GRM fast EMMAX, recomputed; LFMM = lfmm_p).
## Run: Rscript module_sticklebacks/18_manhattan_q.R

suppressMessages({ library(LDscnR); library(data.table); library(ggplot2) })
source("module_sim/R/_config.R")
mod <- "/Users/petrikem/gitlab/LDscnR-paper/module_sticklebacks"
R2LINK <- 0.1; DC <- 5e5; LMIN <- 10L
gcta_grm <- function(X){ p<-colMeans(X)/2; k<-p>0&p<1; X<-X[,k,drop=F]; p<-p[k]
  Z<-sweep(sweep(X,2,2*p,"-"),2,sqrt(2*p*(1-p)),"/"); tcrossprod(Z)/ncol(Z) }
col_vec <- c("#E31A1C","#1F78B4","#33A02C","#FF7F00","#6A3D9A","#B15928","#A6CEE3","#FB9A99","#FDBF6F",
             "#CAB2D6","#B2DF8A","#FFFF99","#8DD3C7","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69",
             "#FCCDE5","#BC80BD","#CCEBC5","#FFED6F","#E7298A","#66A61E","#E6AB02","#A6761D","#666666")

x <- readRDS(file.path(mod, "sim_machinery_3sp.rds"))
C_emx <- x$C_emx; C_lfmm <- x$C_lfmm; tc <- x$tau_c; tc_lfmm <- x$tau_lfmm
sr <- readRDS(file.path(mod, "snp_stats_aligned.rds")); setDT(sr)
LDW <- readRDS("/Users/petrikem/gitlab/LDscnR-paper/3sp_data/ld_ws_3sp_MAF01.rds")[sr$marker, ]
decs <- as.data.table(readRDS(file.path(mod, "decay_sum_3sp.rds")))
e <- new.env(); load("/Users/petrikem/gitlab/LDscnR-paper/LFMM_3sp/data/3sp_data.RData", envir = e)
GTs <- e$GTs_3sp; colnames(GTs) <- e$map_3sp$marker; GTs <- GTs[, sr$marker]
eco <- as.integer(e$pheno_3sp$ecotype == "Marine")

## engine q-values: EMMAX = low-ld_w GRM fast EMMAX (same as the run); LFMM = lfmm_p
K <- gcta_grm(GTs[, which(LDW[, "0.95"] < 0.05 & is.finite(LDW[, "0.95"]))])
p_emx <- fast_emmax_p(fast_emmax_setup(GTs, K), eco)
sr[, q_emx  := p.adjust(p_emx, "fdr")]
sr[, q_lfmm := p.adjust(lfmm_p, "fdr")]

uni <- unique(c(names(C_emx)[C_emx > 0], names(C_lfmm)[C_lfmm > 0]))
edges <- build_edge_cache(uni, sr[, .(marker, Chr, Pos)], GTs, decs, r2_link = R2LINK, dcap = DC); rm(GTs); gc()
region_map <- function(Cv, tau) { mk <- names(Cv)[Cv >= tau]; reg <- cluster_from_cache(mk, edges); reg <- reg[lengths(reg) >= LMIN]
  if (!length(reg)) return(data.table(marker = character(), region = integer()))
  rbindlist(lapply(seq_along(reg), function(i) data.table(marker = reg[[i]], region = i))) }
r_emx <- region_map(C_emx, tc); r_lf <- region_map(C_lfmm, tc_lfmm)

chr_lev <- paste0("Chr", c(1:18,20,21)); sr[, Chr := factor(Chr, levels = chr_lev)]
manh <- function(rmap, qcol, ttl, fn) {
  d <- copy(sr[, .(marker, Chr, Pos, q = get(qcol))]); d <- rmap[d, on = "marker"]
  d[, region := ifelse(is.na(region), "ns", as.character(region))]
  d[, y := -log10(pmax(q, 1e-30))]
  p <- ggplot() +
    geom_point(data = d[region == "ns"], aes(Pos, y), color = "grey72", size = 0.35, alpha = 0.4) +
    geom_point(data = d[region != "ns"], aes(Pos, y, color = region), size = 1.6) +
    geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.5) +
    facet_wrap(~Chr, nrow = 1, scales = "free_x") +
    scale_color_manual(values = rep(col_vec, 6)) +
    labs(x = "Genomic position", y = expression(-log[10](italic(q))), title = ttl) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none", panel.grid = element_blank(), strip.text = element_text(face = "bold", size = 7),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.spacing.x = unit(0.04, "lines"))
  ggsave(file.path(mod, fn), p, width = 18, height = 4.2, dpi = 200); cat("wrote", fn, "\n") }
manh(r_emx, "q_emx",  sprintf("EMMAX (low-ld_w GRM), null-calibrated tau_C=%.2f, l_min=10 -> %d outlier regions (incl Eda/Chr4)", tc, uniqueN(r_emx$region)),
     "fig_manhattan_q_EMX.png")
manh(r_lf,  "q_lfmm", sprintf("LFMM @ genomic-control-mapped tau_C=%.2f, l_min=10 -> %d outlier regions across 14 chr", tc_lfmm, uniqueN(r_lf$region)),
     "fig_manhattan_q_LFMM.png")
