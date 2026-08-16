## module_sticklebacks/15_sim_vs_poster_manhattan.R
## Can the SIM machinery (per-SNP C) reproduce the POSTER Manhattan? Everything held
## identical to 14 (poster's p-values emx_p_GC/lfmm_P, C>0.05 selection, final r2=0.4
## re-cluster >=10, -log10(q) facet-by-Chr) EXCEPT the C-score:
##   poster (14): folded C = frac of cells where SNP is in a >=l_min cluster
##   sim    (15): per-SNP C = frac of (rho,q*) cells where SNP is candidate & BH-FDR<0.05
## If the two Manhattans give similar clusters, the sim machinery achieves poster-like
## results with the per-SNP C (+ the final >=10 re-cluster as the post-filter).
## Run: Rscript module_sticklebacks/15_sim_vs_poster_manhattan.R

suppressMessages({ library(data.table); library(igraph); library(ggplot2) })
source("module_sim/R/_config.R")                       # cscore_count
mod <- "/Users/petrikem/gitlab/LDscnR-paper/module_sticklebacks"; CORES <- 4L
QSTAR <- seq(0, 0.95, by = 0.05); ALPHA <- 0.05

precompute_LD_edges <- function(GTs, map, r2_min = 0.1, max_bp = Inf, cores = 1) {
  mk <- intersect(colnames(GTs), map$marker); ms <- copy(map[marker %in% mk]); setkey(ms, Chr, Pos)
  rbindlist(parallel::mclapply(unique(ms$Chr), function(ch) { cm <- ms[Chr == ch]; cmk <- cm$marker
    if (length(cmk) < 2) return(NULL); g <- as.matrix(GTs[, cmk, drop = FALSE]); storage.mode(g) <- "double"
    R <- cor(g, use = "pairwise.complete.obs")^2; R[is.na(R)] <- 0; diag(R) <- 0
    idx <- which(R >= r2_min, arr.ind = TRUE); idx <- idx[idx[, 1] < idx[, 2], , drop = FALSE]; if (!nrow(idx)) return(NULL)
    pos <- cm$Pos; d <- data.table(marker1 = cmk[idx[, 1]], marker2 = cmk[idx[, 2]], r2 = R[idx], dist_bp = abs(pos[idx[, 1]] - pos[idx[, 2]]))
    if (is.finite(max_bp)) d <- d[dist_bp <= max_bp]; d }, mc.cores = cores), fill = TRUE) }
LD_igraph_components <- function(el, markers, r2_th = 0.4, bp_th = Inf) {
  markers <- unique(markers); if (length(markers) < 2) return(data.table(marker = markers, CL_id = seq_along(markers), n_loci = 1L))
  ed <- el[marker1 %in% markers & marker2 %in% markers & r2 >= r2_th]; if (is.finite(bp_th)) ed <- ed[dist_bp <= bp_th]
  if (!nrow(ed)) return(data.table(marker = markers, CL_id = seq_along(markers), n_loci = 1L))
  g <- graph_from_data_frame(ed[, .(marker1, marker2)], directed = FALSE, vertices = data.table(name = markers))
  comp <- components(g); data.table(marker = names(comp$membership), CL_id = as.integer(comp$membership), n_loci = comp$csize[as.integer(comp$membership)]) }

sr <- as.data.table(readRDS("/Users/petrikem/gitlab/LD-scaling-genome-scans/empirical_data/3sp/SNP_res_3sp.rds"))
map <- sr[, .(Chr, Pos, marker, emx_p_GC, lfmm_P, lfmm_q, emx_q)]
e <- new.env(); load("/Users/petrikem/gitlab/LDscnR-paper/LFMM_3sp/data/3sp_data.RData", envir = e)
GTs <- e$GTs_3sp; colnames(GTs) <- e$map_3sp$marker
LDW <- readRDS("/Users/petrikem/gitlab/LDscnR-paper/3sp_data/ld_ws_3sp_MAF01.rds"); LDW[is.na(LDW)] <- 0
map <- map[marker %in% rownames(LDW)]; LDW <- LDW[map$marker, ]; RHO <- colnames(LDW)

col_vec <- rep(c("#B2DF8A","#FFD92F","firebrick","#33A02C","#7FC97F","#CAB2D6","#FB8072","#E6AB02","steelblue",
  "#FB9A99","#1B9E77","#BC80BD","#E31A1C","#7570B3","#A6761D","#A6CEE3","salmon","forestgreen","#BF5B17"), 30)
run <- function(pcol, qcol, tag) {
  C <- cscore_count(map[[pcol]], LDW, RHO, QSTAR, ALPHA)      # SIM per-SNP C on the poster's p-values
  outl <- names(C)[C > 0.05]
  cat(sprintf("[%s] per-SNP C>0.05: %d markers\n", tag, length(outl)))
  if (!length(outl)) return(invisible())
  el <- precompute_LD_edges(GTs[, outl, drop = FALSE], map[marker %in% outl], r2_min = 0.1, max_bp = 1e6, cores = CORES)
  cls <- LD_igraph_components(el, outl, r2_th = 0.4, bp_th = Inf)
  mh <- cls[n_loci >= 10, .(marker, CL_id)][map, on = "marker"]; mh[, CL_id := as.character(CL_id)][is.na(CL_id), CL_id := "ns"]
  cat(sprintf("[%s] SIM per-SNP C -> %d >=10-SNP clusters across %d chr\n", tag, mh[CL_id!="ns", uniqueN(CL_id)], mh[CL_id!="ns", uniqueN(Chr)]))
  p <- ggplot(mh, aes(Pos, -log10(get(qcol)))) +
    geom_point(data = mh[CL_id == "ns"], color = "grey35", size = 0.5, alpha = 0.3) +
    geom_point(data = mh[CL_id != "ns"], aes(color = CL_id), size = 1.8) +
    geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.6) +
    facet_wrap(~ factor(Chr, levels = paste0("Chr", c(1:18,20,21))), nrow = 1, scales = "free_x") +
    scale_color_manual(values = rep(col_vec, 10)) +
    labs(x = "Genomic position", y = bquote(-log[10](italic(q))~"("*.(tag)*")"),
         title = sprintf("SIM machinery (per-SNP C, alpha=0.05) Manhattan: %s, C>0.05, >=10-SNP clusters (r2=0.4)", tag)) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none", panel.grid = element_blank(), strip.text = element_text(face = "bold", size = 7),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.spacing.x = unit(0.05, "lines"))
  ggsave(file.path(mod, sprintf("fig_sim_manhattan_%s.png", tag)), p, width = 18, height = 4.5, dpi = 200)
  cat("wrote fig_sim_manhattan_", tag, ".png\n", sep = "") }
run("emx_p_GC", "emx_q",  "EMX")
run("lfmm_P",   "lfmm_q", "LFMM")
