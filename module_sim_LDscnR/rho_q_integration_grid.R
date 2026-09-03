## =====================================================================
## module_sim_LDscnR / rho_q_integration_grid.R
##
## The C-score's OWN integration grid -- the 400 analyses C averages over, as
## opposed to the tau x l_min OPERATING grid (null_operating_grid.R) which is
## what happens downstream once C exists.
##
## Each cell is one complete analysis: keep the markers whose ld_w at that rho
## sits above the q* quantile, BH-adjust WITHIN that candidate set, take
## alpha < 0.05, cluster, count regions. C is the fraction of cells in which a
## marker is a hit, so tau_C = 0.05 means "a hit in at least 5% of them" -- and
## this figure is what that fraction is a fraction OF.
##
## Why the grid is canonical, and the tau x l_min one is not: rho and q* are
## both naturally bounded on [0,1], so "fraction of the grid" is a well-defined
## quantity. l_min has no natural upper bound, so a C-squared computed over
## l_min 2-10 and one over 40-100 are not comparable. That asymmetry is the
## reason tau_C = 0.05 can be principled where a C-squared cut cannot.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/rho_q_integration_grid.R
## Env: SIM_INPUTS, NULL_DIR, OUT, ENV (default 1)
## =====================================================================
## The C-score's OWN grid: for each (rho, q*), take the markers whose ld_w at
## that rho sits above the q* quantile, BH-adjust WITHIN that candidate set, keep
## alpha < 0.05, cluster, and count regions. C is the fraction of these 400 cells
## in which a marker is a hit -- so this figure shows what C is averaging over.
suppressMessages({library(data.table); library(ggplot2); library(patchwork); library(LDscnR)})
A <- "/Volumes/Nemo/Nemo_sim/analysis_inputs"
D <- "module_sim_LDscnR/results/nulls_V2_c1"; E <- 1; LMIN <- 3; ALPHA <- 0.05
panel <- readRDS(file.path(A, sprintf("panel_V2_c1_env%d.rds", E)))
map <- flag_true_qtns(as.data.table(panel$map))
qtn <- map[true_pos_QTN %in% TRUE, .(Chr=as.character(Chr), Pos)]
ed  <- readRDS(file.path(D, sprintf("scan_V2_c1_env%d_emmax_genetic_B100.rds", E)))$edges
QS  <- seq(0, 0.95, by = 0.05)
RH  <- colnames(panel$ld_ws)

cell <- function(p, rc, q) {
  lw <- panel$ld_ws[, rc]
  cand <- which(lw >= stats::quantile(lw, q, na.rm = TRUE))
  if (!length(cand)) return(c(0L, 0L, 0L))
  qv <- p.adjust(p[cand], "BH")
  hit <- cand[which(qv < ALPHA)]
  if (!length(hit)) return(c(0L, 0L, 0L))
  ra <- ld_regions(map$marker[hit], ed); ra <- ra[lengths(ra) >= LMIN]
  if (!length(ra)) return(c(length(hit), 0L, 0L))
  co <- rbindlist(lapply(ra, function(m) { mm <- map[marker %in% m]
    data.table(chr=as.character(mm$Chr[1]), lo=min(mm$Pos), hi=max(mm$Pos)) }))
  tp <- sum(vapply(seq_len(nrow(co)), function(k)
    any(qtn$Chr==co$chr[k] & qtn$Pos>=co$lo[k] & qtn$Pos<=co$hi[k]), logical(1)))
  c(length(hit), nrow(co), tp)
}
g <- rbindlist(lapply(c("emmax","lfmm"), function(eng) {
  p <- readRDS(file.path(A, sprintf("pvals_V2_c1_env%d_%s_%s_B100.rds", E, eng,
        if (eng=="emmax") "genetic" else "env_orth")))$p_obs
  rbindlist(lapply(RH, function(rc) rbindlist(lapply(QS, function(q) {
    v <- cell(p, rc, q)
    data.table(engine=eng, rho=as.numeric(sub("rho_","",rc)), qstar=q,
               markers=v[1], regions=v[2], tp=v[3]) })))) }))
fwrite(g, file.path(D, sprintf("rho_q_grid_env%d.csv", E)))

hm <- function(v, ttl, pal) ggplot(g, aes(factor(qstar), factor(rho), fill=.data[[v]])) +
  geom_tile(colour="white", linewidth=.25) +
  facet_wrap(~engine, nrow=1) +
  scale_fill_distiller(palette=pal, direction=1, name=NULL, trans="sqrt") +
  labs(x=expression(q^"*"~"(ld_w quantile kept)"), y=expression(rho~"(ld_w window)"), title=ttl) +
  theme_bw(base_size=10) +
  theme(strip.background=element_blank(), panel.grid=element_blank(),
        axis.text.x=element_text(size=6, angle=90, vjust=.5), axis.text.y=element_text(size=6))
p <- (hm("regions", sprintf("regions at l_min = %d, per (rho, q*) cell", LMIN), "Blues") /
      hm("tp", "of which contain a detectable QTN", "Greens")) +
  plot_annotation(
    title=sprintf("V2_c1 env%d -- the C-score's integration grid, %d cells per engine", E, length(RH)*length(QS)),
    subtitle="Each cell is one analysis: keep the top (1-q*) of ld_w at that rho, BH within it, alpha 0.05.\nC counts the fraction of these cells in which a marker is a hit.",
    theme=theme(plot.subtitle=element_text(size=8, colour="grey35")))
ggsave(sprintf("module_sim_LDscnR/figures/rho_q_grid_env%d.png", E), p, width=13, height=10, dpi=170)
cat(sprintf("  wrote rho_q_grid_env%d.png\n", E))
print(g[, .(cells=.N, cells_with_regions=sum(regions>0), max_regions=max(regions),
            cells_with_TP=sum(tp>0), max_TP=max(tp)), by=engine])
