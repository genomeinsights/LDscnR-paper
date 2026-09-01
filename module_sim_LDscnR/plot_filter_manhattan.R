## Manhattan showing WHY filter-then-test gains discoveries: it does not change
## any p-value, it changes the BH threshold. Genome-wide BH over ~86k units is
## severe; BH over the k selected units is not. Same data, different line.
suppressMessages({library(data.table); library(LDscnR); library(ggplot2)})
SIM  <- Sys.getenv("SIM_DATA", "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5")
OUT  <- Sys.getenv("OUT", "module_sim_LDscnR/figures")
CELL <- Sys.getenv("CELL", "V0.5_c1"); TAG <- Sys.getenv("TAG", "nobgs")
ENV  <- as.integer(Sys.getenv("ENV", "3")); FILES <- 1:10
KSEL <- as.integer(Sys.getenv("KSEL", "5000")); ALPHA <- 0.05
WIN  <- as.numeric(Sys.getenv("WINDOW", "50")) * 1000
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

units_for <- function(i) {
  f <- sprintf("%s/adapt_%s_chr%d_%s_env%d.rds", SIM, TAG, i, CELL, ENV)
  if (!file.exists(f)) return(NULL)
  x <- readRDS(f); m <- as.data.table(x$map)
  pr <- ld_prune_and_eMLG(GTs = x$GTs, stage1 = x$complexity_reduction$stage1,
          LD_decay = x$LD_decay, ld_w_col = "ld_w_095", ld_w_threshold = 0.025,
          score_threshold = 0.80, min_r2_rho = 0.5, distance_threshold = 1e5,
          compute_unflagged_eMLG = FALSE, cores = 1)
  stopifnot(identical(sort(pr$pruned), sort(x$grm_markers)))
  g  <- as.data.table(pr$groups)
  ms <- rbindlist(lapply(seq_len(nrow(g)), function(k)
          data.table(marker = g$members[[k]], CL_id = g$group_id[k])))
  m <- merge(m, ms, by = "marker", all.x = TRUE)[!is.na(CL_id)]
  drv <- m[true_QTN %in% TRUE & MAF > 0.1 & p_Va > 0.05]
  m[, d_qtn := Inf]
  if (nrow(drv)) for (ch in unique(m$Chr)) {
    dd <- drv[Chr == ch]; if (!nrow(dd)) next
    ii <- which(m$Chr == ch)
    D <- abs(outer(m$Pos[ii], dd$Pos, "-")); m[ii, d_qtn := D[cbind(seq_along(ii), max.col(-D))]]
  }
  m <- m[!(true_QTN %in% TRUE)]
  rep_mk <- intersect(pr$pruned, m$marker)
  u <- m[, .(ld_w = median(ld_w_095, na.rm=TRUE), n_loci = .N,
             Chr = Chr[1], Pos = median(Pos), d_qtn = min(d_qtn)), by = CL_id]
  rp <- m[marker %in% rep_mk, .(CL_id, p = emx_p)][, .SD[1], by = CL_id]
  u  <- merge(u, rp, by = "CL_id")[is.finite(p)]
  u[, `:=`(set = i, qtn_pos = list(drv$Pos), chr_lab = paste0("set", i, "_", Chr))]
  u[]
}
u <- rbindlist(lapply(FILES, units_for), fill = TRUE)
cat(sprintf("  panel: %s %s env%d -- %d units over %d chromosome arms\n",
            CELL, TAG, ENV, nrow(u), uniqueN(u$chr_lab)))

## cumulative genome axis
setorder(u, set, Chr, Pos)
off <- u[, .(len = max(Pos)), by = chr_lab][, .(chr_lab, off = cumsum(c(0, head(len, -1))) )]
u <- merge(u, off, by = "chr_lab"); u[, gpos := (Pos + off)/1e6]

## BH cutoffs: genome-wide vs inside the ld_w-selected set
bh_cut <- function(p) { q <- p.adjust(p, "BH"); if (any(q < ALPHA)) max(p[q < ALPHA]) else ALPHA/length(p) }
sel <- head(order(-u$ld_w), KSEL)
u[, selected := FALSE]; u[sel, selected := TRUE]
cut_gw  <- bh_cut(u$p)
cut_sel <- bh_cut(u$p[sel])
u[, `:=`(sig_gw = p <= cut_gw, sig_sel = selected & p <= cut_sel, tp = d_qtn < WIN)]
newly <- u[sig_sel & !sig_gw]
lost  <- u[sig_gw & !sig_sel]          ## genome-wide hits OUTSIDE the selected set
cat(sprintf("  genome-wide BH cutoff p = %.3g  -> %d discoveries (%d true)\n",
            cut_gw, sum(u$sig_gw), sum(u$sig_gw & u$tp)))
cat(sprintf("  filtered   BH cutoff p = %.3g  -> %d discoveries (%d true)\n",
            cut_sel, sum(u$sig_sel), sum(u$sig_sel & u$tp)))
cat(sprintf("  NEWLY found by filtering: %d, of which %d within %.0f kb of a driving QTN (%.0f%%)\n",
            nrow(newly), sum(newly$tp), WIN/1000, 100*mean(newly$tp)))
cat(sprintf("  LOST by filtering (significant genome-wide, not in the selected set): %d, of which %d QTN-proximal\n",
            nrow(lost), sum(lost$tp)))
cat(sprintf("  accounting: %d - %d + %d = %d\n", sum(u$sig_gw), nrow(lost), nrow(newly), sum(u$sig_sel)))

u[, cls := fifelse(sig_sel & !sig_gw, "gained by filtering",
            fifelse(sig_gw, "found by genome-wide BH",
             fifelse(selected, "selected, not significant", "not selected")))]
u[, cls := factor(cls, levels = c("not selected","selected, not significant",
                                  "found by genome-wide BH","gained by filtering"))]
pal <- c("not selected" = "grey82", "selected, not significant" = "#9FB8C4",
         "found by genome-wide BH" = "#1F3F51", "gained by filtering" = "#C1622F")
qtn <- u[d_qtn == 0 | d_qtn < 1, .(gpos)]
qmark <- u[, .(d = min(d_qtn)), by = chr_lab][d < WIN]
qpos <- u[d_qtn < 1e4, .(gpos = median(gpos)), by = chr_lab]

p <- ggplot(u[order(cls)], aes(gpos, -log10(p), colour = cls)) +
  geom_vline(data = qpos, aes(xintercept = gpos), colour = "grey45", linetype = "22", linewidth = .35) +
  geom_point(size = .5, alpha = .85) +
  geom_hline(yintercept = -log10(cut_gw),  colour = "#1F3F51", linewidth = .45) +
  geom_hline(yintercept = -log10(cut_sel), colour = "#C1622F", linewidth = .45, linetype = "42") +
  annotate("text", x = max(u$gpos)*.995, y = -log10(cut_gw)  + .28, hjust = 1, size = 2.6,
           colour = "#1F3F51", label = sprintf("genome-wide BH over %s units", format(nrow(u), big.mark=","))) +
  annotate("text", x = max(u$gpos)*.995, y = -log10(cut_sel) - .34, hjust = 1, size = 2.6,
           colour = "#C1622F", label = sprintf("BH inside the %s selected units", format(KSEL, big.mark=","))) +
  scale_colour_manual(values = pal, name = NULL) +
  labs(x = "genome position (Mb, 20 chromosome arms concatenated)", y = expression(-log[10](p)~"(EMMAX)"),
       title = sprintf("Filtering changes the threshold, not the p-values  -  %s %s env%d", CELL, TAG, ENV),
       subtitle = sprintf(paste0("Selecting the %s clusters with highest ld_w relaxes the BH cutoff from p = %.2g to p = %.2g. ",
                          "%d units are gained; %.0f%% lie within %.0f kb of a driving QTN (dashed grey)."),
                          format(KSEL, big.mark=","), cut_gw, cut_sel, nrow(newly), 100*mean(newly$tp), WIN/1000)) +
  guides(colour = guide_legend(override.aes = list(size = 2.4))) +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.subtitle = element_text(colour = "grey30", size = 7.6))
fn <- file.path(OUT, sprintf("filter_manhattan_%s_%s_env%d.png", CELL, TAG, ENV))
ggsave(fn, p, width = 11, height = 5.2, dpi = 190)
cat(sprintf("  written: %s\n", fn))
