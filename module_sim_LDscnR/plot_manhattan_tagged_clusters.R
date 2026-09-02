## Single-SNP manhattan. Every SNP is a point at its own -log10(p). The stage-2
## clusters that contain at least one SIGNIFICANT SNP are coloured, one colour
## per cluster (recycled), and ALL members of such a cluster take that colour --
## so the horizontal extent of each flagged cluster is visible, not just the SNP
## that triggered it. Driving QTN keep the triangle shape and take the colour of
## the cluster they belong to (grey if that cluster was not flagged).
suppressMessages({library(data.table); library(LDscnR); library(ggplot2)})
SIM  <- Sys.getenv("SIM_DATA", "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5")
OUT  <- Sys.getenv("OUT", "module_sim_LDscnR/figures")
CELL <- Sys.getenv("CELL", "V0.5_c1"); TAG <- Sys.getenv("TAG", "nobgs")
ENV  <- as.integer(Sys.getenv("ENV", "3")); FILES <- 1:10
ALPHA <- as.numeric(Sys.getenv("ALPHA", "0.05"))
## SEL = "ld_w" reproduces the mechanism figure at SNP level: significance is BH
## inside the top-KSEL SNPs by ld_w, so only clusters that survive FILTERING get
## a colour. SEL = "none" is the conventional scan.
SEL   <- Sys.getenv("SEL", "none")
KSEL  <- as.integer(Sys.getenv("KSEL", "5000"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

per_file <- function(i) {
  f <- sprintf("%s/adapt_%s_chr%d_%s_env%d.rds", SIM, TAG, i, CELL, ENV)
  if (!file.exists(f)) return(NULL)
  x <- readRDS(f); m <- as.data.table(x$map)
  pr <- ld_prune_and_eMLG(GTs = x$GTs, stage1 = x$complexity_reduction$stage1,
          LD_decay = x$LD_decay, ld_w_col = "ld_w_095", ld_w_threshold = 0.025,
          score_threshold = 0.80, min_r2_rho = 0.5, distance_threshold = 1e5,
          compute_unflagged_eMLG = FALSE, cores = 1)
  stopifnot(identical(sort(pr$pruned), sort(x$grm_markers)))
  g  <- as.data.table(pr$groups)
  ms <- ld_group_map(g, prefix = i)[, .(marker, CL = group_id)]
  m <- merge(m, ms, by = "marker", all.x = TRUE)
  m[, `:=`(set = i, chr_lab = paste0("s", i, "_", Chr),
           driving = true_QTN %in% TRUE & MAF > 0.1 & p_Va > 0.05)]
  m[, .(marker, Chr, Pos, p = emx_p, ld_w = ld_w_095, CL, set, chr_lab, driving, true_QTN)]
}
m <- rbindlist(lapply(FILES, per_file), fill = TRUE)
m <- m[is.finite(p)]
cat(sprintf("  %d SNPs, %d in a stage-2 cluster, %d driving QTN\n",
            nrow(m), sum(!is.na(m$CL)), sum(m$driving, na.rm = TRUE)))

## significance. Conventional scan, or BH inside the ld_w-selected SNPs.
m[, q_all := p.adjust(p, "BH")]
cut_all <- if (any(m$q_all < ALPHA)) max(m$p[m$q_all < ALPHA]) else ALPHA/nrow(m)
m[, sel := TRUE]
if (SEL == "ld_w") {
  m[, sel := FALSE]
  m[head(order(-ld_w), KSEL), sel := TRUE]
  qq <- rep(NA_real_, nrow(m)); qq[m$sel] <- p.adjust(m$p[m$sel], "BH")
  m[, q := qq]
} else {
  m[, q := q_all]
}
m[, sig := !is.na(q) & q < ALPHA]
cut_used <- if (any(m$sig)) max(m$p[m$sig]) else ALPHA/sum(m$sel)
cat(sprintf("  selection %s (k = %s): BH cutoff p = %.3g; conventional cutoff %.3g\n",
            SEL, if (SEL == "none") "all" else format(KSEL, big.mark=","), cut_used, cut_all))
flagged <- unique(m[sig == TRUE & !is.na(CL)]$CL)
cat(sprintf("  %d significant SNPs, tagging %d distinct stage-2 clusters\n",
            sum(m$sig), length(flagged)))
cat(sprintf("  driving QTN inside a flagged cluster: %d of %d\n",
            sum(m$driving & m$CL %in% flagged, na.rm = TRUE), sum(m$driving, na.rm = TRUE)))

setorder(m, set, Chr, Pos)
off <- m[, .(len = max(Pos)), by = chr_lab][, .(chr_lab, off = cumsum(c(0, head(len, -1))))]
m <- merge(m, off, by = "chr_lab"); m[, gpos := (Pos + off)/1e6]


pal <- c("#1F6F8B","#C1622F","#2E7156","#8E5AA8","#B0392B","#3D7EA6","#7A6A1F",
         "#456A8C","#9E4630","#2F7F6F","#6B4E8C","#8A5A2B","#4E7A3C","#A2553F")
fl  <- sort(flagged)
cm  <- setNames(pal[(seq_along(fl) - 1) %% length(pal) + 1], fl)
m[, col := "grey84"]
m[CL %in% fl, col := cm[CL]]
m[, in_flagged := CL %in% fl]

sub_txt <- if (SEL == "none") {
  sprintf(paste0("%s significant SNPs tag %d stage-2 clusters; every member of a tagged cluster takes its colour. ",
                 "Triangles: %d driving QTN, coloured by their own cluster (grey if untagged)."),
          format(sum(m$sig), big.mark=","), length(fl), sum(m$driving, na.rm=TRUE))
} else {
  sprintf(paste0("Only clusters surviving the filter are coloured: %s SNPs pass BH inside the %s selected, tagging %d clusters. ",
                 "Solid line = conventional cutoff; dashed = the relaxed cutoff filtering buys. Triangles: %d driving QTN, %d in a coloured cluster."),
          format(sum(m$sig), big.mark=","), format(KSEL, big.mark=","), length(fl),
          sum(m$driving, na.rm=TRUE), sum(m$driving & m$CL %in% fl, na.rm=TRUE))
}

p <- ggplot() +
  geom_hline(yintercept = -log10(cut_all), colour = "#1F3F51", linewidth = .4) +
  {if (SEL != "none") geom_hline(yintercept = -log10(cut_used), colour = "#C1622F",
                                 linewidth = .4, linetype = "42") else NULL} +
  geom_point(data = m[in_flagged == FALSE & driving == FALSE],
             aes(gpos, -log10(p)), colour = "grey84", size = .32, alpha = .7, shape = 16) +
  geom_point(data = m[in_flagged == TRUE & driving == FALSE],
             aes(gpos, -log10(p), colour = col), size = .60, alpha = .95, shape = 16) +
  geom_point(data = m[driving == TRUE],
             aes(gpos, -log10(p), fill = col), shape = 24, size = 2.9,
             colour = "black", stroke = .45) +
  scale_colour_identity() + scale_fill_identity() +
  annotate("text", x = max(m$gpos) * .995, y = -log10(cut_all) + .35, hjust = 1, size = 2.6,
           colour = "#1F3F51", label = sprintf("BH 0.05 over all %s SNPs", format(nrow(m), big.mark = ","))) +
  {if (SEL != "none") annotate("text", x = max(m$gpos) * .995, y = -log10(cut_used) - .45, hjust = 1,
       size = 2.6, colour = "#C1622F",
       label = sprintf("BH 0.05 inside the %s ld_w-selected SNPs", format(KSEL, big.mark = ","))) else NULL} +
  labs(x = "genome position (Mb, 20 chromosome arms concatenated)",
       y = expression(-log[10](p)~"(EMMAX, single SNP)"),
       title = sprintf("%s  -  %s %s env%d",
                       if (SEL == "none") "Single-SNP scan, coloured by the stage-2 cluster each significant SNP tags"
                       else "Single-SNP scan: only the stage-2 clusters surviving the ld_w filter are coloured",
                       CELL, TAG, ENV),
       subtitle = sub_txt) +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "none",
        plot.subtitle = element_text(colour = "grey30", size = 7.4))
fn <- file.path(OUT, sprintf("manhattan_tagged_clusters_%s%s_%s_env%d.png",
                             if (SEL == "none") "" else paste0(SEL, KSEL, "_"), CELL, TAG, ENV))
ggsave(fn, p, width = 11.5, height = 5.2, dpi = 190)
cat(sprintf("  written: %s\n", fn))
