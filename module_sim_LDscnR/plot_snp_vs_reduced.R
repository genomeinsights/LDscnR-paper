## Single-SNP EMMAX against the complexity-reduced view of the SAME p-values.
##
## There is only ONE EMMAX run: every marker has its own p-value, computed
## against a GRM built from the pruned (LD-central) markers. "Complexity reduced"
## is not a second analysis -- it is the same p-values restricted to one
## representative per stage-2 cluster. What differs is therefore not the test but
## the number of tests, and hence the BH threshold:
##
##   all SNPs                 ~295k tests
##   cluster representatives  ~100k tests   (about 3x fewer)
##
## Both panels show the same driving QTN, so the comparison is what the reduction
## costs in coverage against what it buys in threshold.
suppressMessages({library(data.table); library(LDscnR); library(ggplot2); library(patchwork)})
SIM  <- Sys.getenv("SIM_DATA", "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5")
OUT  <- Sys.getenv("OUT", "module_sim_LDscnR/figures")
CELL <- Sys.getenv("CELL", "V0.5_c1"); TAG <- Sys.getenv("TAG", "nobgs")
ENV  <- as.integer(Sys.getenv("ENV", "3")); FILES <- 1:10; ALPHA <- 0.05
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
  m[, `:=`(is_rep = marker %in% pr$pruned, set = i,
           chr_lab = paste0("s", i, "_", Chr),
           driving = true_QTN %in% TRUE & MAF > 0.1 & p_Va > 0.05)]
  m[, .(marker, Chr, Pos, p = emx_p, is_rep, set, chr_lab, driving)]
}
m <- rbindlist(lapply(FILES, per_file), fill = TRUE)[is.finite(p)]
setorder(m, set, Chr, Pos)
off <- m[, .(len = max(Pos)), by = chr_lab][, .(chr_lab, off = cumsum(c(0, head(len, -1))))]
m <- merge(m, off, by = "chr_lab"); m[, gpos := (Pos + off)/1e6]

bh_cut <- function(p) { q <- p.adjust(p, "BH"); if (any(q < ALPHA)) max(p[q < ALPHA]) else ALPHA/length(p) }
cut_all <- bh_cut(m$p)
cut_rep <- bh_cut(m[is_rep == TRUE]$p)
qtn_all <- m[driving == TRUE]
cat(sprintf("  all SNPs        : %s tests, BH cutoff p = %.3g, %s significant, %d/%d driving QTN significant\n",
    format(nrow(m), big.mark=","), cut_all, format(sum(m$p <= cut_all), big.mark=","),
    sum(qtn_all$p <= cut_all), nrow(qtn_all)))
r <- m[is_rep == TRUE]
cat(sprintf("  representatives : %s tests, BH cutoff p = %.3g, %s significant, %d/%d driving QTN significant\n",
    format(nrow(r), big.mark=","), cut_rep, format(sum(r$p <= cut_rep), big.mark=","),
    sum(r[driving == TRUE]$p <= cut_rep), sum(m$driving & m$is_rep)))
cat(sprintf("  driving QTN that ARE their cluster's representative: %d of %d\n",
    sum(m$driving & m$is_rep), nrow(qtn_all)))

mk <- function(d, cutoff, ttl, sub) {
  ggplot(d, aes(gpos, -log10(p))) +
    geom_point(colour = "grey78", size = .33, alpha = .75, shape = 16) +
    geom_hline(yintercept = -log10(cutoff), colour = "#1F3F51", linewidth = .42) +
    geom_point(data = d[driving == TRUE], aes(gpos, -log10(p)),
               shape = 24, size = 2.5, fill = "#C1622F", colour = "black", stroke = .4) +
    annotate("text", x = max(m$gpos)*.995, y = -log10(cutoff) + .4, hjust = 1, size = 2.5,
             colour = "#1F3F51", label = sprintf("BH 0.05 (p = %.2g)", cutoff)) +
    coord_cartesian(ylim = c(0, max(-log10(m$p)) * 1.03)) +
    labs(x = NULL, y = expression(-log[10](p)), title = ttl, subtitle = sub) +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(),
          plot.subtitle = element_text(colour = "grey30", size = 7.4))
}
pA <- mk(m, cut_all, "A  Every SNP tested",
         sprintf("%s tests. %s significant; %d of %d driving QTN pass.",
                 format(nrow(m), big.mark=","), format(sum(m$p <= cut_all), big.mark=","),
                 sum(qtn_all$p <= cut_all), nrow(qtn_all)))
pB <- mk(r, cut_rep, "B  One LD-central representative per stage-2 cluster",
         sprintf("%s tests (%.1fx fewer), so the BH cutoff is %.1fx less severe. %s significant; %d of %d driving QTN are themselves representatives.",
                 format(nrow(r), big.mark=","), nrow(m)/nrow(r), cut_rep/cut_all,
                 format(sum(r$p <= cut_rep), big.mark=","),
                 sum(m$driving & m$is_rep), nrow(qtn_all))) +
      labs(x = "genome position (Mb, 20 chromosome arms concatenated)")
ggsave(file.path(OUT, sprintf("snp_vs_reduced_%s_%s_env%d.png", CELL, TAG, ENV)),
       pA / pB, width = 11, height = 7, dpi = 190)
cat(sprintf("  written: %s\n", file.path(OUT, sprintf("snp_vs_reduced_%s_%s_env%d.png", CELL, TAG, ENV))))
