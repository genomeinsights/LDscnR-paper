## module_sticklebacks/04_manhattan.R
## Genome-wide SNP Manhattan (moduleB style) for the three cluster-level
## analyses: all SNPs grey, member SNPs of each method's SIGNIFICANT clusters
## coloured one colour per cluster. Three stacked panels.
## Run from LDscnR-paper/:  Rscript module_sticklebacks/04_manhattan.R

suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
mod <- "/Users/petrikem/gitlab/LDscnR-paper/module_sticklebacks"
args <- commandArgs(trailingOnly = TRUE)
res_file <- if (length(args) >= 1) args[1] else file.path(mod, "results_three_nulls.rds")
out_png  <- if (length(args) >= 2) args[2] else file.path(mod, "fig_manhattan_three.png")
R  <- readRDS(res_file)
for (k in names(R)) { R[[k]]$clusters <- as.data.table(R[[k]]$clusters)
                      R[[k]]$candidates <- as.data.table(R[[k]]$candidates) }
sr <- readRDS(file.path(mod, "snp_stats_aligned.rds")); setDT(sr)

## ---- shared genome coordinate transform (moduleB style) -------------------
mapx <- copy(sr[, .(marker, Chr, Pos)])
mapx[, chr_num := suppressWarnings(as.integer(gsub("[^0-9]", "", Chr)))]
clen <- mapx[, .(len = max(Pos)), by = .(Chr, chr_num)][order(chr_num)]
spg  <- 0.005 * sum(clen$len)
clen[, start := c(0, head(cumsum(len + spg), -1))]
mapx <- clen[, .(Chr, start)][mapx, on = "Chr"][, x := Pos + start]
xkey <- setNames(mapx$x, mapx$marker)
cmid <- clen[, .(Chr, mid = start + len / 2)]
shade <- clen[chr_num %% 2 == 0, .(xmin = start, xmax = start + len)]
chr_lab <- function() scale_x_continuous(breaks = cmid$mid, labels = gsub("^Chr", "", cmid$Chr),
                                         expand = c(0.01, 0.01))

## per-SNP y (-log10 p, clamped) + genome x
sr[, `:=`(x = xkey[marker],
          nlp_emx  = -log10(pmax(emx_p,  1e-30)),
          nlp_lfmm = -log10(pmax(lfmm_p, 1e-30)))]
PAL <- rep(c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#A65628","#F781BF",
             "#1B9E77","#D95F02","#7570B3","#66A61E","#E7298A","#1F78B4","#33A02C"), 8)

## ALL_MEMBERS = TRUE colours every candidate member of a significant cluster
## (significant SNPs + the LD-recruited non-significant ones, which fall below the
## line); FALSE colours only the significant members.
ALL_MEMBERS <- TRUE
panel <- function(res, ycol, title) {
  sig_ids <- res$clusters[significant == TRUE, cluster]
  memb <- if (ALL_MEMBERS) res$candidates[cluster %in% sig_ids, .(marker, cluster)]
          else             res$candidates[cluster %in% sig_ids & sig == TRUE, .(marker, cluster)]
  memb[, x := xkey[marker]][, y := sr[[ycol]][match(marker, sr$marker)]]
  clo  <- memb[, .(mx = min(x)), by = cluster][order(mx)][, col := PAL[seq_len(.N)]]
  memb[, col := clo$col[match(cluster, clo$cluster)]]
  ## thin the grey background: keep all suggestive SNPs + a systematic sample
  bg <- sr[sr[[ycol]] > 1 | (seq_len(nrow(sr)) %% 25 == 0)]
  ggplot() +
    geom_rect(data = shade, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), fill = "grey93") +
    geom_point(data = bg, aes(x, .data[[ycol]]), color = "grey80", size = 0.25, alpha = 0.5) +
    geom_point(data = memb, aes(x, y, color = col), size = 1.1) +
    geom_hline(yintercept = -log10(res$pcut), linetype = 2, color = "red", linewidth = 0.4) +
    scale_color_identity() + chr_lab() +
    labs(x = "Chromosome", y = expression(-log[10](p)),
         title = sprintf("%s  (%d significant clusters)", title, length(sig_ids))) +
    theme_bw(base_size = 10) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
}

P <- panel(R$emx_perm, "nlp_emx",  "EMMAX + permutation") /
     panel(R$emx_bg,   "nlp_emx",  "EMMAX + background") /
     panel(R$lfmm_bg,  "nlp_lfmm", "LFMM + background")
ggsave(out_png, P, width = 16, height = 11, units = "in", dpi = 150)
cat("wrote fig_manhattan_three.png\n")
