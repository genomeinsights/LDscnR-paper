## module_sim/07_consensus_manhattan.R
## Four-panel Manhattan (EMMAX single, LFMM single, EMMAX ld_w, LFMM ld_w) over
## the SHARED region frame from 06_consensus_regions.R. The key point: colour is
## tied to the REGION, not the method, so a region caught by several methods is
## the SAME colour in every panel -- you can read across panels to see which
## method misses which region. TP regions get distinct colours; FP regions are
## grey; true QTNs marked with a triangle. y = genome-wide -log10 FDR(q); the
## dashed line is q=0.05 (the single-SNP outlier cut). ld_w panels colour the
## SNPs the cluster-null KEPT, which need not coincide with q<0.05.
## Reads consensus_V{V}_c{c}_env{env}.rds (sets/regions/lab) + light-pools the map.
## Run from LDscnR-paper/:  Rscript module_sim/07_consensus_manhattan.R [V c env]
## Output (git-ignored): module_sim/consensus_manhattan_V{V}_c{c}_env{env}.png

suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
mod <- "/Users/petrikem/gitlab/LDscnR-paper/module_sim"
dir_data <- file.path(mod, "data"); dir_fig <- file.path(mod, "figures")
for (d in c(dir_data, dir_fig)) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
a   <- commandArgs(trailingOnly = TRUE)
V   <- if (length(a) >= 1) a[1] else "2"
cc  <- if (length(a) >= 2) a[2] else "1"
env <- if (length(a) >= 3) a[3] else "1"
MIN_SNP <- if (length(a) >= 4) as.integer(a[4]) else 2L      # >=k-SNP filter for fig 2
C   <- readRDS(file.path(dir_data, sprintf("consensus_V%s_c%s_env%s.rds", V, cc, env)))
sets <- C$sets; regions <- C$regions; lab <- C$lab
map <- as.data.table(C$map)                        # compact: marker,Chr,Pos,emx_p,lfmm_p,true_pos_QTN

## genome coordinate transform (moduleB style)
map[, chr_num := as.integer(factor(Chr, levels = unique(Chr)))]
clen <- map[, .(len = max(Pos)), by = .(Chr, chr_num)][order(chr_num)]
spg  <- 0.005 * sum(clen$len); clen[, start := c(0, head(cumsum(len + spg), -1))]
map  <- clen[, .(Chr, start)][map, on = "Chr"][, x := Pos + start]
cmid <- clen[, .(Chr, mid = start + len / 2)]
shade <- clen[chr_num %% 2 == 0, .(xmin = start, xmax = start + len)]
xkey <- setNames(map$x, map$marker)
map[, `:=`(q_emx = p.adjust(emx_p, "BH"), q_lfmm = p.adjust(lfmm_p, "BH"))]
map[, `:=`(nlq_emx = -log10(pmax(q_emx, 1e-30)), nlq_lfmm = -log10(pmax(q_lfmm, 1e-30)))]

## Two visual channels:
##   COLOUR = shared (>=2-method) region identity, same across panels, so a region
##            several methods agree on (TP or FP) is the same colour everywhere;
##            method-unique regions are grey. A coloured point NOT over a QTN
##            triangle is a SHARED false positive.
##   SHAPE  = within-method SNP support: + = single-SNP outlier (that method tags
##            the region with just 1 SNP), filled circle = multi-SNP region.
region_dt <- as.data.table(C$region_dt)
setorder(region_dt, region)
m2r <- integer(0); for (i in seq_along(regions)) m2r[regions[[i]]] <- i
reg_x   <- vapply(regions, function(mk) stats::median(xkey[mk], na.rm = TRUE), numeric(1))
support <- region_dt$support; is_tp <- region_dt$is_TP
multi   <- which(support >= 2L); multi_ord <- multi[order(reg_x[multi])]
PAL <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFD400","#A65628","#F781BF",
         "#1B9E77","#D95F02","#7570B3","#66A61E","#E7298A","#1F78B4","#00CED1","#8B4513")
reg_col <- rep("grey55", length(regions))
reg_col[multi_ord] <- PAL[(seq_along(multi_ord) - 1) %% length(PAL) + 1]

qtn <- map[true_pos_QTN == TRUE]
panel <- function(nm, ycol, drop_single_snp = FALSE) {
  nsn <- region_dt[[paste0("nsnp_", nm)]]                 # per-region SNP count for this method
  out <- map[marker %in% sets[[nm]]]
  out[, region := m2r[marker]][, `:=`(col = reg_col[region], nsnp = nsn[region])]
  out[, shp := ifelse(nsnp == 1L, "single", "multi")]
  if (drop_single_snp) out <- out[nsnp >= MIN_SNP]
  n_hit <- uniqueN(out$region)
  ggplot() +
    geom_rect(data = shade, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), fill = "grey93") +
    geom_point(data = map, aes(x, .data[[ycol]]), color = "grey80", size = 0.2, alpha = 0.3) +
    geom_hline(yintercept = -log10(0.05), linetype = 2, color = "red", linewidth = 0.35) +
    geom_point(data = qtn, aes(x, y = -0.4), shape = 25, fill = "black", color = "black", size = 1.3) +
    geom_point(data = out, aes(x, .data[[ycol]], color = col, shape = shp), size = 1.4, stroke = 0.7) +
    scale_color_identity() +
    scale_shape_manual(values = c(multi = 16, single = 3),
                       labels = c(multi = "multi-SNP region", single = "single-SNP outlier"),
                       name = NULL, drop = FALSE) +
    scale_x_continuous(breaks = cmid$mid, labels = gsub("R([0-9]+)_Chr", "\\1.", cmid$Chr),
                       expand = c(0.01, 0.01)) +
    labs(x = NULL, y = expression(-log[10](q)),
         title = sprintf("%s  (regions hit: %d)", nm, n_hit)) +
    theme_bw(base_size = 9) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
          axis.text.x = element_text(size = 6), legend.position = "right")
}

make_fig <- function(drop_single_snp, suffix, subtitle) {
  P <- panel("EMMAX single", "nlq_emx",  drop_single_snp) /
       panel("LFMM single",  "nlq_lfmm", drop_single_snp) /
       panel("EMMAX ld_w",   "nlq_emx",  drop_single_snp) /
       panel("LFMM ld_w",    "nlq_lfmm", drop_single_snp) +
       patchwork::plot_layout(guides = "collect") +
       patchwork::plot_annotation(
         title = sprintf("V%s_c%s_env%s: legacy single-SNP vs ld_w+null on a shared region frame%s",
                         V, cc, env, if (drop_single_snp) " (single-SNP outliers removed)" else ""),
         subtitle = subtitle)
  fn <- sprintf("consensus_manhattan%s_V%s_c%s_env%s.png", suffix, V, cc, env)
  ggsave(file.path(dir_fig, fn), P, width = 16, height = 10, dpi = 150)
  cat("wrote", fn, "\n")
}

make_fig(FALSE, "",           "colour = shared (>=2-method) region identity; grey = method-unique; + = single-SNP outlier; triangle = true QTN; dashed = q0.05")
make_fig(TRUE,  "_snpfilter", sprintf("single-SNP outliers removed: only regions a method tags with >=%d SNPs remain", MIN_SNP))
cat(sprintf("shared (>=2-method) regions: %d (%d TP, %d FP) | method-unique: %d\n",
    length(multi), sum(is_tp[multi]), sum(!is_tp[multi]), sum(support == 1L)))
