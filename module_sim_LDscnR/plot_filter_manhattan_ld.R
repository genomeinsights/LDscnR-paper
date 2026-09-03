## Manhattan of the pooled panel, coloured by r2 with the nearest driving QTN and
## with QTN-carrying clusters drawn as triangles.
##
## r2 is computed from genotypes for every marker within R2WIN of a driving QTN,
## then summarised per cluster as the MAXIMUM over its members -- a cluster tags
## a QTN if any member does. Causal markers are still EXCLUDED from the tested
## set (they never contribute a p-value); the triangle marks the cluster that
## carried one, which is what an analyst would actually recover.
suppressMessages({library(data.table); library(LDscnR); library(ggplot2)})
SIM  <- Sys.getenv("SIM_DATA", "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5")
OUT  <- Sys.getenv("OUT", "module_sim_LDscnR/figures")
CELL <- Sys.getenv("CELL", "V0.5_c1"); TAG <- Sys.getenv("TAG", "nobgs")
ENV  <- as.integer(Sys.getenv("ENV", "3")); FILES <- 1:10
KSEL <- as.integer(Sys.getenv("KSEL", "5000")); ALPHA <- 0.05
R2WIN <- as.numeric(Sys.getenv("R2WIN", "2e6"))
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
  ms <- data.table(marker = unlist(g$members, use.names = FALSE), CL_id = rep.int(g$group_id, lengths(g$members)))
  m <- merge(m, ms, by = "marker", all.x = TRUE)[!is.na(CL_id)]

  drv <- m[true_QTN %in% TRUE & MAF > 0.1 & p_Va > 0.05]
  m[, r2_qtn := 0]
  if (nrow(drv)) for (j in seq_len(nrow(drv))) {
    ch <- as.character(drv$Chr[j])
    ii <- which(as.character(m$Chr) == ch & abs(m$Pos - drv$Pos[j]) <= R2WIN)
    if (!length(ii)) next
    r2 <- suppressWarnings(cor(x$GTs[, drv$marker[j]], x$GTs[, m$marker[ii]],
                               use = "pairwise.complete.obs")^2)
    r2[!is.finite(r2)] <- 0
    m[ii, r2_qtn := pmax(r2_qtn, as.numeric(r2))]
  }
  ## which clusters carried a causal variant, before we drop them
  qtn_cl <- unique(m[true_QTN %in% TRUE]$CL_id)
  m <- m[!(true_QTN %in% TRUE)]
  rep_mk <- intersect(pr$pruned, m$marker)
  u <- m[, .(r2_qtn = max(r2_qtn), Chr = Chr[1], Pos = median(Pos), n_loci = .N,
             ld_w = median(ld_w_095, na.rm = TRUE)), by = CL_id]
  rp <- m[marker %in% rep_mk, .(CL_id, p = emx_p)][, .SD[1], by = CL_id]
  u  <- merge(u, rp, by = "CL_id")[is.finite(p)]
  u[, `:=`(has_qtn = CL_id %in% qtn_cl, set = i, chr_lab = paste0("set", i, "_", Chr))]
  u[]
}
u <- rbindlist(lapply(FILES, units_for), fill = TRUE)
cat(sprintf("  %d units; %d carry a driving QTN; %d have r2 >= 0.2 with one\n",
            nrow(u), sum(u$has_qtn), sum(u$r2_qtn >= 0.2)))

setorder(u, set, Chr, Pos)
off <- u[, .(len = max(Pos)), by = chr_lab][, .(chr_lab, off = cumsum(c(0, head(len, -1))))]
u <- merge(u, off, by = "chr_lab"); u[, gpos := (Pos + off)/1e6]

bh_cut <- function(p) { q <- p.adjust(p, "BH"); if (any(q < ALPHA)) max(p[q < ALPHA]) else ALPHA/length(p) }
sel <- head(order(-u$ld_w), KSEL)
cut_gw <- bh_cut(u$p); cut_sel <- bh_cut(u$p[sel])
cat(sprintf("  BH cutoffs: genome-wide %.3g, inside %d selected %.3g\n", cut_gw, KSEL, cut_sel))

setorder(u, r2_qtn)                      # high-LD points drawn last, on top
p <- ggplot(u, aes(gpos, -log10(p))) +
  geom_hline(yintercept = -log10(cut_gw),  colour = "#1F3F51", linewidth = .45) +
  geom_hline(yintercept = -log10(cut_sel), colour = "#C1622F", linewidth = .45, linetype = "42") +
  geom_point(data = u[has_qtn == FALSE], aes(colour = r2_qtn), size = .55, alpha = .9, shape = 16) +
  geom_point(data = u[has_qtn == TRUE], aes(fill = r2_qtn), shape = 24, size = 2.6,
             colour = "black", stroke = .45) +
  scale_colour_gradientn(colours = c("grey84","#BBD3DE","#6FA8C0","#1F6F8B","#8C2F1E"),
                         values = c(0,.08,.25,.55,1), limits = c(0,1), name = expression(r^2~"with driving QTN")) +
  scale_fill_gradientn(colours = c("grey84","#BBD3DE","#6FA8C0","#1F6F8B","#8C2F1E"),
                       values = c(0,.08,.25,.55,1), limits = c(0,1), guide = "none") +
  annotate("text", x = max(u$gpos)*.995, y = -log10(cut_gw) + .3, hjust = 1, size = 2.6,
           colour = "#1F3F51", label = sprintf("genome-wide BH over %s units", format(nrow(u), big.mark=","))) +
  annotate("text", x = max(u$gpos)*.995, y = -log10(cut_sel) - .34, hjust = 1, size = 2.6,
           colour = "#C1622F", label = sprintf("BH inside the %s ld_w-selected units", format(KSEL, big.mark=","))) +
  labs(x = "genome position (Mb, 20 chromosome arms concatenated)",
       y = expression(-log[10](p)~"(EMMAX)"),
       title = sprintf("Clusters coloured by LD with the nearest driving QTN  -  %s %s env%d", CELL, TAG, ENV),
       subtitle = sprintf(paste0("Triangles: the %d clusters that carry a causal variant (excluded from testing). ",
                          "%d of %s clusters reach r2 >= 0.2 with a QTN."),
                          sum(u$has_qtn), sum(u$r2_qtn >= 0.2), format(nrow(u), big.mark=","))) +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        legend.key.width = unit(1.5, "cm"),
        plot.subtitle = element_text(colour = "grey30", size = 7.6))
fn <- file.path(OUT, sprintf("filter_manhattan_ldQTN_%s_%s_env%d.png", CELL, TAG, ENV))
ggsave(fn, p, width = 11, height = 5.2, dpi = 190)
cat(sprintf("  written: %s\n", fn))
