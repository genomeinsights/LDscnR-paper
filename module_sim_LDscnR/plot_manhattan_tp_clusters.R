## Manhattan with each TRUE-POSITIVE stage-2 cluster in its own colour.
##
## A cluster is a true positive under the OR convention used throughout this
## module: it tags a driving QTN (some member at r2 >= r2min within dmax) AND it
## is the BEST such cluster for that QTN. Any other cluster tagging the same QTN
## is a false positive by dedup -- one region per QTN, the one in highest LD.
## Colours are recycled across TP clusters; everything else is grey. Triangles
## mark clusters carrying a causal variant (excluded from testing).
suppressMessages({library(data.table); library(LDscnR); library(ggplot2)})
SIM  <- Sys.getenv("SIM_DATA", "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5")
OUT  <- Sys.getenv("OUT", "module_sim_LDscnR/figures")
CELL <- Sys.getenv("CELL", "V0.5_c1"); TAG <- Sys.getenv("TAG", "nobgs")
ENV  <- as.integer(Sys.getenv("ENV", "3")); FILES <- 1:10
KSEL <- as.integer(Sys.getenv("KSEL", "5000")); ALPHA <- 0.05
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
  m <- merge(m, ms, by = "marker", all.x = TRUE)[!is.na(CL)]
  th  <- score_thresholds(as.data.table(x$LD_decay$decay_sum),
                          rho_r2 = 0.75, rho_d = 0.95, dmax_cap = 1e5)
  drv <- m[true_QTN %in% TRUE & MAF > 0.1 & p_Va > 0.05]

  ## per (cluster, QTN) maximum r2, restricted to the matching window
  link <- NULL
  if (nrow(drv)) link <- rbindlist(lapply(seq_len(nrow(drv)), function(j) {
    ch <- as.character(drv$Chr[j])
    near <- m[as.character(Chr) == ch & abs(Pos - drv$Pos[j]) < th$dmax & !(true_QTN %in% TRUE)]
    if (!nrow(near)) return(NULL)
    r2 <- suppressWarnings(cor(x$GTs[, drv$marker[j]], x$GTs[, near$marker],
                               use = "pairwise.complete.obs")^2)
    r2[!is.finite(r2)] <- 0
    dt <- data.table(CL = near$CL, r2 = as.numeric(r2))[r2 >= th$r2min]
    if (!nrow(dt)) return(NULL)
    dt[, .(r2 = max(r2)), by = CL][, qtn := paste0(i, "_", ch, "_", drv$Pos[j])][]
  }))
  qtn_cl <- unique(m[true_QTN %in% TRUE]$CL)
  m <- m[!(true_QTN %in% TRUE)]
  u <- m[, .(Chr = Chr[1], Pos = median(Pos), p = NA_real_, n = .N), by = CL]
  rp <- m[marker %in% pr$pruned, .(CL, pp = emx_p)][, .SD[1], by = CL]
  u <- merge(u[, !"p"], rp, by = "CL"); setnames(u, "pp", "p")
  u <- u[is.finite(p)]
  u[, `:=`(has_qtn = CL %in% qtn_cl, set = i, chr_lab = paste0("s", i, "_", Chr))]
  list(u = u, link = link, r2min = th$r2min, dmax = th$dmax)
}
parts <- lapply(FILES, per_file); parts <- Filter(Negate(is.null), parts)
u    <- rbindlist(lapply(parts, `[[`, "u"), fill = TRUE)
link <- rbindlist(lapply(parts, `[[`, "link"), fill = TRUE)

## DEDUP: for each QTN keep only the highest-r2 cluster. That cluster is the TP.
tp <- if (nrow(link)) link[order(-r2)][, .SD[1], by = qtn] else data.table(CL=character(), qtn=character())
u[, tp_qtn := NA_character_]
if (nrow(tp)) u[tp, on = "CL", tp_qtn := i.qtn]
cat(sprintf("  %d clusters; %d tag a QTN; %d are TRUE POSITIVES after dedup (%d QTN reachable)\n",
            nrow(u), if (nrow(link)) uniqueN(link$CL) else 0L, sum(!is.na(u$tp_qtn)),
            if (nrow(tp)) nrow(tp) else 0L))

setorder(u, set, Chr, Pos)
off <- u[, .(len = max(Pos)), by = chr_lab][, .(chr_lab, off = cumsum(c(0, head(len,-1))))]
u <- merge(u, off, by = "chr_lab"); u[, gpos := (Pos + off)/1e6]
bh <- function(p) { q <- p.adjust(p, "BH"); if (any(q < ALPHA)) max(p[q < ALPHA]) else ALPHA/length(p) }
sel <- head(order(-u$n), 0)  # placeholder, ld_w selection not needed for this figure
cut_gw <- bh(u$p)

pal <- c("#1F6F8B","#C1622F","#2E7156","#8E5AA8","#B0392B","#3D7EA6","#7A6A1F",
         "#456A8C","#9E4630","#2F7F6F","#6B4E8C","#8A5A2B")
tpl <- sort(unique(na.omit(u$tp_qtn)))
u[, col := "grey80"]
if (length(tpl)) {
  cm <- setNames(pal[(seq_along(tpl)-1) %% length(pal) + 1], tpl)
  u[!is.na(tp_qtn), col := cm[tp_qtn]]
}
u[, is_tp := !is.na(tp_qtn)]

p <- ggplot(u[order(is_tp)], aes(gpos, -log10(p))) +
  geom_hline(yintercept = -log10(cut_gw), colour = "#1F3F51", linewidth = .45) +
  geom_point(data = u[is_tp == FALSE], colour = "grey80", size = .45, alpha = .8, shape = 16) +
  geom_point(data = u[is_tp == TRUE & has_qtn == FALSE], aes(colour = col), size = 1.9, shape = 16) +
  geom_point(data = u[is_tp == TRUE & has_qtn == TRUE],  aes(fill = col), size = 3.1, shape = 24,
             colour = "black", stroke = .45) +
  scale_colour_identity() + scale_fill_identity() +
  annotate("text", x = max(u$gpos)*.995, y = -log10(cut_gw) + .3, hjust = 1, size = 2.6,
           colour = "#1F3F51", label = sprintf("genome-wide BH over %s clusters", format(nrow(u), big.mark=","))) +
  labs(x = "genome position (Mb, 20 chromosome arms concatenated)",
       y = expression(-log[10](p)~"(EMMAX)"),
       title = sprintf("True-positive stage-2 clusters, one colour each  -  %s %s env%d", CELL, TAG, ENV),
       subtitle = sprintf(paste0("%d of %s clusters are true positives: they tag a driving QTN (r2 >= %.2f within %.0f kb) ",
                          "AND are the best cluster for it after dedup. Grey: everything else. Triangles carry a causal variant."),
                          sum(u$is_tp), format(nrow(u), big.mark=","), parts[[1]]$r2min, parts[[1]]$dmax/1000)) +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "none",
        plot.subtitle = element_text(colour = "grey30", size = 7.6))
fn <- file.path(OUT, sprintf("manhattan_tp_clusters_%s_%s_env%d.png", CELL, TAG, ENV))
ggsave(fn, p, width = 11, height = 5.2, dpi = 190)
cat(sprintf("  written: %s\n", fn))
