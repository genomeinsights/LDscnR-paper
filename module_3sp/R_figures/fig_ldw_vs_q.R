## =============================================================================
## module_3sp/R_figures/fig_ldw_vs_q.R
##
## PK: plot marker-level ld_w against -log10(q), coloured by outlier region, for
## Chr1 and Chr4, using the SAME region colours as figure_manhattan.R's EMMAX
## panel -- not a fresh palette, the LITERAL SAME region-to-colour assignment,
## so a region that reads e.g. orange in the Manhattan figure reads orange here
## too. That assignment depends on processing order across ALL 20 chromosomes
## (see figure_manhattan.R's header: a rotating colour pointer `cur` persists
## ACROSS chromosomes while `used` resets per chromosome), so it is recomputed
## here over the FULL region table exactly as figure_manhattan.R does, then
## subset down to Chr1/Chr4 -- computing it only on two chromosomes would give
## the same regions different colours.
##
## PK follow-up: also produce an ALL-CHROMOSOMES-POOLED version -- same data,
## same colours, no per-chromosome faceting, all 39 regions' points on one
## panel. Rasterised (ggrastr), same reason as figure_manhattan.R's point
## layer: ~790k points is too many for a legible/portable vector PDF.
##
## EMMAX consensus only (not LFMM): ld_w and Stage-1/Stage-2 complexity
## reduction are properties of the EMMAX-side pipeline throughout this module
## (00_config.R, 09_ld_recombination.R, figureS5/S6's Chr1/Chr4 diagnostics);
## "outlier region" with no engine specified means the primary consensus
## regions used everywhere else, not a second, LFMM-specific colour scheme.
##
## AXES: ld_w on x (the predictor -- local LD support, independent of any
## test), -log10(q) on y (the outcome -- EMMAX consensus single-marker
## significance, same q as figure_manhattan.R's point cloud: BH over the
## marker-wise EMMAX scan, not the cluster-level q used for testing). Markers
## outside any significant region are grey, matching the Manhattan figure's
## background-point convention (not separated by odd/even chromosome here,
## since there is only ever one chromosome's background per panel).
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(ggplot2); library(patchwork); library(ggrastr)})
devtools::load_all("~/gitlab/LDscnR")
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_3sp"), "R", "00_config.R"))
STAGE <- "fig_ldw_vs_q"
say("=== %s ===\n\n", STAGE)

b  <- readRDS(file.path(PATHS$out, "02_bundle", "bundle.rds"))
sc <- readRDS(file.path(PATHS$out, "03_EMMAX", "scan.rds"))
map <- b$map
CHRS <- c("Chr1", "Chr4")

## ---- 1. marker-level q, EXACTLY figure_manhattan.R's point-cloud quantity -------------
say("[1] marker-wise EMMAX scan (BH q), for the point cloud\n")
pm_emmax <- emmax_fast(emmax_setup(b$GTs, b$GRM), b$eco)
qv <- p.adjust(pm_emmax, "BH")

## ---- 2. region colours, recomputed EXACTLY as figure_manhattan.R's mk_panel() ---------
## (same PAL source, same per-chromosome "used", same cross-chromosome rotating `cur`) so
## the mapping from region -> colour is identical to what is already printed in the
## Manhattan figure, not merely drawn from the same palette.
say("[2] region colours (EMMAX consensus, all 20 chromosomes, for a matching assignment)\n")
R <- copy(sc$consensus$test$regions)
PAL <- LDscnR:::default_cluster_colours(); R[, col := NA_character_]; cur <- 0L
for (ch in unique(R$Chr)) { used <- character(0)
  for (i in R[Chr == ch, which = TRUE]) {
    repeat { cur <- cur %% length(PAL) + 1L; if (!(PAL[cur] %in% used)) break }
    R$col[i] <- PAL[cur]; used <- c(used, PAL[cur]) } }
stopifnot(R[, .(ok = uniqueN(col) == .N), by = Chr][, all(ok)])
R[, rid := .I]

## ---- 3. per-marker table, GENOME-WIDE (subset to Chr1/Chr4 for that figure below) -----
say("[3] assembling per-marker ld_w/q/region, genome-wide\n")
M <- data.table(marker = map$marker, Chr = map$Chr, pos = map$Pos,
               ld_w = b$ld_ws[map$marker, "rho_0.95"], y = -log10(pmax(qv, .Machine$double.xmin)))
setkey(R, Chr, from, to)
M[, rid := foverlaps(M[, .(Chr, from = pos, to = pos)], R,
      by.x = c("Chr","from","to"), type = "within", mult = "first", nomatch = NA)$rid]
M <- M[is.finite(ld_w) & is.finite(y)]
CV <- setNames(R$col, as.character(R$rid))
M[, cc := ifelse(!is.na(rid), as.character(rid), "..bg")]
M[, ord := !is.na(rid)]; setorder(M, ord)
say("    %s markers genome-wide ; %d in a significant region\n", format(nrow(M), big.mark=","), sum(M$ord))

## ---- 4a. Chr1 + Chr4, one panel each ---------------------------------------------------
mk_panel <- function(ch, Msrc = M, rasterised = FALSE) {
  Mc <- Msrc[Chr == ch]
  CVc <- c(CV[intersect(names(CV), unique(Mc$rid))], "..bg" = "grey75")
  n_reg <- R[Chr == ch, .N]
  pt <- geom_point(aes(ld_w, y, colour = cc, size = ord, alpha = ord))
  if (rasterised) pt <- rasterise(pt, dpi = 200)
  ggplot(Mc) + pt +
    scale_colour_manual(values = CVc, guide = "none") +
    scale_size_manual(values = c("TRUE" = 1.3, "FALSE" = 0.5), guide = "none") +
    scale_alpha_manual(values = c("TRUE" = 0.9, "FALSE" = 0.35), guide = "none") +
    geom_hline(yintercept = -log10(ALPHA), linetype = "dashed", linewidth = 0.4, colour = "grey20") +
    labs(x = expression(ld[w]), y = expression(-log[10](q)*"  (single marker)"),
        title = sprintf("%s (%d regions)", ch, n_reg)) +
    theme_bw(12) + theme(panel.grid = element_blank())
}
FIG_CHR <- mk_panel("Chr1") | mk_panel("Chr4")

OUT_DIR <- file.path(PATHS$out, STAGE); dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_PDF_CHR <- file.path(PATHS$figures, "ldw_vs_q_chr1_chr4.pdf")
ggsave(OUT_PDF_CHR, FIG_CHR, width = 11, height = 4.6, device = cairo_pdf)
ggsave(file.path(PATHS$figures, "ldw_vs_q_chr1_chr4.png"), FIG_CHR, width = 11, height = 4.6, dpi = 200)
say("\n[4] wrote %s\n", OUT_PDF_CHR)

## ---- 4b. all chromosomes pooled, one panel -----------------------------------------------
say("[5] all-chromosomes-pooled panel\n")
CV_ALL <- c(CV, "..bg" = "grey75")
FIG_ALL <- ggplot(M) +
  rasterise(geom_point(aes(ld_w, y, colour = cc, size = ord, alpha = ord)), dpi = 200) +
  scale_colour_manual(values = CV_ALL, guide = "none") +
  scale_size_manual(values = c("TRUE" = 1.0, "FALSE" = 0.3), guide = "none") +
  scale_alpha_manual(values = c("TRUE" = 0.85, "FALSE" = 0.25), guide = "none") +
  geom_hline(yintercept = -log10(ALPHA), linetype = "dashed", linewidth = 0.4, colour = "grey20") +
  labs(x = expression(ld[w]), y = expression(-log[10](q)*"  (single marker)"),
      title = sprintf("All chromosomes (%d regions)", nrow(R))) +
  theme_bw(12) + theme(panel.grid = element_blank())
OUT_PDF_ALL <- file.path(PATHS$figures, "ldw_vs_q_genomewide.pdf")
ggsave(OUT_PDF_ALL, FIG_ALL, width = 7.5, height = 5.5, device = cairo_pdf)
ggsave(file.path(PATHS$figures, "ldw_vs_q_genomewide.png"), FIG_ALL, width = 7.5, height = 5.5, dpi = 200)
say("[6] wrote %s\n", OUT_PDF_ALL)

## ---- 4c. Manhattan OF ld_w (physical position on x), Chr1 + Chr4 -----------------------
## Same M/R/CV as above, same region colours -- just position on x and ld_w on y instead of
## ld_w on x and -log10(q) on y. No significance threshold line: ld_w has no such threshold.
mk_ldw_manhattan <- function(ch) {
  Mc <- M[Chr == ch]
  CVc <- c(CV[intersect(names(CV), unique(Mc$rid))], "..bg" = "grey75")
  n_reg <- R[Chr == ch, .N]
  ggplot(Mc) +
    geom_point(aes(pos / 1e6, ld_w, colour = cc, size = ord, alpha = ord)) +
    scale_colour_manual(values = CVc, guide = "none") +
    scale_size_manual(values = c("TRUE" = 1.3, "FALSE" = 0.5), guide = "none") +
    scale_alpha_manual(values = c("TRUE" = 0.9, "FALSE" = 0.35), guide = "none") +
    labs(x = "position (Mb)", y = expression(ld[w]),
        title = sprintf("%s (%d regions)", ch, n_reg)) +
    theme_bw(12) + theme(panel.grid = element_blank())
}
FIG_LDW_MAN <- mk_ldw_manhattan("Chr1") | mk_ldw_manhattan("Chr4")
OUT_PDF_LDWMAN <- file.path(PATHS$figures, "ldw_manhattan_chr1_chr4.pdf")
ggsave(OUT_PDF_LDWMAN, FIG_LDW_MAN, width = 11, height = 4.6, device = cairo_pdf)
ggsave(file.path(PATHS$figures, "ldw_manhattan_chr1_chr4.png"), FIG_LDW_MAN, width = 11, height = 4.6, dpi = 200)
say("[7] wrote %s\n", OUT_PDF_LDWMAN)

write_receipt(STAGE, inputs = c(file.path(PATHS$out, "02_bundle", "_receipt.rds"),
                                file.path(PATHS$out, "03_EMMAX", "_receipt.rds")),
             params = list(chrs = CHRS), outputs = c(OUT_PDF_CHR, OUT_PDF_ALL, OUT_PDF_LDWMAN))
say("\n[8] receipt: %s\n", receipt_path(STAGE))
