## =============================================================================
## module_9sp/R_figures/figure_manhattan.R
##
## Manhattan for the 9sp EMMAX result, following module_3sp/R_figures/
## figure_manhattan.R's template (single-marker -log10(q), regions coloured,
## one colour never repeated within a chromosome, background greys alternating
## by chromosome parity) with two differences:
##
##   - TWO PANELS ARE CONSENSUS vs SIMES, not EMMAX vs LFMM -- no LFMM run
##     exists for 9sp yet (only 03_EMMAX.R so far), so this compares the two
##     STATISTICS within the one engine we have, not two engines.
##   - NO ECOPEAK RUG -- no external-validation reference exists for
##     nine-spined stickleback (PK, 2026-09-07), so there is nothing to mark
##     regions against. Regions are still coloured; there is just no on/off-
##     peak annotation below the axis.
##
## The marker-wise EMMAX p-values are RECOMPUTED here (not saved by
## 03_EMMAX.R), same reason and same cost as 3sp's own figure script: cheap
## (~seconds), and re-deriving beats adding a second copy of a large per-marker
## vector to scan.rds for a figure-only need.
##
## [!] CHROMOSOME ORDERING: sorts by the NUMERIC chromosome index (as.integer
## on the digits after "Chr"), not order(Chr) on the character string --
## lexicographic sort on "Chr1","Chr10","Chr11",...,"Chr2",... would put
## chromosomes in the wrong genomic order.
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(ggplot2)
  library(patchwork); library(ggrastr)})
devtools::load_all("~/gitlab/LDscnR")
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_9sp"), "R", "00_config.R"))
STAGE <- "figure_manhattan"
say("=== %s ===\n\n", STAGE)

b  <- readRDS(file.path(PATHS$out, "02_bundle", "bundle.rds"))
sc <- readRDS(file.path(PATHS$out, "03_EMMAX", "scan.rds"))
map <- b$map; eco_resid <- sc$eco_resid

## ---- per-marker q, EMMAX, on the SAME covariate-residualised phenotype 03_EMMAX.R tested
say("[1] marker-wise EMMAX scan on the lineage-residualised phenotype (BH q)\n")
pm_emmax <- emmax_fast(emmax_setup(b$GTs, b$GRM), eco_resid)
Q <- list(Consensus = NULL, Simes = p.adjust(pm_emmax, "BH"))   # Consensus panel uses its own unit-level q via regions below; point cloud below is marker-wise for BOTH panels (same underlying scan), matching 3sp's design (each panel's point cloud is single-marker, only the coloured/region overlay differs)
Q$Consensus <- Q$Simes   # identical point cloud; the two panels differ in which regions are overlaid, not in the background scan

REG <- list(Consensus = copy(sc$consensus$test$regions), Simes = copy(sc$simes$test$regions))
for (eng in names(REG)) REG[[eng]][, rid := .I]
say("[2] Consensus: %d regions ; Simes: %d regions\n", nrow(REG$Consensus), nrow(REG$Simes))

## numeric chromosome order -- NOT order(Chr) on the character string (see header)
chr_num_of <- function(ch) as.integer(sub("Chr", "", ch))
OFF <- map[, .(mx = max(Pos)), by = Chr]
OFF[, chr_n := chr_num_of(Chr)]
setorder(OFF, chr_n)
OFF[, off := cumsum(c(0, head(mx, -1)) + 2e6)]
BG <- c("grey80", "grey62")

top_expand <- 0.08
mk_panel <- function(eng, lab) {
  R <- REG[[eng]]; qv <- Q[[eng]]
  PAL <- LDscnR:::default_cluster_colours(); R[, col := NA_character_]; cur <- 0L
  for (ch in unique(R$Chr)) { used <- character(0)
    for (i in R[Chr == ch, which = TRUE]) {
      repeat { cur <- cur %% length(PAL) + 1L; if (!(PAL[cur] %in% used)) break }
      R$col[i] <- PAL[cur]; used <- c(used, PAL[cur]) } }
  stopifnot(R[, .(ok = uniqueN(col) == .N), by = Chr][, all(ok)])

  M <- data.table(Chr = map$Chr, pos = map$Pos, y = -log10(pmax(qv, .Machine$double.xmin)))
  setkey(R, Chr, from, to)
  M[, rid := foverlaps(M[, .(Chr, from = pos, to = pos)], R,
        by.x = c("Chr","from","to"), type = "within", mult = "first", nomatch = NA)$rid]
  M <- merge(M, OFF[, .(Chr, off, chr_n)], by = "Chr")[, gx := pos + off]
  CV <- setNames(R$col, as.character(R$rid))
  CV_all <- c(CV, setNames(BG, c("..bg0", "..bg1")))
  M[, cc := ifelse(!is.na(rid), as.character(rid), ifelse(chr_n %% 2 == 0, "..bg0", "..bg1"))]
  M[, ord := !is.na(rid)]; setorder(M, ord)

  ## rug: one tick per region at its genomic midpoint, below the axis -- same purpose as
  ## 3sp's EcoPeak rug (a quick visual census of where regions fall along the genome), but
  ## ONE ROW/no colour split, since there is no external reference to mark on/off of. Sized
  ## as a fraction of the panel's FULL RENDERED SPAN (data range + top/bottom expansion),
  ## not of y_max alone -- see figure_manhattan.R (3sp)'s header for why that distinction
  ## matters when panels have different y_max/top_expand.
  BOTTOM_MULT <- 0.10
  y_max <- max(M$y, na.rm = TRUE)   ## 2 markers have NA q (BH input NA) -> NA y; without
                                     ## na.rm this silently NA's y_max and every rug position
  full_span <- y_max * (1 + top_expand + BOTTOM_MULT)
  RB <- merge(R, OFF[, .(Chr, off)], by = "Chr")[, gx := (from + to)/2 + off]
  rug_h <- 0.018 * full_span
  rug_y <- -0.045 * full_span

  ggplot(M) +
    rasterise(geom_point(aes(gx, y, colour = cc, size = ord, alpha = ord)), dpi = 200) +
    scale_colour_manual(values = CV_all, guide = "none") +
    scale_size_manual(values = c("TRUE" = 1.1, "FALSE" = 0.3), guide = "none") +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.4), guide = "none") +
    geom_hline(yintercept = -log10(ALPHA), linetype = "dashed", linewidth = 0.45, colour = "grey20") +
    geom_segment(data = RB, aes(x = gx, xend = gx, y = rug_y - rug_h, yend = rug_y + rug_h),
                 colour = "grey30", linewidth = 0.6) +
    scale_x_continuous(breaks = OFF$off + OFF$mx/2, labels = OFF$chr_n, expand = c(0.005, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0.14, top_expand))) +
    labs(x = NULL, y = expression(-log[10](q)*"  (single marker)"), title = lab) +
    theme_bw(12) +
    theme(panel.grid = element_blank(), plot.title = element_text(size = 12))
}

pA <- mk_panel("Consensus", sprintf("A  Consensus test -- %d regions", nrow(REG$Consensus)))
pB <- mk_panel("Simes", sprintf("B  Simes -- %d regions", nrow(REG$Simes))) + labs(x = "chromosome")
FIG <- pA / pB

OUT_DIR <- file.path(PATHS$out, STAGE); dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_PDF <- file.path(PATHS$figures, "manhattan_consensus_simes.pdf")
ggsave(OUT_PDF, FIG, width = 15, height = 8.4, device = cairo_pdf)
ggsave(file.path(PATHS$figures, "manhattan_consensus_simes.png"), FIG, width = 15, height = 8.4, dpi = 200)
say("\n[3] wrote %s\n", OUT_PDF)
write_receipt(STAGE, inputs = c(file.path(PATHS$out, "02_bundle", "_receipt.rds"),
                                file.path(PATHS$out, "03_EMMAX", "_receipt.rds")),
             params = list(), outputs = OUT_PDF)
say("    receipt: %s\n", receipt_path(STAGE))
