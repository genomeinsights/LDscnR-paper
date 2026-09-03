## =============================================================================
## module_3sp/R/05_manhattan.R
##
## THE JOINT EMMAX/LFMM MANHATTAN, reproducing R_3sp/128_joint_manhattan.R's
## exact visual convention -- two stacked panels, one colour per JOINT region
## identical across both, colour rotation never repeated within a chromosome,
## EcoPeaks shaded, a region drawn only in the panel of the engine that found
## it -- but sourced from the new package's saved 03_scan.R/04_lfmm.R output
## rather than the ad hoc IDX-based logic that script used.
##
## PER-ENGINE REGIONS ARE ALREADY DEFINED (stage2_discovered, from 03/04);
## THE UNION STEP HERE IS A SEPARATE, SIMPLER PHYSICAL MERGE, matching what
## 128_joint_manhattan.R did: this is a VISUALISATION grouping so a shared
## locus reads as one colour across panels, not a new scientific claim about
## where the regions are -- the per-engine regions themselves are unchanged.
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(ggplot2)
  library(patchwork); library(ggrastr)})
devtools::load_all("~/gitlab/LDscnR")
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_3sp"), "R", "00_config.R"))
STAGE <- "05_manhattan"
say("=== %s ===\n\n", STAGE)

BUNDLE_PATH <- file.path(PATHS$out, "02_bundle", "bundle.rds")
b   <- readRDS(BUNDLE_PATH)
sc  <- readRDS(file.path(PATHS$out, "03_scan", "scan.rds"))
lf  <- readRDS(file.path(PATHS$out, "04_lfmm", "lfmm.rds"))
map <- b$map

## per-marker p-values for the point clouds. EMMAX's single-marker scan is not
## saved by 03_scan.R (only cluster-level results are), so it is recomputed here
## -- cheap (under a second, per today's own measurements). LFMM's is saved.
Pm <- emmax_setup(b$GTs, b$GRM)
PV <- list(EMMAX = emmax_fast(Pm, b$eco), LFMM = lf$lfmm_p)

## per-engine regions, ALREADY ASSEMBLED (stage2_discovered) -- not recomputed.
RG <- list(EMMAX = sc$consensus$test$regions[, .(Chr, from, to)],
          LFMM  = lf$test$regions[, .(Chr, from, to)])
say("[1] EMMAX %d regions (consensus) ; LFMM %d regions (Simes)\n", nrow(RG$EMMAX), nrow(RG$LFMM))

## ---- joint region set: physical union + merge across engines ---------------
## LDscnR:::.physical_merge(), not a hand-rolled cumsum -- a hand-rolled version here
## originally reproduced a real bug in the package itself (cummax(to) computed globally
## across the sorted table rather than per chromosome, so a chromosome sorting after one
## with a large `to` silently over-merged: 5 real LFMM regions on Chr7 with pairwise gaps
## of 1.0-9.6 Mb collapsed into one 18.9 Mb block). Fixed in the package and reused here
## rather than fixed twice.
J <- LDscnR:::.physical_merge(rbindlist(RG), REGION_GAP_CHECK)
setorder(J, Chr, from); J[, rid := .I]
setkey(J, Chr, from, to)
for (eng in names(RG)) { ov <- foverlaps(RG[[eng]], J, by.x = c("Chr","from","to"),
    type = "any", mult = "first", nomatch = NA)
  J[[paste0("in_", eng)]] <- J$rid %in% ov$rid }
J[, shared := in_EMMAX & in_LFMM]
say("[2] %d JOINT regions: %d EMMAX-only, %d LFMM-only, %d SHARED\n", nrow(J),
    sum(J$in_EMMAX & !J$in_LFMM), sum(!J$in_EMMAX & J$in_LFMM), sum(J$shared))

## ---- colour rotation, never repeated within a chromosome --------------------
PAL <- LDscnR:::default_cluster_colours()
J[, col := NA_character_]; cur <- 0L
for (ch in unique(J$Chr)) { used <- character(0)
  for (i in J[Chr == ch, which = TRUE]) {
    repeat { cur <- cur %% length(PAL) + 1L; if (!(PAL[cur] %in% used)) break }
    J$col[i] <- PAL[cur]; used <- c(used, PAL[cur]) } }
stopifnot(J[, .(ok = uniqueN(col) == .N), by = Chr][, all(ok)])
say("[3] %d of %d palette colours used; no repeat within any chromosome: TRUE\n",
    uniqueN(J$col), length(PAL))

## ---- EcoPeak shading ---------------------------------------------------------
LIFT <- path.expand("~/gitlab/LDscnR-paper/kingman2021/data/liftover")
ROMAN <- c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV",
           "XV","XVI","XVII","XVIII","XIX","XX","XXI")
.rd <- function(f) { x <- fread(file.path(LIFT, f), header=FALSE,
                                col.names=c("chr","start","end","pv"))
  x[, chr_num := match(sub("^chr","",chr), ROMAN)]
  x[!is.na(chr_num), .(Chr = paste0("Chr", chr_num), start, end)] }
KG <- unique(rbindlist(lapply(ECOPEAK_BEDS, .rd)))
setorder(KG, Chr, start)
KG[, grp := cumsum(c(TRUE, Chr[-1] != Chr[-.N] | start[-1] > cummax(end)[-.N])), by = Chr]
SHADE <- KG[, .(from = min(start), to = max(end)), by = .(Chr, grp)]

## ---- plot --------------------------------------------------------------------
M <- data.table(Chr = map$Chr, pos = map$Pos)[, i := .I]
setkey(J, Chr, from, to)
M[, rid := foverlaps(M[, .(Chr, from = pos, to = pos)], J,
      by.x = c("Chr","from","to"), type = "within", mult = "first", nomatch = NA)$rid]
OFF <- M[, .(mx = max(pos)), by = Chr][order(Chr)][, off := cumsum(c(0, head(mx,-1)) + 1e7)][]
M <- merge(M, OFF[, .(Chr, off)], by = "Chr"); M[, gx := pos + off]
SH <- merge(SHADE, OFF[, .(Chr, off)], by = "Chr")[, `:=`(x1 = from + off, x2 = to + off)]
CV <- setNames(J$col, as.character(J$rid))

## ONE combined colour scale per panel, not two stacked via ggnewscale::new_scale_colour().
## The two-scale-per-panel version (grey background on one scale, region colour on a second,
## stacked scale) produced a ggplot2 aesthetics-length error under patchwork's combination of
## the two panels -- every aesthetic length was verified internally consistent at each step
## (data.table diagnostics all matched), so this reads as ggnewscale/patchwork state leaking
## across the two independently-built plots when combined via `/`, not a data bug. A single
## combined categorical colour column sidesteps the two-scales-per-panel pattern entirely and
## gives the identical visual result.
BG <- c("grey80", "grey65")
mk <- function(eng, lab) {
  d <- copy(M)[, y := -log10(pmax(PV[[eng]][i], .Machine$double.xmin))]
  keep <- J$rid[J[[paste0("in_", eng)]]]
  chr_n <- as.integer(gsub("Chr", "", d$Chr))
  d[, cc := ifelse(!is.na(rid) & rid %in% keep, as.character(rid),
                   ifelse(chr_n %% 2 == 0, "..bg0", "..bg1"))]
  CV_all <- c(CV, setNames(BG, c("..bg0", "..bg1")))
  d[, ord := cc %in% names(CV)]   # draw highlighted points on top
  setorder(d, ord)
  ggplot(d) +
    geom_rect(data = SH, aes(xmin = x1, xmax = x2, ymin = -Inf, ymax = Inf),
              fill = "grey60", alpha = 0.22) +
    rasterise(geom_point(aes(gx, y, colour = cc, size = ord, alpha = ord)), dpi = 200) +
    scale_colour_manual(values = CV_all, guide = "none") +
    scale_size_manual(values = c("TRUE" = 1.1, "FALSE" = 0.3), guide = "none") +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.4), guide = "none") +
    scale_x_continuous(breaks = OFF$off + OFF$mx/2, labels = gsub("Chr","",OFF$Chr), expand = c(0.01,0)) +
    labs(x = NULL, y = expression(-log[10](p)), title = lab) +
    theme_bw(12) + theme(strip.background = element_blank(), panel.grid.minor = element_blank())
}
pA <- mk("EMMAX", sprintf("A  EMMAX, consensus test -- %d regions (%d shared with LFMM)",
                          nrow(RG$EMMAX), sum(J$shared)))
pB <- mk("LFMM",  sprintf("B  LFMM, Simes -- %d regions (%d shared with EMMAX)",
                          nrow(RG$LFMM), sum(J$shared))) + labs(x = "chromosome")
FIG <- pA / pB + plot_annotation(
  title = sprintf("Stage-1 clusters (>= %d markers) scored jointly across engines: %d joint regions",
                  SIZE_FLOOR, nrow(J)),
  subtitle = "One colour per joint region, identical across panels, never repeated within a chromosome. Grey shading = published EcoPeaks. A region is drawn only in the panel of the engine that found it.")

dir.create(stage_dir(STAGE), recursive = TRUE, showWarnings = FALSE)
OUT_PDF <- file.path(stage_dir(STAGE), "joint_manhattan.pdf")
OUT_PNG <- file.path(stage_dir(STAGE), "joint_manhattan.png")
ggsave(OUT_PDF, FIG, width = 16, height = 9, device = cairo_pdf)
ggsave(OUT_PNG, FIG, width = 16, height = 9, dpi = 200)
fwrite(J, file.path(stage_dir(STAGE), "joint_regions.csv"))
say("\n[4] wrote %s\n         %s\n", OUT_PDF, OUT_PNG)

write_receipt(STAGE,
  inputs = c(file.path(PATHS$out, "03_scan", "_receipt.rds"),
             file.path(PATHS$out, "04_lfmm", "_receipt.rds")),
  params = list(region_gap_check = REGION_GAP_CHECK, size_floor = SIZE_FLOOR),
  outputs = c(OUT_PDF, OUT_PNG))
say("    receipt: %s\n", receipt_path(STAGE))
