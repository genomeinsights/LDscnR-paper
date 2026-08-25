## =====================================================================
## kingman2021 / R/17_overlap_sweep.R
##
## Corrected successor to module_sticklebacks_LDscnR/kingman_overlap_sweep.R. Two fixes:
##
##   (a) REGION-FLOOR MASKING. prop_hit and fold are degenerate when a cell contains only
##       a handful of regions: a single region that happens to hit gives prop_hit = 1.00,
##       and a near-empty rotation null gives fold in the hundreds. In the original sweep
##       this produced a bright band across the HIGH-tau / HIGH-l_min corner -- cells
##       holding exactly ONE region -- which is the opposite of where the real signal is.
##       Cells with n_regions < MIN_REGIONS are now NA in prop_hit / fold / F1, so the
##       F1-optimal cell is chosen only among cells that can support the statistic.
##
##   (b) fold no longer divides by a 1e-6 floor (which manufactured folds up to 1e6 when
##       the null was empty); an empty null gives NA.
##
## B is also raised to match the headline rotation null in R/12 rather than B=200.
##
## Use the Kingman EcoPeaks as an external positive control to calibrate the
## (tau_C x l_min) operating point: sweep the grid, cluster the 3sp scan at each
## cell, and score its overlap with the Kingman global-specific EcoPeaks. Four
## truth-anchored views over the same grid as the structure-null landscape:
##   n_hit    : # regions overlapping an EcoPeak      (raw overlap; recall-ish)
##   prop_hit : n_hit / n_regions                     (precision-ish; size-confounded)
##   fold     : n_hit / expected under a within-chromosome rotation null (honest
##              enrichment, corrects for region number + size)
##   recall   : fraction of the EcoPeaks that are recovered
## The F1 of (prop_hit, recall) locates the sweet spot; the paper's operating point
## (tau_C=0.05, l_min=3) and the F1-optimal cell are marked.
##
## Reads kingman2021's lifted peak BED READ-ONLY; writes only to this module.
## Run from the LDscnR-paper root:
##   Rscript module_sticklebacks_LDscnR/kingman_overlap_sweep.R [EMMAX|LFMM] [c155.specific|c150.specific|c155.sensitive|c150.sensitive]
## =====================================================================
suppressMessages({ library(data.table); library(LDscnR); library(ggplot2); library(patchwork) })

a       <- commandArgs(trailingOnly = TRUE)
METHOD  <- if (length(a) >= 1) toupper(a[1]) else "EMMAX"       # EMMAX | LFMM
PEAKSET <- if (length(a) >= 2) a[2] else "c155.specific"        # c155.specific | c150.specific | c155.sensitive | c150.sensitive
TAG     <- sprintf("%s_%s", METHOD, gsub("[.]", "", PEAKSET))
BND     <- "module_sticklebacks_LDscnR/data/3sp_LDscnR_data.rds"   # read-only
PEAK    <- sprintf("kingman2021/data/liftover/pv_%s.bed", PEAKSET)   # lifted EcoPeaks (gasAcu1), read-only
OUTFIG  <- "kingman2021/figures"; OUTRES <- "kingman2021/data"
for (p in c(OUTFIG, OUTRES)) if (!dir.exists(p)) dir.create(p, recursive = TRUE)
ROMAN  <- c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV",
            "XVI","XVII","XVIII","XIX","XX","XXI")
ZISSOU <- c("#3B9AB2","#78B7C5","#EBCC2A","#E1AF00","#F21A00")   # Wes Anderson (Zissou1)
TAUS   <- seq(0.05, 0.90, by = 0.05); LMINS <- 1:15; B <- 2000L; set.seed(1)
OP_TAU <- 0.05; OP_LMIN <- 3L
MIN_REGIONS <- 5L        # below this, prop_hit / fold / F1 are not estimable -> NA

## ---- data + C-score + edges --------------------------------------------------------
d  <- readRDS(BND); m3 <- as.data.table(d$map)
C  <- if (METHOD == "LFMM") {
  ld_cscore(m3$lfmm_p, d$ld_ws, alpha = 0.05, qstar = seq(0,.95,.05))
} else {
  ld_cscore(emmax_fast(emmax_setup(d$GTs, d$GRM), d$eco), d$ld_ws, alpha = 0.05, qstar = seq(0,.95,.05))
}
names(C) <- m3$marker
edges <- ld_edges(names(C)[C > 0], d$GTs, m3[, .(marker, Chr, Pos)],
                  as.data.table(d$LD_decay$decay_sum), rho_ld = 0.60, dcap = 5e5)
mpos <- stats::setNames(m3$Pos, m3$marker); mchr <- stats::setNames(as.integer(gsub("Chr","",m3$Chr)), m3$marker)
rng  <- m3[, .(lo = min(Pos), hi = max(Pos)), by = .(chr_num = as.integer(gsub("Chr","",Chr)))]; setkey(rng, chr_num)

Pk <- fread(PEAK, header = FALSE, col.names = c("chr","start","end","pv"))
Pk[, chr_num := match(sub("^chr","",chr), ROMAN)]; Pk <- Pk[!is.na(chr_num), .(chr_num,start,end)]
setkey(Pk, chr_num, start, end); n_peaks <- nrow(Pk)
hit_ids <- function(R) { if (!nrow(R)) return(integer(0))
  o <- foverlaps(R[, .(chr_num,start,end,id=.I)], Pk, by.x=c("chr_num","start","end"), type="any", nomatch=NULL); unique(o$id) }
peaks_hit <- function(R) { if (!nrow(R)) return(0L)
  Rk <- R[, .(chr_num,start,end)]; setkey(Rk, chr_num, start, end)
  o <- foverlaps(Pk[, .(chr_num,start,end,pid=.I)], Rk, by.x=c("chr_num","start","end"), type="any", nomatch=NULL); uniqueN(o$pid) }

## ---- sweep: cluster once per tau, filter by l_min ----------------------------------
## REPLOT=1 regenerates the figure from the saved grid (no permutations) -- use after a
## cosmetic change to the panels.
CSV <- file.path(OUTRES, sprintf("overlap_sweep_%s.csv", TAG))
REPLOT <- nzchar(Sys.getenv("REPLOT")) && file.exists(CSV)
grid <- if (REPLOT) fread(CSV) else rbindlist(lapply(TAUS, function(tau) {
  mk <- names(C)[C >= tau]
  regs <- if (length(mk)) ld_regions(mk, edges) else list()
  rbindlist(lapply(LMINS, function(lm) {
    R <- regs[lengths(regs) >= lm]
    if (!length(R)) return(data.table(tau=tau, lmin=lm, n_regions=0L, n_hit=0L, prop_hit=NA_real_,
                                       fold=NA_real_, recall=0))
    RD <- rbindlist(lapply(R, function(x) data.table(chr_num=unname(mchr[x[1]]), start=min(mpos[x]), end=max(mpos[x]))))
    RD[, span := end - start]
    nh <- length(hit_ids(RD)); nr <- nrow(RD)
    r  <- rng[J(RD$chr_num)]; nn <- integer(B)
    for (b in seq_len(B)) { st <- r$lo + floor(runif(nr) * pmax(1, r$hi - r$lo - RD$span))
      nn[b] <- length(hit_ids(data.table(chr_num=RD$chr_num, start=st, end=st+RD$span))) }
    mn <- mean(nn)
    data.table(tau=tau, lmin=lm, n_regions=nr, n_hit=nh, prop_hit=nh/nr,
               fold=if (mn > 0) nh/mn else NA_real_, recall=peaks_hit(RD)/n_peaks) }))
}))
if (REPLOT) cat(sprintf("[%s] REPLOT: reusing %s\n", METHOD, CSV))
## (a) region-floor mask: statistics that divide by n_regions are not estimable below it
grid[, estimable := n_regions >= MIN_REGIONS]
grid[!(estimable), `:=`(prop_hit = NA_real_, fold = NA_real_)]
grid[, F1 := ifelse(estimable & prop_hit>0 & recall>0, 2*prop_hit*recall/(prop_hit+recall), NA_real_)]
best <- grid[which.max(fifelse(is.na(F1), -Inf, F1))]
cat(sprintf("[%s] %d/%d cells estimable at n_regions >= %d\n",
            METHOD, sum(grid$estimable), nrow(grid), MIN_REGIONS))
cat(sprintf("[%s] peaks=%d ; F1-optimal cell: tau=%.2f l_min=%d (n_reg=%d, prop=%.2f, recall=%.2f, fold=%.1f)\n",
            METHOD, n_peaks, best$tau, best$lmin, best$n_regions, best$prop_hit, best$recall, best$fold))
grid[, fold_cap := fold]     # no cap needed once the region floor removes the degenerate cells
fwrite(grid, file.path(OUTRES, sprintf("overlap_sweep_%s.csv", TAG)))

## ---- heatmaps ----------------------------------------------------------------------
mark <- data.table(tau=c(OP_TAU, best$tau), lmin=c(OP_LMIN, best$lmin), what=c("paper op.","F1 max"))
hm <- function(fill, title, sub) {
  ggplot(grid, aes(factor(tau), factor(lmin), fill = .data[[fill]])) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_point(data = mark, aes(factor(tau), factor(lmin), shape = what), inherit.aes = FALSE, size = 2.6, stroke = 1.1) +
    scale_shape_manual(values = c("paper op." = 4, "F1 max" = 1), name = NULL) +
    scale_fill_gradientn(colours = ZISSOU, na.value = "grey92", name = NULL) +
    labs(x = expression(tau[C]), y = expression(l[min]), title = title, subtitle = sub) +
    theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, size = 5.5),
          axis.text.y = element_text(size = 5.5), plot.title = element_text(face = "bold", size = 9),
          plot.subtitle = element_text(size = 7.5), legend.key.height = unit(0.7, "lines"))
}
g <- (hm("n_hit", "Regions overlapping an EcoPeak", "raw overlap (recall-like)") |
      hm("prop_hit", "Proportion of regions overlapping", sprintf("precision-like (grey: n_regions < %d)", MIN_REGIONS)) |
      hm("recall", "Fraction of EcoPeaks recovered", "recall")) /
     (hm("fold_cap", "Enrichment vs rotation null", sprintf("fold vs rotation null (grey: n_regions < %d)", MIN_REGIONS)) |
      hm("F1", "F1 of precision and recall", "the sweet spot") |
      hm("n_regions", "Number of regions", "context")) +
  plot_annotation(title = sprintf("3sp %s regions vs Kingman %s EcoPeaks over the tau_C x l_min grid", METHOD, PEAKSET),
                  subtitle = sprintf("x = paper operating point (tau=%.2f, l_min=%d);  o = F1-optimal cell (tau=%.2f, l_min=%d)",
                                     OP_TAU, OP_LMIN, best$tau, best$lmin))
ggsave(file.path(OUTFIG, sprintf("kingman_sweep_%s.png", TAG)), g, width = 15, height = 8.5, dpi = 160)
cat(sprintf("[%s] wrote %s/kingman_sweep_%s.png\n", METHOD, OUTFIG, TAG))
cat(sprintf("[%s] paper op (tau=%.2f,l_min=%d): F1=%.4f ; rank %d of %d estimable cells\n",
            METHOD, OP_TAU, OP_LMIN, grid[tau==OP_TAU & lmin==OP_LMIN, F1],
            grid[!is.na(F1)][order(-F1), which(tau==OP_TAU & lmin==OP_LMIN)][1], sum(!is.na(grid$F1))))
