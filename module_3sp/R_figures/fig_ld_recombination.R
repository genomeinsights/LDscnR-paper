## =============================================================================
## module_3sp/R_figures/fig_ld_recombination.R
##
## FIGURES ONLY (PK: kept out of R/09_ld_recombination.R, which is analysis).
## Loads that script's saved numbers (out/09_ld_recombination/ld_recombination.rds)
## plus the bundle directly, and renders every figure this supplementary material
## needs: figS:ROC, figS:lddecay (marker-resolution overlay), an illustrative
## n_win_decay=100 track (Chr1+Chr4 only, NOT part of the analysis), the Stage 1
## vs Stage 2 diagnostic (Chr1+Chr4, canonical n_win_decay=20), and a reproduction
## of LDscnR_manuscript/figures/figure5_ld_recombination.pdf on module_3sp's own
## canonical fit. All written to PATHS$figures.
##
## Two of these (the n_win=100 refit, the Chr1/Chr4 Stage-2 prune) run real
## computation that exists ONLY to feed a figure -- not analysis in the sense
## R/09_ld_recombination.R's numbers are (nothing here is cited as a result),
## so they live here rather than splitting hairs over where "compute" ends and
## "plot" begins.
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(ggplot2); library(patchwork)})
devtools::load_all("~/gitlab/LDscnR")
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_3sp"), "R", "00_config.R"))
say("=== fig_ld_recombination ===\n\n")
invisible(check_ldscnr())

A <- readRDS(file.path(PATHS$out, "09_ld_recombination", "ld_recombination.rds"))
W <- A$windows; MK <- A$MK; RESULTS <- A$roc; pooled <- A$cor$pooled; bych <- A$cor$bych
st <- A$cor$sign_p; chr_best <- A$chr_best; chr_med <- A$chr_med

b <- readRDS(file.path(PATHS$out, "02_bundle", "bundle.rds"))
map <- b$map

## ---- 1. figS:ROC -------------------------------------------------------------
say("[1] figS:ROC\n")
roc_df <- rbindlist(lapply(RESULTS, function(r) rbindlist(list(
  data.table(q = r$q, stat = "a",    spec = r$roc_a$specificities,   sens = r$roc_a$sensitivities,
             auc = as.numeric(pROC::auc(r$roc_a))),
  data.table(q = r$q, stat = "ld_w", spec = r$roc_ldw$specificities, sens = r$roc_ldw$sensitivities,
             auc = as.numeric(pROC::auc(r$roc_ldw)))
))))
roc_df[, panel := sprintf("bottom %.0f%% (AUC: a=%.2f, ld[w]=%.2f)",
                          100*q, auc[stat=="a"][1], auc[stat=="ld_w"][1]), by = q]
roc_df[, panel := factor(panel, levels = unique(panel[order(q)]))]

p_roc <- ggplot(roc_df, aes(1 - spec, sens, colour = stat)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey70") +
  geom_path(linewidth = 0.8) +
  scale_colour_manual(values = c(a = "#1b7837", ld_w = "#762a83"),
                      labels = c(a = "decay rate a", ld_w = expression(ld[w])),
                      name = NULL) +
  facet_wrap(~panel, nrow = 1) +
  coord_equal() +
  labs(x = "1 - specificity", y = "sensitivity") +
  theme_bw(12) + theme(panel.grid.minor = element_blank(), legend.position = "top")

OUT_ROC <- file.path(PATHS$figures, "figS_roc_low_recombination.pdf")
ggsave(OUT_ROC, p_roc, width = 10, height = 4.0, device = cairo_pdf)
say("    wrote %s\n", OUT_ROC)

## ---- 2. figS:lddecay ---------------------------------------------------------
say("\n[2] figS:lddecay -- %s (rho=%.2f, strongest), %s (rho=%.2f, near-median)\n",
    chr_best, bych[Chr==chr_best]$rho, chr_med, bych[Chr==chr_med]$rho)

## PK: overlay the two window-level tracks ON TOP of the marker-level ld_w scatter, in
## different colours, one panel per chromosome -- not three stacked facet rows. `a` and the
## pedigree rate live on genuinely different scales from ld_w (decay-rate units ~1e-3-1e-2 and
## cM/Mb ~0-10, against ld_w's r^2 range ~0-0.8), so a shared y-axis would flatten both lines
## to nothing; each window-level track is min-max rescaled onto ld_w's OWN range for THIS
## chromosome so the three are visually comparable in shape/phase, not in absolute value.
## Said explicitly in the y-axis label rather than left implicit.
rescale_to <- function(x, target_max, target_min = 0)
  target_min + (x - min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)) * (target_max - target_min)

mk_track <- function(ch) {
  mk <- MK[Chr == ch]
  wc <- W[Chr == ch][order(start)]
  wc[, mid := (start+end)/2]
  ld_max <- max(mk$ld_w, na.rm = TRUE)
  L <- rbindlist(list(
    data.table(Pos = mk$Pos, y = mk$ld_w,                              track = "ld_w (per marker)"),
    data.table(Pos = wc$mid, y = rescale_to(wc$a,    ld_max),           track = "decay rate a (per window, rescaled)"),
    data.table(Pos = wc$mid, y = rescale_to(wc$rate, ld_max),           track = "pedigree rate (per window, rescaled)")
  ))
  L[, track := factor(track, levels = c("ld_w (per marker)", "decay rate a (per window, rescaled)",
                                        "pedigree rate (per window, rescaled)"))]
  COL <- c("ld_w (per marker)" = "#762a83", "decay rate a (per window, rescaled)" = "#1b7837",
          "pedigree rate (per window, rescaled)" = "#d95f02")
  ggplot(L, aes(Pos/1e6, y, colour = track)) +
    geom_point(data = L[track == "ld_w (per marker)"], size = 0.25, alpha = 0.35) +
    geom_line(data = L[track != "ld_w (per marker)"], linewidth = 0.7) +
    scale_colour_manual(values = COL, name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1, linetype = c(NA,1,1),
                                                      shape = c(16,NA,NA)))) +
    labs(x = "Mb", y = "ld_w", title = ch) +
    theme_bw(11) + theme(panel.grid.minor = element_blank(), legend.position = "top")
}
p_tracks <- (mk_track(chr_best) | mk_track(chr_med)) +
  plot_layout(guides = "collect") & theme(legend.position = "top")

OUT_LDW <- file.path(PATHS$figures, "figS_lddecay_tracks.pdf")
ggsave(OUT_LDW, p_tracks, width = 11, height = 4.2, device = cairo_pdf)
say("    wrote %s\n", OUT_LDW)
say("    suggested caption note: the y-axis is ld_w's own r^2 scale; the decay-rate `a`\n")
say("    and pedigree-rate tracks are min-max rescaled onto that same axis per chromosome\n")
say("    so their SHAPE is comparable to ld_w's, not their absolute value.\n")

## ---- 3. illustrative decay track at n_win=100, Chr1 & Chr4 ONLY (NOT the analysis) -----
## Ported from LDscnR's own vignette (rendered/LDscnR_complexity_reduction.pdf): re-fitting
## decay at a finer window count only changes the TRACK'S SMOOTHNESS, not what it says --
## the vignette's own point, reproduced here on the real panel rather than its small bundled
## demo subset. PK: run this SEPARATELY from the canonical n_win=20 fit everything else in
## this pipeline uses (never touches b$LD_decay), and restrict it to two chromosomes -- a
## genome-wide refit at n_win=100 is the same multi-hour cost as 02_bundle.R's canonical fit,
## for a figure that exists to illustrate one point, not to feed any analysis.
say("\n[3] illustrative track figure: separate compute_LD_decay(n_win_decay=100), Chr1+Chr4 only\n")
## ld_w itself is NOT recomputed here -- matching the vignette exactly (its track_panel()
## reuses the ONE ld_w computed at the top from the keep_el=TRUE fit, and only the `a` track
## varies with n_win_decay). compute_ld_w() needs pairwise edges, which a keep_el=FALSE decay
## object does not carry ("No edge list" -- the exact failure 02_bundle.R's own comments
## already document); ld_w is also the real, canonical, already-computed value everything
## else in this pipeline uses (MK, from the analysis script), so reusing it is correct.
TRACK_CHRS <- c("Chr1", "Chr4")
keep_tc <- map$Chr %in% TRACK_CHRS
gts_tc  <- b$GTs[, keep_tc, drop = FALSE]
map_tc  <- copy(map[keep_tc])
gds_tc_path <- file.path(PATHS$cache, "3sp_track_chr1_chr4.gds")
if (file.exists(gds_tc_path)) unlink(gds_tc_path)
gds_tc <- create_gds_from_geno(geno = gts_tc, map = map_tc, gds_tc_path)
on.exit(try(SNPRelate::snpgdsClose(gds_tc), silent = TRUE), add = TRUE)

set.seed(SEEDS[["bundle"]])
ld_decay_100 <- compute_LD_decay(gds_tc, keep_el = FALSE, slide = DECAY_ARGS$slide,
                                 ld_method = DECAY_ARGS$ld_method, n_win_decay = 100, cores = 1)
say("    %s markers, %s chromosomes, decay refit at n_win_decay = 100 (ld_w unchanged, canonical)\n",
    format(nrow(map_tc), big.mark=","), uniqueN(map_tc$Chr))

sc01 <- function(v) { v <- as.numeric(v); (v - min(v, na.rm=TRUE)) / diff(range(v, na.rm=TRUE)) }
track_panel <- function(ch) {
  lw <- MK[Chr == ch & is.finite(ld_w)]
  lw <- data.table(x = lw$Pos/1e6, y = sc01(lw$ld_w))
  d <- as.data.table(ld_decay_100$by_chr[[ch]]$decay)[regime == "structured" & is.finite(a)]
  dw <- data.table(x = ((d$start+d$end)/2)/1e6, y = sc01(d$a))[order(x)]
  ggplot() +
    geom_point(data = lw, aes(x, y), colour = "grey55", size = 0.25, alpha = 0.35) +
    geom_line(data = dw, aes(x, y), colour = "#D55E00", linewidth = 0.5) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    labs(x = "position (Mbp)", y = "scaled (0-1)",
        title = sprintf("%s  (%d a-windows, n_win_decay = 100)", ch, nrow(dw))) +
    theme_classic(base_size = 9) + theme(plot.title = element_text(size = 9))
}
p_track100 <- (track_panel("Chr1") / track_panel("Chr4")) + plot_annotation(tag_levels = "a")
OUT_TRACK100 <- file.path(PATHS$figures, "figS_lddecay_nwin100_illustrative.pdf")
ggsave(OUT_TRACK100, p_track100, width = 8, height = 6, device = cairo_pdf)
say("    wrote %s\n", OUT_TRACK100)

## ---- 4. Stage 1 vs Stage 2 (CANONICAL n_win=20, the real analysis), Chr1 & Chr4 ------
## Uses the bundle's ACTUAL stage1 (n_win_decay=20, CR_RHO) untouched -- plot_pruning_
## comparison() filters pruned_stage1 by chr internally, so passing the whole-genome object
## is correct, not wasteful. `result` (Stage 2) is computed only for Chr1/Chr4's own clusters
## -- same real config LDW_FLAG/SCORE_THRESHOLD/DISTANCE_THRESHOLD/rho as ld_outlier_test()'s
## own stage2_discovered branch uses, WITH ONE DELIBERATE DIFFERENCE: min_n_loci_flag is left
## at its default (Inf), not forced to 1. ld_outlier_test() sets min_n_loci_flag=1 safely only
## because it is always called on an already-tiny, pre-filtered set (just the significant
## outlier clusters, tens to low hundreds). Copying that onto a WHOLE CHROMOSOME'S clusters
## (42,321 on Chr1 alone) flagged literally every cluster into ONE merge run and hit R's
## vector memory limit (36 GB) -- caught by actually running it, not assumed safe by analogy.
## With min_n_loci_flag at its default, flagging is driven by ld_w_threshold alone, which is
## what "elevated local LD" is supposed to mean at this scale.
say("\n[4] Stage 1 vs Stage 2 (canonical n_win_decay=20), Chr1 & Chr4\n")
cl_pc <- as.data.table(b$stage1$clusters)[Chr %in% TRACK_CHRS]
mk_pc <- unlist(cl_pc$members, use.names = FALSE)
ms_pc <- as.data.table(b$stage1$map_snp)[marker %chin% mk_pc]
sub_pc <- structure(list(map_snp = ms_pc, clusters = cl_pc, pruned = cl_pc$core_snp),
                    class = "ld_complexity_reduction")
result_pc <- ld_prune_and_eMLG(
  GTs = b$GTs[, mk_pc, drop = FALSE], stage1 = sub_pc, ld_w_col = "ld_w_095",
  ld_w_threshold = LDW_FLAG, LD_decay = b$LD_decay, rho = 0.95,
  score_threshold = SCORE_THRESHOLD, min_r2_rho = b$stage1$params$rho,
  distance_threshold = DISTANCE_THRESHOLD, compute_unflagged_eMLG = FALSE,
  min_n_loci_eMLG = 1, cores = 1)

for (ch in TRACK_CHRS) {
  p_pc <- plot_pruning_comparison(chr = ch, pruned_stage1 = b$stage1, result = result_pc,
                                  map = map, out_folder = PATHS$figures)
  say("    wrote %s/%s_stage1_vs_combined_high.png\n", PATHS$figures, ch)
}

## ---- 5. figure5 reproduction: local decay rate `a` vs the pedigree map -----------------
## Reproduces LDscnR_manuscript/figures/figure5_ld_recombination.pdf, but from module_3sp's
## OWN canonical n_win_decay=20 fit (`W`, from the analysis script) rather than R_3sp_blocks'
## n_win=5 "canonical" -- that script's own canonical differs from this pipeline's (00_
## config.R fixed n_win_decay=20), so this is not merely a re-plot of the same numbers.
## Genuine cross-check: the pooled Spearman rho below IS the analysis script's own `pooled`/
## `bych` values, not recomputed -- one set of numbers, shown two ways (scatter here, ROC above).
say("\n[5] figure5 reproduction: local decay rate a vs pedigree map (n_win_decay=20)\n")
## PK: the two panel titles ran together illegibly at theme_bw()'s default title size when
## combined side by side via patchwork (each panel is only ~half the figure width) -- shrunk
## and wrapped onto two lines rather than shortened to the point of losing the numbers.
p_fig5a <- ggplot(W, aes(rate, a)) +
  geom_point(alpha = 0.35, size = 1, colour = "grey30") +
  geom_smooth(method = "lm", se = FALSE, colour = "#1565C0", linewidth = 0.6) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "pedigree map rate (cM/Mb, log scale)", y = "local LD-decay rate a (log scale)",
      title = sprintf("Canonical fit, n_win_decay = 20\npooled Spearman %+.3f over %s windows",
                      pooled, format(nrow(W), big.mark=","))) +
  theme_bw(11) + theme(plot.title = element_text(size = 9))

bych_plot <- copy(bych)[order(-rho)]
bych_plot[, Chr_lab := factor(gsub("Chr","",Chr), levels = gsub("Chr","",Chr))]
p_fig5b <- ggplot(bych_plot, aes(Chr_lab, rho)) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_point(colour = "#1565C0", size = 2) +
  labs(x = "chromosome (ordered by within-chromosome rho)", y = "within-chromosome Spearman",
      title = sprintf("Positive on %d/%d chromosomes at n_win_decay = 20\nsign p = %.2g",
                      sum(bych$rho > 0), nrow(bych), st)) +
  theme_bw(11) + theme(plot.title = element_text(size = 9))

p_fig5 <- p_fig5a | p_fig5b
OUT_FIG5 <- file.path(PATHS$figures, "figure5_ld_recombination_reproduced.pdf")
ggsave(OUT_FIG5, p_fig5, width = 11, height = 4.2, device = cairo_pdf)
say("    wrote %s\n", OUT_FIG5)
say("    NOTE: no multi-n_win (5/10/20/50) sweep here, unlike the old figure's right panel --\n")
say("    module_3sp only ever fits the canonical n_win_decay=20; a genome-wide refit at three\n")
say("    more window counts was not requested and would be the same multi-hour cost as\n")
say("    02_bundle.R's own decay fit, three times over. Flagged rather than silently matched.\n")

say("\n[6] all figures written to %s\n", PATHS$figures)
