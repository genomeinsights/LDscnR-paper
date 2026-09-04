## =============================================================================
## module_3sp/R/09_ld_recombination.R
##
## ANALYSIS ONLY (PK: figure code lives separately, in R_figures/, so it does not
## clutter this file). This script answers the question and saves the numbers;
## R_figures/fig_ld_recombination.R renders them plus the illustrative material.
##
## SUPPLEMENTARY: does the window-level LD-decay rate `a`, or the per-marker
## local LD support `ld_w`, track a REAL pedigree recombination map -- and can
## either serve as a classifier of low-recombination windows? materials_and_
## methods.tex already states specific AUC numbers (0.91-0.98 for `a`, 0.73-
## 0.86 for `ld_w` at the bottom 5/10/25% of windows) citing Fig.~figS:ROC and
## Fig.~figS:lddecay, but those figures did not exist -- the numbers, it turns
## out, were computed on Formica ant hybrid data (formica_hybrid/module0_ld_
## pruning/), not stickleback (found by cross-repo search, 2026-09-04). PK's
## call: rebuild this on the actual stickleback panel rather than reuse or
## reference the ant analysis, so the supplement matches the dataset the rest
## of the methods describe.
##
## DATA, ALL FROM MODULE_3SP'S OWN BUNDLE, NOT REDERIVED:
##   - per-window decay rate `a`: LD_decay$by_chr[[ch]]$decay (start, end, a),
##     the SAME windowed fits underlying decay_sum's per-chromosome summary --
##     n_win_decay = 20 (00_config.R), matching the n_win=20 sweep point
##     R_3sp_blocks/05_decay_recombination.R already validated on this exact
##     pedigree map (pooled/within-chromosome Spearman rho), so the pooled
##     correlation below is a direct cross-check against that prior result,
##     not a fresh, unverified number.
##   - per-marker ld_w: bundle$ld_ws[, "rho_0.95"], aggregated to the SAME
##     windows by taking the median over markers falling inside each window
##     (foverlaps on physical position) -- this is the marker-resolution
##     statistic 00_config.R and ld_outlier_test() both use downstream.
##   - the pedigree map itself: PATHS$recmap (stickleback_recomb_3crosses_
##     gasAcu1.tsv) -- crossovers counted in three real F2 crosses, no
##     population LD in its own estimation, so it is a non-circular reference
##     (R_3sp_blocks' own header makes this point and it still holds here).
##
## THE CLASSIFICATION QUESTION IS DIFFERENT FROM THE CORRELATION QUESTION.
## R_3sp_blocks/05_decay_recombination.R asks whether `a` correlates with the
## map, genome-wide (Figure 4/5's territory). This script asks a narrower,
## more actionable question: restricted to the LOW-recombination tail (the
## regions complexity reduction most needs to compress), how well do `a` and
## `ld_w` separate low-recombination windows from the rest? That is a ROC/AUC
## question, not a correlation, and only `a` was ever evaluated this way on
## stickleback before this script.
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(pROC)})
devtools::load_all("~/gitlab/LDscnR")
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_3sp"), "R", "00_config.R"))
STAGE <- "09_ld_recombination"
say("=== %s ===\n\n", STAGE)
invisible(check_ldscnr())

b <- readRDS(file.path(PATHS$out, "02_bundle", "bundle.rds"))
map <- b$map

INPUTS <- c(file.path(PATHS$out, "02_bundle", "_receipt.rds"), PATHS$recmap)
PARAMS <- list(q_grid = c(0.05, 0.10, 0.25), ld_w_col = "rho_0.95", sweep_nwin = c(10, 50))
if (!stage_stale(STAGE, INPUTS, PARAMS) && !nzchar(Sys.getenv("FORCE"))) {
  say("\nNothing to do. Set FORCE=1 to rerun anyway.\n"); quit(save = "no")
}
OUT_DIR <- file.path(PATHS$out, STAGE); dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## ---- 1. per-window `a`, from the bundle's own windowed decay fits ----------
say("[1] per-window decay rate `a` (LD_decay$by_chr, n_win_decay = %d)\n", DECAY_ARGS$n_win_decay)
W <- rbindlist(lapply(names(b$LD_decay$by_chr), function(ch) {
  d <- as.data.table(b$LD_decay$by_chr[[ch]]$decay)[, .(start, end, a)]
  d[, Chr := ch][]
}))
W <- W[is.finite(start) & is.finite(end) & is.finite(a) & a > 0]
W[, wid := .I]
say("    %s windows across %d chromosomes ; median width %.0f kb\n",
    format(nrow(W), big.mark=","), uniqueN(W$Chr), median(W$end - W$start)/1e3)

## ---- 2. per-window median ld_w, from the same marker set the panel scans --
say("\n[2] per-window median ld_w (%s)\n", PARAMS$ld_w_col)
MK <- data.table(marker = map$marker, Chr = map$Chr, Pos = map$Pos,
                 ld_w = b$ld_ws[map$marker, PARAMS$ld_w_col])
setkey(MK, Chr, Pos)
setkey(W, Chr, start, end)
JW <- foverlaps(MK[, .(Chr, start = Pos, end = Pos, ld_w)], W,
               by.x = c("Chr","start","end"), type = "any", nomatch = NULL)
LDW <- JW[, .(ld_w_med = median(ld_w), n_markers = .N), by = wid]
W <- merge(W, LDW, by = "wid", all.x = TRUE)
say("    %s / %s windows have >= 1 marker (median %.0f markers/window)\n",
    format(sum(!is.na(W$ld_w_med)), big.mark=","), format(nrow(W), big.mark=","),
    median(W$n_markers, na.rm = TRUE))

## ---- 3. pedigree map rate per window ----------------------------------------
say("\n[3] pedigree map: %s\n", basename(PATHS$recmap))
M0 <- fread(PATHS$recmap)
M0[, Chr := paste0("Chr", as.integer(gsub("[^0-9]", "", as.character(new_chr))))]
M0 <- M0[is.finite(rate_consensus_cMperMb) & is.finite(new_start) & is.finite(new_end)]
say("    %s bins, median width %.0f kb, median rate %.2f cM/Mb\n",
    format(nrow(M0), big.mark=","), median(M0$new_end - M0$new_start)/1e3,
    median(M0$rate_consensus_cMperMb))

MB <- M0[, .(Chr, start = new_start, end = new_end, rate = rate_consensus_cMperMb)]
setkey(MB, Chr, start, end)
setkey(W, Chr, start, end)
JR <- foverlaps(W[, .(Chr, start, end, wid)], MB, by.x = c("Chr","start","end"),
                type = "any", nomatch = NULL)
RW <- JR[, .(rate = median(rate)), by = wid]
W <- merge(W, RW, by = "wid", all.x = TRUE)
W <- W[is.finite(rate) & is.finite(a) & is.finite(ld_w_med)]
say("    %s windows retain a finite rate, a, AND ld_w_med (usable set)\n",
    format(nrow(W), big.mark=","))

## ---- 4. correlation cross-check (against R_3sp_blocks' n_win=20 result) ----
say("\n[4] cross-check: pooled/within-chromosome Spearman rho(a, rate)\n")
pooled <- cor(W$a, W$rate, method = "spearman")
bych <- W[, .(rho = if (.N >= 5) cor(a, rate, method = "spearman") else NA_real_), by = Chr][is.finite(rho)]
st <- binom.test(sum(bych$rho > 0), nrow(bych))$p.value
say("    pooled rho %+.3f | within-chr %d/%d positive | sign p %.3g\n",
    pooled, sum(bych$rho > 0), nrow(bych), st)

## ---- 5. ROC / AUC at the low-recombination tail -----------------------------
say("\n[5] ROC/AUC: `a` and ld_w_med as classifiers of low-recombination windows\n")
roc_one <- function(q) {
  y <- as.integer(W$rate <= quantile(W$rate, q))
  ra <- pROC::roc(y, W$a, quiet = TRUE, direction = "auto")
  rl <- pROC::roc(y, W$ld_w_med, quiet = TRUE, direction = "auto")
  say("    q = %.2f (n_low = %4d/%d) | AUC(a) = %.3f [%.3f-%.3f] | AUC(ld_w) = %.3f [%.3f-%.3f]\n",
      q, sum(y), length(y), pROC::auc(ra), pROC::ci.auc(ra)[1], pROC::ci.auc(ra)[3],
      pROC::auc(rl), pROC::ci.auc(rl)[1], pROC::ci.auc(rl)[3])
  list(q = q, roc_a = ra, roc_ldw = rl)
}
RESULTS <- lapply(PARAMS$q_grid, roc_one)
names(RESULTS) <- as.character(PARAMS$q_grid)

## ---- 6. representative chromosomes for the marker-resolution figure --------
## Picked here (analysis), used by the figures script: the strongest and a near-median
## within-chromosome rho(a, rate), so the figure shows both a clean case and a typical one.
bych_ord <- bych[order(-rho)]
chr_best <- bych_ord$Chr[1]
chr_med  <- bych_ord$Chr[which.min(abs(bych_ord$rho - median(bych_ord$rho)))]
say("\n[6] representative chromosomes: %s (rho=%.2f, strongest), %s (rho=%.2f, near-median)\n",
    chr_best, bych[Chr==chr_best]$rho, chr_med, bych[Chr==chr_med]$rho)

## ---- 7. multi-n_win sensitivity: n_win=10 and n_win=50, genome-wide (PK) ---------------
## Completes figure5's right panel (old R_3sp_blocks version compared within-chromosome
## rho(a, rate) across n_win in {5,10,20,50}). PK: only 10 and 50 needed -- n_win=20 is
## already the canonical fit (section 1/`bych` above), and 5 was not requested.
##
## GENOME-WIDE, unlike the Chr1+Chr4-only illustrative n_win=100 refit in R_figures/: this
## feeds a real reported number (Figure 5), not an illustration, so every chromosome is
## needed, not two. keep_el=FALSE: only the windowed `a` fit is used here, never ld_w or
## clustering, so no edge lists are kept -- a real saving over 02_bundle.R's canonical fit,
## which needs edges afterward for ld_complexity_reduction(). Still genuinely expensive (two
## more genome-wide decay fits) -- 00_config.R's own comment ("the sweep at 10/50 is stage
## 09, not here") always meant this to land in this script, just deferred until asked for.
say("\n[7] multi-n_win sensitivity: genome-wide decay at n_win_decay = 10 and 50\n")
gds_sweep_path <- file.path(PATHS$cache, "3sp_sweep.gds")
if (file.exists(gds_sweep_path)) unlink(gds_sweep_path)
gds_sweep <- create_gds_from_geno(geno = b$GTs, map = map, gds_sweep_path)
on.exit(try(SNPRelate::snpgdsClose(gds_sweep), silent = TRUE), add = TRUE)

sweep_bych <- function(nwin) {
  set.seed(SEEDS[["bundle"]])
  ld <- compute_LD_decay(gds_sweep, keep_el = FALSE, slide = DECAY_ARGS$slide,
                         ld_method = DECAY_ARGS$ld_method, n_win_decay = nwin, cores = 1)
  Wn <- rbindlist(lapply(names(ld$by_chr), function(ch) {
    d <- as.data.table(ld$by_chr[[ch]]$decay)[, .(start, end, a)]
    d[, Chr := ch][]
  }))
  Wn <- Wn[is.finite(start) & is.finite(end) & is.finite(a) & a > 0]
  Wn[, wid := .I]
  setkey(Wn, Chr, start, end)
  Jn <- foverlaps(Wn[, .(Chr, start, end, wid)], MB, by.x = c("Chr","start","end"),
                  type = "any", nomatch = NULL)
  Rn <- Jn[, .(rate = median(rate)), by = wid]
  Wn <- merge(Wn, Rn, by = "wid", all.x = TRUE)
  Wn <- Wn[is.finite(rate) & is.finite(a)]
  bychn <- Wn[, .(rho = if (.N >= 5) cor(a, rate, method = "spearman") else NA_real_), by = Chr][is.finite(rho)]
  say("    n_win_decay = %3d: %s windows ; within-chr %d/%d positive\n",
      nwin, format(nrow(Wn), big.mark=","), sum(bychn$rho > 0), nrow(bychn))
  bychn
}
bych_10 <- sweep_bych(10)
bych_50 <- sweep_bych(50)
BYCH_MULTI <- list(`10` = bych_10, `20` = bych, `50` = bych_50)

## ---- 8. save + receipt ---------------------------------------------------------
## Saves MK too (the per-marker ld_w table) -- the figures script needs it and it is cheap
## (one column of ~790k values), far cheaper than having R_figures/ recompute it from the
## bundle's raw ld_ws matrix.
OUT_RDS <- file.path(OUT_DIR, "ld_recombination.rds")
saveRDS(list(windows = W, MK = MK, roc = RESULTS, cor = list(pooled = pooled, bych = bych, sign_p = st),
            bych_multi = BYCH_MULTI, chr_best = chr_best, chr_med = chr_med), OUT_RDS)
write_receipt(STAGE, inputs = INPUTS, params = PARAMS, outputs = OUT_RDS)
say("\n[8] wrote %s\n    receipt: %s\n", OUT_RDS, receipt_path(STAGE))
say("\n    Figures: R_figures/fig_ld_recombination.R\n")
