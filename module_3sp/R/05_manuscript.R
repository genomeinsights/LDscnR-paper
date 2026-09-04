## =============================================================================
## module_3sp/R/05_manuscript.R
##
## Deliver today's EMMAX + LFMM scan to the manuscript: figure4's slot (its
## caption -- aligned Manhattans, consistent region colours, rotation and
## permutation results -- is a two-panel EMMAX/LFMM figure), plus a values
## file for the macros this module can honestly source.
##
## RENAMED from 06_manuscript.R (PK), wired to figure_manhattan.R's two-panel
## figure -- the design that survived review (Eda-style arrows for Eda and the
## chr1 inversion, two-row EcoPeak rug, span-normalised rug height across
## panels) -- superseding both 07_figure_emmax.R's single EMMAX panel and
## 05_manhattan.R's joint-colour engine comparison, which this delivered
## before. Both of those scripts are now DELETED, not just unused: PK's call,
## since they were different versions of the same figure and keeping dead
## versions around invites someone shipping the wrong one. The StickJoint*
## macros 05_manhattan.R used to source are dropped for the same reason --
## checked first that no manuscript prose cites them (only their own
## definitions in values_panel_module3sp.tex did).
##
## [!] PARTIAL, DELIBERATELY. values_panel.tex (R_3sp_blocks-generated, the
## OLDER pipeline) still defines macros this module has NOT rerun: the
## GCTA-vs-centred estimator comparison, the kinship-basis sweep, the
## decay-vs-map sweep. Overwriting those with a guess, or silently leaving
## them stale under a filename that claims to be current, would both be worse
## than what this does instead -- write ONLY the macros sourced from today's
## 02/03/04 outputs to a SEPARATE file, loaded AFTER values_panel.tex so it
## wins for exactly those macros and nothing else. Same layering values.tex
## already uses for its own hand-set descriptors.
## =============================================================================
suppressMessages({library(data.table); library(LDscnR)})
devtools::load_all("~/gitlab/LDscnR")
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_3sp"), "R", "00_config.R"))
STAGE <- "05_manuscript"
say("=== %s ===\n\n", STAGE)

sc <- readRDS(file.path(PATHS$out, "03_EMMAX", "scan.rds"))
lf <- readRDS(file.path(PATHS$out, "04_lfmm", "lfmm.rds"))
b  <- readRDS(file.path(PATHS$out, "02_bundle", "bundle.rds"))

OUT_DIR <- file.path(PATHS$out, STAGE); dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## ---- figure ------------------------------------------------------------------
## Both source and destination are in PATHS$figures (PK): every module_3sp figure lands in
## one directory, browsed from there to pick what moves to LDscnR_manuscript/figures/ later.
invisible(file.copy(file.path(PATHS$figures, "both_engines_manhattan.pdf"),
         file.path(PATHS$figures, "figure4_stickleback_validation.pdf"), overwrite = TRUE))
say("[1] figures/figure4_stickleback_validation.pdf <- figures/both_engines_manhattan.pdf\n")

## ---- values (module_3sp-sourced subset only) ---------------------------------
con <- sc$consensus; sim <- sc$simes
V <- list()
add <- function(name, value, src, note = "") V[[length(V)+1]] <<-
  list(name = name, value = as.character(value), src = src, note = note)

add("StickTestedUnits", format(nrow(con$test$units), big.mark=","), "03_EMMAX.R")
add("StickSignificantClusters", sum(con$test$units$significant), "03_EMMAX.R", "consensus, EMMAX")
add("StickSimesClusters", sum(sim$test$units$significant), "03_EMMAX.R", "Simes, EMMAX")
sig_con <- con$test$units[significant==TRUE]$unit_id; sig_sim <- sim$test$units[significant==TRUE]$unit_id
add("StickJaccardStatistic", sprintf("%.3f", length(intersect(sig_con,sig_sim))/length(union(sig_con,sig_sim))),
    "03_EMMAX.R", "consensus vs Simes, EMMAX, same units")
add("StickJaccardEngine", sprintf("%.3f",
    { sig_lfmm <- lf$test$units[significant==TRUE]$unit_id
      length(intersect(sig_sim,sig_lfmm))/length(union(sig_sim,sig_lfmm)) }),
    "04_lfmm.R", "EMMAX-Simes vs LFMM-Simes, same statistic")
add("StickReportedRegions", nrow(con$test$regions), "03_EMMAX.R", "consensus, EMMAX")
add("StickEcoPeakRegions", con$rotation$observed, "03_EMMAX.R")
add("StickEcoPeakPct", sprintf("%.1f", 100*con$rotation$observed/nrow(con$test$regions)), "03_EMMAX.R")
add("StickMedianOccupancy", sprintf("%.3f", median(con$test$regions$occupancy)), "03_EMMAX.R")
add("EcoPeakRotationP", sprintf("%.4f", con$rotation$p), "03_EMMAX.R", "consensus, EMMAX")
add("StickConsensusRotFold", sprintf("%.2f", con$rotation$fold), "03_EMMAX.R", "consensus, EMMAX")
add("StickRegionalPermP", sprintf("%.4f", con$null$p), "03_EMMAX.R", "consensus, EMMAX, 1000 surrogates")
## RENAMED from StickRealisedFDR (PK, 2026-09-05): con$null$realised_fdr is
## mean(surrogates)/observed -- an empirical null-to-observed DISCOVERY RATIO, not a
## realised false-discovery proportion (the truth status of observed candidates is
## unknown, so no FDR/FDP can be "realised" from this quantity alone). The LDscnR
## package field name (ld_outlier_perm.R) is unchanged -- renaming that is a separate,
## larger API decision -- but nothing manuscript-facing should call this an FDR.
add("StickConsensusNullObsRatio", sprintf("%.1f", 100*con$null$realised_fdr), "03_EMMAX.R", "consensus, EMMAX")
add("StickSimesNullObsRatio", sprintf("%.1f", 100*sim$null$realised_fdr), "03_EMMAX.R", "Simes, EMMAX")

## marker-wise EMMAX BH count: 03_EMMAX.R does not save its marker-level p-vector (it is
## used once to build the Simes test then discarded), so recompute it here -- cheap
## (~10s for 790,578 markers with the fitted GRM already in the bundle).
p_marker_emmax <- emmax_fast(emmax_setup(b$GTs, b$GRM), b$eco)
add("StickMarkerWiseEMMAX", sum(p.adjust(p_marker_emmax, method="BH") < 0.05, na.rm=TRUE), "03_EMMAX.R (recomputed)",
    "marker-wise BH, q<0.05")
add("StickMarkerWiseLFMM", sum(p.adjust(lf$lfmm_p, method="BH") < 0.05, na.rm=TRUE), "04_lfmm.R", "marker-wise BH, q<0.05")

## EMMAX Simes: region/EcoPeak counts were not previously emitted (only the significant-
## unit count, StickSimesClusters, was).
add("StickSimesRegions", nrow(sim$test$regions), "03_EMMAX.R", "Simes, EMMAX")
add("StickSimesEcoPeakRegions", sim$rotation$observed, "03_EMMAX.R", "Simes, EMMAX")
add("StickSimesEcoPeakPct", sprintf("%.1f", 100*sim$rotation$observed/nrow(sim$test$regions)), "03_EMMAX.R", "Simes, EMMAX")
add("StickSimesRotFold", sprintf("%.2f", sim$rotation$fold), "03_EMMAX.R", "Simes, EMMAX")
add("StickSimesRotP", sprintf("%.4f", sim$rotation$p), "03_EMMAX.R", "Simes, EMMAX")

## Stage-1 cluster inventory (02_bundle.R output, not 03_EMMAX.R -- listed here anyway so
## every module_3sp-sourced macro stays in this one file).
cl_sizes <- b$stage1$clusters$n_snps
add("StickAllClusters", format(length(cl_sizes), big.mark=","), "02_bundle.R")
add("StickSingletonClusters", format(sum(cl_sizes==1), big.mark=","), "02_bundle.R")
add("StickSingletonPct", sprintf("%.1f", 100*mean(cl_sizes==1)), "02_bundle.R")
add("StickMultiClusters", format(sum(cl_sizes>=2), big.mark=","), "02_bundle.R", "clusters with >=2 markers")
add("StickTestedMedianSize", median(cl_sizes[cl_sizes>=SIZE_FLOOR]), "02_bundle.R", "median size among tested (>=size floor) units")
add("StickTestedMaxSize", max(cl_sizes[cl_sizes>=SIZE_FLOOR]), "02_bundle.R", "largest tested unit")

## Ecotype split (02_bundle.R phenotype vector).
add("StickFreshwaterN", sum(b$eco==0), "02_bundle.R")
add("StickMarineN", sum(b$eco==1), "02_bundle.R")

## new: LFMM, not present in the old values_panel.tex at all
add("StickLFMMSignificantClusters", sum(lf$test$units$significant), "04_lfmm.R", "Simes only -- LFMM p-values are precomputed")
add("StickLFMMRegions", nrow(lf$test$regions), "04_lfmm.R")
add("StickLFMMEcoPeakRegions", lf$rotation$observed, "04_lfmm.R")
add("StickLFMMRotFold", sprintf("%.2f", lf$rotation$fold), "04_lfmm.R")
add("StickLFMMRotP", sprintf("%.4f", lf$rotation$p), "04_lfmm.R", "not significant -- state this, do not omit it")
## StickJoint* (from 05_manhattan.R, since deleted) dropped, not carried forward: that
## script's joint-colour region merge across engines is no longer computed anywhere in
## this module, and nothing in the manuscript prose cited these four macros (checked --
## only their own definitions in values_panel_module3sp.tex did). If a joint EMMAX/LFMM
## region count is wanted again later, it needs a real replacement, not a stale number.

OUT_TEX <- file.path(OUT_DIR, "values_panel_module3sp.tex")
con_f <- file(OUT_TEX, "w")
writeLines(c(
"% ============================================================================",
"% values_panel_module3sp.tex -- macros sourced from module_3sp (this repo's",
"% pipeline), loaded AFTER values_panel.tex so these override it for exactly the",
"% macros listed here and no others. values_panel.tex's remaining macros --",
"% the GCTA-vs-centred estimator comparison, the kinship-basis sweep, the",
"% decay-vs-map sweep -- are NOT rerun by module_3sp yet and are left as is.",
"%",
"% GENERATED by module_3sp/R/05_manuscript.R. DO NOT EDIT BY HAND.",
"% ============================================================================",""), con_f)
for (v in V) writeLines(sprintf("\\providecommand{\\%s}{}\\renewcommand{\\%s}{%s}%s",
  v$name, v$name, v$value,
  sprintf("  %% %s%s", v$src, if (nzchar(v$note)) paste0(" -- ", v$note) else "")), con_f)
close(con_f)
say("[2] wrote %s -- %d macros\n", OUT_TEX, length(V))

write_receipt(STAGE,
  inputs = c(file.path(PATHS$out, "02_bundle", "_receipt.rds"),
             file.path(PATHS$out, "03_EMMAX", "_receipt.rds"),
             file.path(PATHS$out, "04_lfmm", "_receipt.rds"),
             file.path(PATHS$out, "figure_manhattan", "_receipt.rds")),
  params = list(), outputs = c(file.path(PATHS$figures, "figure4_stickleback_validation.pdf"), OUT_TEX))
say("\n[3] receipt: %s\n", receipt_path(STAGE))
