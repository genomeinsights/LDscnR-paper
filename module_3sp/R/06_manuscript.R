## =============================================================================
## module_3sp/R/06_manuscript.R
##
## Deliver today's LFMM scan and joint Manhattan to the manuscript: figure4's
## slot (its caption -- aligned Manhattans, consistent region colours, rotation
## and permutation results -- fits this two-panel figure far better than the
## single EMMAX panel that occupied it before), plus a values file for the
## macros this module can honestly source.
##
## [!] PARTIAL, DELIBERATELY. values_panel.tex (R_3sp_blocks-generated, the
## OLDER pipeline) still defines macros this module has NOT rerun: the
## GCTA-vs-centred estimator comparison, the kinship-basis sweep, the
## decay-vs-map sweep. Overwriting those with a guess, or silently leaving
## them stale under a filename that claims to be current, would both be worse
## than what this does instead -- write ONLY the macros sourced from today's
## 02/03/04/05 outputs to a SEPARATE file, loaded AFTER values_panel.tex so it
## wins for exactly those macros and nothing else. Same layering values.tex
## already uses for its own hand-set descriptors.
## =============================================================================
suppressMessages({library(data.table)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_3sp"), "R", "00_config.R"))
STAGE <- "06_manuscript"
say("=== %s ===\n\n", STAGE)

sc <- readRDS(file.path(PATHS$out, "03_scan", "scan.rds"))
lf <- readRDS(file.path(PATHS$out, "04_lfmm", "lfmm.rds"))
J  <- fread(file.path(PATHS$out, "05_manhattan", "joint_regions.csv"))

OUT_DIR <- file.path(PATHS$out, STAGE); dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## ---- figure ------------------------------------------------------------------
## PK: use 07_figure_emmax.R's template (single-marker q vs cluster-level regions, EcoPeak
## triangles, Eda arrow) for figure4, not the two-panel EMMAX/LFMM comparison from
## 05_manhattan.R. This also resolves a mismatch flagged earlier: the slot's caption says
## "SNP- and cluster-level Manhattan plots", which this template IS, and the two-engine
## comparison was not. 05_manhattan.R's joint figure remains available as supplementary/
## engine-comparison material -- not deleted, just not what figure4 delivers.
file.copy(file.path(PATHS$out, "07_figure_emmax", "emmax_consensus_manhattan.pdf"),
         file.path(OUT_DIR, "figure4_stickleback_validation.pdf"), overwrite = TRUE)
say("[1] figure4_stickleback_validation.pdf <- 07_figure_emmax/emmax_consensus_manhattan.pdf\n")

## ---- values (module_3sp-sourced subset only) ---------------------------------
con <- sc$consensus; sim <- sc$simes
V <- list()
add <- function(name, value, src, note = "") V[[length(V)+1]] <<-
  list(name = name, value = as.character(value), src = src, note = note)

add("StickTestedUnits", format(nrow(con$test$units), big.mark=","), "03_scan.R")
add("StickSignificantClusters", sum(con$test$units$significant), "03_scan.R", "consensus, EMMAX")
add("StickSimesClusters", sum(sim$test$units$significant), "03_scan.R", "Simes, EMMAX")
sig_con <- con$test$units[significant==TRUE]$unit_id; sig_sim <- sim$test$units[significant==TRUE]$unit_id
add("StickJaccardStatistic", sprintf("%.3f", length(intersect(sig_con,sig_sim))/length(union(sig_con,sig_sim))),
    "03_scan.R", "consensus vs Simes, EMMAX, same units")
add("StickJaccardEngine", sprintf("%.3f",
    { sig_lfmm <- lf$test$units[significant==TRUE]$unit_id
      length(intersect(sig_sim,sig_lfmm))/length(union(sig_sim,sig_lfmm)) }),
    "04_lfmm.R", "EMMAX-Simes vs LFMM-Simes, same statistic")
add("StickReportedRegions", nrow(con$test$regions), "03_scan.R", "consensus, EMMAX")
add("StickEcoPeakRegions", con$rotation$observed, "03_scan.R")
add("StickEcoPeakPct", sprintf("%.1f", 100*con$rotation$observed/nrow(con$test$regions)), "03_scan.R")
add("StickMedianOccupancy", sprintf("%.3f", median(con$test$regions$occupancy)), "03_scan.R")
add("EcoPeakRotationP", sprintf("%.4f", con$rotation$p), "03_scan.R", "consensus, EMMAX")
add("StickRegionalPermP", sprintf("%.4f", con$null$p), "03_scan.R", "consensus, EMMAX, 1000 surrogates")
add("StickRealisedFDR", sprintf("%.1f", 100*con$null$realised_fdr), "03_scan.R")
## StickSingleMarkerHits deliberately NOT emitted here: 03_scan.R does not save the
## single-marker BH count, and emitting a placeholder string would render as literal
## garbled text if ever cited -- worse than leaving values_panel.tex's existing (real,
## if from the older pipeline) value in place. Add it once 03_scan.R saves the real number.

## new: LFMM, not present in the old values_panel.tex at all
add("StickLFMMSignificantClusters", sum(lf$test$units$significant), "04_lfmm.R", "Simes only -- LFMM p-values are precomputed")
add("StickLFMMRegions", nrow(lf$test$regions), "04_lfmm.R")
add("StickLFMMRotFold", sprintf("%.2f", lf$rotation$fold), "04_lfmm.R")
add("StickLFMMRotP", sprintf("%.4f", lf$rotation$p), "04_lfmm.R", "not significant -- state this, do not omit it")
add("StickJointRegions", nrow(J), "05_manhattan.R")
add("StickJointShared", sum(J$shared), "05_manhattan.R")
add("StickJointEMMAXOnly", sum(J$in_EMMAX & !J$in_LFMM), "05_manhattan.R")
add("StickJointLFMMOnly", sum(!J$in_EMMAX & J$in_LFMM), "05_manhattan.R")

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
"% GENERATED by module_3sp/R/06_manuscript.R. DO NOT EDIT BY HAND.",
"% ============================================================================",""), con_f)
for (v in V) writeLines(sprintf("\\providecommand{\\%s}{}\\renewcommand{\\%s}{%s}%s",
  v$name, v$name, v$value,
  sprintf("  %% %s%s", v$src, if (nzchar(v$note)) paste0(" -- ", v$note) else "")), con_f)
close(con_f)
say("[2] wrote %s -- %d macros\n", OUT_TEX, length(V))

write_receipt(STAGE,
  inputs = c(file.path(PATHS$out, "03_scan", "_receipt.rds"),
             file.path(PATHS$out, "04_lfmm", "_receipt.rds"),
             file.path(PATHS$out, "05_manhattan", "_receipt.rds")),
  params = list(), outputs = c(file.path(OUT_DIR, "figure4_stickleback_validation.pdf"), OUT_TEX))
say("\n[3] receipt: %s\n", receipt_path(STAGE))
