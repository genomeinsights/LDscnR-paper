## =============================================================================
## module_3sp/R/11_table_sensitivity.R
##
## Renders 10_sensitivity.R's output as a LaTeX tabular -- same style as the
## manuscript's existing tables/table4_reported_regions.tex (plain tabular,
## \hline rules, no booktabs/siunitx dependency), so \TableOrPlaceholder can
## \input{} it directly once wired into main.tex. Written to PATHS$tables (the
## same "one place to browse before hand-picking" convention PATHS$figures
## already uses) -- NOT auto-copied into LDscnR_manuscript/tables/, matching
## how figures are moved over by hand, not automatically.
##
## perm_p and realised_fdr both come from the SAME permutation run
## (ld_outlier_perm()), not from BH's nominal alpha. realised_fdr =
## mean(surrogate significant count) / observed significant count.
##
## TERMINOLOGY (PK, 2026-09-05): this is an empirical NULL-TO-OBSERVED
## DISCOVERY RATIO, not a realised false-discovery proportion/rate -- the
## truth status of the observed candidates is unknown, so no FDR/FDP can be
## "realised" from this quantity alone. It answers "how many discoveries
## would this null typically produce, relative to what we actually saw?",
## not "what fraction of what we saw is false." Reported here as "Null/obs."
## The LDscnR package field name (ld_outlier_perm.R's realised_fdr) is
## unchanged -- renaming the package API is a separate, larger decision --
## but nothing manuscript-facing should call this an FDR.
##
## Neither perm_p nor this ratio exists for LFMM (no permutation null -- no
## fast per-surrogate path for its precomputed p-values, same reason stated
## throughout 04_lfmm.R).
## =============================================================================
suppressMessages({library(data.table)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_3sp"), "R", "00_config.R"))
say("=== 11_table_sensitivity ===\n\n")

SENS_PATH <- file.path(PATHS$out, "10_sensitivity", "sensitivity.rds")
s <- readRDS(SENS_PATH)
S <- copy(s$summary)

## ---- formatting ----------------------------------------------------------------
ANALYSIS_LAB <- c(consensus = "EMMAX (consensus)", simes = "EMMAX (Simes)", lfmm = "LFMM (Simes)")
S[, analysis_lab := ANALYSIS_LAB[analysis]]

## grid label: parse "rho=X.XX_sf=N" back into a formatted, math-mode string, tagging the
## canonical point explicitly rather than relying on row position to convey it.
S[, is_canonical := grid == "rho=0.50_sf=8"]
S[, c("rho_val","sf_val") := {
  m <- regmatches(grid, regexec("rho=([0-9.]+)_sf=([0-9]+)", grid))
  list(vapply(m, `[`, "", 2), vapply(m, `[`, "", 3))
}]
S[, grid_lab := sprintf("$\\rho=%s$, floor$=%s$%s", rho_val, sf_val,
                        ifelse(is_canonical, " (canonical)", ""))]

## on-EcoPeak as "n/N (pp%)"
S[, on_ecopeak_lab := sprintf("%d/%d (%.0f\\%%)", on_ecopeak, regions, 100*on_ecopeak/regions)]

## p-values: both ld_region_rotation() and ld_outlier_perm() compute p as
## (1 + #surrogates >= observed) / (B + 1), which has a hard floor of 1/(B+1) -- p can
## never be SMALLER than that, only equal to it (when zero surrogates matched or beat the
## observed count). At the floor, report "$<10^{-4}$"-style (not a false-precision decimal),
## using <= with a small tolerance since the floor and the observed p are computed the same
## way but may differ in the last bit. perm_p's floor depends on the ARM (NPERM_CONSENSUS =
## 1000 vs NPERM_SIMES = 200 -- LFMM has no permutation null at all, already "--").
fmt_p <- function(p, floor_n) {
  floor_val <- 1/(floor_n+1)
  ifelse(is.na(p), "--",
    ifelse(p <= floor_val * (1 + 1e-9), sprintf("$<%.4f$", floor_val), sprintf("%.4f", p)))
}
S[, rotation_p_lab := fmt_p(rotation_p, N_ROTATIONS)]
S[, perm_floor_n := fifelse(analysis == "consensus", NPERM_CONSENSUS,
                            fifelse(analysis == "simes", NPERM_SIMES, NA_integer_))]
S[, perm_p_lab := fmt_p(perm_p, perm_floor_n)]
S[, null_obs_lab := ifelse(is.na(realised_fdr), "--", sprintf("%.1f\\%%", 100*realised_fdr))]
S[, fold_lab := sprintf("%.2f", rotation_fold)]

## ---- assemble rows, one \hline between grid points (not between every analysis row) ----
GROUP_ORDER <- c("rho=0.50_sf=8", "rho=0.50_sf=4", "rho=0.50_sf=12", "rho=0.50_sf=16",
                 "rho=0.50_sf=24", "rho=0.50_sf=48", "rho=0.30_sf=8", "rho=0.70_sf=8")
S[, grid := factor(grid, levels = GROUP_ORDER)]
setorder(S, grid)
ARM_ORDER <- c("consensus","simes","lfmm")
S[, analysis := factor(analysis, levels = ARM_ORDER)]
setorder(S, grid, analysis)

## longtable, not resizebox+tabular (PK, 2026-09-05): 10 columns still overflows a normal
## text page even at \scriptsize, and a non-breaking tabular forced 24 rows onto one page
## via \resizebox made the text hard to read. longtable breaks across pages at whatever size
## the surrounding landscape+\scriptsize environment sets (see Supplementary.tex), and carries
## its own \caption+\label so this file is \input directly, NOT wrapped in \TableOrPlaceholder
## (which supplies its own caption/label and would duplicate both).
lines <- c(
  "% GENERATED by module_3sp/R/11_table_sensitivity.R -- do not edit by hand.",
  "\\begin{longtable}{llrrrlrrrr}",
  paste0("\\caption{One-factor-at-a-time sensitivity of the stickleback outlier scan to the ",
         "complexity-reduction $\\rho$ and the Stage-1 size floor. Null/obs.\\ is the mean number ",
         "of significant units across regional phenotype permutations divided by the observed ",
         "number, i.e.\\ an empirical null-to-observed discovery ratio -- not a realised ",
         "false-discovery proportion, because the truth status of observed candidates is ",
         "unknown.}\\label{tab:stickleback-sensitivity}\\\\"),
  "\\hline",
  "Grid point & Analysis & Tested & Sig. & Regions & On EcoPeak & Fold & Rotation $p$ & Perm. $p$ & Null/obs. \\\\",
  "\\hline",
  "\\endfirsthead",
  "\\multicolumn{10}{l}{\\tablename~\\thetable{} continued}\\\\",
  "\\hline",
  "Grid point & Analysis & Tested & Sig. & Regions & On EcoPeak & Fold & Rotation $p$ & Perm. $p$ & Null/obs. \\\\",
  "\\hline",
  "\\endhead")
for (g in levels(S$grid)) {
  rows <- S[grid == g]
  for (i in seq_len(nrow(rows))) {
    r <- rows[i]
    grid_cell <- if (i == 1) r$grid_lab else ""
    lines <- c(lines, sprintf("%s & %s & %d & %d & %d & %s & %s & %s & %s & %s \\\\",
      grid_cell, r$analysis_lab, r$tested, r$significant, r$regions, r$on_ecopeak_lab,
      r$fold_lab, r$rotation_p_lab, r$perm_p_lab, r$null_obs_lab))
  }
  lines <- c(lines, "\\hline")
}
lines <- c(lines, "\\end{longtable}")

OUT_TEX <- file.path(PATHS$tables, "table5_sensitivity.tex")
writeLines(lines, OUT_TEX)
say("[1] wrote %s -- %d grid points x %d analyses = %d rows\n",
    OUT_TEX, uniqueN(S$grid), uniqueN(S$analysis), nrow(S))
say("\n    Placement: LDscnR_manuscript/Supplementary.tex, stickleback validation subsection,\n")
say("    landscape + \\scriptsize, \\input directly (longtable carries its own caption+label).\n")
