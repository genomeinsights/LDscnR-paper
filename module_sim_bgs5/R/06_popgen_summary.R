## module_sim/R/06_popgen_summary.R
##
## Characterize the SIMULATION ITSELF (not outlier-detection performance) for
## one (tag,cell,rep,env): Fst, number of detectable QTN, total Va, and a
## local-adaptation proxy. PK: summarize the simulations -- Fst, mean number
## of detectable QTN, Va, level of local adaptation, some estimate of the
## effect of bgs.
##
## Methodology checked against the superseded module_sim_LDscnR/
## parse_and_regen_sim_data.R before choosing what to reuse vs simplify:
##   - detectable QTN / per-QTN Va: EXACT match to flag_true_qtns() and
##     get_va() (2*p*(1-p)*allelic_value^2) already used in R/04_score.R --
##     nothing new, just carried forward here.
##   - Fst: the old code used SNPRelate::snpgdsFst(method="W&C84") on NEMO's
##     own patch id (here: env$pop, confirmed present in every bundle) --
##     reused as-is rather than Nei's Fst (the other option there), since
##     W&C84 is the more standard estimator and was the code's own primary
##     choice.
##   - local adaptation: the old local_adaptation() (Kawecki & Ebert
##     sympatric-allopatric contrasts, Jensen-corrected home/away Gaussian
##     fitness) is considerably more elaborate than what this summary needs
##     and requires per-patch local optima not carried in these bundles.
##     Simplified to its most basic, well-established proxy instead: the R^2
##     between each individual's additive breeding value at the QTN loci
##     (GTs %*% allelic_values) and its local environment (env$env) -- a
##     population is "locally adapted" to the extent its genotypes track the
##     environment, which this measures directly with data already on hand.
##   - BGS effect: no bgs-vs-nobgs collation existed in the old code (its own
##     bgs_windows() flagged the matched-ratio estimator as correct but never
##     implemented it). Not computed per-combo here at all -- R/07_popgen_pool.R
##     folds bgs into these same per-cell summaries as a colour, and the BGS
##     effect is then whatever visual/numeric contrast that produces between
##     the two tags at each cell, rather than a separate derived quantity.
suppressMessages({library(data.table); library(LDscnR); library(SNPRelate)})
devtools::load_all("~/gitlab/LDscnR")
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim"), "R", "00_config.R"))
STAGE <- "06_popgen_summary"
say("=== %s ===\n\n", STAGE)
invisible(check_ldscnr())

TARGET_TAG <- Sys.getenv("SIM_TAG", TAGS[1]); TARGET_CELL <- Sys.getenv("SIM_CELL", CELLS[1]); TARGET_ENV <- as.integer(Sys.getenv("SIM_ENV", ENVS[1])); TARGET_REP <- as.integer(Sys.getenv("SIM_REP", REPS[1]))
combo_id <- sprintf("%s_%s_rep%d_env%d", TARGET_TAG, TARGET_CELL, TARGET_REP, TARGET_ENV)

BUNDLE_PATH <- file.path(PATHS$out, "02_bundle",
  sprintf("bundle_%s_rep%d_%s_env%d.rds", TARGET_TAG, TARGET_REP, TARGET_CELL, TARGET_ENV))
if (!file.exists(BUNDLE_PATH)) stop("R/02_bundle.R has not produced: ", basename(BUNDLE_PATH))
b <- readRDS(BUNDLE_PATH)
GTs <- b$GTs; map <- b$map; env <- b$env

INPUTS <- BUNDLE_PATH
PARAMS <- list(maf_min = 0.1, p_va_min = 0.05, fst_method = "W&C84")
if (!stage_stale(STAGE, INPUTS, PARAMS, target = combo_id) && !nzchar(Sys.getenv("FORCE"))) {
  say("\nNothing to do. Set FORCE=1 to rerun anyway.\n"); quit(save = "no")
}

## ---- 1. Va, detectable QTN (identical to R/04_score.R) -------------------------
say("[1] Va, detectable QTN\n")
qtn_rows <- which(map$type == "QTN")
map[, Va := NA_real_]
map$Va[qtn_rows] <- vapply(qtn_rows, function(i) {
  p <- mean(GTs[, i]) / 2
  2 * p * (1 - p) * map$allelic_values[i]^2
}, 0)
map <- flag_true_qtns(map, va_col = "Va", maf_col = "MAF", maf_min = 0.1, p_va_min = 0.05)
n_qtn_total <- length(qtn_rows)
n_qtn_detectable <- sum(map$true_pos_QTN, na.rm = TRUE)
va_total <- sum(map$Va[qtn_rows], na.rm = TRUE)
va_detectable <- sum(map$Va[map$true_pos_QTN], na.rm = TRUE)
say("    %d QTN total (%d detectable) ; Va total = %.5f (detectable = %.5f)\n",
    n_qtn_total, n_qtn_detectable, va_total, va_detectable)

## ---- 2. Fst (W&C84), population = NEMO's own patch id (env$pop) ---------------
say("[2] Fst (W&C84), %d populations\n", uniqueN(env$pop))
gds_path <- file.path(PATHS$cache, paste0("popgen_", combo_id, ".gds"))
dir.create(dirname(gds_path), recursive = TRUE, showWarnings = FALSE)
gds <- create_gds_from_geno(geno = GTs, map = map, gds_path)
on.exit({ try(snpgdsClose(gds), silent = TRUE); unlink(gds_path) }, add = TRUE)
fst <- snpgdsFst(gds, population = factor(env$pop), method = "W&C84",
                 autosome.only = FALSE, verbose = FALSE)
say("    Fst = %.4f\n", fst$MeanFst)

## ---- 3. local adaptation proxy: breeding value ~ env, R^2 ----------------------
say("[3] local adaptation: breeding value (all QTN) ~ environment, R^2\n")
bv <- as.numeric(GTs[, qtn_rows, drop = FALSE] %*% map$allelic_values[qtn_rows])
local_adapt_r2 <- if (n_qtn_total > 0 && sd(bv) > 0) cor(bv, env$env)^2 else NA_real_
say("    R^2 = %s\n", if (is.na(local_adapt_r2)) "NA" else sprintf("%.4f", local_adapt_r2))

summary_row <- data.table(tag = TARGET_TAG, cell = TARGET_CELL, rep = TARGET_REP, env = TARGET_ENV,
  n_qtn_total = n_qtn_total, n_qtn_detectable = n_qtn_detectable,
  Va_total = va_total, Va_detectable = va_detectable,
  Fst = fst$MeanFst, local_adapt_r2 = local_adapt_r2)

OUT <- file.path(stage_dir(STAGE), sprintf("popgen_%s_rep%d_%s_env%d.rds", TARGET_TAG, TARGET_REP, TARGET_CELL, TARGET_ENV))
dir.create(stage_dir(STAGE), recursive = TRUE, showWarnings = FALSE)
saveRDS(summary_row, OUT)
write_receipt(STAGE, inputs = INPUTS, params = PARAMS, outputs = OUT, target = combo_id)
say("\n[4] wrote %s\n    receipt: %s\n", OUT, receipt_path(STAGE, combo_id))
