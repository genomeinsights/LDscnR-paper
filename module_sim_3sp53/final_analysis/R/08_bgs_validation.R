## =============================================================================
## final_analysis/R/08_bgs_validation.R
##
## Compact BGS quality-control validation (CLAUDE_REANALYSIS_INSTRUCTIONS.md,
## "One compact BGS validation"): does BGS produce its intended genomic
## effect -- reduced nucleotide diversity, more so in low-recombination
## regions (background selection's classic signature: linked selection's
## footprint is larger where recombination is rarer)? Paired bgs-vs-nobgs
## comparison at matched (cell, rep, env), overall and by recombination-rate
## stratum, computed BOTH pre- and post-MAF-filter (parse_nemo_run(...,
## keep_prefilter=TRUE), 01_parse_nemo.R) -- a fresh reparse of raw NEMO
## output per combo, not a read of the analysis-ready out_final_v1/
## 02_build_ld_units.R cache, so this script runs standalone off raw NEMO +
## recmaps only.
##
## CONFIRMED, not just a disclosed risk: computing on the pipeline's
## standard MAF>0.10-filtered marker set alone gives a BACKWARDS, spurious
## result -- among markers that SURVIVE the filter, mean He is HIGHER under
## bgs (e.g. +0.020 overall, 95% CI excludes 0), the opposite of the
## intended signature. The pre-filter marker-RETENTION ratio tells the real
## story: bgs keeps ~17% fewer segregating markers overall (ratio 0.833),
## and the loss is steepest at low recombination (ratio 0.738) vs high
## (0.931) -- the correct classic BGS signature, exactly backwards from
## what post-filter He alone would suggest. MAF filtering removes most of
## what BGS actually did (drives rare variants to loss) before per-marker
## He ever sees it; the residual positive He among survivors is consistent
## with associative overdominance (a heterozygote at a surviving neutral
## site is less likely to also be homozygous for a linked, partially
## recessive deleterious allele) -- a real, secondary effect, not the
## primary validated result. This reproduces, in this exact dataset, a
## measurement trap an earlier BGS investigation in this project had
## already flagged (MAF-filtered heterozygosity vs raw SNP density).
## PRIMARY reportable metric: pre-filter ratio_n_markers by recombination
## bin (figureS_simulation_bgs_qc.R's left panel). Post-filter He is kept
## and reported (right panel) as the trap itself, not silently dropped.
##
## Recombination-rate bins: same scheme as the old exploratory
## module_sim_3sp53/R/08_bgs_recomb.R (reused as implementation material,
## not re-derived, per the instructions' "Sources of truth" ranking). Bin
## edges are cut from the FULL rec_map<rep>.rds (all 1e6 reference
## positions), not from whichever markers survive MAF filtering -- keeps
## bin boundaries identical across cells/tags/envs at the same rep.
## rec_rate is ~33% exactly zero (recombination coldspots), so zero is its
## own bin rather than folded into a quantile spanning 0 to a positive
## value: {zero, then tertiles of the nonzero distribution}.
##
## Uncertainty: a rep-cluster bootstrap on the paired (bgs - nobgs) He
## difference, matching the primary pipeline's own convention (instructions,
## "Pooling and uncertainty": "For BGS versus no-BGS contrasts, retain the
## matched map, environment and seed identifiers in the same resample").
## =============================================================================
suppressMessages({library(data.table); library(parallel)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "R", "00_config.R"))
source(file.path(PATHS$module, "R", "01_parse_nemo.R"))
say("=== bgs_validation ===\n\n")

BIN_LABELS <- c("zero", "low", "medium", "high")

.recmap_bins <- function(rep) {
  recmap <- readRDS(file.path(PATHS$raw_recmap_dir, sprintf("rec_map%d.rds", rep)))
  nz <- recmap$rec_rate[!is.na(recmap$rec_rate) & recmap$rec_rate > 0]
  breaks <- c(-Inf, 0, stats::quantile(nz, c(1 / 3, 2 / 3), na.rm = TRUE), Inf)
  recmap[, Chr := paste0("Chr", Chr)]
  data.table::setkey(recmap, Chr, bp)
  list(recmap = recmap, breaks = breaks)
}

## [!] FOUND running this on nobgs_V0.5_c1_rep1_env1 vs bgs_V0.5_c1_rep1_env1:
## bgs combos carry an extra locus type, "delet" (the deleterious loci that
## implement background selection itself -- 137 markers, present in NO
## nobgs combo, none), that MAF-passes and would otherwise be folded into
## the "bgs" He average with no nobgs counterpart at all. Left in, this
## produced a spurious, backwards result (bgs LOOKED more diverse than
## nobgs -- the classic BGS signature is REDUCED diversity, especially at
## low recombination) purely from analysis-set composition, not biology.
## Restricted to marker types present in BOTH designs (ntrl, QTN) so the
## comparison is apples-to-apples -- the classic BGS signature is measured
## at linked NEUTRAL/QTN sites, not at the deleterious loci themselves,
## which are under direct selection by construction and were never meant
## to represent standing background diversity.
.he_by_bin <- function(GTs, map, rb) {
  keep <- map$type %in% c("ntrl", "QTN")
  GTs <- GTs[, keep, drop = FALSE]; map <- map[keep]
  p <- colSums(GTs) / nrow(GTs) / 2
  he <- 2 * p * (1 - p)
  rec_rate <- rb$recmap[.(map$Chr, map$Pos), rec_rate, on = c("Chr", "bp"), mult = "first"]
  bin <- cut(rec_rate, breaks = rb$breaks, labels = BIN_LABELS, include.lowest = TRUE)
  dt <- data.table(he = he, bin = bin)
  overall <- dt[, .(bin = "overall", he_mean = mean(he), n_markers = .N)]
  by_bin <- dt[!is.na(bin), .(he_mean = mean(he), n_markers = .N), by = bin]
  rbind(overall, by_bin, use.names = TRUE)
}

summarise_bgs_validation <- function(B = N_BOOTSTRAP) {
  combos <- CJ(cell = CELLS_ALL, rep = REPS_ALL, env = ENVS_ALL)
  say("[1] He (overall + by recomb bin) for %d (cell,rep,env) triples x 2 tags\n", nrow(combos))

  rows <- mclapply(seq_len(nrow(combos)), function(i) {
    cell <- combos$cell[i]; rep <- combos$rep[i]; env <- combos$env[i]
    rb <- .recmap_bins(rep)
    out <- list()
    for (tag in TAGS_ALL) {
      res <- tryCatch(parse_nemo_run(tag, cell, rep, env, offset = 0L, keep_prefilter = TRUE),
                      error = function(e) { message(sprintf("[bgs_validation] %s_%s_rep%d_env%d: %s", tag, cell, rep, env, conditionMessage(e))); NULL })
      if (is.null(res)) next
      hb_post <- .he_by_bin(res$GTs, res$map, rb)[, filt := "post_maf"]
      hb_pre  <- .he_by_bin(res$GTs_prefilter, res$map_prefilter, rb)[, filt := "pre_maf"]
      hb <- rbind(hb_post, hb_pre)
      hb[, `:=`(tag = tag, cell = cell, rep = rep, env = env)]
      out[[tag]] <- hb
    }
    if (!length(out)) return(NULL)
    rbindlist(out)
  }, mc.cores = 7)
  dt <- rbindlist(rows)
  say("    %d rows\n", nrow(dt))

  say("\n[2] pairing bgs vs nobgs at matched (cell, rep, env, bin, filt)\n")
  wide <- dcast(dt, cell + rep + env + bin + filt ~ tag, value.var = c("he_mean", "n_markers"))
  wide <- wide[!is.na(he_mean_nobgs) & !is.na(he_mean_bgs)]
  wide[, diff := he_mean_bgs - he_mean_nobgs]
  wide[, diff_n_markers := n_markers_bgs - n_markers_nobgs]
  say("    %d paired rows\n", nrow(wide))

  ## PK, after auditing the first version of this script: n_markers was
  ## computed per combo but never carried through to the pooled output --
  ## added here as its own paired diff/ratio (bootstrapped the same way as
  ## the He diff), not just reported as two unpaired means, so retained-
  ## marker COUNT (not just per-marker He among survivors) is a directly
  ## checkable part of the compact validation, matching this script's own
  ## stated rationale for reporting n_markers in the first place.
  say("\n[3] rep-cluster bootstrap (B=%d), grand-pooled across cells, pre- vs post-MAF-filter\n", B)
  ci <- function(x, probs = c(0.025, 0.975)) stats::quantile(x, probs, na.rm = TRUE, names = FALSE)
  n_rep <- length(REPS_ALL)
  set.seed(SEEDS[["bootstrap"]])
  draws <- matrix(sample.int(n_rep, size = n_rep * B, replace = TRUE), nrow = n_rep, ncol = B)
  mult <- apply(draws, 2, tabulate, nbins = n_rep)

  pool_bin <- function(bin_name, filt_name) {
    sub <- wide[bin == bin_name & filt == filt_name]
    rep_means <- sub[, .(diff = mean(diff), he_nobgs = mean(he_mean_nobgs), he_bgs = mean(he_mean_bgs),
                         n_markers_nobgs = mean(n_markers_nobgs), n_markers_bgs = mean(n_markers_bgs),
                         diff_n_markers = mean(diff_n_markers)), by = rep]
    rep_means <- rep_means[match(REPS_ALL, rep)]
    boot <- as.numeric(crossprod(mult, rep_means$diff)) / colSums(mult)
    boot_nm <- as.numeric(crossprod(mult, rep_means$diff_n_markers)) / colSums(mult)
    data.table(filt = filt_name, bin = bin_name, n_pairs = nrow(sub),
              n_markers_nobgs = mean(rep_means$n_markers_nobgs), n_markers_bgs = mean(rep_means$n_markers_bgs),
              mean_diff_n_markers = mean(rep_means$diff_n_markers),
              diff_n_markers_ci_lo = ci(boot_nm)[1], diff_n_markers_ci_hi = ci(boot_nm)[2],
              ratio_n_markers = mean(rep_means$n_markers_bgs) / mean(rep_means$n_markers_nobgs),
              he_nobgs = mean(rep_means$he_nobgs), he_bgs = mean(rep_means$he_bgs),
              mean_diff_he = mean(rep_means$diff), diff_ci_lo = ci(boot)[1], diff_ci_hi = ci(boot)[2])
  }
  pooled <- rbindlist(lapply(c("pre_maf", "post_maf"), function(ft)
    rbindlist(lapply(c("overall", BIN_LABELS), pool_bin, filt_name = ft))))
  pooled[, `:=`(bin = factor(bin, levels = c("overall", BIN_LABELS)), filt = factor(filt, levels = c("pre_maf", "post_maf")))]
  setorder(pooled, filt, bin)
  print(pooled)

  list(raw = dt, paired = wide, pooled = pooled)
}

if (sys.nframe() == 0L) {
  res <- summarise_bgs_validation()
  dir.create("results", showWarnings = FALSE)
  fwrite(res$pooled, "results/simulation_bgs_validation.tsv", sep = "\t")
  say("\nwrote results/simulation_bgs_validation.tsv\n")
}
