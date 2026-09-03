## =====================================================================
## module_sim_LDscnR / analyse_one_dataset.R
##
## Analyse ONE dataset from p-values. Method-agnostic by construction: the only
## things this script reads are a panel and a set of p-values.
##
##   Rscript module_sim_LDscnR/analyse_one_dataset.R <panel.rds> <pvals.rds> [outdir]
##
## INPUT 1 -- panel.rds, the LD description of the study. A list with:
##   GTs        genotype dosages, individuals x markers, colnames = markers
##   map        data.table with marker, Chr, Pos (one row per marker, same order)
##   ld_ws      local-LD weights, markers x rho windows, rownames = markers
##   decay_sum  per-chromosome LD-decay summary
##   map$true_QTN   OPTIONAL. Present in simulations, absent on real data. If
##                  present, the evaluation section at the end runs; if not, it
##                  is skipped and nothing else changes.
##
## INPUT 2 -- pvals.rds, the output of ONE association method. A list with:
##   p_obs      named numeric, one p-value per marker, names = map$marker
##   p_perm     the permuted datasets, as a markers-by-B matrix, a list of
##              vectors, or anything ld_null_from_p() accepts
##   basis, engine   free-text labels, carried into every printed output
##
## THIS IS THE WHOLE POINT OF THE SPLIT. The panel depends on the panel; the
## p-values depend on the method. To analyse the same data with a different
## method -- LFMM instead of EMMAX, BayPass, a GLM, anything -- you produce a new
## pvals.rds and rerun this script against the same panel. Nothing here knows or
## cares which scanner produced the numbers, only that the observed data and
## every permuted dataset went through the SAME one.
##
## Producing p-values is your job and is deliberately outside this script. For
## the simulations, module_sim_LDscnR/prepare_sim_inputs.R does it.
##
## Env vars (all optional): TAU, LMIN, RHO_LD, DCAP, ALPHA, FDR, RUN_C2.
## =====================================================================
suppressMessages({ library(data.table); library(LDscnR) })
## LDscnR keeps %||% internal, and base R only gained it in 4.4 -- define it
## locally so this script runs on older R too.
`%||%` <- function(a, b) if (is.null(a)) b else a

a <- commandArgs(trailingOnly = TRUE)
if (length(a) < 2) stop("usage: analyse_one_dataset.R <panel.rds> <pvals.rds> [outdir]")
PANEL_F <- a[1]; PVALS_F <- a[2]
OUT <- if (length(a) >= 3) a[3] else Sys.getenv("OUT", "module_sim_LDscnR/results")
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)

## Operating point. alpha is FIXED at 0.05 -- the alpha sweep belongs to the
## PR-AUC comparison of methods, not to region discovery. tau and l_min are
## nominal: ld_region_c2() integrates over them rather than trusting one cell.
##
## l_min = 3 IS A CHOICE, NOT A CALIBRATION, and it is provisional (PK, 2026-08-26).
## calibrate_lmin() was meant to remove the need for a choice, but on these
## simulations it returns 13-510 across ten genomes (only 2 of 9 at or below 20),
## because it takes 1 + the 99th percentile of the null's largest-region size and
## the surrogate region-size distribution is heavy-tailed. It is not usable as an
## operating rule.
##
## The choice is between two objectives that disagree, measured over 9 genomes of
## V2_c1 (module_sim_LDscnR/pr_cscore_vs_alpha.R, operating_point_vs_truth.R):
##
##   l_min      C-score advantage over BH alpha      PR at (tau 0.05, l_min) vs
##              (paired PR-AUC difference)           the per-genome optimum
##      1              +0.100                                 --
##      3              +0.069  (C wins 7/9)                  33%
##     10              +0.018  (C wins 5/9)                  60%
##     20              +0.005  (C wins 4/9)                   --
##
## The C-score and the region-size filter are SUBSTITUTES: both suppress isolated
## spurious markers, so tightening l_min raises absolute performance while eroding
## the C-score's margin over plain BH. l_min = 3 keeps the margin the method is
## claimed on; l_min = 10 -- which is what Fang et al. (2021) published -- gives
## better absolute PR. Neither is wrong; they answer different questions.
##
## Nothing here is significant at n = 9, so this should be revisited once more
## (V, c) cells are analysed.
PAR <- list(alpha  = as.numeric(Sys.getenv("ALPHA",  "0.05")),
            qstar  = seq(0, 0.95, by = 0.05),
            tau    = as.numeric(Sys.getenv("TAU",    "0.05")),
            l_min  = as.integer(Sys.getenv("LMIN",   "3")),
            fdr    = as.numeric(Sys.getenv("FDR",    "0.05")),
            rho_ld = as.numeric(Sys.getenv("RHO_LD", "0.75")),
            dcap   = as.numeric(Sys.getenv("DCAP",   "1e5")),   # matches the bundles' stage-2 cap
            tau_grid  = seq(0.05, 0.30, by = 0.05),
            lmin_grid = c(2L, 3L, 5L, 10L, 20L))
RUN_C2 <- Sys.getenv("RUN_C2", "1") == "1"

## ---- load and check the two inputs -----------------------------------------
panel <- readRDS(PANEL_F); P <- readRDS(PVALS_F)
map <- as.data.table(panel$map)
for (f in c("GTs", "map", "ld_ws", "decay_sum")) if (is.null(panel[[f]])) stop("panel is missing `", f, "`")
for (f in c("p_obs", "p_perm")) if (is.null(P[[f]])) stop("pvals is missing `", f, "`")
if (!all(c("marker", "Chr", "Pos") %in% names(map))) stop("panel$map needs marker, Chr and Pos")
if (length(P$p_obs) != nrow(map))
  stop(sprintf("p_obs has %d values but the panel has %d markers", length(P$p_obs), nrow(map)))
## Names are checked, never used to reorder: a p-vector of the right length in
## the wrong order is silently wrong, and the caller's pipeline is where that
## originates. ld_cscore() enforces this too; failing here just names the file.
if (!is.null(names(P$p_obs)) && !identical(names(P$p_obs), map$marker))
  stop("p_obs names do not match panel$map$marker -- reorder before calling")

nB <- if (is.matrix(P$p_perm)) ncol(P$p_perm) else length(P$p_perm)
cat(sprintf("[1] panel %s: %d markers x %d individuals, %d chromosomes\n",
            basename(PANEL_F), nrow(map), nrow(panel$GTs), uniqueN(map$Chr)))
cat(sprintf("    p-values %s: engine %s, basis %s, B = %d, p-floor = %.4f\n",
            basename(PVALS_F), P$engine %||% "?", P$basis %||% "?", nB, 1 / (1 + nB)))
cat(sprintf("    truth in panel: %s\n",
            if ("true_QTN" %in% names(map)) sprintf("yes, %d QTN", sum(map$true_QTN %in% TRUE)) else "no (real data)"))
flush.console()

## ---- the analysis: p-values in, regions out --------------------------------
cat(sprintf("\n[2] ld_scan: tau=%.2f l_min=%d rho_ld=%.2f alpha=%.2f fdr=%.2f c2=%s\n",
            PAR$tau, PAR$l_min, PAR$rho_ld, PAR$alpha, PAR$fdr, RUN_C2)); flush.console()
t0 <- Sys.time()
fit <- ld_scan(P$p_obs, P$p_perm, panel$ld_ws, map[, .(marker, Chr, Pos)], panel$GTs, panel$decay_sum,
               tau = PAR$tau, l_min = PAR$l_min, fdr = PAR$fdr,
               alpha = PAR$alpha, qstar = PAR$qstar,
               rho_ld = PAR$rho_ld, dcap = PAR$dcap,
               c2 = RUN_C2, tau_grid = PAR$tau_grid, lmin_grid = PAR$lmin_grid,
               basis = P$basis %||% "user-supplied", engine = P$engine %||% "user-supplied",
               verbose = TRUE)
cat(sprintf("    done in %.1f min\n\n", as.numeric(Sys.time() - t0, units = "mins")))

## The gate first, deliberately: a p-value drawn from a null whose surrogates
## reproduce the observed background is not interpretable, and no amount of care
## downstream recovers it.
print(fit$gate)
if (!all(fit$gate$pass))
  cat("\n!! This null FAILED its gate. The regions below are NOT interpretable.\n")
print(fit$regions)

regions <- copy(fit$regions$regions)
if (RUN_C2 && !is.null(fit$c2)) {
  regions <- merge(regions, fit$c2$regions[, .(chr, lo, hi, c2)], by = c("chr", "lo", "hi"), all.x = TRUE)
  setorder(regions, -c2, -s_R)
  cat(sprintf("\n[3] C-squared over %d usable grid cells of %d\n", fit$c2$n_usable, fit$c2$n_cells))
  cat("    (q_R says 'not a structure artifact'; when regions pile up at the p-floor it\n")
  cat("     cannot rank them, and C-squared is what does.)\n")
}

## ---- evaluation against truth: simulations only ----------------------------
## Runs only if the panel carries true_QTN. On real data this section is silently
## skipped -- which is exactly why the analysis above never looks at it.
if ("true_QTN" %in% names(map) && nrow(regions)) {
  ## true_pos_QTN, not true_QTN: flag_true_qtns() keeps only QTNs with MAF > 0.1
  ## carrying >= 5% of their chromosome's additive variance. The others are
  ## undetectable in principle and counting them understates recall.
  qtn <- map[true_pos_QTN %in% TRUE, .(Chr = as.character(Chr), Pos)]
  regions[, has_qtn := vapply(seq_len(.N), function(i)
    any(qtn$Chr == chr[i] & qtn$Pos >= lo[i] & qtn$Pos <= hi[i]), logical(1))]
  ## A pooled cell is a GENOME and about half its chromosomes carry no detectable
  ## QTN -- as in any real genome. The QTN-free set is derived from the truth, not
  ## assumed to be "the Chr2s": some Chr1s carry none either, and hardcoding _Chr2
  ## both undercounts and implies Chr2 is a designated control.
  qtn_free_chr <- as.character(map[, .(nq = sum(true_pos_QTN %in% TRUE)), by = Chr][nq == 0, Chr])
  regions[, on_qtn_free := chr %in% qtn_free_chr]
  sig <- regions[sig == TRUE]
  rec <- sum(vapply(seq_len(nrow(qtn)), function(j) nrow(sig) > 0 &&
    any(sig$chr == qtn$Chr[j] & sig$lo <= qtn$Pos[j] & sig$hi >= qtn$Pos[j]), logical(1)))
  q <- stats::p.adjust(P$p_obs, "BH")
  cat("\n[4] TRUTH (simulation only)\n")
  cat(sprintf("    %d true QTN | %d significant regions | %d contain a QTN | %d on a QTN-free chromosome (%.0f%% vs %.0f%% by chance)\n",
              nrow(qtn), nrow(sig), sum(sig$has_qtn), sum(sig$on_qtn_free),
              100 * mean(sig$on_qtn_free), 100 * length(qtn_free_chr) / uniqueN(map$Chr)))
  cat(sprintf("    recall %d/%d | precision %.3f\n", rec, nrow(qtn),
              sum(sig$has_qtn) / max(nrow(sig), 1)))
  cat(sprintf("    single-SNP BH on the SAME p-values: %d markers, %d on neutral chromosomes, %d are QTN\n",
              sum(q < PAR$fdr), sum(q < PAR$fdr & as.character(map$Chr) %in% qtn_free_chr),
              sum(q[map$true_QTN %in% TRUE] < PAR$fdr)))
}

stem <- sub("^pvals_", "", sub("[.]rds$", "", basename(PVALS_F)))
fwrite(regions, file.path(OUT, sprintf("regions_%s.csv", stem)))
saveRDS(fit,    file.path(OUT, sprintf("scan_%s.rds", stem)))
cat(sprintf("\n[5] wrote regions_%s.csv and scan_%s.rds to %s\n", stem, stem, OUT))
