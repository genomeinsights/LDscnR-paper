## =============================================================================
## module_sim_LDscnR / decay_window_test.R
##
## Does a wider LD-decay window improve the fit? PK asked for 20 -> 50; the bgs5
## bundles are built at n_win_decay = 10 (itself raised from 5), so 10/20/50.
##
## THE CRITERION IS DETERMINACY, NOT THE VALUE OF a. Per ?compute_LD_decay, the
## fitted a/b/c are conditional on the settings, so a from two window sizes are
## not two estimates of one quantity and "which is right" is not answerable that
## way. What IS answerable: the two chromosomes in a bundle come from the SAME
## simulation, so a real difference in their decay rate should be small and
## their disagreement is a noise measure.
##
## THAT DISAGREEMENT IS WHAT |a - a_pred|/a_pred MEASURES HERE, AND ONLY HERE.
## The chromosomes are equal-size, so the size-corrected a_pred is exactly their
## MEAN and the deviation is identical for both by arithmetic -- it is half
## their relative difference, not agreement with an independent prediction. Two
## consequences: the independent unit is the FILE (n = 10), not the chromosome
## (n = 20), and this quantity does not carry the meaning it would on a real
## genome with many differently-sized chromosomes.
##
## RESULT: 10 -> 20 HELPS, 20 -> 50 DOES NOT. Paired over files, w = 20 beats
## w = 10 on 9 of 10 (sign p = 0.021) and median disagreement falls 0.273 ->
## 0.095; w = 50 is 5 of 10 against w = 20 (p = 1.0) at ~2.5x the windows.
## A subset of 5 files looked monotone through w = 50 and was misleading.
##
## THE COST IS THAT a REORDERS. Spearman(a_w10, a_w50) across chromosomes is
## only 0.60 and the median |change| in a is 72%, so ANY BETWEEN-CHROMOSOME
## RESULT THAT RANKS CHROMOSOMES BY a IS NOT ROBUST TO THIS SETTING and has to
## be recomputed at whatever window it will be reported at.
##
## The installed LDscnR predates compute_LD_decay()'s `seed` argument, so the
## RNG is pinned in the caller immediately before each call. b is estimated from
## an n_sub_bg subsample and moves between unseeded runs enough to shift derived
## thresholds; with the pin it is identical across all three windows (0.02996),
## which is the check that the pin worked.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/decay_window_test.R
## Env: SIM_DATA, CELL, TAG, ENV, FILES, WINDOWS, CORES, OUT
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(parallel)})
SIM   <- Sys.getenv("SIM_DATA", "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5")
OUT   <- Sys.getenv("OUT", "module_sim_LDscnR/results/operating_points")
CELL  <- Sys.getenv("CELL", "V0.5_c1")
TAG   <- Sys.getenv("TAG", "nobgs")
ENV   <- as.integer(Sys.getenv("ENV", "1"))
FILES <- as.integer(strsplit(Sys.getenv("FILES", "1,2,3,4,5,6,7,8,9,10"), ",")[[1]])
WIN   <- as.integer(strsplit(Sys.getenv("WINDOWS", "10,20,50"), ",")[[1]])
CORES <- as.integer(Sys.getenv("CORES", "10"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
## Everything except n_win_decay is held at the bundles' own DECAY_ARGS, so the
## only thing that differs between the three fits is the window count.
ARGS <- list(el_data_folder = NULL, min_maf_decay = 0.1, q = 0.95, n_sub_bg = 5000,
             overlap = 0.5, max_SNPs_decay = Inf, prob_robust = 0.95, max_pairs = 5000,
             ld_method = "corr", n_strata = 20, keep_el = FALSE, slide = 1000,
             rho_targets = c(0.99), cores = 1)

one <- function(i) {
  f <- sprintf("%s/adapt_%s_chr%d_%s_env%d.rds", SIM, TAG, i, CELL, ENV)
  if (!file.exists(f)) return(NULL)
  x  <- readRDS(f)
  gp <- tempfile(fileext = ".gds"); on.exit(unlink(gp), add = TRUE)
  g  <- create_gds_from_geno(geno = x$GTs, map = as.data.table(x$map), gp)
  rbindlist(lapply(WIN, function(w) {
    set.seed(1L)
    d <- do.call(compute_LD_decay, c(list(gds = g), ARGS, list(n_win_decay = w)))
    as.data.table(d$decay_sum)[, .(w = w, file = i, Chr, a, b, c, a_pred, n_w_used)] }))
}
r <- rbindlist(Filter(Negate(is.null), mclapply(FILES, one, mc.cores = CORES)))
r[, d95 := d_from_rho(a_pred, 0.95)]
fwrite(r, file.path(OUT, "decay_window_test_raw.csv"))

S <- r[, .(n_chr = .N, n_w_used = mean(n_w_used), a_med = median(a),
           a_spread = max(a) / min(a), b_med = median(b), c_med = median(c),
           d95_kb = median(d95) / 1e3,
           disagreement = median(abs(a - a_pred) / a_pred)), by = w]
fwrite(S, file.path(OUT, "decay_window_test_summary.csv"))
print(S)

## Paired over FILES, for the reason in the header: the two chromosomes in a
## file share one value of this quantity, so treating them as 20 independent
## points doubles n and halves every p-value.
r[, dev := abs(a - a_pred) / a_pred]
F <- dcast(r[, .(dev = dev[1]), by = .(file, w)], file ~ w, value.var = "dev")
ws <- as.character(WIN)
cat(sprintf("\npaired over %d files:\n", nrow(F)))
for (j in seq_along(ws)[-1]) {
  x <- F[[ws[j - 1]]]; y <- F[[ws[j]]]
  cat(sprintf("  %s -> %s : %d/%d better, sign p = %.4f  (median %.3f -> %.3f)\n",
    ws[j - 1], ws[j], sum(y < x), length(x),
    stats::binom.test(sum(y < x), length(x))$p.value, median(x), median(y)))
}
W <- dcast(r, file + Chr ~ w, value.var = "a")
cat("\nrank stability of a across windows (Spearman):\n")
for (j in seq_along(ws)[-1]) cat(sprintf("  %s vs %s : %.3f   median |change| %.1f%%\n",
  ws[1], ws[j], stats::cor(W[[ws[1]]], W[[ws[j]]], method = "spearman"),
  100 * median(abs(W[[ws[j]]] - W[[ws[1]]]) / W[[ws[1]]])))
