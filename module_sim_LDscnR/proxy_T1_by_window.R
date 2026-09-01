## =============================================================================
## module_sim_LDscnR / proxy_T1_by_window.R
##
## T1 scored at four window counts, against TRUE recombination. Three questions
## in one run:
##
##   (i)  WHICH WINDOW COUNT, ON AN EXTERNAL CRITERION. decay_window_test.R
##        chose 20 on an internal consistency measure that turned out to be
##        arithmetic (two equal-size chromosomes make a_pred their mean). This
##        scores each window count against the simulation's own rec_rate, which
##        is external and needs no such assumption.
##   (ii) 2c's PREDICTION. They see 0 of 180 windows with a == 0 on the panel
##        against my 17 of 363, but their windows are ~6.3 Mb wide against my
##        ~2 Mb, wide enough to average a cold block together with its warm
##        surroundings. If a == 0 is a RESOLUTION effect rather than a dataset
##        difference, my own a == 0 fraction must RISE as windows narrow.
##  (iii) WHETHER THE CAP IS A CONSTANT. If the divergent-window fraction moves
##        with window count, a cap dimensioned once does not transfer between
##        window settings and the two have to be set together.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/proxy_T1_by_window.R
## Env: SIM_DATA, CELL, TAG, ENV, FILES, WINDOWS, CORES, OUT
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(parallel)})
SIM   <- Sys.getenv("SIM_DATA", "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5")
OUT   <- Sys.getenv("OUT", "module_sim_LDscnR/results/operating_points")
CELL  <- Sys.getenv("CELL", "V0.5_c1")
TAG   <- Sys.getenv("TAG", "nobgs")
ENV   <- as.integer(Sys.getenv("ENV", "1"))
FILES <- as.integer(strsplit(Sys.getenv("FILES", "1,2,3,4,5,6,7,8,9,10"), ",")[[1]])
WIN   <- as.integer(strsplit(Sys.getenv("WINDOWS", "5,10,20,50"), ",")[[1]])
CORES <- as.integer(Sys.getenv("CORES", "10"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
ARGS <- list(el_data_folder = NULL, min_maf_decay = 0.1, q = 0.95, n_sub_bg = 5000,
             overlap = 0.5, max_SNPs_decay = Inf, prob_robust = 0.95, max_pairs = 5000,
             ld_method = "corr", n_strata = 20, keep_el = FALSE, slide = 1000,
             rho_targets = c(0.99), cores = 1)

one <- function(i) {
  f <- sprintf("%s/adapt_%s_chr%d_%s_env%d.rds", SIM, TAG, i, CELL, ENV)
  if (!file.exists(f)) return(NULL)
  x  <- readRDS(f); m <- as.data.table(x$map)
  gp <- tempfile(fileext = ".gds"); on.exit(unlink(gp), add = TRUE)
  g  <- create_gds_from_geno(geno = x$GTs, map = m, gp)
  rbindlist(lapply(WIN, function(w) {
    set.seed(1L)
    D <- do.call(compute_LD_decay, c(list(gds = g), ARGS, list(n_win_decay = w)))
    rbindlist(lapply(names(D$by_chr), function(ch) {
      d  <- as.data.table(D$by_chr[[ch]]$decay)
      mm <- m[Chr == ch][order(Pos)]
      tru <- vapply(seq_len(nrow(d)), function(k) {
        s <- mm[Pos >= d$start[k] & Pos <= d$end[k]]
        if (nrow(s) < 2) return(NA_real_)
        (max(s$cM) - min(s$cM)) / ((d$end[k] - d$start[k]) / 1e6) }, numeric(1))
      data.table(w = w, file = i, Chr = ch, a = d$a, c = d$c, contrast = d$contrast,
                 regime = d$regime, cMMb = tru, span = d$end - d$start) })) })) }

W <- rbindlist(Filter(Negate(is.null), mclapply(FILES, one, mc.cores = CORES)))
W <- W[is.finite(a) & is.finite(cMMb)]
fwrite(W, file.path(OUT, "proxy_T1_by_window_raw.csv"))

## Per-chromosome Spearman, then summarised across chromosomes -- windows overlap,
## so a pooled p-value would be optimistic; the sign test across chromosomes is not.
R <- W[, .(rho = suppressWarnings(stats::cor(a, cMMb, method = "spearman")), n = .N),
       by = .(w, file, Chr)][is.finite(rho)]
S <- R[, .(n_chr = .N, n_win_per_chr = as.integer(median(n)),
           rho_med = median(rho), rho_q25 = quantile(rho, .25), rho_q75 = quantile(rho, .75),
           n_pos = sum(rho > 0),
           sign_p = stats::binom.test(sum(rho > 0), .N)$p.value), by = w][order(w)]
Z <- W[, .(span_kb_med = median(span) / 1e3, pct_a_zero = 100 * mean(a == 0),
           n_a_zero = sum(a == 0), n_win = .N,
           d95_max_kb = suppressWarnings(max(d_from_rho(a[a > 0], 0.95)) / 1e3)), by = w][order(w)]
S <- merge(S, Z, by = "w")
fwrite(S, file.path(OUT, "proxy_T1_by_window_summary.csv"))
print(S)
cat("\ncovered chromosomes per window count (a drop would break the comparison):\n")
print(W[, .(n_chr = uniqueN(paste0(file, Chr))), by = w][order(w)])
