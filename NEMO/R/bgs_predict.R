## ---------------------------------------------------------------------------
## bgs_predict.R -- analytic B = Ne/N0 landscape for the Nemo chromosome maps
##
## Decide the deleterious-trait parameters that give a DETECTABLE background-
## selection effect before spending simulation time.
##
## Model (Hudson & Kaplan 1995; Nordborg, Charlesworth & Charlesworth 1996):
##
##     B(x) = exp( - sum_j  U_j * t_j / (t_j + r_jx)^2 )
##
##   U_j   diploid mutation rate at deleterious locus j (= 2 * delet_mutation_rate;
##         Nemo: _genomic_mut_rate = 2 * nloc * delet_mutation_rate)
##   t_j   heterozygous effect h*s_j  (Nemo _effects[0][j], h = delet_dominance_mean)
##   r_jx  Morgans between locus j and neutral site x
##
## For uniform density mu (diploid del. mutations per Morgan) this integrates to
## E = 2*mu: BGS strength is set by mutational INPUT PER MORGAN, not by s. That
## is why low-recombination regions are hit hardest, and why raising s alone
## changes nothing.
##
## Genetic distance is taken from pos_nemo * genetic_map_resolution -- what Nemo
## actually recombines over -- NOT the cM column of the rds map. The map builder
## (recombination_map_with_LRR.R) sets pos_nemo = trunc(cM * max(bp)/2000), and
## the .ini declares 3.3e-5 cM per unit, so the simulated chromosome is 53.1 cM
## where the real Formica map is 107.7 cM. Using the rds cM would halve E.
## ---------------------------------------------------------------------------

suppressMessages(library(data.table))

## Resolve paths relative to this script, so the folder can be moved or handed on.
script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) normalizePath(dirname(f)) else getwd()
}
NEMO <- Sys.getenv("BGS_NEMO_DIR", dirname(script_dir()))   # R/ -> NEMO/

MAP_DIR  <- Sys.getenv("BGS_MAP_DIR", file.path(NEMO, "params_V3_rds"))
MAP_RES  <- 3.3e-5      # ntrl/delet/quanti_genetic_map_resolution, cM per pos_nemo unit
WIN      <- 5e5         # analysis window, matches the popgen pipeline
N_E      <- 18432       # patch_capacity 8 * patch_number 2304
THIN     <- 100L        # evaluate B at every THIN'th neutral site

scenario <- function(label, n_delet, mu_locus, meanlog, sdlog, h = 0.5)
  list(label = label, n_delet = n_delet, mu_locus = mu_locus,
       meanlog = meanlog, sdlog = sdlog, h = h)

## --- B landscape for one chromosome ----------------------------------------
b_landscape <- function(chr, sc, seed = 1L) {
  set.seed(seed)
  ntrl <- chr[type == "ntrl"]

  ## Deleterious positions: taken from the map itself when it already carries the
  ## right number of loci (params_V4), otherwise idealised as uniform in PHYSICAL
  ## space -- which is how both map generations place them (KS vs uniform-in-bp
  ## p = 0.37, vs uniform-in-cM p = 1e-4).
  d_map <- chr[type == "delet", pos_nemo]
  if (length(d_map) == sc$n_delet) {
    d_M <- d_map * MAP_RES / 100
  } else {
    d_bp <- seq(min(chr$bp), max(chr$bp), length.out = sc$n_delet)
    d_M  <- approx(chr$bp, chr$pos_nemo, xout = d_bp, rule = 2)$y * MAP_RES / 100
  }

  s <- pmin(rlnorm(sc$n_delet, sc$meanlog, sc$sdlog), 1)   # Nemo truncates s > 1
  t <- sc$h * s
  U <- 2 * sc$mu_locus

  ## Drift barrier: the Hudson-Kaplan form diverges as t -> 0, but mutations with
  ## 2*Ne*t < 1 are effectively neutral and cause no BGS. Drop them.
  keep <- t > 1 / (2 * N_E)
  d_M <- d_M[keep]; t <- t[keep]

  idx <- seq(1L, nrow(ntrl), by = THIN)
  x_M <- ntrl$pos_nemo[idx] * MAP_RES / 100

  E <- numeric(length(x_M))
  for (j in seq_along(d_M)) E <- E + U * t[j] / (t[j] + abs(x_M - d_M[j]))^2

  data.table(bp = ntrl$bp[idx], M = x_M, B = exp(-E))
}

window_summary <- function(land, chr) {
  land[, win := floor(bp / WIN)]
  rec <- chr[, .(M_lo = min(pos_nemo), M_hi = max(pos_nemo),
                 bp_lo = min(bp), bp_hi = max(bp)), by = .(win = floor(bp / WIN))]
  rec[, rec_cM_Mb := (M_hi - M_lo) * MAP_RES / ((bp_hi - bp_lo) / 1e6)]
  w <- merge(land[, .(B = mean(B), n = .N), by = win], rec[, .(win, rec_cM_Mb)], by = "win")
  w[is.finite(rec_cM_Mb) & n >= 10]
}

## --- non-equilibrium correction --------------------------------------------
## The burn-in starts monomorphic (ntrl_init_model 0 = "all monomorphic", see
## ttneutralgenes_bitstring.cc:74) and runs BURNIN generations, which is only
## ~0.5 Ne. Diversity is still building, and in that mutation-limited regime
## H(t) ~ 2*u*t is nearly independent of Ne -- so equilibrium B OVERSTATES the
## reduction that will actually be observed. Sanity check on the regime: 2*u*t
## = 2 * 5e-7 * 1e4 = 0.010 genome-wide, and the production runs report
## hs_ntrl = 0.107 over the ~8% of loci passing MAF 0.01, i.e. ~0.0086
## genome-wide. The sims are indeed mutation-limited.
##
##   H(t) = theta * (1 - exp(-t / (2*Ne)))     with theta = 4*Ne*u
##   under BGS, Ne -> Ne*B, so
##   H_bgs/H_nobgs = B*(1 - exp(-t/(2*Ne*B))) / (1 - exp(-t/(2*Ne)))
##
## Ne is taken as the census size, which is conservative: a smaller true Ne puts
## the run closer to equilibrium and makes the contrast LARGER, not smaller.
BURNIN <- as.numeric(Sys.getenv("BGS_BURNIN_GENS", "10000"))

pi_ratio <- function(B, t = BURNIN, Ne = N_E)
  B * (1 - exp(-t / (2 * Ne * B))) / (1 - exp(-t / (2 * Ne)))

report <- function(sc, chrs) {
  ws <- rbindlist(lapply(chrs, function(chr) window_summary(b_landscape(chr, sc), chr)),
                  idcol = "chr")
  ws[, rec_q := cut(rec_cM_Mb, quantile(rec_cM_Mb, seq(0, 1, .2)),
                    include.lowest = TRUE, labels = paste0("Q", 1:5))]
  ws[, pi_rel := pi_ratio(B)]
  q <- ws[, .(rec_cM_Mb = round(median(rec_cM_Mb), 2), B = round(mean(B), 3),
              pi_rel = round(mean(pi_rel), 3), n_win = .N), by = rec_q][order(rec_q)]
  U_chr <- 2 * sc$n_delet * sc$mu_locus
  s_med <- exp(sc$meanlog)
  n_eff <- mean(sc$h * pmin(rlnorm(2e5, sc$meanlog, sc$sdlog), 1) > 1 / (2 * N_E))
  cat(sprintf("\n== %s\n   delet_loci %d/chr  delet_mutation_rate %.3g  lognormal(%.2f, %.2f)\n",
              sc$label, sc$n_delet, sc$mu_locus, sc$meanlog, sc$sdlog))
  cat(sprintf("   U %.4f per chr (diploid); 2 chr -> Haldane load %.1f%%\n",
              U_chr, 100 * (1 - exp(-2 * U_chr))))
  cat(sprintf("   median s %.2e -> h*s %.2e = %.0f x 1/(2N); %.0f%% of loci effectively selected; reach h*s = %.2f cM\n",
              s_med, sc$h * s_med, sc$h * s_med * 2 * N_E, 100 * n_eff, sc$h * s_med * 100))
  print(q)
  cat(sprintf("   equilibrium: mean B %.3f, Q1/Q5 contrast %.1f%%\n",
              mean(ws$B), 100 * (1 - q$B[1] / q$B[5])))
  cat(sprintf("   at %d gens:  mean pi %.3f of nobgs, Q1/Q5 contrast %.1f%%  <- what will be observed\n",
              BURNIN, mean(ws$pi_rel), 100 * (1 - q$pi_rel[1] / q$pi_rel[5])))
  invisible(ws)
}

## --- parameter scan --------------------------------------------------------
## The question this answers: which (delet_loci, delet_mutation_rate, DFE) give a
## pi contrast a single run can actually resolve? The bar is ~5%, set by the 4.5%
## per-window CV of pi. Anything under ~10% is not worth cluster time.

load_maps <- function(dir) lapply(1:10, function(i)
  readRDS(file.path(dir, sprintf("rec_map%d.rds", i)))[Chr == 1][order(bp)])

chrs <- load_maps(file.path(NEMO, "params_V4", "rds"))
cat(sprintf("%d maps | %.1f Mb | Nemo map %.1f cM | mean %.2f cM/Mb | burn-in %d gens\n",
            length(chrs), diff(range(chrs[[1]]$bp)) / 1e6,
            max(chrs[[1]]$pos_nemo) * MAP_RES,
            max(chrs[[1]]$pos_nemo) * MAP_RES / (diff(range(chrs[[1]]$bp)) / 1e6), BURNIN))

## The params_V3 baseline is only needed to show why the old settings failed; those
## rds maps are not shipped with this folder (they live in Nemo_v3/chromosome_maps_500kb_rds).
V3 <- Sys.getenv("BGS_V3_RDS", "")
if (nzchar(V3) && dir.exists(V3))
  report(scenario("CURRENT (bgs2, params_V3)", 100, 2e-7, -9.425, 3.36), load_maps(V3))

## Candidates. Deleterious positions come from the params_V4 maps for the 1000-locus
## rows; the 100-locus rows are evaluated on idealised uniform-in-bp spacing, which is
## how those maps place them anyway.
scan <- list(
  scenario("100 loci, u 3e-5",       100,  3.0e-5, -5.30, 2.00),
  scenario("1000 loci, u 5e-6",     1000,  5.0e-6, -4.60, 1.50),
  scenario("1000 loci, u 1.5e-5 *", 1000,  1.5e-5, -4.60, 1.50),   # chosen
  scenario("1000 loci, u 3e-5",     1000,  3.0e-5, -4.60, 1.50),
  scenario("1000 loci, tight DFE",  1000,  1.5e-5, -4.60, 0.75),
  scenario("1000 loci, broad DFE",  1000,  1.5e-5, -4.60, 2.50)
)

res <- rbindlist(lapply(scan, function(sc) {
  ws <- report(sc, chrs)
  ws[, q := cut(rec_cM_Mb, quantile(rec_cM_Mb, seq(0, 1, .2)),
                include.lowest = TRUE, labels = paste0("Q", 1:5))]
  qq <- ws[, .(p = mean(pi_rel)), by = q][order(q)]
  data.table(setting = sc$label,
             load_pct = round(100 * (1 - exp(-4 * sc$n_delet * sc$mu_locus)), 1),
             B = round(mean(ws$B), 3),
             pi_rel = round(mean(ws$pi_rel), 3),
             contrast_pct = round(100 * (1 - qq$p[1] / qq$p[5]), 1))
}))
res[, detectable := ifelse(contrast_pct > 10, "yes", "MARGINAL")]

cat("\n\n=== scan summary (pi contrast realised at", BURNIN, "generations) ===\n")
print(res)
cat("\n* = the setting written into ini/ by make_run.R\n")
cat("Detection bar ~5% (pi window CV 4.5%); 'yes' requires >10%.\n")
