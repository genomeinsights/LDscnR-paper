## module_sim/05_distance_cap.R
## Is the whole-chromosome clustering (structure regime) fixable by a tighter
## distance cap? Compare the decay-derived cap (rho_d=0.95 ~ 2.2Mb here) with a
## fixed 500kb cap on the same pooled genome, reporting TP/FP AND the built-in
## structure diagnostics.
## FINDING: the fixed cap DOES fragment the whole-chromosome clusters (median
## span 1.00 -> ~0.5) but precision gets WORSE (one big FP per chromosome becomes
## several small FP sub-clusters), and ldw_tracks_structure / frac_chr_ldw_positive
## are IDENTICAL at both caps -- they don't depend on clustering geometry. The
## distance cap is not the lever; the structure is. Confirms the diagnostics, not
## the cap, are what flags this regime.
## Run from LDscnR-paper/:  Rscript module_sim/05_distance_cap.R [V c env]
## Output (git-ignored): module_sim/distance_cap_V{V}_c{c}_env{env}.rds

source("module_sim/_config.R")
a   <- commandArgs(trailingOnly = TRUE)
V   <- if (length(a) >= 1) a[1] else "1"
cc  <- if (length(a) >= 2) a[2] else "2"
env <- if (length(a) >= 3) a[3] else "3"

P   <- pool_group(sprintf("^adapt_bgs_chr[0-9]+_V%s_c%s_env%s\\.rds$", V, cc, env))
map <- flag_true_positive_QTNs(as.data.table(P$map))
th  <- score_thresholds(P$decay$decay_sum)
cat(sprintf("V%s_c%s_env%s pooled: %d SNPs, %d chr | r2min=%.2f dmax(score)=%.0fkb\n",
    V, cc, env, nrow(map), uniqueN(map$Chr), th$r2min, th$dmax / 1e3))

## cap = NULL uses rho_d (decay-derived, wide); cap = 5e5 is the fixed 500kb test
run <- function(pcol, lab, cap) {
  args <- list(setNames(map[[pcol]], map$marker), P$ldw[map$marker],
               map[, .(marker, Chr, Pos)], P$GTs, null = "background",
               r2_link = 0.5, B = 1000, cores = 4, verbose = FALSE)
  if (is.null(cap)) { args$rho_ld <- 0.95; args$rho_d <- 0.95; args$LD_decay <- P$decay }
  else             { args$distance_threshold <- cap }
  r    <- do.call(ld_outlier_clusters, args)
  sig  <- as.data.table(r$clusters)[significant == TRUE]
  qtab <- precompute_QTN_LD(P$GTs, map, r$candidates$marker, 2e6, cores = 4)
  ev   <- if (nrow(sig)) evaluate_ORs(lapply(sig$members, identity), map, qtab, th$r2min, th$dmax)
          else list(TP = 0, FP = 0, Precision = NA_real_, Recall = 0)
  dg   <- r$diagnostics
  cat(sprintf("[%-5s | %-8s] q*=%.3f sig=%d | TP=%d FP=%d Prec=%.2f Rec=%.2f | span=%.2f whole_chr=%s | ldw_struct=%s(%.0f%%)\n",
      lab, if (is.null(cap)) "rho_d.95" else "500kb", r$ld_w_threshold, nrow(sig),
      ev$TP, ev$FP, ev$Precision, ev$Recall, dg$median_span_fraction, dg$whole_chr_clusters,
      dg$ldw_tracks_structure, 100 * dg$frac_chr_ldw_positive))
  list(ev = ev, diagnostics = dg)
}

cat("=== distance-cap comparison (structure diagnostics are cap-invariant) ===\n")
res <- list(
  emx_rho  = run("emx_p",  "EMMAX", NULL),
  emx_500  = run("emx_p",  "EMMAX", 5e5),
  lfmm_rho = run("lfmm_p", "LFMM",  NULL),
  lfmm_500 = run("lfmm_p", "LFMM",  5e5))
saveRDS(res, file.path(mod, sprintf("distance_cap_V%s_c%s_env%s.rds", V, cc, env)))
