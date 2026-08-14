## module_sim/01_score_pooled.R
## Candidate-first + background-null caller scored against the QTN ground truth
## on a pooled 20-chromosome genome. This is THE canonical sim benchmark.
## Args (optional): V c env   (default V1 c2 env3 -- the hard, structure-limited
## case; use "2 1 1" for the easy, well-behaved case).
## Run from LDscnR-paper/:  Rscript module_sim/01_score_pooled.R [V c env]
## Output (git-ignored): module_sim/score_V{V}_c{c}_env{env}.rds

source("module_sim/_config.R")
a   <- commandArgs(trailingOnly = TRUE)
V   <- if (length(a) >= 1) a[1] else "1"
cc  <- if (length(a) >= 2) a[2] else "2"
env <- if (length(a) >= 3) a[3] else "3"
pat <- sprintf("^adapt_bgs_chr[0-9]+_V%s_c%s_env%s\\.rds$", V, cc, env)

P   <- pool_group(pat)
map <- flag_true_positive_QTNs(as.data.table(P$map))
th  <- score_thresholds(P$decay$decay_sum)
cat(sprintf("V%s_c%s_env%s pooled: %d SNPs, %d chr (%d QTN-bearing) | true_pos_QTN=%d | r2min=%.2f dmax=%.0fkb\n",
    V, cc, env, nrow(map), uniqueN(map$Chr), uniqueN(map[type == "QTN", Chr]),
    map[true_pos_QTN == TRUE, .N], th$r2min, th$dmax / 1e3))

score <- function(pcol, lab) {
  r <- ld_outlier_clusters(setNames(map[[pcol]], map$marker), P$ldw[map$marker],
       map[, .(marker, Chr, Pos)], P$GTs, null = "background",
       rho_ld = 0.95, rho_d = 0.95, LD_decay = P$decay,
       rmsc_grid = seq(0, 0.99, 0.005), B = 1000, cores = 4, verbose = FALSE)
  sig  <- as.data.table(r$clusters)[significant == TRUE]
  qtab <- precompute_QTN_LD(P$GTs, map, r$candidates$marker, 2e6, cores = 4)
  ev   <- if (nrow(sig)) evaluate_ORs(lapply(sig$members, identity), map, qtab, th$r2min, th$dmax)
          else list(TP = 0, FP = 0, Precision = NA_real_, Recall = 0)
  ## FP location: neutral chromosome (R*_Chr2) vs QTN-chromosome non-QTN cluster
  nfp_neutral <- if (nrow(sig)) {
    ass <- classify_ORs(lapply(sig$members, identity), map, qtab, th$r2min, th$dmax)
    sig[, `:=`(isTP = ass$is_TP, neu = grepl("_Chr2$", Chr))]
    sig[isTP == FALSE & neu == TRUE, .N]
  } else 0L
  dg <- r$diagnostics
  cat(sprintf("[%-5s] q*=%.3f no_elbow=%s cand=%d sig=%d | TP=%d FP=%d (neutral FP=%d) Prec=%.2f Rec=%.2f\n",
      lab, r$ld_w_threshold, dg$no_elbow, nrow(r$candidates), nrow(sig),
      ev$TP, ev$FP, nfp_neutral, ev$Precision, ev$Recall))
  cat(sprintf("        DIAG ldw_tracks_structure=%s (%.0f%% chr) | whole_chr_clusters=%s (median span=%.2f)\n",
      dg$ldw_tracks_structure, 100 * dg$frac_chr_ldw_positive,
      dg$whole_chr_clusters, dg$median_span_fraction))
  saveRDS(list(clusters = as.data.table(r$clusters), candidates = r$candidates,
               qstar = r$ld_w_threshold, ev = ev, diagnostics = dg),
          file.path(mod, sprintf("score_V%s_c%s_env%s_%s.rds", V, cc, env, lab)))
  ev
}

cat("=== candidate-first + background null ===\n")
invisible(score("emx_p",  "EMMAX"))
invisible(score("lfmm_p", "LFMM"))
saveRDS(list(map = map, decay = P$decay$decay_sum),
        file.path(mod, sprintf("meta_V%s_c%s_env%s.rds", V, cc, env)))
