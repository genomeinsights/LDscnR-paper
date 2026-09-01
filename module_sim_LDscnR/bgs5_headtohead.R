## =====================================================================
## module_sim_LDscnR / bgs5_headtohead.R
##
## alpha = 0.05 against tau_C = 0.05 on the bgs5 bundles: 4 cells x 10
## environments x 2 arms x 2 engines. NO SURROGATES NEEDED -- the comparison
## uses only p_obs, ld_ws, the edge list and truth, all of which the bundles
## already carry. Only the gate and q_R need a null, and those are answered
## separately.
##
## This supersedes the V2_c1-only null run for the C-vs-alpha question: that set
## is ONE cell, ONE arm (nobgs -- it has zero deleterious loci), and a different
## generation from the PR-AUC benchmark. This runs on the same bundles as the
## benchmark, so the two compose.
##
## l_min is held IDENTICAL between the arms at each row, so the region machinery
## -- pooling, edges, clustering, truth rule -- is the same and the ONLY
## difference is q < 0.05 against C >= 0.05. Paired within genome.
##
## Truth is containment of a flag_true_qtns() QTN in the region span, matching
## analyse_one_dataset.R.
##
## Run from the LDscnR-paper root (long: edge-building dominates):
##   Rscript module_sim_LDscnR/bgs5_headtohead.R
## Env: SIM_DATA, OUT, CELLS, TAU, LMINS
## =====================================================================
suppressMessages({library(data.table); library(LDscnR)})

SIM   <- Sys.getenv("SIM_DATA", "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5")
OUT   <- Sys.getenv("OUT", "module_sim_LDscnR/results/bgs5_headtohead")
CELLS <- strsplit(Sys.getenv("CELLS", "V0.5_c1,V0.5_c2,V1_c1.5,V2_c1"), ",")[[1]]
TAU   <- as.numeric(Sys.getenv("TAU", "0.05"))
LMINS <- as.integer(strsplit(Sys.getenv("LMINS", "2,3,5,10,20"), ",")[[1]])
QSTAR <- seq(0, 0.95, by = 0.05)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

pool <- function(cell, tag, env) {
  ff <- list.files(SIM, full.names = TRUE, pattern = sprintf(
    "^adapt_%s_chr[0-9]+_%s_env%d[.]rds$", tag, gsub("\\.", "[.]", cell), env))
  if (!length(ff)) return(NULL)
  ff <- ff[order(as.integer(sub(".*_chr([0-9]+)_.*", "\\1", basename(ff))))]
  maps <- gts <- lws <- dss <- vector("list", length(ff))
  for (i in seq_along(ff)) {
    d <- readRDS(ff[i]); m <- as.data.table(d$map)
    m[, `:=`(Chr = paste0("R", i, "_", Chr), marker = paste0("R", i, "_", marker))]
    G <- d$GTs; colnames(G) <- m$marker
    lw <- d$ld_ws; rownames(lw) <- m$marker
    ds <- as.data.table(d$LD_decay$decay_sum); ds[, Chr := paste0("R", i, "_", Chr)]
    maps[[i]] <- m; gts[[i]] <- G; lws[[i]] <- lw; dss[[i]] <- ds
  }
  map <- flag_true_qtns(rbindlist(maps, fill = TRUE))
  list(map = map, GTs = do.call(cbind, gts)[, map$marker],
       ld_ws = do.call(rbind, lws)[map$marker, ], decay_sum = rbindlist(dss, fill = TRUE))
}

## INCREMENTAL AND RESUMABLE. At ~7 min a genome this is a multi-hour run, and
## accumulating in memory means a crash at hour 9 loses everything. Each genome's
## rows are appended as soon as they exist, and a genome already present in the
## csv is skipped, so the job can be killed and restarted without losing work.
PART <- file.path(OUT, "bgs5_headtohead.csv")
done <- if (file.exists(PART)) unique(fread(PART)[, .(cell, tag, env)]) else NULL
if (!is.null(done)) cat(sprintf("  resuming: %d genomes already done\n", nrow(done)))
appendrow <- function(dt) fwrite(dt, PART, append = file.exists(PART))

res <- list(); k <- 0L; skipped <- list()
for (cell in CELLS) for (tag in c("bgs", "nobgs")) for (env in 1:10) {
  if (!is.null(done) && nrow(done[cell == ..cell & tag == ..tag & env == ..env])) {
    cat(sprintf("  [skip] %s %s env%d -- already in csv\n", cell, tag, env)); next }
  P <- pool(cell, tag, env)
  if (is.null(P)) { skipped[[length(skipped)+1L]] <-
    data.table(cell, tag, env, why = "no files"); next }
  qtn <- P$map[true_pos_QTN %in% TRUE, .(Chr = as.character(Chr), Pos)]
  if (!nrow(qtn)) { skipped[[length(skipped)+1L]] <-
    data.table(cell, tag, env, why = "no detectable QTN"); next }
  t0 <- Sys.time()
  ## Edges are built on the UNION OF THE TWO CANDIDATE SETS, not the whole
  ## genome. ld_regions() only ever uses edges among the markers handed to it,
  ## and r2 between two markers does not depend on what else is in the set, so
  ## the restricted edge list gives byte-identical regions. It is also 5000x
  ## faster: the union is ~0.3% of markers (726 of 228,643 on one genome) and
  ## ld_edges drops from 501 s to 0.1 s, which is 97% of the per-genome runtime.
  ## Verified identical on both rules before adopting.
  ## NOTE: the union must cover EVERY selection this script makes. If a looser
  ## tau or alpha is ever added, widen it or the regions will be wrong.
  Cvs <- list(); qs <- list()
  for (eng in c("emmax", "lfmm")) {
    pe <- P$map[[if (eng == "emmax") "emx_p" else "lfmm_p"]]
    if (is.null(pe) || all(is.na(pe))) next
    Ce <- ld_cscore(pe, P$ld_ws, alpha = 0.05, rho = colnames(P$ld_ws), qstar = QSTAR)
    cv <- rep(0, nrow(P$map)); names(cv) <- P$map$marker; cv[names(Ce)] <- Ce
    Cvs[[eng]] <- cv; qs[[eng]] <- p.adjust(pe, "BH")
  }
  if (!length(Cvs)) { skipped[[length(skipped)+1L]] <-
    data.table(cell, tag, env, why = "no p-values for either engine"); next }
  uni <- unique(unlist(c(
    lapply(Cvs, function(v) names(v)[which(v >= TAU)]),
    lapply(qs,  function(v) P$map$marker[which(v < 0.05)]))))
  if (!length(uni)) { skipped[[length(skipped)+1L]] <-
    data.table(cell, tag, env, why = "no candidate markers"); next }
  ed <- ld_edges(uni, P$GTs, P$map[, .(marker, Chr, Pos)],
                 P$decay_sum, rho_ld = 0.75, dcap = 5e5)
  ## CANONICAL SCORING: qtn_ld_table() + evaluate_ors() from the package, called
  ## exactly as bgs_vs_nobgs_prauc.R and grm_comparison.R call them. Focal-QTN
  ## assignment, the r2/distance match and the dedup all live there.
  ##
  ## ONE DIFFERENCE FROM THE MANUSCRIPT, reported rather than silently resolved.
  ## evaluate_ors() scores dedup losers that point at a real QTN as NEITHER TP
  ## nor FP -- it calls them `extras` and drops them from the denominator:
  ##     FP <- sum(!is_TP & !extra)
  ## The manuscript figure states FP = total ORs - TP, i.e. those satellite
  ## regions ARE false positives. The two differ whenever a locus fragments, and
  ## the difference inflates precision. Both are recorded: *_prec / *_PR follow
  ## evaluate_ors, *_precS / *_PRS use FP = n - TP.
  th_   <- score_thresholds(P$decay_sum, rho_r2 = 0.75, rho_d = 0.95, dmax_cap = 5e5)
  qtab  <- qtn_ld_table(P$GTs, P$map, uni, 2e6, cores = 1)
  n_true <- sum(P$map$true_pos_QTN %in% TRUE)
  sc <- function(mk, L) {
    if (!length(mk)) return(c(0L, 0L, 0L, 0L))
    ra <- ld_regions(mk, ed); ra <- ra[lengths(ra) >= L]
    if (!length(ra)) return(c(0L, 0L, 0L, 0L))
    ev <- evaluate_ors(ra, P$map, qtab, th_$r2min, th_$dmax)
    c(length(ra), ev$TP, ev$FP, ev$extras)
  }
  for (eng in names(Cvs)) {
    Cv <- Cvs[[eng]]; q <- qs[[eng]]
    rows <- rbindlist(lapply(LMINS, function(L) {
      a <- sc(P$map$marker[which(q < 0.05)], L); c_ <- sc(names(Cv)[which(Cv >= TAU)], L)
      data.table(cell, tag, env, engine = eng, l_min = L,
                 n_qtn = nrow(qtn), n_markers = nrow(P$map),
                 alpha_n = a[1], alpha_tp = a[2], alpha_fp = a[3], alpha_extra = a[4],
                 C_n = c_[1], C_tp = c_[2], C_fp = c_[3], C_extra = c_[4]) }))
    appendrow(rows); k <- k + 1L; res[[k]] <- rows
  }
  cat(sprintf("  %-9s %-5s env%-2d | %d QTN, %s markers | %.1f min\n", cell, tag, env,
      nrow(qtn), format(nrow(P$map), big.mark=","),
      as.numeric(Sys.time() - t0, units = "mins"))); flush.console()
}
out <- if (file.exists(PART)) fread(PART) else rbindlist(res)
if (length(skipped)) { fwrite(rbindlist(skipped), file.path(OUT, "bgs5_headtohead_SKIPPED.csv"))
  cat(sprintf("\n*** %d genome/engine combinations skipped, see the SKIPPED csv ***\n",
              length(skipped))) }
cat(sprintf("\ncoverage: %d rows from %d genomes (expected %d)\n",
            nrow(out), uniqueN(out[, .(cell, tag, env)]),
            length(CELLS) * 2L * 10L * 2L * length(LMINS)))
## evaluate_ors convention (extras excluded from FP)
out[, `:=`(a_prec = fifelse(alpha_tp + alpha_fp > 0, alpha_tp/(alpha_tp + alpha_fp), NA_real_),
           C_prec = fifelse(C_tp + C_fp > 0, C_tp/(C_tp + C_fp), NA_real_),
           a_rec  = alpha_tp/n_qtn, C_rec = C_tp/n_qtn)]
## manuscript convention: FP = total ORs - TP, so satellites count against you
out[, `:=`(a_precS = fifelse(alpha_n > 0, alpha_tp/alpha_n, NA_real_),
           C_precS = fifelse(C_n > 0, C_tp/C_n, NA_real_))]
out[, `:=`(a_PR = a_prec*a_rec, C_PR = C_prec*C_rec,
           a_PRS = a_precS*a_rec, C_PRS = C_precS*C_rec)]
cat("\n=== ABSOLUTE performance, manuscript OR scoring (l_min 5) ===\n")
print(out[l_min == 5, .(n = .N, QTN = round(mean(n_qtn), 1),
      a_reg = round(mean(alpha_n), 1), a_TP = round(mean(alpha_tp), 2),
      a_rec = paste0(round(100*mean(a_rec)), "%"), a_prec = round(mean(a_prec, na.rm=TRUE), 2),
      a_precS = round(mean(a_precS, na.rm=TRUE), 2), a_PR = round(mean(a_PR, na.rm=TRUE), 3),
      C_reg = round(mean(C_n), 1), C_TP = round(mean(C_tp), 2),
      C_rec = paste0(round(100*mean(C_rec)), "%"), C_prec = round(mean(C_prec, na.rm=TRUE), 2),
      C_precS = round(mean(C_precS, na.rm=TRUE), 2), C_PR = round(mean(C_PR, na.rm=TRUE), 3),
      C_extra = round(mean(C_extra), 1)), by = cell][order(cell)])
cat("\n=== head to head under the strict rule ===\n")
print(out[, .(n = .N, a_PR = round(mean(a_PR, na.rm=TRUE), 3),
              C_PR = round(mean(C_PR, na.rm=TRUE), 3),
              dPR = round(mean(C_PR - a_PR, na.rm=TRUE), 3),
              C_wins = sum(C_PR > a_PR, na.rm=TRUE), ties = sum(C_PR == a_PR, na.rm=TRUE),
              alpha_wins = sum(C_PR < a_PR, na.rm=TRUE)), by = l_min][order(l_min)])
