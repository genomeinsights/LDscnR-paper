## =====================================================================
## module_sim_LDscnR / tau_C_sweep.R
##
## Sweeps tau_C against a FIXED alpha = 0.05, scored by the manuscript OR rule.
##
## The comparison is deliberately ASYMMETRIC and that is the point: in practice
## nobody sweeps alpha, they run BH at 0.05 and report what survives. So the
## question is not "which curve is higher" but "is there a tau_C at which the
## C-score beats the conventional rule, and is it a value anyone could pick in
## advance". tau_C = 0.05 has the principled reading -- significant at 0.05 in
## at least 5% of the (rho, q*) analyses, both axes bounded on [0,1] -- so if
## the optimum sits far from it, the principled value is not the useful one.
##
## The alpha arm is recomputed at every tau so the pairing is exact, though it
## does not vary; the repetition is what makes each row a matched comparison.
##
## rho_ld = 0.9 (the sweep optimum) and the manuscript matching thresholds
## (r2_{rho=0.75}, d_{rho=0.95}) are held fixed throughout.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/tau_C_sweep.R
## Env: SIM_DATA, OUT, CELLS, TAUS, LMINS, ENVS, RHO_LD
## =====================================================================
suppressMessages({library(data.table); library(LDscnR)})
## scoring distance cap. 1e5 matches the bundles' clustering distance_threshold
## (commit 8dbb09a harmonised these); earlier runs used a stale 5e5, which scored
## regions built at 100 kb against a 500 kb truth window.
DCAP <- as.numeric(Sys.getenv("DCAP", "1e5"))
SIM   <- Sys.getenv("SIM_DATA", "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5")
OUT   <- Sys.getenv("OUT", "module_sim_LDscnR/results/bgs5_headtohead")
CELLS <- strsplit(Sys.getenv("CELLS", "V0.5_c1,V0.5_c2,V1_c1.5,V2_c1"), ",")[[1]]
TAUS  <- as.numeric(strsplit(Sys.getenv("TAUS",
           "0.01,0.02,0.05,0.10,0.15,0.20,0.30,0.40,0.50,0.70"), ",")[[1]])
LMINS <- as.integer(strsplit(Sys.getenv("LMINS", "2,3,5,10"), ",")[[1]])
ENVS  <- as.integer(strsplit(Sys.getenv("ENVS", "1,2,3,4,5,6,7,8,9,10"), ",")[[1]])
RHO_LD <- as.numeric(Sys.getenv("RHO_LD", "0.9"))
QSTAR <- seq(0, 0.95, by = 0.05)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
OUTF <- file.path(OUT, "tau_C_sweep.csv")

pool <- function(cell, tag, env) {
  ff <- list.files(SIM, full.names = TRUE, pattern = sprintf(
    "^adapt_%s_chr[0-9]+_%s_env%d[.]rds$", tag, gsub("\\.", "[.]", cell), env))
  if (!length(ff)) return(NULL)
  ff <- ff[order(as.integer(sub(".*_chr([0-9]+)_.*", "\\1", basename(ff))))]
  maps <- gts <- lws <- dss <- vector("list", length(ff))
  for (i in seq_along(ff)) { d <- readRDS(ff[i]); m <- as.data.table(d$map)
    m[, `:=`(Chr = paste0("R", i, "_", Chr), marker = paste0("R", i, "_", marker))]
    G <- d$GTs; colnames(G) <- m$marker; lw <- d$ld_ws; rownames(lw) <- m$marker
    ds <- as.data.table(d$LD_decay$decay_sum); ds[, Chr := paste0("R", i, "_", Chr)]
    maps[[i]] <- m; gts[[i]] <- G; lws[[i]] <- lw; dss[[i]] <- ds }
  map <- flag_true_qtns(rbindlist(maps, fill = TRUE))
  list(map = map, GTs = do.call(cbind, gts)[, map$marker],
       ld_ws = do.call(rbind, lws)[map$marker, ], decay_sum = rbindlist(dss, fill = TRUE)) }

done <- if (file.exists(OUTF)) unique(fread(OUTF)[, .(cell, tag, env)]) else NULL
for (cell in CELLS) for (tag in c("bgs", "nobgs")) for (env in ENVS) {
  if (!is.null(done) && nrow(done[cell == ..cell & tag == ..tag & env == ..env])) next
  P <- pool(cell, tag, env); if (is.null(P)) next
  nq <- sum(P$map$true_pos_QTN %in% TRUE); if (!nq) next
  th <- score_thresholds(P$decay_sum, rho_r2 = 0.75, rho_d = 0.95, dmax_cap = DCAP)
  Cvs <- qs <- list()
  for (eng in c("emmax", "lfmm")) {
    pe <- P$map[[if (eng == "emmax") "emx_p" else "lfmm_p"]]
    if (is.null(pe) || all(is.na(pe))) next
    Ce <- ld_cscore(pe, P$ld_ws, alpha = 0.05, rho = colnames(P$ld_ws), qstar = QSTAR)
    cv <- rep(0, nrow(P$map)); names(cv) <- P$map$marker; cv[names(Ce)] <- Ce
    Cvs[[eng]] <- cv; qs[[eng]] <- p.adjust(pe, "BH") }
  if (!length(Cvs)) next
  ## union spans the WIDEST tau (lowest) so every threshold's selection is inside it
  uni <- unique(unlist(c(lapply(Cvs, function(v) names(v)[which(v >= min(TAUS))]),
                         lapply(qs,  function(v) P$map$marker[which(v < 0.05)]))))
  if (!length(uni)) next
  ed   <- ld_edges(uni, P$GTs, P$map[, .(marker, Chr, Pos)], P$decay_sum,
                   rho_ld = RHO_LD, dcap = 5e5)
  qtab <- qtn_ld_table(P$GTs, P$map, uni, 2e6, cores = 1)
  ev <- function(mk, L) { if (!length(mk)) return(NULL)
    ra <- ld_regions(mk, ed); ra <- ra[lengths(ra) >= L]
    if (!length(ra)) return(NULL)
    evaluate_ors(ra, P$map, qtab, th$r2min, th$dmax) }
  t0 <- Sys.time(); rows <- list(); k <- 0L
  for (eng in names(Cvs)) for (L in LMINS) {
    a <- ev(P$map$marker[which(qs[[eng]] < 0.05)], L)
    for (tau in TAUS) {
      cc <- ev(names(Cvs[[eng]])[which(Cvs[[eng]] >= tau)], L)
      k <- k + 1L
      rows[[k]] <- data.table(cell, tag, env, engine = eng, l_min = L, tau = tau, n_qtn = nq,
        a_n = if (is.null(a)) 0L else a$TP + a$FP, a_TP = if (is.null(a)) 0L else a$TP,
        a_PR = if (is.null(a)) 0 else a$PR,
        C_n = if (is.null(cc)) 0L else cc$TP + cc$FP, C_TP = if (is.null(cc)) 0L else cc$TP,
        C_PR = if (is.null(cc)) 0 else cc$PR) } }
  fwrite(rbindlist(rows), OUTF, append = file.exists(OUTF))
  cat(sprintf("  %-9s %-5s env%-2d | %d tau x %d l_min | %.1f min\n", cell, tag, env,
      length(TAUS), length(LMINS), as.numeric(Sys.time() - t0, units = "mins")))
  flush.console() }
d <- fread(OUTF)
cat(sprintf("\n=== tau_C sweep vs FIXED alpha 0.05 | rho_ld %.2f | %d genomes ===\n",
            RHO_LD, uniqueN(d[, .(cell, tag, env)])))
print(d[, .(n = .N, C_PR = round(mean(C_PR, na.rm = TRUE), 3),
            a_PR = round(mean(a_PR, na.rm = TRUE), 3),
            dPR = round(mean(C_PR - a_PR, na.rm = TRUE), 3),
            C_wins = sum(C_PR > a_PR, na.rm = TRUE)), by = .(l_min, tau)][order(l_min, tau)])
cat("\n=== best tau per cell (l_min 5), and what tau = 0.05 gives there ===\n")
b <- d[l_min == 5, .(C_PR = mean(C_PR, na.rm = TRUE)), by = .(cell, tau)]
print(merge(b[, .SD[which.max(C_PR)], by = cell][, .(cell, best_tau = tau, best_PR = round(C_PR, 3))],
            b[tau == 0.05, .(cell, PR_at_0.05 = round(C_PR, 3))], by = "cell"))
