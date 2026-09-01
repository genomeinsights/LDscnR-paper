## =====================================================================
## module_sim_LDscnR / rho_ld_sweep.R
##
## Is the C-score's disadvantage at low l_min an OVER-SPLITTING artefact of the
## clustering threshold rather than a property of the method?
##
## tau_LD = ld_from_rho(b, c, rho_ld) = b + (c - b)(1 - rho_ld), so HIGHER
## rho_ld gives a LOWER r2 threshold, looser linking, and more merging. The
## manuscript's over-clustering failure mode -- satellite ORs around one QTN,
## each scored a false positive after dedup -- is produced by tau_LD too HIGH,
## i.e. rho_ld too LOW. If the C-score is being penalised for fragmenting the
## same loci BH finds whole, raising rho_ld should close the gap.
##
## The MATCHING thresholds stay fixed at the manuscript values
## (r2_{rho=0.75}, d_{rho=0.95}); only the CLUSTERING threshold sweeps. Sweeping
## both at once would confound "how regions are built" with "how they are
## scored against truth".
##
## Cheap because edges are built on the candidate union (a few hundred markers),
## not the whole genome: pool once per genome, then vary rho_ld inside.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/rho_ld_sweep.R
## Env: SIM_DATA, OUT, CELLS, RHOS, LMINS, ENVS
## =====================================================================
suppressMessages({library(data.table); library(LDscnR)})
## scoring distance cap. 1e5 matches the bundles' clustering distance_threshold
## (commit 8dbb09a harmonised these); earlier runs used a stale 5e5, which scored
## regions built at 100 kb against a 500 kb truth window.
DCAP <- as.numeric(Sys.getenv("DCAP", "1e5"))
SIM   <- Sys.getenv("SIM_DATA", "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5")
OUT   <- Sys.getenv("OUT", "module_sim_LDscnR/results/bgs5_headtohead")
CELLS <- strsplit(Sys.getenv("CELLS", "V0.5_c1,V0.5_c2,V1_c1.5,V2_c1"), ",")[[1]]
RHOS  <- as.numeric(strsplit(Sys.getenv("RHOS", "0.5,0.6,0.75,0.85,0.9,0.95,0.99"), ",")[[1]])
LMINS <- as.integer(strsplit(Sys.getenv("LMINS", "2,3,5,10"), ",")[[1]])
ENVS  <- as.integer(strsplit(Sys.getenv("ENVS", "1,2,3,4,5,6,7,8,9,10"), ",")[[1]])
TAU <- 0.05; QSTAR <- seq(0, 0.95, by = 0.05)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
OUTF <- file.path(OUT, "rho_ld_sweep.csv")

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
  if (!sum(P$map$true_pos_QTN %in% TRUE)) next
  th <- score_thresholds(P$decay_sum, rho_r2 = 0.75, rho_d = 0.95, dmax_cap = DCAP)
  Cvs <- qs <- list()
  for (eng in c("emmax", "lfmm")) {
    pe <- P$map[[if (eng == "emmax") "emx_p" else "lfmm_p"]]
    if (is.null(pe) || all(is.na(pe))) next
    Ce <- ld_cscore(pe, P$ld_ws, alpha = 0.05, rho = colnames(P$ld_ws), qstar = QSTAR)
    cv <- rep(0, nrow(P$map)); names(cv) <- P$map$marker; cv[names(Ce)] <- Ce
    Cvs[[eng]] <- cv; qs[[eng]] <- p.adjust(pe, "BH") }
  if (!length(Cvs)) next
  uni <- unique(unlist(c(lapply(Cvs, function(v) names(v)[which(v >= TAU)]),
                         lapply(qs, function(v) P$map$marker[which(v < 0.05)]))))
  if (!length(uni)) next
  qtab <- qtn_ld_table(P$GTs, P$map, uni, 2e6, cores = 1)
  t0 <- Sys.time(); rows <- list(); k <- 0L
  for (rl in RHOS) {
    ed <- ld_edges(uni, P$GTs, P$map[, .(marker, Chr, Pos)], P$decay_sum,
                   rho_ld = rl, dcap = 5e5)
    for (eng in names(Cvs)) {
      sets <- list(alpha = P$map$marker[which(qs[[eng]] < 0.05)],
                   C     = names(Cvs[[eng]])[which(Cvs[[eng]] >= TAU)])
      for (L in LMINS) {
        v <- lapply(sets, function(mk) { if (!length(mk)) return(NULL)
          ra <- ld_regions(mk, ed); ra <- ra[lengths(ra) >= L]
          if (!length(ra)) return(NULL)
          evaluate_ors(ra, P$map, qtab, th$r2min, th$dmax) })
        k <- k + 1L
        rows[[k]] <- data.table(cell, tag, env, engine = eng, rho_ld = rl, l_min = L,
          n_qtn = sum(P$map$true_pos_QTN %in% TRUE),
          a_n = if (is.null(v$alpha)) 0L else v$alpha$TP + v$alpha$FP,
          a_TP = if (is.null(v$alpha)) 0L else v$alpha$TP,
          a_PR = if (is.null(v$alpha)) 0 else v$alpha$PR,
          C_n = if (is.null(v$C)) 0L else v$C$TP + v$C$FP,
          C_TP = if (is.null(v$C)) 0L else v$C$TP,
          C_PR = if (is.null(v$C)) 0 else v$C$PR) } } }
  fwrite(rbindlist(rows), OUTF, append = file.exists(OUTF))
  cat(sprintf("  %-9s %-5s env%-2d | %d rho x %d l_min | %.1f min\n", cell, tag, env,
      length(RHOS), length(LMINS), as.numeric(Sys.time() - t0, units = "mins")))
  flush.console()
}
d <- fread(OUTF)
cat("\n=== PR by rho_ld, pooled over cells and engines ===\n")
print(d[, .(n = .N, a_PR = round(mean(a_PR, na.rm = TRUE), 3),
            C_PR = round(mean(C_PR, na.rm = TRUE), 3),
            dPR = round(mean(C_PR - a_PR, na.rm = TRUE), 3),
            C_wins = sum(C_PR > a_PR, na.rm = TRUE),
            alpha_wins = sum(C_PR < a_PR, na.rm = TRUE)),
        by = .(rho_ld, l_min)][order(l_min, rho_ld)])
cat("\n=== region counts: does raising rho_ld merge the satellites? ===\n")
print(d[l_min == 2, .(a_regions = round(mean(a_n), 1), C_regions = round(mean(C_n), 1),
                      ratio = round(mean(C_n)/mean(a_n), 2)), by = rho_ld][order(rho_ld)])
