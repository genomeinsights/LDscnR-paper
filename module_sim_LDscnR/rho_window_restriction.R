## =====================================================================
## module_sim_LDscnR / rho_window_restriction.R
##
## Does restricting the C-score's rho window to [0.5, 1] help?
##
## rho sets the width of the ld_w window: d_from_rho(a, rho) = rho / (a(1-rho)),
## so rho = 0.05 measures LD over a very NARROW span and rho = 0.99 over a wide
## one. A narrow window reports immediate physical linkage -- the background LD
## architecture -- whereas a selective sweep leaves a WIDER footprint. So the
## low-rho columns plausibly contribute linkage signal rather than selection
## signal, and dropping them should sharpen C rather than merely shrink it.
##
## The C grid stays bounded on [0,1] either way, so tau_C keeps its reading as
## "significant at 0.05 in at least tau of the analyses" -- it is a fraction of
## a smaller grid, not a different kind of quantity.
##
## V0.5_c2 is excluded: r2_link there is 0.607 (its background LD b = 0.42 drags
## the decay-relative threshold up), PR is 0.014-0.018 for both rules, and it is
## not a regime any conclusion should rest on.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/rho_window_restriction.R
## Env: SIM_DATA, OUT, CELLS, RHO_MIN, LMINS, TAU
## =====================================================================
suppressMessages({library(data.table); library(LDscnR)})
## scoring distance cap. 1e5 matches the bundles' clustering distance_threshold
## (commit 8dbb09a harmonised these); earlier runs used a stale 5e5, which scored
## regions built at 100 kb against a 500 kb truth window.
DCAP <- as.numeric(Sys.getenv("DCAP", "1e5"))
SIM   <- Sys.getenv("SIM_DATA", "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5")
OUT   <- Sys.getenv("OUT", "module_sim_LDscnR/results/bgs5_headtohead")
CELLS <- strsplit(Sys.getenv("CELLS", "V0.5_c1,V1_c1.5,V2_c1"), ",")[[1]]
RMIN  <- as.numeric(Sys.getenv("RHO_MIN", "0.5"))
LMINS <- as.integer(strsplit(Sys.getenv("LMINS", "2,3,5,10"), ",")[[1]])
TAU   <- as.numeric(Sys.getenv("TAU", "0.05"))
QSTAR <- seq(0, 0.95, by = 0.05); RHO_LD <- 0.75; DCAP <- 1e5
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
OUTF <- file.path(OUT, "rho_window_restriction.csv")

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
for (cell in CELLS) for (tag in c("bgs","nobgs")) for (env in 1:10) {
  if (!is.null(done) && nrow(done[cell == ..cell & tag == ..tag & env == ..env])) next
  P <- pool(cell, tag, env); if (is.null(P)) next
  nq <- sum(P$map$true_pos_QTN %in% TRUE); if (!nq) next
  th <- score_thresholds(P$decay_sum, rho_r2 = 0.75, rho_d = 0.95, dmax_cap = DCAP)
  rho_all  <- colnames(P$ld_ws)
  rho_high <- rho_all[as.numeric(sub("rho_", "", rho_all)) >= RMIN]
  rows <- list(); k <- 0L; t0 <- Sys.time()
  for (eng in c("emmax","lfmm")) {
    p <- P$map[[if (eng == "emmax") "emx_p" else "lfmm_p"]]
    if (is.null(p) || all(is.na(p))) next
    q <- p.adjust(p, "BH")
    Cs <- lapply(list(all = rho_all, high = rho_high), function(rr) {
      Ce <- ld_cscore(p, P$ld_ws, alpha = 0.05, rho = rr, qstar = QSTAR)
      cv <- rep(0, nrow(P$map)); names(cv) <- P$map$marker; cv[names(Ce)] <- Ce; cv })
    uni <- unique(c(unlist(lapply(Cs, function(v) names(v)[which(v >= TAU)])),
                    P$map$marker[which(q < 0.05)]))
    if (!length(uni)) next
    ed   <- ld_edges(uni, P$GTs, P$map[, .(marker, Chr, Pos)], P$decay_sum,
                     rho_ld = RHO_LD, dcap = DCAP)
    qtab <- qtn_ld_table(P$GTs, P$map, uni, 2e6, cores = 1)
    ev <- function(mk, L) { if (!length(mk)) return(NULL)
      ra <- ld_regions(mk, ed); ra <- ra[lengths(ra) >= L]
      if (!length(ra)) return(NULL)
      evaluate_ors(ra, P$map, qtab, th$r2min, th$dmax) }
    for (L in LMINS) {
      a  <- ev(P$map$marker[which(q < 0.05)], L)
      cA <- ev(names(Cs$all)[which(Cs$all >= TAU)], L)
      cH <- ev(names(Cs$high)[which(Cs$high >= TAU)], L)
      k <- k + 1L
      rows[[k]] <- data.table(cell, tag, env, engine = eng, l_min = L, n_qtn = nq,
        n_rho_all = length(rho_all), n_rho_high = length(rho_high),
        a_PR  = if (is.null(a))  0 else a$PR,
        Call_PR = if (is.null(cA)) 0 else cA$PR, Call_n = if (is.null(cA)) 0L else cA$TP + cA$FP,
        Chigh_PR = if (is.null(cH)) 0 else cH$PR, Chigh_n = if (is.null(cH)) 0L else cH$TP + cH$FP,
        Call_TP = if (is.null(cA)) 0L else cA$TP, Chigh_TP = if (is.null(cH)) 0L else cH$TP) } }
  if (length(rows)) {
    fwrite(rbindlist(rows), OUTF, append = file.exists(OUTF))
    cat(sprintf("  %-9s %-5s env%-2d | %.1f min\n", cell, tag, env,
        as.numeric(Sys.time() - t0, units = "mins"))); flush.console() } }
d <- fread(OUTF)
cat(sprintf("\n=== C on ALL rho vs rho >= %.2f | %d genomes | tau %.2f ===\n",
            RMIN, uniqueN(d[, .(cell, tag, env)]), TAU))
print(d[, .(n = .N, alpha = round(mean(a_PR), 3),
            C_all = round(mean(Call_PR), 3), C_high = round(mean(Chigh_PR), 3),
            gain = round(mean(Chigh_PR - Call_PR), 3),
            high_wins = sum(Chigh_PR > Call_PR), ties = sum(Chigh_PR == Call_PR)),
        by = l_min][order(l_min)])
cat("\n=== per cell, l_min 5 ===\n")
print(d[l_min == 5, .(n = .N, alpha = round(mean(a_PR), 3),
                      C_all = round(mean(Call_PR), 3), C_high = round(mean(Chigh_PR), 3),
                      reg_all = round(mean(Call_n), 1), reg_high = round(mean(Chigh_n), 1),
                      TP_all = round(mean(Call_TP), 2), TP_high = round(mean(Chigh_TP), 2)),
        by = cell][order(cell)])
