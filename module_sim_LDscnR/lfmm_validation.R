## =====================================================================
## module_sim_LDscnR / lfmm_validation.R
##
## Validate the LFMM-without-a-structure-null recipe on the sims (where truth is
## known). LFMM has no cheap structure null, so its C-score threshold tau_C is
## transferred from EMMAX (which does have one, via the spatial surrogate null) by
## genomic control on the C-score: gc_map_tauc matches LFMM's C-score TAIL MASS to
## EMMAX's null-calibrated tau_C. That assumes LFMM is as powerful as EMMAX (same
## number of true outliers); if LFMM is genuinely more powerful, tail-matching
## OVER-corrects and discards real hits. The sims settle it.
##
## Reports, replicate-averaged over env1-5 (truth = true_pos_QTN):
##   * PR-AUC of the C-score under EMMAX vs LFMM p-values (equal power => similar);
##   * whether LFMM's extra tail markers are TRUE QTN (power) or NEUTRAL (inflation);
##   * the gc_map_tauc tau_C for LFMM vs LFMM's TRUTH-OPTIMAL tau_C -- does the
##     transfer land at a good operating point, over-correct, or under-correct?
##
## Run from the LDscnR-paper root (use the local cache + cores):
##   SIM_DATA=module_sim_LDscnR/data/cache SIM_CORES=5 \
##     Rscript module_sim_LDscnR/lfmm_validation.R [V c]
## =====================================================================

suppressMessages({ library(data.table); library(LDscnR) })

a  <- commandArgs(trailingOnly = TRUE)
V  <- if (length(a) >= 1) a[1] else "2"
CC <- if (length(a) >= 2) a[2] else "1"
SIM_DATA <- Sys.getenv("SIM_DATA", "/Volumes/Nemo/Nemo_sim/regen_sim_data_nobgs")
TAG <- "nobgs"; OUTRES <- "module_sim_LDscnR/results"

## ENVS is derived from what is on disk, not hardcoded, so adding replicates does
## not mean editing this file. Only COMPLETE cells are used: pooling needs all
## CHR_N chromosomes of an env, and a partial cell would silently pool a smaller
## genome. Override with SIM_ENVS=1,2,3.
CHR_N <- 10L
discover_envs <- function(dir, tag, V, CC, chr_n = CHR_N) {
  fs <- list.files(dir, pattern = sprintf("^adapt_%s_chr[0-9]+_V%s_c%s_env[0-9]+[.]rds$", tag, V, CC))
  if (!length(fs)) stop("no bundles matching V", V, "_c", CC, " (", tag, ") in ", dir)
  e <- as.integer(sub(".*_env([0-9]+)[.]rds$", "\\1", fs)); tab <- table(e)
  full <- as.integer(names(tab)[tab == chr_n]); short <- setdiff(as.integer(names(tab)), full)
  if (length(short)) message(sprintf("  [envs] incomplete cells skipped: %s",
    paste(sprintf("env%d (%d/%d chr)", short, as.integer(tab[as.character(short)]), chr_n), collapse = ", ")))
  sort(full)
}
ENVS <- { .e <- Sys.getenv("SIM_ENVS", "")
          if (nzchar(.e)) as.integer(strsplit(.e, ",")[[1]]) else discover_envs(SIM_DATA, TAG, V, CC) }
if (!length(ENVS)) stop("no complete env cells in ", SIM_DATA)
message(sprintf("  [envs] using %d env cell(s): %s", length(ENVS), paste(ENVS, collapse = ",")))
CORES  <- as.integer(Sys.getenv("SIM_CORES", "1")); QTAB_C <- if (CORES > 1L) 1L else 4L
if (CORES > 1L) { data.table::setDTthreads(1L); Sys.setenv(OMP_NUM_THREADS = "1") }  # avoid fork thread oversubscription
PAR <- list(qstar = seq(0, 0.95, by = 0.05), alpha = c(0.001, 0.01, 0.05, 0.1),
            rho_ld = 0.75, rho_d = 0.95, dcap = 1e5, B = 100L, seed = 1L,
            lmin = c(1L, 2L, 4L), fdr = 0.05, lmin_q = 0.99, lmin_tau = 0.05,
            tau_grid = seq(0.02, 1, by = 0.02))
gcta_grm <- function(X) { p <- colMeans(X)/2; k <- p>0 & p<1; X <- X[,k,drop=FALSE]; p <- p[k]
  Z <- sweep(sweep(X,2,2*p,"-"),2,sqrt(2*p*(1-p)),"/"); tcrossprod(Z)/ncol(Z) }

## pool a (V,c,env): per-file EMMAX engine on the saved GRM, pooled ld_ws/map/lfmm_p
pool_cell <- function(env) {
  files <- list.files(SIM_DATA, full.names = TRUE,
                      pattern = sprintf("^adapt_%s_chr[0-9]+_V%s_c%s_env%s\\.rds$", TAG, V, CC, env))
  files <- files[order(as.integer(sub(".*_chr([0-9]+)_.*", "\\1", basename(files))))]
  maps <- gts <- ldws <- decs <- prep <- mk_i <- vector("list", length(files)); Yobs <- coords <- NULL
  for (i in seq_along(files)) {
    d <- readRDS(files[i]); m <- as.data.table(d$map)
    m[, `:=`(Chr = paste0("R", i, "_", Chr), marker = paste0("R", i, "_", marker))]
    G <- d$GTs; colnames(G) <- m$marker; lw <- d$ld_ws; rownames(lw) <- m$marker
    ds <- as.data.table(d$LD_decay$decay_sum); ds[, Chr := paste0("R", i, "_", Chr)]
    prep[[i]] <- emmax_setup(G, d$GRM); mk_i[[i]] <- m$marker
    maps[[i]] <- m; gts[[i]] <- G; ldws[[i]] <- lw; decs[[i]] <- ds
    if (is.null(Yobs)) { Yobs <- d$env$env; coords <- cbind(d$env$x, d$env$y) }
  }
  map <- flag_true_qtns(rbindlist(maps, fill = TRUE))
  list(map = map, GTs = do.call(cbind, gts)[, map$marker], LDW = do.call(rbind, ldws)[map$marker, ],
       decay_sum = rbindlist(decs, fill = TRUE), prep = prep, mk_i = mk_i, Yobs = Yobs, coords = coords)
}
emmax_pooled <- function(P, Yv) unlist(lapply(seq_along(P$prep), function(i)
  stats::setNames(emmax_fast(P$prep[[i]], Yv), P$mk_i[[i]])))[P$map$marker]
cscore <- function(pv, LDW) ld_cscore(pv, LDW, alpha = PAR$alpha, rho = colnames(LDW), qstar = PAR$qstar)

## EMMAX structured null (pooled spatial surrogate), cached
null_bundle <- function(P, env) {
  cache <- file.path(OUTRES, sprintf("nullE_V%s_c%s_env%s.rds", V, CC, env))
  if (file.exists(cache)) return(readRDS(cache))
  n <- length(P$Yobs); Dm <- as.matrix(stats::dist(P$coords)); l <- stats::median(Dm[lower.tri(Dm)])
  eK <- eigen(exp(-0.5*(Dm/l)^2), symmetric = TRUE); Lv <- pmax(eK$values, 0); Vk <- eK$vectors
  surr <- function() as.numeric(stats::resid(stats::lm((Vk %*% (sqrt(Lv)*stats::rnorm(n))) ~ P$Yobs)))
  spC <- function(pv) { C <- cscore(pv, P$LDW); C[C > 0] }
  C_obs <- cscore(emmax_pooled(P, P$Yobs), P$LDW)
  set.seed(PAR$seed); C_surr <- vector("list", PAR$B)
  for (b in seq_len(PAR$B)) { C_surr[[b]] <- spC(emmax_pooled(P, surr())); if (b %% 20 == 0) { cat(sprintf(" env%d null surrogate %d/%d\n", env, b, PAR$B)); utils::flush.console() } }
  bundle <- structure(list(C_obs = C_obs, C_surr = C_surr,
    universe = unique(c(names(C_obs)[C_obs > 0], unlist(lapply(C_surr, names)))),
    basis = "spatial", B = PAR$B), class = "ld_null")
  saveRDS(bundle, cache); bundle
}

per_env <- function(env) {
  P <- pool_cell(env)
  th <- score_thresholds(P$decay_sum, rho_r2 = PAR$rho_ld, rho_d = PAR$rho_d, dmax_cap = PAR$dcap)
  nb <- null_bundle(P, env)
  C_emx  <- nb$C_obs
  C_lfmm <- cscore(P$map$lfmm_p, P$LDW); names(C_lfmm) <- P$map$marker
  uni <- unique(c(nb$universe, names(C_lfmm)[C_lfmm > 0]))
  edges <- ld_edges(uni, P$GTs, P$map[, .(marker, Chr, Pos)], P$decay_sum, rho_ld = PAR$rho_ld, dcap = PAR$dcap)
  qtab <- qtn_ld_table(P$GTs, P$map, uni, 2e6, cores = QTAB_C)
  ## EMMAX null-calibrated tau_C, transferred to LFMM by C-score genomic control
  op_l  <- calibrate_lmin(nb, edges, tau = PAR$lmin_tau, q = PAR$lmin_q)
  tau_e <- calibrate_tauc(nb, edges, l_min = op_l, fdr = PAR$fdr, tau_grid = PAR$tau_grid)
  tau_l <- if (is.na(tau_e)) NA_real_ else gc_map_tauc(tau_e, C_emx, C_lfmm)
  ## evaluate a C-vector at a fixed l_min across the tau sweep -> PR curve + an operating point
  evalC <- function(C, lmin, tau_op) {
    sc <- function(t) { mk <- names(C)[C >= t]; reg <- if (length(mk)) ld_regions(mk, edges) else list()
      evaluate_ors(reg[lengths(reg) >= lmin], P$map, qtab, th$r2min, th$dmax) }
    curve <- rbindlist(lapply(PAR$tau_grid, function(t) { ev <- sc(t)
      data.table(tau = t, recall = ev$Recall, precision = ev$Precision, TP = ev$TP, FP = ev$FP) }))
    auc <- pr_auc(curve$recall, curve$precision)
    ## truth-optimal tau on this curve (max F1), and the operating point at tau_op
    f1 <- with(curve, ifelse(is.finite(precision) & (precision + recall) > 0,
                             2 * precision * recall / (precision + recall), 0))
    opt <- curve[which.max(f1)]
    at  <- if (is.finite(tau_op)) { e <- sc(tau_op)
      data.table(tau = tau_op, recall = e$Recall, precision = e$Precision, TP = e$TP, FP = e$FP) } else curve[NA][1]
    list(auc = auc, opt = opt, at = at)
  }
  e_emx  <- evalC(C_emx,  op_l, tau_e)
  e_lf   <- evalC(C_lfmm, op_l, tau_l)                 # LFMM at the GC-transferred tau
  ## is LFMM's extra tail true or neutral? markers with C >= tau_e (a common low bar)
  lf_tail <- names(C_lfmm)[C_lfmm >= tau_e]
  tp_frac <- if (length(lf_tail)) mean(lf_tail %in% P$map[true_pos_QTN == TRUE | max_LD_with_QTN > th$r2min, marker]) else NA_real_
  cat(sprintf("env%d: tau_e=%.3f l_min=%d ; tau_lfmm(GC)=%.3f ; LFMM opt tau=%.3f | PR-AUC emx=%.3f lfmm=%.3f\n",
      env, tau_e, op_l, tau_l, e_lf$opt$tau, e_emx$auc, e_lf$auc)); utils::flush.console()
  data.table(env = env, l_min = op_l,
             prauc_emx = e_emx$auc, prauc_lfmm = e_lf$auc,
             tau_emx = tau_e, tau_lfmm_gc = tau_l, tau_lfmm_opt = e_lf$opt$tau,
             lfmm_gc_prec = e_lf$at$precision, lfmm_gc_rec = e_lf$at$recall, lfmm_gc_TP = e_lf$at$TP, lfmm_gc_FP = e_lf$at$FP,
             lfmm_opt_prec = e_lf$opt$precision, lfmm_opt_rec = e_lf$opt$recall,
             lfmm_tail_true_frac = tp_frac)
}

cat(sprintf("V%s_c%s: LFMM GC validation over env %s\n", V, CC, paste(ENVS, collapse = ",")))
res <- rbindlist(if (CORES > 1L) parallel::mclapply(ENVS, per_env, mc.cores = CORES) else lapply(ENVS, per_env))
fwrite(res, file.path(OUTRES, "lfmm_validation.csv"))
cat("\n=== per-env ===\n"); print(res)
cat("\n=== means (env1-5) ===\n")
print(res[, lapply(.SD, function(x) round(mean(x, na.rm = TRUE), 3)),
          .SDcols = c("prauc_emx","prauc_lfmm","tau_emx","tau_lfmm_gc","tau_lfmm_opt",
                      "lfmm_gc_prec","lfmm_gc_rec","lfmm_opt_prec","lfmm_opt_rec","lfmm_tail_true_frac")])
