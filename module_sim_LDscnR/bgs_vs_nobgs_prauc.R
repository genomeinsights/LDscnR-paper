## =====================================================================
## module_sim_LDscnR / bgs_vs_nobgs_prauc.R
##
## Does background selection change the relative performance of the two engines,
## and of the C-score against plain BH? And how much does it cost overall?
##
## FOUR COMBINATIONS per genome: {EMMAX, LFMM} x {C-score, alpha}, each swept over
## its own knob (tau_C / BH alpha), clustered and scored identically, so the only
## differences are the association engine and the selection rule.
##
## No permutations are needed. The C-score comes from observed p-values and ld_w;
## the sweep is scored against truth. Surrogates are only needed for q_R, which
## this does not use.
##
## regen_sim_data_bgs4 holds MATCHED bgs and nobgs bundles for V0.5_c2 and
## V1_c1.5 -- same pipeline, same run, same environment vector -- which is the
## controlled contrast. They are NOT paired at the marker level (independent
## simulations under the two regimes: ~1.2k of ~36k markers shared, and different
## QTN counts), so this compares regimes, not the same genome with BGS added.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/bgs_vs_nobgs_prauc.R [outdir]
## Env: SIM_DATA (default regen_sim_data_bgs4), CELLS, LMINS (default 1,3),
##      CORES (default 1) -- forks the tau/alpha sweeps and qtn_ld_table.
##      The sweeps are the wall-clock: each threshold is an independent
##      ld_regions + evaluate_ors over one marker subset. Forking shares the
##      pooled GTs copy-on-write, so memory is roughly flat in CORES, but a
##      pooled genome is ~2 GB resident before forking -- watch it alongside
##      other jobs rather than setting this to every core you have.
## =====================================================================
suppressMessages({ library(data.table); library(LDscnR) })
`%||%` <- function(a, b) if (is.null(a)) b else a

a <- commandArgs(trailingOnly = TRUE)
OUT <- if (length(a)) a[1] else "module_sim_LDscnR/results"
SIM_DATA <- Sys.getenv("SIM_DATA", "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs4")
CELLS <- strsplit(Sys.getenv("CELLS", "V0.5_c2_env1,V1_c1.5_env1"), ",")[[1]]
LMINS <- as.integer(strsplit(Sys.getenv("LMINS", "1,3"), ",")[[1]])
CORES <- max(1L, as.integer(Sys.getenv("CORES", "1")))

## mclapply returns try-error OBJECTS in place of failed elements rather than
## raising, so a partial failure would flow into rbindlist and silently corrupt
## a PR-AUC. Check explicitly and stop instead.
.lapply <- function(X, FUN) {
  if (CORES <= 1L) return(lapply(X, FUN))
  r <- parallel::mclapply(X, FUN, mc.cores = CORES)
  bad <- vapply(r, inherits, logical(1), "try-error")
  if (any(bad)) stop("mclapply failed on ", sum(bad), " of ", length(r),
                     " elements; first: ", conditionMessage(attr(r[[which(bad)[1]]], "condition")))
  r
}
ALPHA_C <- 0.05; QSTAR <- seq(0, 0.95, by = 0.05)
RHO_LD <- 0.75; RHO_D <- 0.95; DCAP <- 5e5; MAX_TAU <- 40L
ALPHAS <- sort(unique(c(10^seq(-6, log10(0.5), length.out = 30), 0.05)))
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)

pool <- function(tag, cell) {
  fs <- list.files(SIM_DATA, full.names = TRUE,
                   pattern = sprintf("^adapt_%s_chr[0-9]+_%s[.]rds$", tag, gsub("\\.", "[.]", cell)))
  if (!length(fs)) return(NULL)
  fs <- fs[order(as.integer(sub(".*_chr([0-9]+)_.*", "\\1", basename(fs))))]
  maps <- gts <- ldws <- decs <- vector("list", length(fs))
  for (i in seq_along(fs)) {
    d <- readRDS(fs[i]); m <- as.data.table(d$map)
    m[, `:=`(Chr = paste0("R", i, "_", Chr), marker = paste0("R", i, "_", marker))]
    G <- d$GTs; colnames(G) <- m$marker
    lw <- d$ld_ws; rownames(lw) <- m$marker
    ds <- as.data.table(d$LD_decay$decay_sum); ds[, Chr := paste0("R", i, "_", Chr)]
    maps[[i]] <- m; gts[[i]] <- G; ldws[[i]] <- lw; decs[[i]] <- ds
  }
  map <- flag_true_qtns(rbindlist(maps, fill = TRUE))
  list(map = map, GTs = do.call(cbind, gts)[, map$marker],
       ld_ws = do.call(rbind, ldws)[map$marker, ], decay_sum = rbindlist(decs, fill = TRUE),
       n_files = length(fs))
}

res <- list(); k <- 0L
for (cell in CELLS) for (tag in c("bgs", "nobgs")) {
  P <- pool(tag, cell); if (is.null(P)) { cat(sprintf("  [skip] %s / %s: no files\n", tag, cell)); next }
  map <- P$map; n_true <- sum(map$true_pos_QTN %in% TRUE)
  if (!n_true) { cat(sprintf("  [skip] %s / %s: no detectable QTN\n", tag, cell)); next }
  th <- score_thresholds(P$decay_sum, rho_r2 = RHO_LD, rho_d = RHO_D, dmax_cap = DCAP)

  for (eng in c("emmax", "lfmm")) {
    pcol <- if (eng == "emmax") "emx_p" else "lfmm_p"
    p <- map[[pcol]]; if (is.null(p) || all(is.na(p))) next
    t0 <- Sys.time()
    C <- ld_cscore(p, P$ld_ws, alpha = ALPHA_C, rho = colnames(P$ld_ws), qstar = QSTAR)
    q <- stats::p.adjust(p, "BH")
    uni <- unique(c(names(C)[which(C > 0)], map$marker[which(q < max(ALPHAS))]))
    if (!length(uni)) next
    edges <- ld_edges(uni, P$GTs, map[, .(marker, Chr, Pos)], P$decay_sum, rho_ld = RHO_LD, dcap = DCAP)
    qtab  <- qtn_ld_table(P$GTs, map, uni, 2e6, cores = CORES)

    sc <- function(mk) { if (!length(mk)) return(NULL)
      ra <- ld_regions(mk, edges)
      rbindlist(lapply(LMINS, function(L) { r <- ra[lengths(ra) >= L]
        ev <- if (length(r)) evaluate_ors(r, map, qtab, th$r2min, th$dmax)
              else list(Precision = NA_real_, Recall = 0)
        data.table(l_min = L, precision = ev$Precision, recall = ev$Recall) })) }

    tv <- sort(unique(C[C > 0]))
    if (length(tv) > MAX_TAU) tv <- as.numeric(stats::quantile(tv, seq(0, 1, length.out = MAX_TAU), type = 1))
    prC <- rbindlist(.lapply(unique(tv), function(t) { s <- sc(names(C)[which(C >= t)]); if (!is.null(s)) s[, knob := t]; s }))
    prA <- rbindlist(.lapply(ALPHAS,     function(al){ s <- sc(map$marker[which(q < al)]);  if (!is.null(s)) s[, knob := al]; s }))
    for (L in LMINS) {
      k <- k + 1L
      res[[k]] <- data.table(cell, tag, engine = eng, l_min = L, n_true = n_true,
        n_markers = nrow(map), n_Cgt0 = sum(C > 0), n_BH05 = sum(q < 0.05, na.rm = TRUE),
        PR_AUC_C     = tryCatch(pr_auc(prC[l_min == L & !is.na(precision)]$recall,
                                       prC[l_min == L & !is.na(precision)]$precision), error = function(e) NA_real_),
        PR_AUC_alpha = tryCatch(pr_auc(prA[l_min == L & !is.na(precision)]$recall,
                                       prA[l_min == L & !is.na(precision)]$precision), error = function(e) NA_real_))
    }
    cat(sprintf("  %-14s %-5s %-5s | %d QTN, %d markers, C>0 %d | %.1f min\n",
                cell, tag, eng, n_true, nrow(map), sum(C > 0),
                as.numeric(Sys.time() - t0, units = "mins"))); flush.console()
  }
}
out <- rbindlist(res)
fwrite(out, file.path(OUT, "bgs_vs_nobgs_prauc.csv"))
out[, C_minus_alpha := round(PR_AUC_C - PR_AUC_alpha, 3)]
cat(sprintf("\n=== PR-AUC: four combinations, bgs vs nobgs (CORES=%d) ===\n", CORES))
print(out[, .(cell, tag, engine, l_min, n_true,
              C = round(PR_AUC_C, 3), alpha = round(PR_AUC_alpha, 3), C_minus_alpha)])
cat("\n=== effect of BGS on overall performance (mean over cells) ===\n")
print(dcast(out, engine + l_min ~ tag, value.var = "PR_AUC_C", fun.aggregate = function(z) round(mean(z, na.rm=TRUE),3)))
