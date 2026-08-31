## =====================================================================
## module_sim_LDscnR / gc_diagnostic_bgs.R
##
## Is the BGS performance loss genuine, or an artefact of genomic control?
##
## Two tests, because the obvious one cannot answer it alone.
##
## TEST A -- sweep the C-score's internal alpha (0.05 .. 0.5). ld_cscore() fixes
## it at 0.05 within every candidate set regardless of the scan's actual
## inflation. If BGS's loss is threshold stringency, raising alpha should recover
## more under bgs than under nobgs.
##
## TEST B -- recompute the p-values WITHOUT genomic control and rerun both arms.
## This is the decisive one, and the reason is a subtlety in test A:
##
##   PR-AUC over a swept threshold is RANK-BASED. GC divides F by lambda and
##   recomputes p -- a monotone transform. A single GLOBAL lambda would leave
##   marker ordering untouched and could not change the alpha arm's PR-AUC at
##   all. But GC is applied PER CHROMOSOME FILE with its own lambda, so across a
##   pooled genome it rescales chromosomes differently and RE-RANKS markers
##   between them. That can move PR-AUC, and should move it more under BGS, where
##   6 of 10 files exceed the 1.1 threshold rather than 2 of 10 and the lambdas
##   are more dispersed.
##
## So: if removing GC recovers performance under bgs, the loss is per-file
## re-ranking. If it does not, BGS is genuinely destroying detectability.
##
## Raw p-values come from emmax_setup() + emmax_fast() on the bundle's own saved
## GRM -- the same path run_sim_nulls.R now uses for both observed and surrogates.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/gc_diagnostic_bgs.R [outdir]
## =====================================================================
suppressMessages({ library(data.table); library(LDscnR) })

a <- commandArgs(trailingOnly = TRUE)
OUT <- if (length(a)) a[1] else "module_sim_LDscnR/results"
SIM_DATA <- Sys.getenv("SIM_DATA", "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs4")
CELLS <- strsplit(Sys.getenv("CELLS", "V0.5_c2_env1,V1_c1.5_env1"), ",")[[1]]
ALPHA_GRID <- as.numeric(strsplit(Sys.getenv("ALPHA_GRID", "0.05,0.1,0.2,0.5"), ",")[[1]])
LMIN <- 3L; QSTAR <- seq(0, 0.95, by = 0.05)
## 1e5, not 5e5: the stage-2 partition in the bundles moved to
## distance_threshold = 1e5 on 2026-08-29, and the scoring geometry has to match
## it or regions are formed on one distance scale and the partition on another.
## It is also load-bearing rather than nominal: d(rho=0.95) is 636-845 kb on these
## cells, so a 5e5 cap BINDS and is what actually sets the window -- region
## formation is cap-governed, not decay-governed, at either value.
RHO_LD <- 0.75; RHO_D <- 0.95; DCAP <- 1e5; MAX_TAU <- 30L
ALPHAS <- sort(unique(c(10^seq(-6, log10(0.5), length.out = 25), 0.05)))
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)

pool <- function(tag, cell) {
  fs <- list.files(SIM_DATA, full.names = TRUE,
                   pattern = sprintf("^adapt_%s_chr[0-9]+_%s[.]rds$", tag, gsub("\\.", "[.]", cell)))
  if (!length(fs)) return(NULL)
  fs <- fs[order(as.integer(sub(".*_chr([0-9]+)_.*", "\\1", basename(fs))))]
  maps <- gts <- ldws <- decs <- prep <- mk_i <- vector("list", length(fs)); Y <- NULL; gifs <- numeric(0)
  for (i in seq_along(fs)) {
    d <- readRDS(fs[i]); m <- as.data.table(d$map)
    m[, `:=`(Chr = paste0("R", i, "_", Chr), marker = paste0("R", i, "_", marker))]
    G <- d$GTs; colnames(G) <- m$marker
    lw <- d$ld_ws; rownames(lw) <- m$marker
    ds <- as.data.table(d$LD_decay$decay_sum); ds[, Chr := paste0("R", i, "_", Chr)]
    prep[[i]] <- emmax_setup(G, d$GRM); mk_i[[i]] <- m$marker
    gifs <- c(gifs, d$emx_gif %||% NA_real_)
    maps[[i]] <- m; gts[[i]] <- G; ldws[[i]] <- lw; decs[[i]] <- ds
    if (is.null(Y)) Y <- d$env$env
  }
  map <- flag_true_qtns(rbindlist(maps, fill = TRUE))
  list(map = map, GTs = do.call(cbind, gts)[, map$marker],
       ld_ws = do.call(rbind, ldws)[map$marker, ], decay_sum = rbindlist(decs, fill = TRUE),
       prep = prep, mk_i = mk_i, Yobs = Y, gifs = gifs)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

res <- list(); k <- 0L
for (cell in CELLS) for (tag in c("bgs", "nobgs")) {
  P <- pool(tag, cell); if (is.null(P)) next
  map <- P$map; n_true <- sum(map$true_pos_QTN %in% TRUE); if (!n_true) next
  th <- score_thresholds(P$decay_sum, rho_r2 = RHO_LD, rho_d = RHO_D, dmax_cap = DCAP)

  ## the two p-value versions: stored (per-file GC applied when gif > 1.1) and raw
  p_gc  <- map$emx_p
  p_raw <- unlist(lapply(seq_along(P$prep), function(i)
    stats::setNames(emmax_fast(P$prep[[i]], P$Yobs), P$mk_i[[i]])))[map$marker]
  n_gc_fired <- sum(P$gifs > 1.1, na.rm = TRUE)
  cat(sprintf("\n== %s / %-5s | %d QTN | GC fired on %d of %d chromosomes | lambda %.3f-%.3f\n",
              cell, tag, n_true, n_gc_fired, length(P$gifs), min(P$gifs), max(P$gifs))); flush.console()

  for (pv_name in c("gc", "raw")) {
    p <- if (pv_name == "gc") p_gc else p_raw
    q <- stats::p.adjust(p, "BH")
    for (ac in ALPHA_GRID) {
      C <- ld_cscore(p, P$ld_ws, alpha = ac, rho = colnames(P$ld_ws), qstar = QSTAR)
      if (!any(C > 0)) next
      uni <- unique(c(names(C)[which(C > 0)], map$marker[which(q < max(ALPHAS))]))
      edges <- ld_edges(uni, P$GTs, map[, .(marker, Chr, Pos)], P$decay_sum, rho_ld = RHO_LD, dcap = DCAP)
      qtab  <- qtn_ld_table(P$GTs, map, uni, 2e6, cores = 1)
      sc <- function(mk) { if (!length(mk)) return(NULL)
        r <- ld_regions(mk, edges); r <- r[lengths(r) >= LMIN]
        ev <- if (length(r)) evaluate_ors(r, map, qtab, th$r2min, th$dmax) else list(Precision=NA_real_, Recall=0)
        data.table(precision = ev$Precision, recall = ev$Recall) }
      tv <- sort(unique(C[C > 0])); if (length(tv) > MAX_TAU) tv <- as.numeric(stats::quantile(tv, seq(0,1,length.out=MAX_TAU), type=1))
      prC <- rbindlist(lapply(unique(tv), function(t) sc(names(C)[which(C >= t)])))
      prA <- rbindlist(lapply(ALPHAS, function(al) sc(map$marker[which(q < al)])))
      k <- k + 1L
      res[[k]] <- data.table(cell, tag, pvals = pv_name, cscore_alpha = ac, n_true,
        n_gc_fired, n_Cgt0 = sum(C > 0),
        PR_AUC_C = tryCatch(pr_auc(prC[!is.na(precision)]$recall, prC[!is.na(precision)]$precision), error=function(e) NA_real_),
        PR_AUC_alpha = tryCatch(pr_auc(prA[!is.na(precision)]$recall, prA[!is.na(precision)]$precision), error=function(e) NA_real_))
      cat(sprintf("   %-3s alpha_C %.2f : C %.3f | alpha-arm %.3f | C>0 %d\n", pv_name, ac,
                  res[[k]]$PR_AUC_C, res[[k]]$PR_AUC_alpha, sum(C > 0))); flush.console()
    }
  }
}
out <- rbindlist(res); fwrite(out, file.path(OUT, "gc_diagnostic_bgs.csv"))

cat("\n=== TEST A: does raising the C-score's alpha recover BGS losses? (GC'd p-values) ===\n")
print(dcast(out[pvals == "gc"], cell + tag ~ cscore_alpha, value.var = "PR_AUC_C"))
cat("\n=== TEST B: does REMOVING genomic control recover them? (alpha_C = 0.05) ===\n")
print(dcast(out[cscore_alpha == 0.05], cell + tag ~ pvals, value.var = c("PR_AUC_C", "PR_AUC_alpha")))
cat("\n  alpha-arm PR-AUC is RANK-based: a global GC could not change it.\n")
cat("  Any gc-vs-raw difference there is per-file re-ranking across chromosomes.\n")
