## =============================================================================
## filter_then_test_clusters.R -- filter-then-test at STAGE-2 CLUSTER level, on
## a pooled genome-wide panel. The cluster-level counterpart of
## filter_then_test.R (which works on 2-chromosome marker tracks).
##
## WHY CLUSTER LEVEL AND WHY POOLED. The marker-level result showed ld_w beats
## every same-size filter but only beats a genome-wide scan from k ~ 1000 up.
## The suspicion is that the gain scales with the multiplicity burden: the
## tracks carry ~30k markers over 2 chromosomes, whereas the 3sp panel session
## saw a large gain over 108,857 clusters genome-wide. Pooling the 10 map sets
## per (cell, tag, env) gives ~130k clusters -- the same order -- so this tests
## that explanation rather than assuming it.
##
## UNITS are stage-2 clusters, each represented by its pruned marker; the
## cluster's p-value is that marker's. Stage-2 membership is recomputed per
## bundle and REQUIRED to reproduce the stored grm_markers exactly.
##
## MEMORY. Genotypes are needed only for the stage-2 recomputation, so each
## bundle is reduced to a small per-cluster table and then discarded. Peak
## memory is one bundle per worker, not ten.
##
## TRUTH is bp distance to the nearest driving QTN (MAF > 0.1, p_Va > 0.05) with
## the causal markers removed from the tested set -- no shared machinery with ld_w.
##
## Env: SIM_DATA, CELLS, ENVS, FILES, KS, WINDOWS, ALPHA, ENGINE, CORES, OUT, LDAGG
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(parallel)})

SIM   <- Sys.getenv("SIM_DATA", "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5")
OUT   <- Sys.getenv("OUT", "module_sim_LDscnR/results/filter_then_test")
CELLS <- strsplit(Sys.getenv("CELLS", "V0.5_c1,V0.5_c2,V1_c1.5,V2_c1"), ",")[[1]]
ENVS  <- as.integer(strsplit(Sys.getenv("ENVS", "1,2,3,4,5,6,7,8,9,10"), ",")[[1]])
FILES <- as.integer(strsplit(Sys.getenv("FILES", "1,2,3,4,5,6,7,8,9,10"), ",")[[1]])
KS    <- as.integer(strsplit(Sys.getenv("KS", "500,1000,2000,5000,10000,20000"), ",")[[1]])
WIN   <- as.numeric(strsplit(Sys.getenv("WINDOWS", "10,50,100"), ",")[[1]]) * 1000
ALPHA <- as.numeric(Sys.getenv("ALPHA", "0.05"))
ENG   <- Sys.getenv("ENGINE", "emx")
CORES <- as.integer(Sys.getenv("CORES", "1"))
LDAGG <- Sys.getenv("LDAGG", "median")            # how a cluster's ld_w is summarised
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
pcol <- paste0(ENG, "_p")

## ---- one bundle -> one small per-cluster table; genotypes dropped on exit ----
units_for <- function(cell, tag, env, i) {
  f <- sprintf("%s/adapt_%s_chr%d_%s_env%d.rds", SIM, tag, i, cell, env)
  if (!file.exists(f)) return(NULL)
  x <- readRDS(f)
  m <- as.data.table(x$map)
  if (!pcol %in% names(m)) return(NULL)
  pr <- ld_prune_and_eMLG(GTs = x$GTs, stage1 = x$complexity_reduction$stage1,
          LD_decay = x$LD_decay, ld_w_col = "ld_w_095", ld_w_threshold = 0.025,
          score_threshold = 0.80, min_r2_rho = 0.5, distance_threshold = 1e5,
          compute_unflagged_eMLG = FALSE, cores = 1)
  if (!identical(sort(pr$pruned), sort(x$grm_markers)))
    stop(sprintf("stage-2 does not reproduce grm_markers: %s", basename(f)))
  g  <- as.data.table(pr$groups)
  ms <- rbindlist(lapply(seq_len(nrow(g)), function(k)
          data.table(marker = g$members[[k]], CL_id = g$group_id[k])))
  m <- merge(m, ms, by = "marker", all.x = TRUE)[!is.na(CL_id)]

  ## driving QTN, then distance from every marker to the nearest one
  drv <- m[true_QTN %in% TRUE & MAF > 0.1 & p_Va > 0.05]
  m[, d_qtn := Inf]; m[, near_qtn := NA_character_]
  if (nrow(drv)) for (ch in unique(m$Chr)) {
    dd <- drv[Chr == ch]; if (!nrow(dd)) next
    ii <- which(m$Chr == ch)
    D  <- abs(outer(m$Pos[ii], dd$Pos, "-")); j <- max.col(-D)
    m[ii, `:=`(d_qtn = D[cbind(seq_along(ii), j)],
               near_qtn = paste0(cell, "_", i, "_", ch, "_", dd$Pos[j]))]
  }
  m <- m[!(true_QTN %in% TRUE)]                        # causal variants excluded
  if (!nrow(m)) return(NULL)

  rep_mk <- intersect(pr$pruned, m$marker)
  u <- m[, .(ld_w   = if (LDAGG == "max") max(ld_w_095, na.rm = TRUE) else median(ld_w_095, na.rm = TRUE),
             MAF    = median(MAF, na.rm = TRUE),
             n_loci = .N,
             d_qtn  = min(d_qtn), near_qtn = near_qtn[which.min(d_qtn)]), by = CL_id]
  rp <- m[marker %in% rep_mk, .(CL_id, p = get(pcol))][, .SD[1], by = CL_id]
  u  <- merge(u, rp, by = "CL_id")[is.finite(p)]
  if (!nrow(u)) return(NULL)
  u[, `:=`(CL_id = paste0(i, "_", CL_id), file = i)]
  u[]
}

score_sel <- function(u, idx, w) {
  s <- u[idx]
  if (!nrow(s)) return(data.table(n_sel = 0L, n_sig = 0L, n_tp = 0L, precision = NA_real_, qtn_hit = 0L))
  q <- p.adjust(s$p, "BH"); sig <- which(q < ALPHA)
  if (!length(sig)) return(data.table(n_sel = nrow(s), n_sig = 0L, n_tp = 0L, precision = NA_real_, qtn_hit = 0L))
  h <- s[sig]; tp <- h$d_qtn < w
  data.table(n_sel = nrow(s), n_sig = length(sig), n_tp = sum(tp),
             precision = mean(tp), qtn_hit = uniqueN(h$near_qtn[tp]))
}

grid <- CJ(cell = CELLS, tag = c("nobgs", "bgs"), env = ENVS, sorted = FALSE)
cat(sprintf("  %d panels (%d cells x 2 arms x %d envs), %d files each, CORES=%d\n",
            nrow(grid), length(CELLS), length(ENVS), length(FILES), CORES))

one_panel <- function(r) {
  cell <- grid$cell[r]; tag <- grid$tag[r]; env <- grid$env[r]
  u <- tryCatch(rbindlist(lapply(FILES, function(i) units_for(cell, tag, env, i)), fill = TRUE),
                error = function(e) { message("  FAIL ", cell, " ", tag, " env", env, ": ", conditionMessage(e)); NULL })
  if (is.null(u) || !nrow(u)) return(NULL)
  set.seed(9000 + r)
  sels <- list(ld_w = order(-u$ld_w), MAF = order(-u$MAF),
               size = order(-u$n_loci), random = sample.int(nrow(u)))
  meta <- data.table(cell, tag, env, n_units = nrow(u),
                     n_qtn = uniqueN(u$near_qtn[is.finite(u$d_qtn)]))
  out <- rbindlist(lapply(WIN, function(w) {
    base <- cbind(meta, method = "genome_wide", k = nrow(u), window_kb = w/1000,
                  score_sel(u, seq_len(nrow(u)), w))
    per_k <- rbindlist(lapply(KS, function(kk) {
      if (kk >= nrow(u)) return(NULL)
      rbindlist(lapply(names(sels), function(nm)
        cbind(meta, method = nm, k = kk, window_kb = w/1000,
              score_sel(u, head(sels[[nm]], kk), w))))
    }))
    rbind(base, per_k, fill = TRUE)
  }))
  cat(sprintf("    %-9s %-5s env%-3d %6d units, %2d QTN\n", cell, tag, env, nrow(u), meta$n_qtn))
  out
}

res <- if (CORES > 1) {
  mclapply(seq_len(nrow(grid)), one_panel, mc.cores = CORES, mc.preschedule = FALSE)
} else {
  lapply(seq_len(nrow(grid)), one_panel)
}
all <- rbindlist(Filter(Negate(is.null), res), fill = TRUE)
stopifnot(nrow(all) > 0)
fwrite(all, file.path(OUT, sprintf("filter_then_test_clusters_%s.csv", ENG)))
cat(sprintf("\n  written: %s  (%d rows)\n", file.path(OUT, sprintf("filter_then_test_clusters_%s.csv", ENG)), nrow(all)))
