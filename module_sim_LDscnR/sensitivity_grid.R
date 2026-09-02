## =============================================================================
## module_sim_LDscnR / sensitivity_grid.R
##
## The factorial sensitivity grid (TODO_sensitivity.md), run on mini2.
##
##   ld_w_threshold   0.0125 / 0.025 / 0.05   stage 2 -- needs re-clustering
##   min_r2_rho       0.35   / 0.5   / 0.65   stage 2 -- needs re-clustering
##   size floor       1      / 2     / 5      post-filter -- free
##
## 27 cells for the price of 9 re-clusterings, since the floor is applied to an
## existing partition. Both engines are scored off the SAME partition (only the
## p-column changes), so the test axis costs nothing on top.
##
## ADMISSIBILITY WAS CHECKED FIRST (admissible_region.R). At distance_threshold
## 1e5 the three ld_w levels are 95-100% admissible on bgs5, so no cell has to be
## dropped -- but ld_w = 0.10 was excluded on that evidence, because 27.5% of
## chromosomes fall below two flagged clusters there and leave the analysis
## without an error. n_runs is recorded per cell so degeneracy is visible in the
## output rather than inferred.
##
## TWO AXES, NEVER COMBINED (TODO_sensitivity.md):
##   SELECTION  C over the 27 parameter cells at a fixed engine -- stability of
##              the units and the filter.
##   TEST       EMMAX vs LFMM at fixed selection -- stability of the association.
## A single C over all 54 is not interpretable: a region can be stable on one
## axis and unstable on the other and the scalar hides which.
##
## C IS DEFINED ON FIXED 50 kb BINS, not on clusters. The partition itself moves
## between cells, so a cluster is not a stable coordinate and matching clusters
## across partitions is exactly where an earlier overlap-scoping bug lived. A bin
## is "called" in a cell if any significant cluster overlaps it.
##
## RESUME KEY INCLUDES EVERY GRID PARAMETER. An earlier family of scripts in this
## module resumed on (cell,tag,env) alone, which would have silently mixed two
## dcap settings in one table. Output files here are named by the full setting.
##
## Run:  Rscript sensitivity_grid.R
## Env: SIM_DATA, OUT, CELLS, TAGS, ENVS, FILES, LDW, RHO, FLOORS, ALPHA, CORES
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(parallel)})
SIM    <- Sys.getenv("SIM_DATA", "~/sens_grid/data")
OUT    <- Sys.getenv("OUT", "~/sens_grid/out")
CELLS  <- strsplit(Sys.getenv("CELLS", "V0.5_c1,V0.5_c2,V1_c1.5,V2_c1"), ",")[[1]]
TAGS   <- strsplit(Sys.getenv("TAGS", "nobgs,bgs"), ",")[[1]]
ENVS   <- as.integer(strsplit(Sys.getenv("ENVS", "1,2,3,4,5,6,7,8,9,10"), ",")[[1]])
FILES  <- as.integer(strsplit(Sys.getenv("FILES", "1,2,3,4,5,6,7,8,9,10"), ",")[[1]])
LDW    <- as.numeric(strsplit(Sys.getenv("LDW", "0.0125,0.025,0.05"), ",")[[1]])
RHO    <- as.numeric(strsplit(Sys.getenv("RHO", "0.35,0.5,0.65"), ",")[[1]])
FLOORS <- as.integer(strsplit(Sys.getenv("FLOORS", "1,2,5"), ",")[[1]])
ALPHA  <- as.numeric(Sys.getenv("ALPHA", "0.05"))
CORES  <- as.integer(Sys.getenv("CORES", "8"))
BIN    <- as.numeric(Sys.getenv("BIN", "5e4"))
SIM <- path.expand(SIM); OUT <- path.expand(OUT)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
simes <- function(p) { p <- sort(p[is.finite(p)]); n <- length(p)
  if (!n) return(NA_real_); min(n * p / seq_len(n)) }

one_panel <- function(CELL, TAG, ENV) {
  tag <- sprintf("%s_%s_env%d", CELL, TAG, ENV)
  fo  <- file.path(OUT, sprintf("panel_%s.rds", tag))
  if (file.exists(fo)) return(invisible(NULL))
  clu <- list(); lnk <- list(); runsrec <- list()
  for (i in FILES) {
    f <- file.path(SIM, sprintf("adapt_%s_chr%d_%s_env%d.rds", TAG, i, CELL, ENV))
    if (!file.exists(f)) next
    x <- readRDS(f); m <- flag_true_qtns(as.data.table(x$map))
    th <- score_thresholds(as.data.table(x$LD_decay$decay_sum),
                           rho_r2 = 0.75, rho_d = 0.95, dmax_cap = 1e5)
    drv <- m[true_pos_QTN %in% TRUE]
    for (lw in LDW) for (rr in RHO) {
      pr <- ld_prune_and_eMLG(GTs = x$GTs, stage1 = x$complexity_reduction$stage1,
              LD_decay = x$LD_decay, ld_w_col = "ld_w_095", ld_w_threshold = lw,
              score_threshold = 0.80, min_r2_rho = rr, distance_threshold = 1e5,
              compute_unflagged_eMLG = TRUE, cores = 1)
      g  <- as.data.table(pr$groups)
      ms <- ld_group_map(g, prefix = i)[, .(marker, CL = group_id)]
      mm <- merge(m, ms, by = "marker", all.x = TRUE)[!is.na(CL)]
      su <- mm[, .(emx = simes(emx_p), lfm = simes(lfmm_p), n_loci = .N,
                   Chr = as.character(Chr)[1], pmin = min(Pos), pmax = max(Pos),
                   has_qtn = any(true_pos_QTN %in% TRUE)), by = CL]
      su[, `:=`(ldw = lw, rho = rr, chrfile = i)]
      clu[[length(clu)+1]] <- su
      runsrec[[length(runsrec)+1]] <- data.table(ldw = lw, rho = rr, chrfile = i,
                                                 n_flag = sum(!grepl("^[0-9]+_U", su$CL)))
      ## QTN tagging links, per setting -- r2 geometry does not change with the
      ## setting but cluster membership does, so this must be recomputed per cell.
      if (nrow(drv)) lnk[[length(lnk)+1]] <- rbindlist(lapply(seq_len(nrow(drv)), function(j) {
        ch <- as.character(drv$Chr[j])
        near <- mm[as.character(Chr) == ch & abs(Pos - drv$Pos[j]) < th$dmax]
        if (!nrow(near)) return(NULL)
        r2 <- suppressWarnings(cor(x$GTs[, drv$marker[j]], x$GTs[, near$marker],
                                   use = "pairwise.complete.obs")^2)
        d <- data.table(CL = near$CL, r2 = as.numeric(r2))[is.finite(r2) & r2 >= th$r2min]
        if (!nrow(d)) return(NULL)
        d[, .(r2 = max(r2)), by = CL][, `:=`(qtn = paste0(i, "_", drv$marker[j]),
                                             ldw = lw, rho = rr)][] }))
    }
  }
  su <- rbindlist(clu, fill = TRUE); lk <- rbindlist(lnk, fill = TRUE)
  if (!nrow(su)) return(invisible(NULL))
  ## QTN truth is a property of the panel, not of a grid cell: fixed once so the
  ## recall denominator is identical across cells.
  ##
  ## TWO DENOMINATORS ARE RECORDED, deliberately. floor_sweep.R and
  ## snp_vs_cluster_dedup.R divide the count of QTN RECOVERED by the count of
  ## CLUSTERS CONTAINING a QTN -- different units, which coincide only when no
  ## cluster holds two QTN. The correct denominator is the number of distinct
  ## taggable QTN; the legacy one is kept so this grid can be checked against
  ## those tables rather than merely disagreeing with them.
  nq      <- uniqueN(lk$qtn)
  nq_lega <- sum(su[ldw == LDW[LDW == 0.025][1] & rho == 0.5]$has_qtn)
  if (!nq) { saveRDS(NULL, fo); return(invisible(NULL)) }

  res <- list(); calls <- list()
  for (lw in LDW) for (rr in RHO) for (fo_ in FLOORS) for (eng in c("emx","lfm")) {
    s <- su[ldw == lw & rho == rr]
    ok <- s$n_loci >= fo_ & is.finite(s[[eng]])
    q  <- rep(NA_real_, nrow(s)); q[ok] <- p.adjust(s[[eng]][ok], "BH")
    sig <- which(!is.na(q) & q < ALPHA)
    fl  <- s$CL[sig]
    sub  <- lk[ldw == lw & rho == rr][CL %in% fl]
    best <- if (nrow(sub)) sub[order(-r2)][, .SD[1], by = qtn] else sub
    kept <- setdiff(fl, setdiff(sub$CL, best$CL))
    res[[length(res)+1]] <- data.table(ldw=lw, rho=rr, floor=fo_, engine=eng,
      n_test = sum(ok), flagged = length(fl), regions = length(kept),
      tp_reg = if (nrow(best)) uniqueN(best$CL) else 0L,
      tp_qtn = if (nrow(best)) uniqueN(best$qtn) else 0L,
      n_qtn = nq, n_qtn_legacy = nq_lega,
      prec = (if (nrow(best)) uniqueN(best$CL) else 0L) / max(length(kept), 1),
      rec  = (if (nrow(best)) uniqueN(best$qtn) else 0L) / nq)
    ## bin-level calls, the partition-independent coordinate for C
    if (length(sig)) {
      z <- s[sig]
      calls[[length(calls)+1]] <- unique(rbindlist(lapply(seq_len(nrow(z)), function(k)
        data.table(chrfile = z$chrfile[k], Chr = z$Chr[k],
                   bin = seq(floor(z$pmin[k]/BIN), floor(z$pmax[k]/BIN))
        )))[, `:=`(ldw=lw, rho=rr, floor=fo_, engine=eng)])
    }
  }
  ## truth on the same bin grid
  ## truth bins from the CANONICAL partition only -- pooling over settings would
  ## make the reference grid itself depend on the grid being evaluated.
  qb <- unique(su[has_qtn == TRUE & ldw == 0.025 & rho == 0.5,
                  .(chrfile, Chr, bin = floor(pmin/BIN))])
  saveRDS(list(cell = CELL, tag = TAG, env = ENV, n_qtn = nq, n_qtn_legacy = nq_lega,
               perf = rbindlist(res), calls = rbindlist(calls), qtn_bins = qb,
               runs = rbindlist(runsrec)), fo)
  invisible(NULL)
}

grid <- CJ(cell = CELLS, tag = TAGS, env = ENVS, sorted = FALSE)
cat(sprintf("%d panels x %d clusterings x %d floors x 2 engines | CORES=%d\n",
            nrow(grid), length(LDW)*length(RHO), length(FLOORS), CORES))
invisible(mclapply(seq_len(nrow(grid)), function(z)
  tryCatch(one_panel(grid$cell[z], grid$tag[z], grid$env[z]),
           error = function(e) cat("FAIL", grid$cell[z], grid$tag[z], grid$env[z],
                                   conditionMessage(e), "\n")),
  mc.cores = CORES, mc.preschedule = FALSE))
cat("ALL PANELS DONE\n")
