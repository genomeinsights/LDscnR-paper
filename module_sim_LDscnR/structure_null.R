## =============================================================================
## module_sim_LDscnR / structure_null.R
##
## A STRUCTURE-AWARE NULL for the cluster-level test, ported from the design PK
## already validated in Formica_hybrid (moduleB_eMLG_null.R / moduleC_null_regen.R,
## after Li et al. 2018 MER 18:809-824).
##
## WHY THIS MATTERS MORE THAN AS A FEATURE. The paper's central claim is that BH
## controls the false discovery proportion among the UNITS TESTED, so testing LD
## clusters makes the reported FDR a statement about loci. The obvious attack on
## it is whether BH is even valid on Simes p-values over correlated clusters.
## A structure-aware null retires that question instead of arguing it, because
## the FDR estimate below assumes NOTHING about dependence between units.
##
## THE CONSTRUCTION. NSIM null covariates are drawn from the among-population
## covariance -- here the GRM aggregated to population level, which is exactly
## BayPass's Omega under a different name, and like it is estimated from the
## LD-PRUNED marker set (`grm_markers`). Each draw has the structure
## autocorrelation of a real environment and no causal link to genotypes. For
## each unit, k = the number of nulls matching or beating its observed statistic;
## SURVIVORS are the units at the FLOOR, k == 0.
##
## Under the pure-structure null an observed statistic is exchangeable with its
## NSIM nulls, so P(observed is the strict max) = 1/(NSIM+1) FOR EACH UNIT
## SEPARATELY -- no independence between units is needed. Hence
##
##     E[survivors | null] = N_tested / (NSIM + 1)
##     estimated FDR       = E[survivors | null] / observed survivors
##
## THE FLOOR DOES NOT SCALE TO THIS UNIT COUNT, WHICH IS A PORTING FINDING RATHER
## THAN A FLAW IN THE ORIGINAL. Its null expectation is N/(NSIM+1), so a useful
## FDR needs NSIM of order N, and the cost is NSIM x (time per covariate), which
## is proportional to N -- so total cost grows as N^2. Measured here: 13,437
## units per chromosome file, 134,370 per panel, 0.619 s per covariate per file.
## At NSIM = 10,000 that is 17.2 h per panel AND STILL LEAVES A NULL EXPECTATION
## OF 13.4, so the FDR would be ~0.7 before any real signal is considered. The
## design is sound at Formica's N = 32,840; it is not affordable at 134,370.
##
## SO THE DEFAULT ESTIMATOR HERE IS THE PERMUTATION FDR, whose cost is INDEPENDENT
## of N. For a threshold t,
##
##     E[V(t)] = mean over nulls of #{units with p_null <= t}
##     FDP(t)  = E[V(t)] / #{units with p_obs <= t}
##
## and t is chosen as the largest threshold holding FDP(t) <= alpha. This makes
## no independence assumption either -- the null covariates carry the dependence
## between units exactly, because every unit is recomputed under the same draw --
## and NSIM = 1000 is ample, since the quantity estimated is a MEAN COUNT rather
## than a per-unit tail probability. The floor remains available (ESTIMATOR=floor)
## for datasets where N is small enough to afford it, and is what should be used
## when comparing directly against the Formica numbers.
##
## THE SIMS CAN DO THE ONE THING FORMICA COULD NOT: check that estimate against
## the truth. QTN positions are known, so the REALISED false discovery proportion
## of the survivor set is measurable and can be compared with the estimate. That
## comparison is a direct test of the paper's central claim, and it is stronger
## evidence than any precision ratio, because it validates the guarantee rather
## than the performance.
##
## TWO ROUTES, matching the two classes of association method:
##   consensus  the test is run on the cluster's eMLG consensus genotype. This is
##              the only route open to a method that emits no per-marker p-value
##              (BayPass BF, XtX), and it is ALSO what makes a large NSIM
##              affordable: one call per covariate over ~N units instead of over
##              ~300k markers. Ranked second to Simes on power (cluster_summary_
##              test.R), so this measures the price of engine-generality.
##   simes      Simes over member p-values. Better powered, far more expensive
##              under a null this size, so it is run at a reduced NSIM.
##
## WHY NOT PER-UNIT PERMUTATION p-VALUES WITH BH: BH over 1.3e5 units at alpha
## 0.05 needs p resolved to ~4e-7, i.e. NSIM >= 2.5e6. The permutation FDR above
## avoids this because it never needs a per-unit tail probability -- only the
## average number of null units below a threshold, which 1000 draws estimate well.
##
## REUSES THE PACKAGE RATHER THAN REIMPLEMENTING IT. emmax_setup() eigendecomposes
## the kinship once and rotates the genotypes into the whitening basis, so each
## surrogate costs one emmax_fast() call -- ~25x faster, and the difference
## between an affordable null and an unaffordable one here. LDscnR already has
## structured_null(), which generates surrogates from the same basis and uses the
## same "mean surrogate / observed" FDR; it is wired to the C-score, which this
## project retired, so only the surrogate machinery is reused and the scoring is
## done at the stage-2 cluster level instead.
##
## TWO SURROGATE BASES, because the package's choice and PK's Formica design
## differ and the difference is not cosmetic:
##   pop    drawn at the POPULATION level from the block-mean GRM. This is the
##          faithful analogue of BayPass's Omega and of the observed covariate,
##          which is CONSTANT WITHIN A POPULATION here.
##   indiv  drawn at the individual level from K, as structured_null() does. Its
##          surrogates vary within a population where the real environment does
##          not, so they are an imperfect match to the observed covariate and
##          could be anti-conservative. Run for comparison, not as the default.
## Both are residualised against the observed environment, following
## structured_null(), so a surrogate carries none of the real signal.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/structure_null.R
## Env: SIM_DATA, OUT, CELL, TAG, ENVS, FILES, NSIM, ROUTE, BASIS, SEED, CORES
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(parallel)})
SIM   <- path.expand(Sys.getenv("SIM_DATA", "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5"))
OUT   <- path.expand(Sys.getenv("OUT", "module_sim_LDscnR/results/structure_null"))
CELLS <- strsplit(Sys.getenv("CELLS", "V0.5_c1,V0.5_c2,V1_c1.5,V2_c1"), ",")[[1]]
TAGS  <- strsplit(Sys.getenv("TAGS", "nobgs,bgs"), ",")[[1]]
ENVS  <- as.integer(strsplit(Sys.getenv("ENVS", "1,2,3,4,5"), ",")[[1]])
FILES <- as.integer(strsplit(Sys.getenv("FILES", "1,2,3,4,5,6,7,8,9,10"), ",")[[1]])
NSIM  <- as.integer(Sys.getenv("NSIM", "1000"))
ROUTE <- Sys.getenv("ROUTE", "consensus")          # consensus | simes
BASIS <- Sys.getenv("BASIS", "pop")                 # pop | indiv
SEED  <- as.integer(Sys.getenv("SEED", "2026"))
CORES <- as.integer(Sys.getenv("CORES", "4"))
## Threshold grid for the permutation FDR. Log-spaced over the range where a
## genome-wide discovery can live; the estimator reports FDP at each and takes
## the largest t holding FDP <= alpha.
TGRID <- as.numeric(strsplit(Sys.getenv(
  "TGRID", paste(signif(10^seq(-9, -2, by = 0.25), 3), collapse = ",")), ",")[[1]])
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
simes <- function(p) { p <- sort(p[is.finite(p)]); n <- length(p)
  if (!n) return(NA_real_); min(n * p / seq_len(n)) }

## Null covariates at the POPULATION level, because the real environment is
## constant within a population -- an individual-level draw would not be an
## analogue of the observed covariate. The population-level covariance is the
## block mean of the GRM, which is the object BayPass calls Omega.
draw_null_env <- function(GRM, pop, y, nsim, seed, basis = "pop") {
  n <- length(pop)
  if (basis == "pop") {
    pl  <- sort(unique(pop)); np <- length(pl)
    idx <- lapply(pl, function(p) which(pop == p))
    Om  <- matrix(0, np, np)
    for (a in seq_len(np)) for (b in a:np) {
      v <- mean(GRM[idx[[a]], idx[[b]]]); Om[a,b] <- v; Om[b,a] <- v }
    ## the GRM is singular (a zero eigenvalue), so draw through the eigen
    ## decomposition with negatives clamped rather than a Cholesky, which fails.
    e <- eigen(Om, symmetric = TRUE)
    L <- e$vectors %*% diag(sqrt(pmax(e$values, 0)), np)
    set.seed(seed)
    Z <- L %*% matrix(stats::rnorm(np * nsim), np, nsim)
    out <- matrix(NA_real_, n, nsim)
    for (a in seq_len(np)) out[idx[[a]], ] <- rep(Z[a, ], each = length(idx[[a]]))
  } else {
    e <- eigen(GRM, symmetric = TRUE)
    L <- e$vectors %*% diag(sqrt(pmax(e$values, 0)), n)
    set.seed(seed)
    out <- L %*% matrix(stats::rnorm(n * nsim), n, nsim)
  }
  ## residualise against the observed environment, as structured_null() does, so
  ## no surrogate carries the real signal.
  apply(out, 2, function(v) as.numeric(stats::resid(stats::lm(v ~ y))))
}

one_panel <- function(CELL, TAG, ENV) {
  ## the resume key carries every parameter that defines the result -- an
  ## earlier family of scripts here resumed on (cell,tag,env) alone.
  fo <- file.path(OUT, sprintf("null_%s_%s_env%d_%s_%s_N%d.rds",
                               CELL, TAG, ENV, ROUTE, BASIS, NSIM))
  if (file.exists(fo)) return(invisible(NULL))
  obs <- list(); nullmax <- list(); nullcnt <- list(); lnk <- list(); meta <- list()
  for (i in FILES) {
    f <- file.path(SIM, sprintf("adapt_%s_chr%d_%s_env%d.rds", TAG, i, CELL, ENV))
    if (!file.exists(f)) next
    x <- readRDS(f); m <- flag_true_qtns(as.data.table(x$map))
    pr <- ld_prune_and_eMLG(GTs = x$GTs, stage1 = x$complexity_reduction$stage1,
            LD_decay = x$LD_decay, ld_w_col = "ld_w_095", ld_w_threshold = 0.025,
            score_threshold = 0.80, min_r2_rho = 0.5, distance_threshold = 1e5,
            compute_unflagged_eMLG = TRUE, cores = 1)
    g  <- as.data.table(pr$groups)
    ms <- rbindlist(lapply(seq_len(nrow(g)), function(k)
            data.table(marker = g$members[[k]], CL = paste0(i, "_", g$group_id[k]))))
    mm <- merge(m, ms, by = "marker", all.x = TRUE)[!is.na(CL)]
    E  <- pr$eMLG
    if (is.null(E) || !ncol(E)) next
    colnames(E) <- paste0(i, "_", colnames(E))
    yy <- as.numeric(x$env$env)
    ## the SAME null covariates for every chromosome of a panel -- one structure,
    ## one draw, so units are calibrated against a common null rather than
    ## against ten independent ones.
    NE <- draw_null_env(x$GRM, x$env$pop, yy, NSIM, SEED + ENV, BASIS)
    ## one eigendecomposition and one rotation for all NSIM+1 scans
    prep <- emmax_setup(E, x$GRM)
    o <- emmax_fast(prep, yy)
    obs[[length(obs)+1]] <- data.table(CL = colnames(E), p_obs = as.numeric(o))
    ## Two accumulators, neither holding the NSIM x N matrix (the memory failure
    ## mode of the Formica version):
    ##   kk    per-unit exceedance count, for the floor estimator
    ##   cnt   per-null COUNT of units below each threshold on a fixed grid, which
    ##         is all the permutation FDR needs. Its size is NSIM x |grid|, not
    ##         NSIM x N, so the grid is what makes the estimator cheap.
    kk  <- integer(ncol(E))
    cnt <- matrix(0L, NSIM, length(TGRID))
    for (b in seq_len(NSIM)) {
      pn <- as.numeric(emmax_fast(prep, NE[, b]))
      kk <- kk + as.integer(pn <= as.numeric(o))
      cnt[b, ] <- as.integer(vapply(TGRID, function(t) sum(pn <= t, na.rm = TRUE), numeric(1)))
    }
    nullmax[[length(nullmax)+1]] <- data.table(CL = colnames(E), k = kk)
    nullcnt[[length(nullcnt)+1]] <- cnt
    ## truth links, same dedup convention as everywhere else in this module
    th  <- score_thresholds(as.data.table(x$LD_decay$decay_sum),
                            rho_r2 = 0.75, rho_d = 0.95, dmax_cap = 1e5)
    drv <- mm[true_pos_QTN %in% TRUE]
    if (nrow(drv)) lnk[[length(lnk)+1]] <- rbindlist(lapply(seq_len(nrow(drv)), function(j) {
      ch <- as.character(drv$Chr[j])
      near <- mm[as.character(Chr) == ch & abs(Pos - drv$Pos[j]) < th$dmax]
      if (!nrow(near)) return(NULL)
      r2 <- suppressWarnings(cor(x$GTs[, drv$marker[j]], x$GTs[, near$marker],
                                 use = "pairwise.complete.obs")^2)
      d <- data.table(CL = near$CL, r2 = as.numeric(r2))[is.finite(r2) & r2 >= th$r2min]
      if (!nrow(d)) return(NULL)
      d[, .(r2 = max(r2)), by = CL][, qtn := paste0(i, "_", drv$marker[j])][] }))
    meta[[length(meta)+1]] <- data.table(chrfile = i, n_units = ncol(E))
  }
  if (!length(obs)) return(invisible(NULL))
  O <- merge(rbindlist(obs), rbindlist(nullmax), by = "CL")
  ## null counts sum ACROSS chromosomes within a null draw -- the panel, not the
  ## chromosome, is the unit the FDR is stated over.
  NC <- Reduce(`+`, nullcnt)
  lk <- rbindlist(lnk, fill = TRUE); nq <- uniqueN(lk$qtn)
  saveRDS(list(cell = CELL, tag = TAG, env = ENV, route = ROUTE, nsim = NSIM,
               units = O, null_counts = NC, tgrid = TGRID,
               links = lk, n_qtn = nq, meta = rbindlist(meta)), fo)
  cat(sprintf("  %s env%d: %d units, %d survivors (k=0), expected %.1f\n",
              TAG, ENV, nrow(O), sum(O$k == 0), nrow(O)/(NSIM+1)))
  invisible(NULL)
}

grid <- CJ(cell = CELLS, tag = TAGS, env = ENVS, sorted = FALSE)
cat(sprintf("%s | route %s | basis %s | NSIM %d | %d panels | CORES %d\n",
            paste(CELLS, collapse=","), ROUTE, BASIS, NSIM, nrow(grid), CORES))
invisible(mclapply(seq_len(nrow(grid)), function(z)
  tryCatch(one_panel(grid$cell[z], grid$tag[z], grid$env[z]),
           error = function(e) cat("FAIL", grid$cell[z], grid$tag[z], grid$env[z],
                                   conditionMessage(e), "\n")),
  mc.cores = CORES, mc.preschedule = FALSE))
cat("DONE\n")
