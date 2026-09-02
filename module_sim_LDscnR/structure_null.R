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
## THREE SUMMARISATION ROUTES, ALL SCORED FROM THE SAME SURROGATE DRAW, so the
## comparison is between rules rather than between nulls:
##   consensus  the test is run on the cluster's eMLG consensus genotype. The only
##              route open to a method that emits no per-marker p-value (BayPass
##              BF, XtX). Ranked second to Simes under BH (cluster_summary_test.R),
##              so this measures the price of engine-generality.
##   best       the smallest member p-value -- what formica_paper used. Under BH
##              this is badly calibrated, being the minimum of correlated
##              p-values, which is why it ranked poorly before. UNDER A
##              PERMUTATION NULL IT IS CALIBRATED AUTOMATICALLY, because the
##              surrogates are summarised by the same statistic. So the null is
##              expected to rescue a rule that BH could not use, and that
##              prediction is worth recording before the numbers arrive.
##   simes      Simes over member p-values. Best under BH; included here to see
##              whether its advantage survives a null that calibrates its rivals.
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
##   vc     THE CALIBRATED ONE, and the default after the other three failed:
##          y_null ~ MVN(0, vg*K + ve*I) with vg and ve estimated by REML from the
##          OBSERVED phenotype. The three bases below all draw from PURE
##          STRUCTURE with no residual term, so a surrogate is more structured
##          than any real phenotype, which has both a heritable and a residual
##          component. Measured over 40 panels, that miscalibrates in BOTH
##          directions depending on which structure is used -- estimated FDP
##          0.001-0.008 under `pop` and 1.4-3.2 under `spatial`, against a
##          realised 0.57-0.76. Neither brackets the truth by accident: the
##          missing ve*I term is what sets how much of the phenotype is
##          structure. This is the form PK asked the panel session for.
##   spatial drawn from a Gaussian kernel on the (x, y) coordinates, the second
##          basis structured_null() offers. THIS EXISTS BECAUSE THE GENETIC BASIS
##          FAILED A CALIBRATION CHECK: under `pop`, the estimated FDP is
##          0.1-0.8% while the REALISED proportion of discoveries tagging no QTN
##          is 32-54% -- two orders of magnitude apart. Widening the truth set
##          from detectable QTN to ALL QTN changes that by nothing, so the gap is
##          not a scoring artefact. The observed environment is a spatial
##          gradient, and a surrogate drawn from genetic covariance need not
##          carry the spatial autocorrelation that makes the real covariate hard;
##          the 3sp work in this project found the spatial null the stringent one
##          for the same reason.
## Both are residualised against the observed environment, following
## structured_null(), so a surrogate carries none of the real signal.
##
## SURROGATE-MAJOR, NOT CHROMOSOME-MAJOR. BH is applied over the whole panel, so
## a surrogate's p-values must exist for all ten chromosomes at once before it can
## be thresholded. The loop therefore holds ten emmax_setup() rotations (~170 MB)
## and iterates surrogates on the outside. A chromosome-major loop can only
## accumulate counts below fixed thresholds, which cannot reproduce a BH rejection
## count -- and the BH count is the quantity the panel session needs.
##
## THE MULTIPLICITY QUESTION IT ANSWERS. 2c's floor 8 escapes a 300x multiplicity
## reduction (790,578 markers -> 2,631 units) against roughly 6x here, and asks
## whether the surrogate discovery rate scales with the multiplicity escaped
## rather than with the data. Recording BH rejections per FLOOR under the same
## surrogates answers it directly: if raising the floor raises observed
## discoveries without raising E[V], the gain is detection; if E[V] rises with
## it, that is relief made visible in a null instead of argued from a contrast.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/structure_null.R
## Env: SIM_DATA, OUT, CELL, TAG, ENVS, FILES, NSIM, ROUTE, BASIS, SEED, CORES
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(parallel)})
emma.REMLE <- utils::getFromNamespace("emma.REMLE", "LDscnR")
SIM   <- path.expand(Sys.getenv("SIM_DATA", "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5"))
OUT   <- path.expand(Sys.getenv("OUT", "module_sim_LDscnR/results/structure_null"))
CELLS <- strsplit(Sys.getenv("CELLS", "V0.5_c1,V0.5_c2,V1_c1.5,V2_c1"), ",")[[1]]
TAGS  <- strsplit(Sys.getenv("TAGS", "nobgs,bgs"), ",")[[1]]
ENVS  <- as.integer(strsplit(Sys.getenv("ENVS", "1,2,3,4,5"), ",")[[1]])
FILES <- as.integer(strsplit(Sys.getenv("FILES", "1,2,3,4,5,6,7,8,9,10"), ",")[[1]])
NSIM  <- as.integer(Sys.getenv("NSIM", "1000"))
ROUTES <- strsplit(Sys.getenv("ROUTES", "consensus,best,simes"), ",")[[1]]
BASIS  <- Sys.getenv("BASIS", "pop")                # pop | indiv
FLOORS <- as.integer(strsplit(Sys.getenv("FLOORS", "1,2,5,8"), ",")[[1]])
ALPHA  <- as.numeric(Sys.getenv("ALPHA", "0.05"))
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
draw_null_env <- function(GRM, pop, y, nsim, seed, basis = "pop", coords = NULL,
                          prep = NULL) {
  n <- length(pop)
  if (basis == "vc") {
    ## REML components from the observed phenotype, so the surrogate carries the
    ## same balance of structure and noise rather than pure structure.
    if (is.null(prep)) stop("basis = \"vc\" needs a prep from emmax_setup()")
    re <- emma.REMLE(y, prep$Xo, prep$Kn, eig.R = prep$eigR)
    Sig <- re$vg * prep$Kn + re$ve * diag(n)
    e <- eigen(Sig, symmetric = TRUE)
    L <- e$vectors %*% diag(sqrt(pmax(e$values, 0)), n)
    set.seed(seed)
    out <- L %*% matrix(stats::rnorm(n * nsim), n, nsim)
  } else if (basis == "spatial") {
    if (is.null(coords)) stop("basis = \"spatial\" needs coords")
    Dm <- as.matrix(stats::dist(coords)); l <- stats::median(Dm[lower.tri(Dm)])
    e  <- eigen(exp(-0.5 * (Dm / l)^2), symmetric = TRUE)
    L  <- e$vectors %*% diag(sqrt(pmax(e$values, 0)), n)
    set.seed(seed)
    out <- L %*% matrix(stats::rnorm(n * nsim), n, nsim)
  } else if (basis == "pop") {
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
  fo <- file.path(OUT, sprintf("null_%s_%s_env%d_%s_%s_N%d.rds",
                               CELL, TAG, ENV, paste(ROUTES, collapse = "-"), BASIS, NSIM))
  if (file.exists(fo)) return(invisible(NULL))

  ## ---- pass 1: build every chromosome's partition, rotation and truth links
  P <- list(); lnk <- list(); yy <- NULL; GRM <- NULL; POP <- NULL; XY <- NULL
  for (i in FILES) {
    f <- file.path(SIM, sprintf("adapt_%s_chr%d_%s_env%d.rds", TAG, i, CELL, ENV))
    if (!file.exists(f)) next
    x <- readRDS(f); m <- flag_true_qtns(as.data.table(x$map))
    pr <- ld_prune_and_eMLG(GTs = x$GTs, stage1 = x$complexity_reduction$stage1,
            LD_decay = x$LD_decay, ld_w_col = "ld_w_095", ld_w_threshold = 0.025,
            score_threshold = 0.80, min_r2_rho = 0.5, distance_threshold = 1e5,
            compute_unflagged_eMLG = TRUE, cores = 1)
    E <- pr$eMLG
    if (is.null(E) || !ncol(E)) next
    g  <- as.data.table(pr$groups)
    nl <- stats::setNames(g$n_loci, g$group_id)[colnames(E)]   # unit size, for the floor
    colnames(E) <- paste0(i, "_", colnames(E))
    if (is.null(yy)) { yy <- as.numeric(x$env$env); GRM <- x$GRM; POP <- x$env$pop
                       XY <- as.matrix(x$env[, .(x, y)]) }
    ms <- rbindlist(lapply(seq_len(nrow(g)), function(k)
            data.table(marker = g$members[[k]], CL = paste0(i, "_", g$group_id[k]))))
    mm <- merge(m, ms, by = "marker", all.x = TRUE)[!is.na(CL)]
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
    ## the marker-level rotation is what best/simes need; consensus needs only E.
    keepm <- mm$marker[mm$marker %in% colnames(x$GTs)]
    gi    <- match(mm$CL[match(keepm, mm$marker)], colnames(E))
    P[[length(P)+1]] <- list(
      prep  = emmax_setup(E, x$GRM), cl = colnames(E), n_loci = nl,
      prepM = if (any(ROUTES %in% c("best","simes")))
                emmax_setup(x$GTs[, keepm, drop = FALSE], x$GRM) else NULL,
      gidx  = gi)                       # member marker -> unit index within this file
    rm(x, E); gc(FALSE)
  }
  if (!length(P)) return(invisible(NULL))
  CLall <- unlist(lapply(P, `[[`, "cl"), use.names = FALSE)
  NL    <- unlist(lapply(P, `[[`, "n_loci"), use.names = FALSE)

  ## ---- one scan -> three per-unit vectors, pooled across chromosomes
  ## Simes and best are aggregations of the SAME marker scan, so they cost one
  ## extra emmax_fast per file rather than one per route.
  score_routes <- function(yv) {
    cons <- if ("consensus" %in% ROUTES)
      unlist(lapply(P, function(z) as.numeric(emmax_fast(z$prep, yv))), use.names = FALSE) else NULL
    bs <- sm <- NULL
    if (any(ROUTES %in% c("best","simes"))) {
      bl <- list(); sl <- list()
      for (z in P) {
        pm <- as.numeric(emmax_fast(z$prepM, yv))
        D  <- data.table::data.table(g = z$gidx, p = pm)[is.finite(p) & !is.na(g)]
        data.table::setorder(D, g, p)
        D[, r := seq_len(.N), by = g]
        D[, ng := .N, by = g]
        agg <- D[, .(best = p[1L], simes = min(ng * p / r)), by = g]
        b <- rep(NA_real_, length(z$cl)); b[agg$g] <- agg$best
        m2 <- rep(NA_real_, length(z$cl)); m2[agg$g] <- agg$simes
        bl[[length(bl)+1]] <- b; sl[[length(sl)+1]] <- m2
      }
      bs <- unlist(bl, use.names = FALSE); sm <- unlist(sl, use.names = FALSE)
    }
    list(consensus = cons, best = bs, simes = sm)
  }
  OBS <- score_routes(yy)
  bh_rej <- function(pv, keep) {
    v <- pv[keep]; v <- v[is.finite(v)]
    if (!length(v)) return(0L)
    sum(p.adjust(v, "BH") < ALPHA)
  }
  bh_set <- function(pv, keep) {
    idx <- which(keep & is.finite(pv))
    if (!length(idx)) return(character(0))
    CLall[idx][p.adjust(pv[idx], "BH") < ALPHA]
  }
  keeps <- lapply(FLOORS, function(f) NL >= f)
  obs_R <- lapply(ROUTES, function(r) vapply(keeps, function(k) bh_rej(OBS[[r]], k), integer(1)))
  obs_S <- lapply(ROUTES, function(r) lapply(keeps, function(k) bh_set(OBS[[r]], k)))
  names(obs_R) <- names(obs_S) <- ROUTES

  ## ---- surrogates, one draw at a time across the whole panel
  NE <- draw_null_env(GRM, POP, yy, NSIM, SEED + ENV, BASIS, XY, P[[1]]$prep)
  Vmat <- lapply(ROUTES, function(r) matrix(0L, NSIM, length(FLOORS)))
  names(Vmat) <- ROUTES
  for (b in seq_len(NSIM)) {
    SB <- score_routes(NE[, b])
    for (r in ROUTES) Vmat[[r]][b, ] <- vapply(keeps, function(k) bh_rej(SB[[r]], k), integer(1))
  }

  lk <- rbindlist(lnk, fill = TRUE); nq <- uniqueN(lk$qtn)
  ## realised FDP of the observed set, by the module's usual dedup convention:
  ## one region per QTN, best tagger, satellites removed.
  realised_for <- function(SS) vapply(seq_along(FLOORS), function(j) {
    fl <- SS[[j]]; if (!length(fl)) return(NA_real_)
    sub  <- lk[CL %in% fl]
    best <- if (nrow(sub)) sub[order(-r2)][, .SD[1], by = qtn] else sub
    kept <- setdiff(fl, setdiff(sub$CL, best$CL))
    tp   <- if (nrow(best)) uniqueN(best$CL) else 0L
    1 - tp / max(length(kept), 1)
  }, numeric(1))

  out <- rbindlist(lapply(ROUTES, function(r) data.table(
    cell = CELL, tag = TAG, env = ENV, basis = BASIS, nsim = NSIM, route = r,
    floor = FLOORS,
    n_test  = vapply(keeps, sum, integer(1)),
    obs_R   = obs_R[[r]],
    EV      = colMeans(Vmat[[r]]),
    V_med   = apply(Vmat[[r]], 2, stats::median),
    V_max   = apply(Vmat[[r]], 2, max),
    pct_any = 100 * colMeans(Vmat[[r]] > 0),
    est_FDP = colMeans(Vmat[[r]]) / pmax(obs_R[[r]], 1),
    realised_FDP = realised_for(obs_S[[r]]),
    n_qtn = nq)))
  saveRDS(list(summary = out, Vmat = Vmat, floors = FLOORS, routes = ROUTES,
               units = data.table(CL = CLall, n_loci = NL,
                                  p_cons = OBS$consensus, p_best = OBS$best, p_simes = OBS$simes),
               links = lk), fo)
  cat(sprintf("  %s %s env%d: %s\n", CELL, TAG, ENV,
      paste(sprintf("%s f%d R=%d EV=%.2f", out$route, out$floor, out$obs_R, out$EV),
            collapse = " | ")))
  invisible(NULL)
}

grid <- CJ(cell = CELLS, tag = TAGS, env = ENVS, sorted = FALSE)
cat(sprintf("%s | routes %s | basis %s | NSIM %d | %d panels | CORES %d\n",
            paste(CELLS, collapse=","), paste(ROUTES, collapse=","), BASIS, NSIM,
            nrow(grid), CORES))
invisible(mclapply(seq_len(nrow(grid)), function(z)
  tryCatch(one_panel(grid$cell[z], grid$tag[z], grid$env[z]),
           error = function(e) cat("FAIL", grid$cell[z], grid$tag[z], grid$env[z],
                                   conditionMessage(e), "\n")),
  mc.cores = CORES, mc.preschedule = FALSE))
cat("DONE\n")
