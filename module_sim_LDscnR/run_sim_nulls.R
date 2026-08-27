## =====================================================================
## module_sim_LDscnR / run_sim_nulls.R
##
## Permutation / surrogate nulls for the SIM benchmark -- the sim counterpart of
## ~/3sp_lfmm_perm/run_paired_nulls_3sp.R, same paired two-engine design, same
## per-(type,b) seeding so EMMAX and LFMM see the IDENTICAL surrogate phenotypes.
##
## THIS SCRIPT PRODUCES RAW P-VALUES ONLY: the observed per-marker p vector and one
## p vector per surrogate. It computes no C-scores, applies no threshold and builds
## no regions.
##
## OUTPUT CONTRACT -- exactly what module_sim_LDscnR/analyse_one_dataset.R consumes,
## so these files are handed over and run as-is:
##
##   panel_<cell>.rds                     the LD description, method-independent
##     GTs        individuals x markers, colnames = markers
##     map        marker, Chr, Pos (+ true_QTN, simulation-only)
##     ld_ws      markers x rho, rownames = markers
##     decay_sum  per-chromosome decay fit
##
##   pvals_<cell>_<engine>_<basis>_B<B>.rds    ONE association method's scans
##     p_obs      named numeric, names == map$marker, same order
##     p_perm     markers x B matrix of surrogate p-values
##     basis, engine, B, cell, seed provenance
##
## The panel is written once per cell and reused by every engine x basis, so GTs
## and ld_ws are stored once rather than repeated. Run:
##   Rscript module_sim_LDscnR/analyse_one_dataset.R panel_<cell>.rds pvals_<...>.rds
##
## THE NULLS. Two kinds, and they are NOT crossed: a home-field null is drawn from
## one engine's own model of structure, so it is only meaningful for that engine --
## it asks "does this method invent signal from the structure it claims to correct?"
## Running EMMAX against LFMM's latent basis (or the reverse) tests neither method's
## specificity, so those pairs are skipped. The method-agnostic nulls run on both.
##
##   genetic     [EMMAX ONLY] MVN(0, K) on the cell's mean kinship -- EMMAX's own
##               model of structure, so a clean surrogate here means EMMAX does not
##               manufacture false positives from kinship. The 10 per-file GRMs are
##               independent estimates of the same sampling structure (same patches,
##               same positions), so their mean is the cell's K -- and it keeps ONE
##               phenotype per draw, which the pooled design needs.
##   latent      [LFMM ONLY] MVN in the top-K subspace of the pooled genotype PCs --
##               LFMM's own model of structure, the exact analogue for that engine.
##   global_perm env shuffled among the 80 patches. Model-agnostic; destroys the
##               spatial structure entirely.
##   env_orth    THE SIM-SPECIFIC ONE, and method-agnostic -- the fair cross-engine
##               comparison. The 48x48 env field is shifted toroidally
##               and rotated/reflected, then read off at the sampled patches. A
##               torus shift is a RELABELLING of positions, so the empirical
##               autocorrelation function is preserved EXACTLY, while the field is
##               moved off the true one. Draws are screened to |cor| < COR_MAX
##               against the observed env and then residualised, so the surrogate
##               is exactly orthogonal to it: "same spatial structure, none of the
##               signal". This is what `region_perm` is for 3sp -- the arbiter --
##               and it has no analogue there because the lattice has no localities.
##   spatial     MVN(0, Gaussian kernel over patch coords). SUPPLEMENTARY -- the
##               framework demotes it after its gate failure on 3sp, and the sim
##               smoke test agrees (median 66 surrogate markers with C>0 against 0
##               for genetic/global_perm). Off by default; request it explicitly.
##   global_perm DROPPED from the defaults: dominated by env_orth, exactly as the
##               framework argues for the 3sp global permutation. env_orth preserves
##               the spatial autocorrelation a clinal signal would exploit, so
##               anything clearing env_orth clears a plain shuffle. Still available.
##
## Every surrogate is residualised against the observed env (resid(lm(s ~ Yobs))),
## as in the 3sp runner, so no null can accidentally carry the real signal.
##
## Engine x type pairs actually run:
##   emmax : genetic, global_perm, env_orth, spatial
##   lfmm  : latent,  global_perm, env_orth, spatial
##
## COST (measured: emmax_fast 0.10 s/scan, ld_cscore 0.83 s per 31k SNPs -- the
## C-score dominates, not the association; LFMM 12.5 s per file). One pooled draw
## over 10 chromosomes: EMMAX ~9.4 s, LFMM ~2.1 min => B=100 for ONE cell and ONE
## type is ~15.6 min (EMMAX) or ~3.5 h (LFMM). Over the reduced scope (3 (V,c)
## cells x 10 env = 30 cells): EMMAX 4 types ~31 core-h (~2.6 h on 12 cores) per
## subsample stage; LFMM 2 types at B=100 ~210 core-h (~17.5 h on 12 cores).
##
## Run (from the LDscnR-paper root):
##   Rscript module_sim_LDscnR/run_sim_nulls.R V c env [tag]
##   Rscript module_sim_LDscnR/run_sim_nulls.R 2 1 all            # all env of a (V,c)
## Env:
##   SIM_DATA       bundle dir (default $SIM_ROOT/regen_sim_data_nobgs)
##   SIM_NULL_OUT   output dir (default module_sim_LDscnR/results/nulls)
##   SIM_NULL_B     draws per type (default 100)
##   SIM_NULL_TYPES comma list (default genetic,latent,global_perm,env_orth,spatial)
##   SIM_NULL_ENGINES  emmax | emmax,lfmm   (default emmax)
##   SIM_NULL_CORES draws run concurrently (default 1)
## =====================================================================

suppressMessages({ library(data.table); library(parallel)
                   LDSCNR_PATH <- Sys.getenv("LDSCNR_PATH", "")
                   if (nzchar(LDSCNR_PATH) && dir.exists(LDSCNR_PATH))
                     devtools::load_all(LDSCNR_PATH, quiet = TRUE) else library(LDscnR) })

SIM_ROOT <- Sys.getenv("SIM_ROOT", "/Volumes/Nemo/Nemo_sim")
SIM_DATA <- Sys.getenv("SIM_DATA", file.path(SIM_ROOT, "regen_sim_data_nobgs"))
ENV_DIR  <- Sys.getenv("SIM_ENV",  file.path(SIM_ROOT, "env"))
OUTDIR   <- Sys.getenv("SIM_NULL_OUT", "module_sim_LDscnR/results/nulls")
B        <- as.integer(Sys.getenv("SIM_NULL_B", "100"))
CORES    <- as.integer(Sys.getenv("SIM_NULL_CORES", "1"))
## defaults are the framework's bases: home field per engine + the arbiter
TYPES    <- strsplit(Sys.getenv("SIM_NULL_TYPES", "genetic,latent,env_orth"), ",")[[1]]
ENGINES  <- strsplit(Sys.getenv("SIM_NULL_ENGINES", "emmax"), ",")[[1]]
if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)

K_LFMM <- 5L; K_LATENT <- 5L; SIDE <- 48L; SEED0 <- 7654321L; COR_MAX <- 0.3
CANON  <- c("genetic", "latent", "global_perm", "env_orth", "spatial")   # fixed seed map
## home-field nulls belong to ONE engine each; the rest are method-agnostic
HOME   <- c(genetic = "emmax", latent = "lfmm")
TYPES  <- intersect(TYPES, CANON); stopifnot(length(TYPES) > 0)
ENGINES <- intersect(ENGINES, c("emmax", "lfmm")); stopifnot(length(ENGINES) > 0)

a   <- commandArgs(trailingOnly = TRUE)
V   <- if (length(a) >= 1) a[1] else "2"
CC  <- if (length(a) >= 2) a[2] else "1"
ENV <- if (length(a) >= 3) a[3] else "1"
TAG <- if (length(a) >= 4) a[4] else "nobgs"
ENVS <- if (ENV == "all") as.character(1:10) else ENV

## ---- pool one cell (mirrors run_sim_LDscnR.R::pool_cell) --------------
pool_cell <- function(env) {
  files <- list.files(SIM_DATA, full.names = TRUE,
                      pattern = sprintf("^adapt_%s_chr[0-9]+_V%s_c%s_env%s\\.rds$", TAG, V, CC, env))
  if (!length(files)) return(NULL)
  files <- files[order(as.integer(sub(".*_chr([0-9]+)_.*", "\\1", basename(files))))]
  ## an incomplete cell would pool a smaller genome and silently change both the
  ## ld_w quantile and the BH denominator -- refuse rather than half-pool
  if (length(files) != 10L)
    stop(sprintf("expected 10 chromosome files for V%s_c%s_env%s, found %d", V, CC, env, length(files)))
  maps <- gts <- ldws <- decs <- prep <- mk_i <- vector("list", length(files))
  Kacc <- NULL; Yobs <- coords <- pops <- NULL
  for (i in seq_along(files)) {
    d <- readRDS(files[i]); m <- as.data.table(d$map)
    m[, `:=`(Chr = paste0("R", i, "_", Chr), marker = paste0("R", i, "_", marker))]
    G <- d$GTs; colnames(G) <- m$marker
    lw <- d$ld_ws; rownames(lw) <- m$marker
    ## decay fit carries the same R<i>_ prefix, so score_thresholds() downstream can
    ## match a region's chromosome to the curve it was fitted on
    ds <- as.data.table(d$LD_decay$decay_sum); ds[, Chr := paste0("R", i, "_", Chr)]
    prep[[i]] <- emmax_setup(G, d$GRM); mk_i[[i]] <- m$marker
    maps[[i]] <- m; gts[[i]] <- G; ldws[[i]] <- lw; decs[[i]] <- ds
    Kacc <- if (is.null(Kacc)) d$GRM else Kacc + d$GRM
    if (is.null(Yobs)) { Yobs <- d$env$env; coords <- cbind(d$env$x, d$env$y); pops <- d$env$pop }
  }
  ## flag_true_qtns() so the truth column matches what the consumer scores against
  map <- flag_true_qtns(rbindlist(maps, fill = TRUE))
  list(map = map, GTs = do.call(cbind, gts)[, map$marker],
       LDW = do.call(rbind, ldws)[map$marker, ], prep = prep, mk_i = mk_i,
       decay_sum = rbindlist(decs, fill = TRUE),
       Kmean = Kacc / length(files), Yobs = Yobs, coords = coords, pops = pops,
       n_files = length(files))
}

emmax_pooled <- function(P, Yv) unlist(lapply(seq_along(P$prep), function(i)
  stats::setNames(emmax_fast(P$prep[[i]], Yv), P$mk_i[[i]])))[P$map$marker]

## ---- surrogate generators -------------------------------------------
mvn_drawer <- function(K, y) {
  eK <- eigen(K, symmetric = TRUE); Lv <- pmax(eK$values, 0); Vk <- eK$vectors; n <- length(y)
  function() as.numeric(stats::resid(stats::lm(as.numeric(Vk %*% (sqrt(Lv) * stats::rnorm(n))) ~ y)))
}
latent_drawer <- function(GTs, y, K = K_LATENT) {
  Z <- scale(GTs); Z[is.na(Z)] <- 0
  ## Top-K genotype PCs (LFMM's basis). With n (~160) << p (~310k), take them from the
  ## n x n cross-product rather than svd() on the full matrix: LAPACK decomposes the
  ## whole thing regardless of `nu`, which here is a 160 x 310k job for five vectors.
  ## eigen(ZZ') gives the identical left singular vectors at a fraction of the cost.
  e  <- eigen(tcrossprod(Z), symmetric = TRUE)
  U  <- e$vectors[, seq_len(K), drop = FALSE]
  dd <- sqrt(pmax(e$values[seq_len(K)], 0))
  function() as.numeric(stats::resid(stats::lm(as.numeric(U %*% (dd * stats::rnorm(K))) ~ y)))
}
spatial_drawer <- function(coords, y) {
  Dm <- as.matrix(stats::dist(coords)); l <- stats::median(Dm[lower.tri(Dm)])
  mvn_drawer(exp(-0.5 * (Dm / l)^2), y)
}
global_perm_drawer <- function(pops, y) {
  up <- unique(pops)
  function() {
    v <- stats::setNames(sample(tapply(y, pops, `[`, 1)), up)     # patch values shuffled
    as.numeric(stats::resid(stats::lm(as.numeric(v[as.character(pops)]) ~ y)))
  }
}
## env_orth: torus shift + rotation/reflection of the 48x48 field. A shift is a
## relabelling of lattice positions, so the autocorrelation function is preserved
## exactly; screening on |cor| then residualising makes it orthogonal to the truth.
env_orth_drawer <- function(env_id, pops, y) {
  raw  <- scan(file.path(ENV_DIR, paste0("env_", env_id, ".txt")), what = character(), quiet = TRUE)
  raw  <- strsplit(raw, c("}{"), fixed = TRUE)[[1]]
  vals <- as.numeric(gsub("}}", "", gsub("{{", "", raw, fixed = TRUE), fixed = TRUE))
  lat  <- data.table(expand.grid(x = 1:SIDE, y = SIDE:1))[, `:=`(pop = seq_len(SIDE * SIDE), env = vals)]
  key  <- lat[, paste(x, y)]
  idx_of_pop <- match(pops, lat$pop)
  function() {
    for (try in 1:20) {
      sx <- sample.int(SIDE, 1L) - 1L; sy <- sample.int(SIDE, 1L) - 1L
      x2 <- ((lat$x - 1L + sx) %% SIDE) + 1L
      y2 <- ((lat$y - 1L + sy) %% SIDE) + 1L
      if (stats::runif(1) < 0.5) { tmp <- x2; x2 <- y2; y2 <- tmp }      # transpose (90 deg)
      if (stats::runif(1) < 0.5) x2 <- SIDE - x2 + 1L                    # reflect
      srg <- lat$env[match(paste(x2, y2), key)]
      v   <- srg[idx_of_pop]
      if (abs(stats::cor(v, y)) < COR_MAX) break
    }
    as.numeric(stats::resid(stats::lm(v ~ y)))
  }
}

## ---- engines ---------------------------------------------------------
## EMMAX returns p directly (no genomic control, by design -- see below). LFMM
## returns raw F, which the caller converts to p with a FIXED lambda.
##
## GENOMIC CONTROL, and why it is handled here rather than per scan:
##   EMMAX  none at all. The GRM is tuned so lambda sits in [1, 1.1], so there is
##          nothing to correct -- and correcting the observed while leaving the
##          surrogates raw (which is what the stored emx_p did, in 19 of 100
##          files) puts the two on different scales.
##   LFMM   residual inflation is expected and fine, but the correction must be
##          the SAME constant everywhere. lfmm2.test(genomic.control = TRUE)
##          re-estimates lambda inside every scan, so each surrogate would be
##          standardised by its own value -- that shrinks the null's spread and
##          flatters the observed. Instead lambda is estimated ONCE on the
##          observed F and applied unchanged to every permuted F.
##
## lambda is computed on the F scale -- median(F) / qf(0.5, df1, df2) -- not via
## the chi-square convention.
scan_engine <- function(engine, P, y) {
  if (engine == "emmax") return(emmax_pooled(P, y))
  ## LFMM: per chromosome file, then pooled; RAW F, genomic.control = FALSE
  unlist(lapply(seq_along(P$mk_i), function(i) {
    mk <- P$mk_i[[i]]; G <- P$GTs[, mk, drop = FALSE]
    tmp <- tempfile("lfmmsim_"); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
    gf <- file.path(tmp, "g.lfmm"); ef <- file.path(tmp, "e.env")
    LEA::write.lfmm(G, gf); LEA::write.env(y, ef)
    pr <- LEA::lfmm2(gf, ef, K = K_LFMM)
    stats::setNames(suppressWarnings(
      LEA::lfmm2.test(pr, gf, ef, genomic.control = FALSE, full = TRUE))$fscores, mk)
  }))[P$map$marker]
}

## lambda on the F scale, and the p-values it implies
lfmm_df2   <- function(n) n - K_LFMM - 1L                    # d = 1 predictor
emmax_df2  <- function(n) n - 2L
lambda_F   <- function(Fv, df2) stats::median(Fv) / stats::qf(0.5, 1, df2, lower.tail = FALSE)
p_from_F   <- function(Fv, lambda, df2) stats::pf(Fv / lambda, 1, df2, lower.tail = FALSE)
## same thing starting from p, for engines that only hand back p-values
lambda_p   <- function(pv, df2) lambda_F(stats::qf(pv, 1, df2, lower.tail = FALSE), df2)

## Diagnostic, never applied. Fixing lambda at the OBSERVED value is right only if
## the surrogates are inflated to a similar degree: if the observed scan's lambda is
## lifted by real polygenic signal, lambda_obs > lambda_surr and dividing the
## surrogates by it over-corrects them, quieting the null. Recording each
## surrogate's own lambda lets the analysis check that rather than assume it.

## ---- main loop: cell -> type -> engine -------------------------------
for (e in ENVS) {
  P <- pool_cell(e)
  if (is.null(P)) { message("no files for V", V, "_c", CC, "_env", e); next }
  message(sprintf("[cell V%s_c%s_env%s] %d files, %d markers, n=%d", V, CC, e,
                  P$n_files, nrow(P$map), length(P$Yobs)))

  GEN <- list(genetic     = mvn_drawer(P$Kmean, P$Yobs),
              latent      = latent_drawer(P$GTs, P$Yobs),
              global_perm = global_perm_drawer(P$pops, P$Yobs),
              env_orth    = env_orth_drawer(e, P$pops, P$Yobs),
              spatial     = spatial_drawer(P$coords, P$Yobs))

  ## Observed scans are recomputed HERE, through the identical code path the
  ## surrogates take, rather than read from the bundles. The stored emx_p was
  ## genomic-controlled whenever gif > 1.1 and the stored lfmm_p was corrected by
  ## LEA's own per-scan lambda, so neither is commensurable with a raw surrogate.
  ## EMMAX: raw p, no correction anywhere. LFMM: one lambda from the observed F,
  ## reused for every permuted F.
  df2 <- lfmm_df2(length(P$Yobs))
  obs <- list(); lam <- list()
  if ("emmax" %in% ENGINES) {
    obs$emmax <- emmax_pooled(P, P$Yobs)
    lam$emmax <- NA_real_                       # by design: no GC on either side
  }
  if ("lfmm" %in% ENGINES) {
    F_obs      <- scan_engine("lfmm", P, P$Yobs)
    lam$lfmm   <- lambda_F(F_obs, df2)
    obs$lfmm   <- p_from_F(F_obs, lam$lfmm, df2)
    message(sprintf("   LFMM observed lambda (F scale, df2=%d) = %.3f -- applied unchanged to every surrogate",
                    df2, lam$lfmm))
  }
  p_obs_eng <- obs

  ## panel, written once per cell and reused by every engine x basis
  cell_id  <- sprintf("V%s_c%s_env%s", V, CC, e)
  panel_f  <- file.path(OUTDIR, sprintf("panel_%s.rds", cell_id))
  if (!file.exists(panel_f)) {
    saveRDS(list(GTs = P$GTs, map = P$map, ld_ws = P$LDW, decay_sum = P$decay_sum,
                 cell = cell_id, n_ind = length(P$Yobs),
                 env_obs = P$Yobs, coords = P$coords),
            panel_f)
    message(sprintf("   panel: %d markers x %d individuals, %d chromosomes -> %s (%.0f MB)",
                    nrow(P$map), length(P$Yobs), uniqueN(P$map$Chr),
                    basename(panel_f), file.size(panel_f) / 1e6))
  } else message("   panel: reusing ", basename(panel_f))

  for (ty in TYPES) {
    ## skip the type entirely if no active engine takes it (e.g. `latent` with
    ## ENGINES=emmax) -- otherwise the surrogate draws are computed for nothing
    engs_for_ty <- if (ty %in% names(HOME)) intersect(ENGINES, HOME[[ty]]) else ENGINES
    if (!length(engs_for_ty)) next

    ## identical draws across engines and run order (seeded per type and draw)
    phen <- lapply(seq_len(B), function(b) {
      set.seed(SEED0 + 100000L * match(ty, CANON) + b); GEN[[ty]]() })

    for (eng in intersect(c("emmax", "lfmm"), ENGINES)) {
      ## a home-field null is only meaningful for its own engine -- skip the cross
      if (ty %in% names(HOME) && HOME[[ty]] != eng) next
      outf <- file.path(OUTDIR, sprintf("pvals_%s_%s_%s_B%d.rds", cell_id, eng, ty, B))
      if (file.exists(outf)) { message("   ", basename(outf), " exists -> skip"); next }
      t0 <- Sys.time()

      pl <- mclapply(phen, function(y) tryCatch({
                       v <- scan_engine(eng, P, y)
                       ## LFMM surrogates are corrected by the OBSERVED lambda, not
                       ## their own -- the whole point of fixing it
                       if (eng == "lfmm") v <- p_from_F(v, lam$lfmm, df2)
                       v
                     }, error = function(err) err),
                     mc.cores = CORES, mc.preschedule = FALSE)
      ok <- vapply(pl, is.numeric, logical(1))
      if (any(!ok)) message(sprintf("   [%s|%s] dropped %d failed draws", eng, ty, sum(!ok)))
      pl <- pl[ok]

      ## markers x B matrix: marker names are stored ONCE, in `markers`, rather than
      ## repeated on every surrogate vector
      P_surr <- matrix(unlist(pl, use.names = FALSE), nrow = nrow(P$map), ncol = length(pl))

      ## per-surrogate lambda: measured, reported, NOT applied (see above)
      df2_eng   <- if (eng == "lfmm") df2 else emmax_df2(length(P$Yobs))
      lam_surr  <- apply(P_surr, 2, lambda_p, df2 = df2_eng)
      lam_obs_e <- if (eng == "lfmm") lam$lfmm else lambda_p(obs$emmax, df2_eng)

      ## names carried on p_obs so analyse_one_dataset.R can VERIFY the ordering
      ## rather than trust it -- a right-length vector in the wrong order is
      ## silently wrong, and this is where that would originate
      p_obs <- p_obs_eng[[eng]]
      stopifnot(identical(names(p_obs), P$map$marker), nrow(P_surr) == nrow(P$map))
      saveRDS(list(p_obs = p_obs, p_perm = P_surr,
                   basis = ty, engine = eng, B = ncol(P_surr), cell = cell_id,
                   panel = basename(panel_f),
                   ## how genomic control was handled, so the analysis need not guess
                   lambda_obs = lam_obs_e, lambda_surr = lam_surr,
                   gc_mode = if (eng == "emmax") "none"
                             else sprintf("fixed lambda_obs = %.4f on F, df=(1,%d)", lam$lfmm, df2),
                   seed0 = SEED0, canon_index = match(ty, CANON)),
              outf)

      ## gate diagnostic, framework 4: median per-surrogate count of markers below a
      ## nominal cutoff, against the observed. Reported, not applied -- no threshold
      ## in this script decides anything.
      med <- stats::median(colSums(P_surr < 1e-4))
      message(sprintf("   [%-5s | %-11s] B=%d ; median surrogate p<1e-4 = %.0f (obs %d) ; %.1f min -> %s",
                      eng, ty, ncol(P_surr), med, sum(p_obs < 1e-4),
                      as.numeric(Sys.time() - t0, units = "mins"), basename(outf)))
      message(sprintf("       lambda: observed %.3f | surrogates median %.3f (%.3f-%.3f)%s",
                      lam_obs_e, stats::median(lam_surr), min(lam_surr), max(lam_surr),
                      if (eng == "lfmm") "  [obs lambda applied to all]" else "  [no GC]"))
      message(sprintf("       Rscript module_sim_LDscnR/analyse_one_dataset.R %s %s",
                      basename(panel_f), basename(outf)))
    }
  }
  rm(P); gc()
}
