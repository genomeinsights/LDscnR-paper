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
## OUTPUT CONTRACT -- deliberately dataset-agnostic, so that ONE downstream analysis
## consumes the simulations and the empirical data identically. Two kinds of file:
##
##   cell_<id>.rds     context, written once per cell:
##                       markers    character, the pooled marker order
##                       map        marker, Chr, Pos, type, true_QTN (truth is
##                                  simulation-only; absent/NA for empirical data)
##                       ld_ws      markers x rho matrix, the local-LD statistic
##                       decay_sum  per-chromosome decay fit (b, c, a_pred ...)
##                       dataset/cell identifiers
##
##   pnull_<engine>_<null>_<id>.rds   the scans:
##                       p_obs      numeric, one per marker, SAME ORDER as markers
##                       P_surr     markers x B matrix of surrogate p-values
##                       engine, null_type, B, cell, seed provenance
##
## The empirical side (~/3sp_lfmm_perm/run_paired_nulls_3sp.R) currently emits
## ld_null bundles with C-scores already computed; to share the downstream it needs
## the same split -- context once, p-values per basis. Everything downstream of the scan -- the (rho, q*, alpha) grid, tau_C,
## clustering, l_min, the region test -- is the analysis's business, and keeping it
## there means those choices can change without re-running a single scan. It also
## removes a real hazard: this runner previously computed its own C-score with
## alpha = 0.05 while run_sim_LDscnR.R sweeps alpha over four values, so observed and
## surrogate C-scores were not the same statistic.
##
## Unit of work is one (V, c, env) CELL: the 10 chromosome files are pooled (as in
## run_sim_LDscnR.R), each keeps its own saved GRM, and one surrogate phenotype is
## pushed through all ten -> one pooled surrogate C. The ten files share their 80
## patches and their env values exactly, so a lattice-level surrogate maps across
## all of them unchanged.
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
  map <- rbindlist(maps, fill = TRUE)
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
## returns one p per marker, in pooled marker order
scan_engine <- function(engine, P, y) {
  if (engine == "emmax") return(emmax_pooled(P, y))
  ## LFMM: per chromosome file, then pooled (mirrors the observed pipeline)
  unlist(lapply(seq_along(P$mk_i), function(i) {
    mk <- P$mk_i[[i]]; G <- P$GTs[, mk, drop = FALSE]
    tmp <- tempfile("lfmmsim_"); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
    gf <- file.path(tmp, "g.lfmm"); ef <- file.path(tmp, "e.env")
    LEA::write.lfmm(G, gf); LEA::write.env(y, ef)
    pr <- LEA::lfmm2(gf, ef, K = K_LFMM)
    stats::setNames(suppressWarnings(
      LEA::lfmm2.test(pr, gf, ef, genomic.control = TRUE, full = TRUE))$pvalues, mk)
  }))[P$map$marker]
}

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

  ## Observed p per ENGINE, from that engine's own saved scan. EMMAX's emx_p and
  ## LFMM's lfmm_p are different scans and must never be mixed.
  p_obs_eng <- lapply(c(emmax = "emx_p", lfmm = "lfmm_p"),
                      function(col) stats::setNames(P$map[[col]], P$map$marker))

  ## context, written once per cell: everything the downstream needs that is NOT a
  ## scan. Kept out of the per-basis files so ld_ws (markers x rho, ~50 MB) is stored
  ## once rather than repeated for every engine x null combination.
  cell_id  <- sprintf("V%s_c%s_env%s", V, CC, e)
  ctx_file <- file.path(OUTDIR, paste0("cell_", cell_id, ".rds"))
  if (!file.exists(ctx_file)) {
    saveRDS(list(markers = P$map$marker,
                 map = P$map[, .(marker, Chr, Pos, type, true_QTN)],
                 ld_ws = P$LDW, decay_sum = P$decay_sum,
                 dataset = "sim", cell_id = cell_id,
                 cell = list(V = V, c = CC, env = e, tag = TAG, n_files = P$n_files),
                 n_ind = length(P$Yobs), env_obs = P$Yobs, coords = P$coords),
            ctx_file)
    message("   wrote context ", basename(ctx_file))
  }

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
      outf <- file.path(OUTDIR, sprintf("pnull_%s_%s_%s.rds", eng, ty, cell_id))
      if (file.exists(outf)) { message("   ", basename(outf), " exists -> skip"); next }
      t0 <- Sys.time()

      pl <- mclapply(phen, function(y) tryCatch(scan_engine(eng, P, y), error = function(err) err),
                     mc.cores = CORES, mc.preschedule = FALSE)
      ok <- vapply(pl, is.numeric, logical(1))
      if (any(!ok)) message(sprintf("   [%s|%s] dropped %d failed draws", eng, ty, sum(!ok)))
      pl <- pl[ok]

      ## markers x B matrix: marker names are stored ONCE, in `markers`, rather than
      ## repeated on every surrogate vector
      P_surr <- matrix(unlist(pl, use.names = FALSE), nrow = nrow(P$map), ncol = length(pl))

      saveRDS(list(p_obs = unname(p_obs_eng[[eng]]), P_surr = P_surr,
                   engine = eng, null_type = ty, B = ncol(P_surr),
                   dataset = "sim", cell_id = cell_id,
                   cell = list(V = V, c = CC, env = e, tag = TAG, n_files = P$n_files),
                   context = paste0("cell_", cell_id, ".rds"),
                   seed0 = SEED0, canon_index = match(ty, CANON)),
              outf)

      ## gate diagnostic, framework 4: median per-surrogate count of markers below a
      ## nominal cutoff, against the observed. Reported, not applied -- no threshold
      ## in this script decides anything.
      med <- stats::median(colSums(P_surr < 1e-4))
      message(sprintf("   [%-5s | %-11s] B=%d ; median surrogate p<1e-4 = %.0f (obs %d) ; %.1f min -> %s",
                      eng, ty, ncol(P_surr), med, sum(p_obs_eng[[eng]] < 1e-4),
                      as.numeric(Sys.time() - t0, units = "mins"), basename(outf)))
    }
  }
  rm(P); gc()
}
