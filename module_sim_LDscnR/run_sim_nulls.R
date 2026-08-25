## =====================================================================
## module_sim_LDscnR / run_sim_nulls.R
##
## Permutation / surrogate nulls for the SIM benchmark -- the sim counterpart of
## ~/3sp_lfmm_perm/run_paired_nulls_3sp.R, same paired two-engine design, same
## ld_null bundle format, same per-(type,b) seeding so EMMAX and LFMM see the
## IDENTICAL surrogate phenotypes.
##
## Unit of work is one (V, c, env) CELL: the 10 chromosome files are pooled (as in
## run_sim_LDscnR.R), each keeps its own saved GRM, and one surrogate phenotype is
## pushed through all ten -> one pooled surrogate C. The ten files share their 80
## patches and their env values exactly, so a lattice-level surrogate maps across
## all of them unchanged.
##
## THE NULLS
##   genetic     MVN(0, K) on the cell's mean kinship. EMMAX's home-field
##               specificity null; a high-rank stress test for LFMM (whose 5
##               factors cannot span full-rank K). The 10 per-file GRMs are
##               independent estimates of the same sampling structure (same
##               patches, same positions), so their mean is the cell's K -- and it
##               keeps ONE phenotype per draw, which the pooled design needs.
##   latent      MVN in the top-K subspace of the pooled genotype PCs. LFMM's
##               home-field null and the exact analogue of K-MVN for EMMAX.
##   global_perm env shuffled among the 80 patches. Model-agnostic; destroys the
##               spatial structure entirely.
##   env_orth    THE SIM-SPECIFIC ONE. The 48x48 env field is shifted toroidally
##               and rotated/reflected, then read off at the sampled patches. A
##               torus shift is a RELABELLING of positions, so the empirical
##               autocorrelation function is preserved EXACTLY, while the field is
##               moved off the true one. Draws are screened to |cor| < COR_MAX
##               against the observed env and then residualised, so the surrogate
##               is exactly orthogonal to it: "same spatial structure, none of the
##               signal". This is what `region_perm` is for 3sp -- the arbiter --
##               and it has no analogue there because the lattice has no localities.
##   spatial     MVN(0, Gaussian kernel over patch coords). Matches the null
##               run_sim_LDscnR.R already uses; kept for continuity.
##
## Every surrogate is residualised against the observed env (resid(lm(s ~ Yobs))),
## as in the 3sp runner, so no null can accidentally carry the real signal.
##
## COST (measured: emmax_fast 0.10 s/scan, ld_cscore 0.83 s per 31k SNPs; the
## C-score dominates, not the association). One pooled draw over 10 chromosomes
## ~9.4 s => B=100 for ONE cell and ONE type ~15.6 min. All 90 cells x 4 types
## ~94 core-hours ~ 9-10 h on 10 workers, per subsample stage. LFMM is ~25 s per
## file-scan, i.e. ~250x more per draw -- restrict it to a few cells (SIM_NULL_CELLS)
## or a small B; it is OFF by default.
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
TYPES    <- strsplit(Sys.getenv("SIM_NULL_TYPES", "genetic,latent,global_perm,env_orth,spatial"), ",")[[1]]
ENGINES  <- strsplit(Sys.getenv("SIM_NULL_ENGINES", "emmax"), ",")[[1]]
if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)

PAR    <- list(alpha = 0.05, qstar = seq(0, 0.95, by = 0.05))
K_LFMM <- 5L; K_LATENT <- 5L; SIDE <- 48L; SEED0 <- 7654321L; COR_MAX <- 0.3
CANON  <- c("genetic", "latent", "global_perm", "env_orth", "spatial")   # fixed seed map
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
  maps <- gts <- ldws <- prep <- mk_i <- vector("list", length(files))
  Kacc <- NULL; Yobs <- coords <- pops <- NULL
  for (i in seq_along(files)) {
    d <- readRDS(files[i]); m <- as.data.table(d$map)
    m[, `:=`(Chr = paste0("R", i, "_", Chr), marker = paste0("R", i, "_", marker))]
    G <- d$GTs; colnames(G) <- m$marker
    lw <- d$ld_ws; rownames(lw) <- m$marker
    prep[[i]] <- emmax_setup(G, d$GRM); mk_i[[i]] <- m$marker
    maps[[i]] <- m; gts[[i]] <- G; ldws[[i]] <- lw
    Kacc <- if (is.null(Kacc)) d$GRM else Kacc + d$GRM
    if (is.null(Yobs)) { Yobs <- d$env$env; coords <- cbind(d$env$x, d$env$y); pops <- d$env$pop }
  }
  map <- rbindlist(maps, fill = TRUE)
  list(map = map, GTs = do.call(cbind, gts)[, map$marker],
       LDW = do.call(rbind, ldws)[map$marker, ], prep = prep, mk_i = mk_i,
       Kmean = Kacc / length(files), Yobs = Yobs, coords = coords, pops = pops,
       n_files = length(files))
}

emmax_pooled <- function(P, Yv) unlist(lapply(seq_along(P$prep), function(i)
  stats::setNames(emmax_fast(P$prep[[i]], Yv), P$mk_i[[i]])))[P$map$marker]

cscore <- function(pv, LDW) ld_cscore(pv, LDW, alpha = PAR$alpha, rho = colnames(LDW), qstar = PAR$qstar)

## ---- surrogate generators -------------------------------------------
mvn_drawer <- function(K, y) {
  eK <- eigen(K, symmetric = TRUE); Lv <- pmax(eK$values, 0); Vk <- eK$vectors; n <- length(y)
  function() as.numeric(stats::resid(stats::lm(as.numeric(Vk %*% (sqrt(Lv) * stats::rnorm(n))) ~ y)))
}
latent_drawer <- function(GTs, y, K = K_LATENT) {
  Z <- scale(GTs); Z[is.na(Z)] <- 0
  sv <- svd(Z, nu = K, nv = 0)                      # top-K genotype PCs, LFMM's basis
  U <- sv$u; dd <- sv$d[seq_len(K)]
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
scan_engine <- function(engine, P, y) {
  if (engine == "emmax") return(cscore(emmax_pooled(P, y), P$LDW))
  ## LFMM: per chromosome file, then pooled (mirrors the observed pipeline)
  pv <- unlist(lapply(seq_along(P$mk_i), function(i) {
    mk <- P$mk_i[[i]]; G <- P$GTs[, mk, drop = FALSE]
    tmp <- tempfile("lfmmsim_"); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
    gf <- file.path(tmp, "g.lfmm"); ef <- file.path(tmp, "e.env")
    LEA::write.lfmm(G, gf); LEA::write.env(y, ef)
    pr <- LEA::lfmm2(gf, ef, K = K_LFMM)
    stats::setNames(suppressWarnings(
      LEA::lfmm2.test(pr, gf, ef, genomic.control = TRUE, full = TRUE))$pvalues, mk)
  }))[P$map$marker]
  cscore(pv, P$LDW)
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

  ## Observed C per ENGINE, from that engine's own saved scan -- EMMAX's emx_p and
  ## LFMM's lfmm_p are different scans, so scoring LFMM surrogates against EMMAX's
  ## observed C would compare across engines and inflate/deflate everything.
  C_obs_eng <- lapply(c(emmax = "emx_p", lfmm = "lfmm_p"), function(col) {
    C <- cscore(P$map[[col]], P$LDW); names(C) <- P$map$marker; C })

  for (ty in TYPES) {
    ## identical draws across engines and run order (seeded per type and draw)
    phen <- lapply(seq_len(B), function(b) {
      set.seed(SEED0 + 100000L * match(ty, CANON) + b); GEN[[ty]]() })

    for (eng in intersect(c("emmax", "lfmm"), ENGINES)) {
      outf <- file.path(OUTDIR, sprintf("null_%s_%s_V%s_c%s_env%s.rds", eng, ty, V, CC, e))
      if (file.exists(outf)) { message("   ", basename(outf), " exists -> skip"); next }
      t0 <- Sys.time()
      Cl <- mclapply(phen, function(y) tryCatch({ C <- scan_engine(eng, P, y); C[C > 0] },
                                               error = function(err) err),
                     mc.cores = CORES, mc.preschedule = FALSE)
      ok <- vapply(Cl, is.numeric, logical(1))
      if (any(!ok)) message(sprintf("   [%s|%s] dropped %d failed draws", eng, ty, sum(!ok)))
      Cl <- Cl[ok]
      C_obs <- C_obs_eng[[eng]]
      null <- structure(list(C_obs = C_obs, C_surr = Cl,
                             universe = unique(c(names(C_obs)[C_obs > 0], unlist(lapply(Cl, names)))),
                             basis = paste0(eng, "_", ty), engine = eng, B = length(Cl),
                             cell = list(V = V, c = CC, env = e, tag = TAG, n_files = P$n_files),
                             params = PAR), class = "ld_null")
      saveRDS(null, outf)
      message(sprintf("   [%-5s | %-11s] B=%d ; median null C>0/surr = %.0f (obs %d) ; %.1f min -> %s",
                      eng, ty, length(Cl),
                      if (length(Cl)) stats::median(vapply(Cl, length, integer(1))) else NA_real_,
                      sum(C_obs > 0), as.numeric(Sys.time() - t0, units = "mins"), basename(outf)))
    }
  }
  rm(P); gc()
}
