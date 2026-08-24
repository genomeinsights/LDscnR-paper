## =====================================================================
## module_sim_LDscnR / parse_and_regen_sim_data.R
##
## RAW Nemo tarballs  ->  regen-format bundles, in ONE pass.
##
## Merges the two historical steps:
##   Stage A  (was /Volumes/Nemo/Nemo_sim/Parse_sim_data.R)
##            unpack .tgz -> Nemo map/geno -> real recombination map ->
##            marker ids, env values, MAF filter, Va, true_QTN / focal_QTN truth
##   Stage B  (was module_sim_LDscnR/regen_sim_data.R)
##            compute_LD_decay (corr, ld_w in place) -> GRM basis -> emmax (+GC)
##            -> LFMM (K=5, GC) -> save bundle
##   Stage A2 (new)
##            population-genetic, local-adaptation and background-selection
##            summaries, written to their OWN files in $SIM_POPGEN (see the
##            Stage A2 section for the metrics and what each one assumes)
##
## The output is byte-compatible with regen_sim_data.R's bundles, i.e. exactly
## what run_sim_LDscnR.R reads: list(GTs, map, env, LD_decay, ld_ws, GRM,
## grm_markers, grm_method, complexity_reduction, emx_gif).
##
## Everything that the OLD parse script computed and the regen step then
## RECOMPUTED (LD_decay, ld_w, GRM, EMMAX, LFMM) is computed ONCE here, with the
## canonical regen settings. The only knock-on: `rho_d` / `ld_rel` (QTN-chromosome
## truth annotations) are now derived from the canonical corr-based decay fit
## rather than the old parse-time "r"/slide=2000 fit. They are not used by
## run_sim_LDscnR.R; they are kept for backward compatibility of the map.
##
## Run from the LDscnR-paper root (heavy; one subprocess per file recommended):
##   Rscript module_sim_LDscnR/parse_and_regen_sim_data.R  V c env chr [tag]
##   Rscript module_sim_LDscnR/parse_and_regen_sim_data.R  0.5 1 all all bgs
## defaults: V=2 c=1 env=1 chr=1 tag=bgs  ("all" allowed for env and chr)
##
## Env vars:
##   SIM_RAW     raw tarball folder        (default /Volumes/Nemo/Nemo_sim/bgs2)
##   SIM_OUT     bundle output folder      (default /Volumes/Nemo/Nemo_sim/regen_sim_data_bgs2)
##   SIM_MAPS    recombination maps        (default .../maps_500kb_with_allelic_values/chromosome_maps_500kb_rds)
##   SIM_ENV     env_<n>.txt folder        (default /Volumes/Nemo/Nemo_sim/env)
##   SIM_TMP     scratch root for untar    (default tempdir(); one subfolder per file)
##   SIM_GRM     "complexity_chain" (default) or "ld_w_threshold"
##   SIM_CORES   cores for compute_LD_decay (default 1)
##   SIM_PARSED  optional folder; if set, the Stage-A parsed bundle is ALSO saved
##               there (old parsed_sim_data2 format, pre-EMMAX/LFMM)
##   SIM_POPGEN  popgen/local-adaptation output (default /Volumes/Nemo/Nemo_sim/popgen_sim_data);
##               per file: <stem>.rds (summary + per-class diversity + per-patch
##               table + BGS window table + ini settings) and <stem>_summary.csv
##               (one row, for rbind)
##   SIM_WIN     BGS window size in bp (default 5e5, the map's own resolution).
##               The grid is fixed and shared across runs, so the matched
##               B_obs = pi_bgs / pi_nobgs per window can be collated later.
## =====================================================================

suppressMessages({ library(data.table); library(parallel); library(SNPRelate); library(LEA)
                   library(igraph)
                   devtools::load_all("/Users/petrikem/gitlab/LDscnR", quiet = TRUE) })

RAW_DIR    <- Sys.getenv("SIM_RAW",  "/Volumes/Nemo/Nemo_sim/bgs2")
OUT_DIR    <- Sys.getenv("SIM_OUT",  "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs2")
MAP_DIR    <- Sys.getenv("SIM_MAPS", "/Volumes/Nemo/Nemo_sim/maps_500kb_with_allelic_values/chromosome_maps_500kb_rds")
ENV_DIR    <- Sys.getenv("SIM_ENV",  "/Volumes/Nemo/Nemo_sim/env")
TMP_ROOT   <- Sys.getenv("SIM_TMP",  tempdir())
PARSED_DIR <- Sys.getenv("SIM_PARSED", "")
POPGEN_DIR <- Sys.getenv("SIM_POPGEN", "/Volumes/Nemo/Nemo_sim/popgen_sim_data")
WIN_BP     <- as.numeric(Sys.getenv("SIM_WIN", "5e5"))   # BGS window grid (map resolution)
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
if (nzchar(PARSED_DIR) && !dir.exists(PARSED_DIR)) dir.create(PARSED_DIR, recursive = TRUE)
if (!dir.exists(POPGEN_DIR)) dir.create(POPGEN_DIR, recursive = TRUE)

## ---- Stage-A settings (from Parse_sim_data.R) ------------------------
## Nemo samples 80 patches x 4 individuals = 320 (files_sample_patch in the .ini),
## so this keeps 2 individuals from each of 80 patches -- NOT one per patch, which
## is what the original Parse_sim_data.R comment claimed. Same in the old nobgs
## tarballs; the behaviour is unchanged here, only the description is corrected.
KEEP_INDS <- seq(1, 320, by = 2)
MIN_MAF   <- 0.05
SIDE      <- 48                    # 48 x 48 population lattice

## ---- Stage-B canonical settings (from regen_sim_data.R) --------------
RHO_GRID   <- c(seq(0.05, 0.95, by = 0.05), 0.99)          # ld_w columns
## ld_w computed IN PLACE (ld_w_rho) from the edge lists built for the decay fit,
## which are then dropped (keep_el = FALSE) -- no edge lists saved or retained. The
## complexity chain rebuilds what it needs per chromosome from `gds` on the fly
## (ld_complexity_reduction(gds = ...)), so nothing large is kept between steps.
DECAY_ARGS <- list(min_maf_decay = 0.1, q = 0.95, n_sub_bg = 5000, n_win_decay = 5,
                   overlap = 0.5, max_SNPs_decay = Inf, prob_robust = 0.95,
                   max_pairs = 5000, ld_method = "corr", n_strata = 20, keep_el = FALSE,
                   slide = 1000, rho_targets = c(0.99),
                   cores = as.integer(Sys.getenv("SIM_CORES", "1")), ld_w_rho = RHO_GRID)
GRM_METHOD <- Sys.getenv("SIM_GRM", "complexity_chain")    # or "ld_w_threshold"
CR_RHO     <- 0.5                                           # ld_complexity_reduction rho
PRUNE_ARGS <- list(ld_w_col = "ld_w_095", ld_w_threshold = 0.025, score_threshold = 0.80,
                   min_r2 = 0.2, distance_threshold = 5e5, compute_unflagged_eMLG = FALSE)
GRM_LDW_THRESHOLD <- 0.02                                  # used only if GRM_METHOD == "ld_w_threshold"
LFMM_K     <- 5L

## =====================================================================
## Stage-A helpers (verbatim logic from Parse_sim_data.R)
## =====================================================================

## max r2 of every marker with the best QTN on its own chromosome, + which QTN
## and how far away it is -- the truth annotation the TP/FP benchmark scores against.
focal_QTN <- function(map, GTs, qtn_col = "true_QTN") {
  rbindlist(lapply(unique(map$Chr), function(ch) {
    map_chr     <- map[Chr == ch]
    markers_chr <- map_chr$marker
    qtn_chr     <- map_chr[get(qtn_col) == TRUE, marker]

    if (length(qtn_chr) == 0L) {
      return(data.table(Chr = ch, marker = markers_chr, max_LD_with_QTN = 0,
                        focal_QTN = NA_character_, bp_to_focal_QTN = NA_real_))
    }

    r2        <- cor(GTs[, qtn_chr], GTs[, markers_chr], use = "pairwise.complete.obs")^2
    max_ld    <- apply(r2, 2, max, na.rm = TRUE)
    focal_idx <- max.col(t(r2), ties.method = "first")
    focal_qtn <- qtn_chr[focal_idx]

    pos_vec <- map_chr$Pos; names(pos_vec) <- map_chr$marker

    data.table(Chr = ch, marker = markers_chr, max_LD_with_QTN = max_ld,
               focal_QTN = focal_qtn,
               bp_to_focal_QTN = abs(pos_vec[markers_chr] - pos_vec[focal_qtn]))
  }), fill = TRUE)
}

## per-QTN additive variance 2 p (1-p) alpha^2
get_va <- function(map, GTs, qtn_rows) {
  if (length(qtn_rows) == 0L) return(numeric(0))
  sapply(qtn_rows, function(i) {
    a  <- map$allelic_values[i]
    gt <- GTs[, i, drop = TRUE]
    p  <- mean(gt) / 2
    2 * p * (1 - p) * a^2
  })
}

extract_params <- function(base_name) {
  params <- data.table(t(strsplit(base_name, "_", fixed = TRUE)[[1]]))[, 1:6]
  setnames(params, c("sim", "bgs", "Chr", "V", "c", "env"))
  params[, env := gsub(".tgz", "", env)][]
}

## Unpack into a private tmp folder. The number of leading path components to
## strip is DETECTED from the archive itself (the old Nemo_out_* tarballs carry
## the full scratch/... path -> 4; the bgs2 tarballs are already flat -> 0),
## instead of the old hard-coded 4-with-fallback.
unpack_sim <- function(file_gz, tmp_dir) {
  if (dir.exists(tmp_dir)) unlink(tmp_dir, recursive = TRUE)
  dir.create(tmp_dir, recursive = TRUE)

  entries <- utils::untar(file_gz, list = TRUE)
  geno_e  <- grep("snp_geno$", entries, value = TRUE)[1]
  if (is.na(geno_e)) stop("no *.snp_geno entry in ", basename(file_gz))
  n_strip <- length(strsplit(geno_e, "/", fixed = TRUE)[[1]]) - 2L   # .../GENO/<file>

  st <- system2("tar", c("-xzf", shQuote(file_gz),
                         paste0("--strip-components=", n_strip),
                         "-C", shQuote(tmp_dir)))
  if (st != 0) stop("tar failed on ", basename(file_gz))
  invisible(tmp_dir)
}

## =====================================================================
## Stage A2: population-genetic + local-adaptation summaries
##
## Computed from the RAW sample (all 320 individuals = 2 per sampled patch,
## before the MAF filter) plus Nemo's own whole-metapopulation stats file and
## the .ini, and written to their own files -- they are diagnostics ABOUT each
## simulation, not inputs to the outlier benchmark, so they stay out of the
## analysis bundle.
##
## Local adaptation: the sim's fitness function is known exactly from the .ini
## (`selection_model (gaussian, direct)`, `selection_variance` = V, and
## `selection_local_optima` = the same env file), i.e.
##     w(z, theta) = exp(-(z - theta)^2 / (2 V))
## with the local optimum theta_p equal to the patch's env value. That lets us
## report local adaptation on the fitness scale rather than only as a
## correlation:
##   la_cor2     cor(theta, z)^2 -- the legacy statistic, kept for continuity
##   la_slope    slope of z_bar ~ theta; perfect tracking = 1. cor2 is
##               scale-free, so a cline of slope 0.1 and one of slope 1 both
##               give cor2 ~ 1 -- the slope is what separates them
##   w_home      mean w(z_p, theta_p)             ("adaptedness")
##   w_away      mean w(z_p, theta_q), q != p     (foreign optima)
##   delta_SA    w_home - w_away -- the sympatric-allopatric contrast, the
##               standard definition of local adaptation (Kawecki & Ebert 2004;
##               Blanquart et al. 2013 local-vs-foreign). In fitness units and
##               scaled by this run's own V, so it is comparable across the
##               V grid in a way cor2 is not (the same phenotypic mismatch is
##               near-lethal at V = 0.5 and nearly free at V = 2)
##   *_vc        the same, corrected for within-patch variance Va_p via
##               E[w] = sqrt(V/(V+Va_p)) exp(-(z-theta)^2 / (2(V+Va_p))) --
##               patch means alone overstate mean fitness (Jensen)
##   lag_mse     mean (z_bar - theta)^2, and load_sel = lag_mse / (2V) = the
##               mean log-fitness deficit
##   qst_fst     Nemo's Qst - Fst, the classic divergent-selection contrast
## =====================================================================

## --- Nemo .ini: the few settings the summaries depend on ---------------
read_ini <- function(ini_file) {
  ln <- readLines(ini_file, warn = FALSE)
  get1 <- function(key) {
    hit <- grep(paste0("^", key, "[ \t]"), ln, value = TRUE)
    if (!length(hit)) return(NA_character_)
    trimws(sub(paste0("^", key, "[ \t]+"), "", hit[1]))
  }
  list(selection_variance = as.numeric(get1("selection_variance")),
       selection_model    = get1("selection_model"),
       patch_capacity     = as.numeric(get1("patch_capacity")),
       patch_number       = as.numeric(get1("patch_number")),
       quanti_loci        = as.numeric(get1("quanti_loci")),
       ntrl_loci          = as.numeric(get1("ntrl_loci")),
       generations        = as.numeric(get1("generations")),
       random_seed        = get1("random_seed"))
}

## --- Nemo stats file: last generation, scalars + per-patch vectors -----
read_nemo_stats <- function(stats_file) {
  st <- fread(stats_file, na.strings = c("nan", "-nan", "inf", "-inf", "NA"))
  st <- st[.N]                                    ## final generation
  patch_cols <- grep("^adlt\\.q1\\.p[0-9]+$",    names(st), value = TRUE)
  va_cols    <- grep("^adlt\\.Va\\.q1\\.p[0-9]+$", names(st), value = TRUE)
  ord   <- order(as.integer(sub("^adlt\\.q1\\.p",     "", patch_cols)))
  ord_v <- order(as.integer(sub("^adlt\\.Va\\.q1\\.p", "", va_cols)))
  list(
    scalars = st[, .SD, .SDcols = setdiff(names(st), c(patch_cols, va_cols))],
    z_bar   = as.numeric(st[, ..patch_cols][, ..ord]),
    Va_p    = if (length(va_cols)) as.numeric(st[, ..va_cols][, ..ord_v]) else NULL
  )
}

## --- local adaptation from patch trait means + local optima -----------
local_adaptation <- function(z_bar, theta, V, Va_p = NULL, prefix = "") {
  ok <- is.finite(z_bar) & is.finite(theta)
  z <- z_bar[ok]; th <- theta[ok]
  out <- list()
  nm  <- function(x) paste0(prefix, x)

  if (length(z) < 3 || !is.finite(V) || V <= 0) {
    for (k in c("n_patch", "la_cor2", "la_slope", "la_intercept", "w_home", "w_away",
                "delta_SA", "w_home_vc", "w_away_vc", "delta_SA_vc", "lag_mse", "load_sel"))
      out[[nm(k)]] <- NA_real_
    return(out)
  }

  fit <- stats::lm(z ~ th)
  out[[nm("n_patch")]]      <- length(z)
  out[[nm("la_cor2")]]      <- stats::cor(th, z)^2          ## legacy statistic
  out[[nm("la_slope")]]     <- unname(stats::coef(fit)[2])  ## 1 = tracks the gradient fully
  out[[nm("la_intercept")]] <- unname(stats::coef(fit)[1])

  ## Gaussian fitness on patch means: home = own optimum, away = every other
  ## patch's optimum. delta_SA is the sympatric-allopatric contrast.
  W <- exp(-outer(z, th, function(a, b) (a - b)^2) / (2 * V))
  w_home <- mean(diag(W))
  w_away <- (sum(W) - sum(diag(W))) / (length(z) * (length(z) - 1))
  out[[nm("w_home")]]   <- w_home
  out[[nm("w_away")]]   <- w_away
  out[[nm("delta_SA")]] <- w_home - w_away

  ## within-patch variance correction (patch means alone overstate mean fitness)
  if (!is.null(Va_p)) {
    va <- Va_p[ok]; va[!is.finite(va)] <- 0
    Vt  <- V + va
    Wv  <- sqrt(V / Vt) * exp(-outer(z, th, function(a, b) (a - b)^2) / (2 * Vt))
    out[[nm("w_home_vc")]]   <- mean(diag(Wv))
    out[[nm("w_away_vc")]]   <- (sum(Wv) - sum(diag(Wv))) / (length(z) * (length(z) - 1))
    out[[nm("delta_SA_vc")]] <- out[[nm("w_home_vc")]] - out[[nm("w_away_vc")]]
  } else {
    out[[nm("w_home_vc")]] <- out[[nm("w_away_vc")]] <- out[[nm("delta_SA_vc")]] <- NA_real_
  }

  out[[nm("lag_mse")]]  <- mean((z - th)^2)
  out[[nm("load_sel")]] <- mean((z - th)^2) / (2 * V)        ## mean -log w
  out
}

## --- per-locus diversity, computed once and reused ---------------------
locus_stats <- function(GTs_all, map_all, pop) {
  pop_f <- factor(pop)
  cnt <- rowsum(GTs_all, group = pop_f, reorder = TRUE)                  ## allele counts per patch
  ## individuals per patch: a vector unless there are missing genotypes, in which
  ## case it has to be counted per locus (a full n_ind x n_loci matrix)
  nobs <- if (anyNA(GTs_all)) {
    rowsum((!is.na(GTs_all)) * 1, group = pop_f, reorder = TRUE)
  } else {
    as.vector(table(pop_f))                                              ## recycles down columns
  }
  p <- cnt / (2 * nobs)                                                  ## patch x locus freq

  n_all      <- 2 * nobs
  het_within <- 2 * p * (1 - p) * n_all / pmax(n_all - 1, 1)             ## unbiased within-patch
  pbar <- colMeans(p, na.rm = TRUE)

  data.table(marker = map_all$marker, Chr = map_all$Chr, Pos = map_all$Pos,
             type = map_all$type,
             Ho   = colMeans(GTs_all == 1, na.rm = TRUE),
             Hs   = colMeans(het_within, na.rm = TRUE),
             Ht   = 2 * pbar * (1 - pbar),
             pbar = pbar)
}

## --- diversity / differentiation per locus class ----------------------
## Nei's Hs/Ht plus Weir & Cockerham's Fst via SNPRelate. Uses ALL sampled
## individuals, not the 1-per-patch analysis subset.
diversity_stats <- function(ls, map_all, pop, gds_all = NULL) {
  pop_f <- factor(pop)

  by_class <- function(idx, label) {
    if (!length(idx)) return(NULL)
    d  <- ls[idx]
    hs <- mean(d$Hs, na.rm = TRUE); ht <- mean(d$Ht, na.rm = TRUE)
    fst_wc <- NA_real_
    if (!is.null(gds_all) && length(idx) > 1) {
      fst_wc <- tryCatch(
        SNPRelate::snpgdsFst(gds_all, population = pop_f, method = "W&C84",
                             snp.id = d$marker, autosome.only = FALSE,
                             verbose = FALSE)$Fst,
        error = function(e) NA_real_)
    }
    data.table(class = label, n_snp = length(idx),
               Ho = mean(d$Ho, na.rm = TRUE), Hs = hs, Ht = ht,
               Fis = 1 - mean(d$Ho, na.rm = TRUE) / hs,
               Fst_nei = 1 - hs / ht, Fst_wc = fst_wc,
               mean_maf = mean(pmin(d$pbar, 1 - d$pbar), na.rm = TRUE),
               prop_poly = mean(d$pbar > 0 & d$pbar < 1, na.rm = TRUE))
  }

  rbindlist(list(
    by_class(which(ls$type == "ntrl"),                        "ntrl"),
    by_class(which(ls$type == "ntrl" & ls$Chr == "Chr1"),      "ntrl_chr1_QTN"),
    by_class(which(ls$type == "ntrl" & ls$Chr == "Chr2"),      "ntrl_chr2"),
    by_class(which(ls$type == "delet"),                        "delet"),
    by_class(which(ls$type == "QTN"),                          "QTN")
  ), fill = TRUE)
}

## --- Tajima's D from a window's segregating sites ---------------------
tajimas_D <- function(pbar, n_alleles) {
  S <- length(pbar)
  if (S < 3 || n_alleles < 4) return(NA_real_)
  n   <- n_alleles
  pi_ <- sum(2 * pbar * (1 - pbar) * n / (n - 1))
  a1  <- sum(1 / seq_len(n - 1)); a2 <- sum(1 / seq_len(n - 1)^2)
  b1  <- (n + 1) / (3 * (n - 1)); b2 <- 2 * (n^2 + n + 3) / (9 * n * (n - 1))
  c1  <- b1 - 1 / a1;             c2 <- b2 - (n + 2) / (a1 * n) + a2 / a1^2
  e1  <- c1 / a1;                 e2 <- c2 / (a1^2 + a2)
  v <- e1 * S + e2 * S * (S - 1)
  if (!is.finite(v) || v <= 0) return(NA_real_)
  (pi_ - S / a1) / sqrt(v)
}

## --- background-selection windows -------------------------------------
## BGS is a LOCAL effect: diversity is reduced where selected sites are dense
## and recombination is low. So the observable is the diversity landscape
## regressed on (a) local recombination and (b) proximity to deleterious sites,
## both known exactly from the recombination map -- no selection coefficients
## or theory parameters needed.
##
## delet_load_cM is a recombination-weighted count of deleterious SITES,
## sum_j exp(-|cM_i - cM_j| / scale): the SHAPE of a classic B predictor
## (Hudson & Kaplan 1995 / Nordborg et al. 1996) without its s/h/u constants,
## which the .ini does not pin down unambiguously. It is a relative predictor,
## fine for correlations and for ranking windows, not an absolute B.
##
## The window grid is fixed (WIN_BP, aligned to the 500 kb map resolution) and
## identical across runs, so the MATCHED contrast -- B_obs = pi_bgs / pi_nobgs
## per window, for the same chr/V/c/env -- can be computed later by collating
## these tables. That paired estimate is the gold standard here, since the
## nobgs grid mirrors bgs2 exactly; everything in this function is the
## within-run estimate that does not need the paired file.
bgs_windows <- function(ls, map_full, n_ind, win_bp = 5e5, cM_scale = 1) {

  ntrl <- ls[type == "ntrl"]
  if (!nrow(ntrl)) return(NULL)

  ## interpolate cM for any position from the full recombination map
  cm_of <- function(ch, pos) {
    mf <- map_full[Chr == ch]
    stats::approx(mf$bp, mf$cM, xout = pos, rule = 2, ties = "ordered")$y
  }

  rbindlist(lapply(unique(ntrl$Chr), function(ch) {
    d  <- ntrl[Chr == ch]
    mf <- map_full[Chr == ch]
    brk <- seq(0, max(mf$bp) + win_bp, by = win_bp)
    d[, win := cut(Pos, brk, labels = FALSE, include.lowest = TRUE)]

    w <- d[, .(n_snp = .N,
               pi    = mean(Ht, na.rm = TRUE),      ## metapopulation heterozygosity
               Hs    = mean(Hs, na.rm = TRUE),      ## within-patch
               Ho    = mean(Ho, na.rm = TRUE),
               Fst   = 1 - mean(Hs, na.rm = TRUE) / mean(Ht, na.rm = TRUE),
               mean_maf  = mean(pmin(pbar, 1 - pbar), na.rm = TRUE),
               tajima_D  = tajimas_D(pbar, 2 * n_ind)),
           by = win]

    w[, `:=`(Chr = ch, start = brk[win], mid = brk[win] + win_bp / 2)]
    w[, snp_density := n_snp / (win_bp / 1e6)]

    ## local recombination: cM per Mb across the window, from the map itself
    w[, cM_start := cm_of(ch, start)]
    w[, cM_end   := cm_of(ch, start + win_bp)]
    w[, cM_per_Mb := (cM_end - cM_start) / (win_bp / 1e6)]

    ## deleterious sites: count in window, and recombination-weighted proximity
    del <- mf[type == "delet"]
    if (nrow(del)) {
      del_cM <- del$cM
      w[, n_delet := vapply(seq_len(.N), function(i)
          sum(del$bp >= start[i] & del$bp < start[i] + win_bp), integer(1))]
      mid_cM <- cm_of(ch, w$mid)
      w[, delet_load_cM := vapply(seq_along(mid_cM), function(i)
          sum(exp(-abs(mid_cM[i] - del_cM) / cM_scale)), numeric(1))]
    } else {
      w[, `:=`(n_delet = 0L, delet_load_cM = 0)]
    }

    ## QTN sites, for the same treatment of the selected-chromosome effect
    qtn <- mf[type == "QTN"]
    w[, n_qtn := if (nrow(qtn)) vapply(seq_len(.N), function(i)
        sum(qtn$bp >= start[i] & qtn$bp < start[i] + win_bp), integer(1)) else 0L]

    w[, .(Chr, win, start, mid, n_snp, snp_density, pi, Hs, Ho, Fst, mean_maf,
          tajima_D, cM_per_Mb, n_delet, delet_load_cM, n_qtn)]
  }), fill = TRUE)
}

## --- headline BGS statistics from the window table --------------------
bgs_summary <- function(w) {
  out <- list()
  if (is.null(w) || nrow(w) < 10) {
    for (k in c("bgs_n_win", "bgs_cor_pi_rec", "bgs_slope_pi_rec", "bgs_cor_pi_delet",
                "bgs_B_rel_rec", "bgs_B_rel_delet", "bgs_tajD_mean", "bgs_tajD_lowrec_diff",
                "bgs_cor_pi_rec_chr2", "bgs_cor_pi_delet_chr2"))
      out[[k]] <- NA_real_
    return(out)
  }

  ok <- w[is.finite(pi) & is.finite(cM_per_Mb) & cM_per_Mb > 0]
  out$bgs_n_win <- nrow(w)

  ## diversity vs local recombination: BGS predicts a POSITIVE correlation
  out$bgs_cor_pi_rec   <- if (nrow(ok) > 5) stats::cor(ok$pi, log(ok$cM_per_Mb)) else NA_real_
  out$bgs_slope_pi_rec <- if (nrow(ok) > 5)
    unname(stats::coef(stats::lm(pi ~ log(cM_per_Mb), data = ok))[2]) else NA_real_

  ## diversity vs recombination-weighted deleterious-site load: NEGATIVE
  okd <- w[is.finite(pi) & is.finite(delet_load_cM)]
  out$bgs_cor_pi_delet <- if (nrow(okd) > 5 && stats::sd(okd$delet_load_cM) > 0)
    stats::cor(okd$pi, okd$delet_load_cM) else NA_real_

  ## relative B: diversity in the lowest vs highest quintile of each predictor
  rel <- function(d, col) {
    if (nrow(d) < 10) return(NA_real_)
    q <- stats::quantile(d[[col]], c(0.2, 0.8), na.rm = TRUE)
    lo <- mean(d$pi[d[[col]] <= q[1]], na.rm = TRUE)
    hi <- mean(d$pi[d[[col]] >= q[2]], na.rm = TRUE)
    lo / hi
  }
  out$bgs_B_rel_rec   <- rel(ok,  "cM_per_Mb")       ## < 1 if low recombination loses diversity
  out$bgs_B_rel_delet <- 1 / rel(okd, "delet_load_cM")

  out$bgs_tajD_mean <- mean(w$tajima_D, na.rm = TRUE)
  if (nrow(ok) >= 10) {
    q <- stats::quantile(ok$cM_per_Mb, c(0.2, 0.8), na.rm = TRUE)
    out$bgs_tajD_lowrec_diff <- mean(ok$tajima_D[ok$cM_per_Mb <= q[1]], na.rm = TRUE) -
                                mean(ok$tajima_D[ok$cM_per_Mb >= q[2]], na.rm = TRUE)
  } else out$bgs_tajD_lowrec_diff <- NA_real_

  ## same, restricted to the neutral chromosome -- free of the QTN chromosome's
  ## own selection, so a cleaner read on background selection alone
  w2 <- w[Chr == "Chr2" & is.finite(pi)]
  out$bgs_cor_pi_rec_chr2 <- if (nrow(w2[cM_per_Mb > 0]) > 5)
    stats::cor(w2[cM_per_Mb > 0]$pi, log(w2[cM_per_Mb > 0]$cM_per_Mb)) else NA_real_
  out$bgs_cor_pi_delet_chr2 <- if (nrow(w2) > 5 && stats::sd(w2$delet_load_cM) > 0)
    stats::cor(w2$pi, w2$delet_load_cM) else NA_real_
  out
}

## --- isolation by distance, from the saved GRM ------------------------
## Cheap structure statistic that needs nothing new: regress off-diagonal GRM
## entries on the lattice distance between the two individuals' patches.
ibd_from_grm <- function(GRM, env_ind) {
  n <- nrow(GRM)
  if (is.null(GRM) || n < 10) return(list(ibd_grm_slope = NA_real_, ibd_grm_cor = NA_real_))
  d  <- as.matrix(stats::dist(as.matrix(env_ind[, .(x, y)])))
  ut <- upper.tri(d)
  dd <- d[ut]; gg <- GRM[ut]
  ok <- is.finite(dd) & is.finite(gg) & dd > 0
  if (sum(ok) < 10) return(list(ibd_grm_slope = NA_real_, ibd_grm_cor = NA_real_))
  fit <- stats::lm(gg[ok] ~ log(dd[ok]))
  list(ibd_grm_slope = unname(stats::coef(fit)[2]),
       ibd_grm_cor   = stats::cor(gg[ok], log(dd[ok])))
}

## --- assemble everything for one simulation ---------------------------
sim_popgen <- function(files, base_name, params, GTs_all, map_all, pop_vec, env_full, map_full) {

  ini_file   <- files[grepl("\\.ini$", files) & grepl(base_name, files, fixed = TRUE)][1]
  stats_file <- files[grepl("stats/", files, fixed = TRUE) & grepl("\\.txt$", files)][1]

  ini <- if (!is.na(ini_file)) read_ini(ini_file) else list(selection_variance = NA_real_)
  V   <- ini$selection_variance
  ns  <- if (!is.na(stats_file)) read_nemo_stats(stats_file) else NULL

  ## ---- per-patch table: local optimum, trait mean, within-patch Va ----
  ## theta_p IS the env value (selection_local_optima points at the same file)
  patch <- data.table(pop = env_full$pop, x = env_full$x, y = env_full$y, theta = env_full$env)
  patch[, z_bar := if (!is.null(ns) && length(ns$z_bar) == .N) ns$z_bar else NA_real_]
  patch[, Va_p  := if (!is.null(ns) && !is.null(ns$Va_p) && length(ns$Va_p) == .N) ns$Va_p else NA_real_]
  patch[, sampled := pop %in% unique(pop_vec)]

  la_all <- local_adaptation(patch$z_bar, patch$theta, V, patch$Va_p, prefix = "")
  la_smp <- local_adaptation(patch[sampled == TRUE, z_bar], patch[sampled == TRUE, theta],
                             V, patch[sampled == TRUE, Va_p], prefix = "smp_")

  ## ---- diversity / differentiation from the raw genotype sample ------
  gds_path <- tempfile(fileext = ".gds")
  gds_all  <- create_gds_from_geno(geno = GTs_all, map = map_all, gds_path)
  on.exit({ SNPRelate::snpgdsClose(gds_all); unlink(gds_path) }, add = TRUE)
  ls_  <- locus_stats(GTs_all, map_all, pop_vec)          ## per-locus, reused below
  div  <- diversity_stats(ls_, map_all, pop_vec, gds_all)

  ## ---- background selection: the diversity landscape vs recombination and
  ## deleterious-site proximity (see bgs_windows() for what is and isn't assumed)
  win <- bgs_windows(ls_, map_full, n_ind = nrow(GTs_all), win_bp = WIN_BP)
  bgs <- bgs_summary(win)

  ## ---- additive variance in the sample (segregating QTNs only) -------
  qtn <- which(map_all$type == "QTN")
  va_sample <- NA_real_; n_qtn_seg <- 0L
  if (length(qtn)) {
    p_q <- colMeans(GTs_all[, qtn, drop = FALSE]) / 2
    a_q <- map_all$allelic_values[qtn]
    seg <- p_q > 0 & p_q < 1
    n_qtn_seg <- sum(seg)
    va_sample <- sum(2 * p_q[seg] * (1 - p_q[seg]) * a_q[seg]^2, na.rm = TRUE)
  }

  ## ---- one-row summary ------------------------------------------------
  summ <- data.table(
    file = base_name, bgs = params$bgs, chr_set = params$Chr,
    V_sel = V, c_par = params$c, env_id = params$env,
    selection_model = ini$selection_model,
    n_ind_sample = nrow(GTs_all), n_patch_sample = uniqueN(pop_vec),
    n_snp_sample = ncol(GTs_all),
    n_qtn_sample = length(qtn), n_qtn_segregating = n_qtn_seg,
    n_delet_sample = sum(map_all$type == "delet"),
    quanti_loci_total = ini$quanti_loci, patch_capacity = ini$patch_capacity,
    patch_number = ini$patch_number, va_sample = va_sample
  )
  summ <- cbind(summ, as.data.table(la_all), as.data.table(la_smp), as.data.table(bgs))

  ## headline diversity columns (the full per-class table is saved alongside)
  pull <- function(cl, col) { v <- div[class == cl][[col]]; if (length(v)) v[1] else NA_real_ }
  summ[, `:=`(
    fst_ntrl_wc   = pull("ntrl", "Fst_wc"),   fst_ntrl_nei = pull("ntrl", "Fst_nei"),
    hs_ntrl       = pull("ntrl", "Hs"),       ht_ntrl      = pull("ntrl", "Ht"),
    ho_ntrl       = pull("ntrl", "Ho"),       fis_ntrl     = pull("ntrl", "Fis"),
    fst_ntrl_chr1 = pull("ntrl_chr1_QTN", "Fst_wc"),
    fst_ntrl_chr2 = pull("ntrl_chr2", "Fst_wc"),
    fst_delet     = pull("delet", "Fst_wc"),  fst_qtn      = pull("QTN", "Fst_wc")
  )]

  ## Nemo's own whole-metapopulation stats (all patches, its own estimators)
  if (!is.null(ns)) {
    keep <- intersect(c("generation", "adlt.nbr", "extrate", "adlt.ho", "adlt.hsnei",
                        "adlt.htnei", "adlt.fis", "adlt.fst", "adlt.fit", "adlt.q1",
                        "adlt.q1.Va", "adlt.q1.Vb", "adlt.q1.Qst", "adlt.delfreq",
                        "adlt.delfst", "adlt.delsegr", "load", "heterosis", "fitness.mean"),
                      names(ns$scalars))
    sc <- ns$scalars[, ..keep]
    setnames(sc, paste0("nemo_", gsub(".", "_", names(sc), fixed = TRUE)))
    summ <- cbind(summ, sc)
    if (all(c("nemo_adlt_q1_Qst", "nemo_adlt_fst") %in% names(summ)))
      summ[, qst_minus_fst := nemo_adlt_q1_Qst - nemo_adlt_fst]
  }

  list(summary = summ, diversity = div, patch = patch, windows = win, ini = ini)
}

## ---------------------------------------------------------------------
## Stage A: raw tarball -> GTs / map (with truth) / env_ind
## ---------------------------------------------------------------------
parse_raw <- function(file_gz, tmp_dir) {

  base_name <- basename(sub("\\.tgz$|\\.tar\\.gz$", "", file_gz))
  params    <- extract_params(base_name)
  message("Working on ", base_name)

  unpack_sim(file_gz, tmp_dir)
  files <- list.files(tmp_dir, recursive = TRUE, full.names = TRUE)

  if (!any(grepl("snp_geno", files))) {
    message("Simulation did not work\n"); return(NULL)
  }

  ## Nemo outputs for THIS sim (the *_1000_1.{map,snp_geno} pair)
  map_nemo <- fread(files[grepl(".map", files, fixed = TRUE) &
                            grepl(paste0(base_name, "_"), files, fixed = TRUE)])
  GTs      <- fread(files[grepl("snp_geno", files, fixed = TRUE) &
                            grepl(paste0(base_name, "_"), files, fixed = TRUE)])

  nemo_map <- data.table(marker = map_nemo$trait.locus,
                         do.call(rbind, strsplit(map_nemo$trait.locus, ".", fixed = TRUE)))
  setnames(nemo_map, c("V1", "V2"), c("type", "idx"))

  sample_info <- GTs[, .(pop, ID)]
  GTs <- as.matrix(GTs[1:.N, 6:ncol(GTs), with = FALSE])

  ## the ORIGINAL map (bp, cM, rec_rate, allelic_values) for this chromosome set
  map <- readRDS(file.path(MAP_DIR, paste0("rec_map", gsub("chr", "", params$Chr), ".rds")))
  map[, indx := .I]
  ## the FULL map (all loci, incl. deleterious sites and the cM scale) -- the BGS
  ## window statistics need the positions of selected sites whether or not they
  ## are segregating in the sample, so keep it before the subset below
  map_full <- data.table::copy(map)[, Chr := paste0("Chr", Chr)]

  ## map Nemo's per-type running indices back onto the original map rows
  nemo_map[, idx := as.numeric(idx) + 1]   ## Nemo counts from 0
  ntrl_idx   <- nemo_map[type == "ntrl",  idx]
  quanti_idx <- nemo_map[type == "quant", idx]
  delet_idx  <- nemo_map[type == "delet", idx]

  indx_ntrl   <- map[type == "ntrl"][ntrl_idx,   indx]
  indx_quanti <- map[type == "QTN"][quanti_idx,  indx]
  indx_delet  <- map[type == "delet"][delet_idx, indx]

  map[indx_ntrl,   nemo_marker := nemo_map[type == "ntrl",  marker]]
  map[indx_quanti, nemo_marker := nemo_map[type == "quant", marker]]
  map[indx_delet,  nemo_marker := nemo_map[type == "delet", marker]]

  map <- map[nemo_marker %in% colnames(GTs)]
  GTs <- GTs[, map$nemo_marker]
  map[, Pos := bp]
  map[, marker := paste(paste0("Chr", Chr), Pos, sep = ":")]
  colnames(GTs) <- map$marker
  setorder(map, Chr, bp)
  map[, chr_type := ifelse(Chr == 1, "QTN", "ntrl")]   ## Chr1 carries the QTNs, Chr2 is the neutral control
  map[, Chr := paste0("Chr", Chr)]

  ## De-duplicate BEFORE splitting off the analysis subset: duplicate markers are a
  ## property of the map, not of which individuals are kept, so this is the same
  ## removal as before, done once for both copies.
  remove <- map[, duplicated(marker), by = Chr][, which(V1)]
  if (length(remove) > 0) map <- map[-remove, ]

  ## GTs_all: every sampled individual (2 per patch), pre-MAF -- used ONLY by the
  ## population-genetic summaries. GTs: the 1-per-patch analysis subset, as before.
  GTs_all <- GTs[, map$marker]                         ## name-based -> realigns after setorder
  map_all <- data.table::copy(map)
  GTs     <- GTs_all[KEEP_INDS, , drop = FALSE]

  ## ---- environmental values on the 48 x 48 lattice ----
  x        <- gsub("env", "", params$env)
  env_raw  <- scan(file.path(ENV_DIR, paste0("env_", x, ".txt")), what = character(), quiet = TRUE)
  env_raw  <- strsplit(env_raw, c("}{"), fixed = TRUE)[[1]]
  env_vals <- as.numeric(gsub("}}", "", gsub("{{", "", env_raw, fixed = TRUE), fixed = TRUE))

  env <- data.table(expand.grid(x = 1:SIDE, y = SIDE:1))
  env[, pop := 1:(SIDE * SIDE)]
  env[, env := env_vals]
  ## NB: compute the match OUTSIDE the [ ] -- inside env[...], `env$pop` would
  ## resolve to the atomic COLUMN `env`, not the table.
  new_order <- match(sample_info$pop, env$pop)
  env_ind   <- env[new_order][KEEP_INDS]
  env_ind[, indx := .I]

  ## ---- MAF filter (on the analysis subset, as before) ----
  maf <- colSums(GTs) / nrow(GTs) / 2
  map[, MAF := ifelse(maf < 0.5, maf, 1 - maf)]
  keep_snps <- map$MAF > MIN_MAF
  GTs <- GTs[, keep_snps]
  map <- map[MAF > MIN_MAF, ]
  stopifnot(identical(colnames(GTs), map$marker))

  message("## ------------------------------------------------")
  message("Data contains ", nrow(map), " SNPs")
  message("## ------------------------------------------------")

  ## ---- Va per QTN (always create the columns, even with zero QTNs) ----
  qtn_rows <- map[, which(type == "QTN")]
  map[, `:=`(Va = NA_real_, sum_Va = NA_real_, p_Va = NA_real_)]
  if (length(qtn_rows) == 0L) {
    message("No QTNs present in this simulation -- skipping Va/p_Va computation.")
  } else {
    map[qtn_rows, Va     := get_va(map, GTs, qtn_rows)]
    map[qtn_rows, sum_Va := sum(Va)]
    map[qtn_rows, p_Va   := Va / sum_Va]
  }

  ## ---- truth annotation ----
  map[type == "QTN", true_QTN := TRUE]
  fq <- focal_QTN(map, GTs, qtn_col = "true_QTN")
  map[fq, on = "marker", `:=`(focal_QTN       = i.focal_QTN,
                              bp_to_focal_QTN = i.bp_to_focal_QTN,
                              max_LD_with_QTN = i.max_LD_with_QTN)]
  map[type == "QTN", `:=`(focal_QTN = marker, max_LD_with_QTN = 1, bp_to_focal_QTN = 0)]
  ## the neutral chromosome has no true signal by construction
  map[chr_type == "ntrl", bp_to_focal_QTN := max(Pos)]
  map[chr_type == "ntrl", max_LD_with_QTN := 0]

  ## ---- population-genetic + local-adaptation summaries (own output files) ----
  pg <- tryCatch(sim_popgen(files, base_name, params, GTs_all, map_all, sample_info$pop, env, map_full),
                 error = function(e) { message("  !! popgen summary failed: ", conditionMessage(e)); NULL })

  list(GTs = GTs, map = map, env = env_ind, params = params, popgen = pg)
}

## ---------------------------------------------------------------------
## Stage B: parsed genotypes -> regen bundle (LD_decay / ld_w / GRM / EMMAX / LFMM)
## ---------------------------------------------------------------------
regen_from_parsed <- function(d, out_path) {

  GTs <- d$GTs; map <- as.data.table(d$map); env_ind <- d$env
  n <- nrow(GTs)

  gds_path <- tempfile(fileext = ".gds"); on.exit(unlink(gds_path), add = TRUE)
  gds <- create_gds_from_geno(geno = GTs, map = map, gds_path)

  ## LD-decay (corr) + per-SNP ld_w computed in place from the same edge lists
  LD_decay <- do.call(compute_LD_decay, c(list(gds = gds, el_data_folder = NULL), DECAY_ARGS))
  ld_ws    <- LD_decay$ld_ws[map$marker, , drop = FALSE]
  ld95_col <- if ("rho_0.95" %in% colnames(ld_ws)) "rho_0.95" else "0.95"
  map[, ld_w_095 := ld_ws[, ld95_col]]

  ## QTN-chromosome distance/LD annotations, relative to the QTN chromosome's own
  ## fitted decay curve (kept for map compatibility; not used by run_sim_LDscnR.R)
  ds1 <- LD_decay$by_chr[["Chr1"]]$decay_sum
  if (!is.null(ds1)) {
    a <- ds1$a; b <- ds1$b; cc_ <- ds1$c
    map[chr_type == "QTN", rho_d  := a * bp_to_focal_QTN / (a * bp_to_focal_QTN + 1)]
    map[chr_type == "QTN", ld_rel := (max_LD_with_QTN - b) / (cc_ - b)]
    map[chr_type == "QTN", ld_rel := pmin(pmax(ld_rel, 0), 1)]
  }

  ## GRM basis (shared kinship for observed + structured-null EMMAX)
  if (GRM_METHOD == "ld_w_threshold") {
    grm_markers <- map$marker[which(map$ld_w_095 < GRM_LDW_THRESHOLD)]
    stage1 <- NULL
  } else {
    ## gds = : rebuilds each chromosome's edge list on the fly, so keep_el stays FALSE
    stage1 <- ld_complexity_reduction(map = map, LD_decay = LD_decay, rho = CR_RHO, gds = gds)
    grm_markers <- do.call(ld_prune_and_eMLG, c(list(GTs = GTs, stage1 = stage1), PRUNE_ARGS))$pruned
  }
  message(sprintf("  GRM (%s): %d / %d markers", GRM_METHOD, length(grm_markers), nrow(map)))
  GRM <- snpgdsGRM(gds, snp.id = grm_markers, method = "GCTA",
                   verbose = FALSE, autosome.only = FALSE)$grm

  ## EMMAX on the GRM + genomic control if gif > 1.1
  emx <- emmax(env_ind$env, GTs, K = GRM)
  map[, emx_p := emx$pval][, emx_F := emx$F]
  gif <- stats::median(map$emx_F) / stats::qf(0.5, 1, n - 2, lower.tail = FALSE)
  if (gif > 1.1) { map[, emx_F := emx_F / gif]
                   map[, emx_p := stats::pf(emx_F, 1, n - 2, lower.tail = FALSE)] }

  ## LFMM (K=5) with genomic control -- unchanged engine
  tmp <- tempfile(); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  write.lfmm(GTs, file.path(tmp, "geno.lfmm")); write.env(env_ind$env, file.path(tmp, "grad.env"))
  proj <- lfmm2(file.path(tmp, "geno.lfmm"), file.path(tmp, "grad.env"), K = LFMM_K)
  pv   <- suppressWarnings(lfmm2.test(proj, file.path(tmp, "geno.lfmm"), file.path(tmp, "grad.env"),
                                      genomic.control = TRUE, full = TRUE))
  map[, lfmm_p := pv$pvalues][, lfmm_F := pv$fscores / pv$gif]

  ## same column order as the regen_sim_data.R bundles (names are what downstream
  ## code uses, but keep it identical so positional habits don't bite)
  canon <- c("Chr", "type", "bp", "cM", "rec_rate", "pos_nemo", "allelic_values", "indx",
             "nemo_marker", "Pos", "marker", "chr_type", "MAF", "Va", "sum_Va", "p_Va",
             "true_QTN", "focal_QTN", "bp_to_focal_QTN", "max_LD_with_QTN",
             "lfmm_p", "lfmm_F", "rho_d", "ld_rel", "ld_w_095", "emx_p", "emx_F")
  setcolorder(map, intersect(canon, names(map)))

  saveRDS(list(GTs = GTs, map = map, env = env_ind, LD_decay = LD_decay, ld_ws = ld_ws,
               GRM = GRM, grm_markers = grm_markers, grm_method = GRM_METHOD,
               complexity_reduction = if (is.null(stage1)) NULL else list(stage1 = stage1),
               emx_gif = gif),
          out_path)
  message("  wrote ", basename(out_path))
  invisible(list(out_path = out_path, GRM = GRM))
}

## ---------------------------------------------------------------------
## One raw file, end to end
## ---------------------------------------------------------------------
process_file <- function(file_gz, out_path) {
  tmp_dir <- file.path(TMP_ROOT, paste0("unpack_", basename(sub("\\.tgz$|\\.tar\\.gz$", "", file_gz))))
  on.exit({ unlink(tmp_dir, recursive = TRUE); gc() }, add = TRUE)

  d <- parse_raw(file_gz, tmp_dir)
  if (is.null(d)) return(invisible(NULL))

  if (nzchar(PARSED_DIR)) {
    saveRDS(list(GTs = d$GTs, map = d$map, env = d$env),
            file.path(PARSED_DIR, paste0(basename(sub("\\.tgz$|\\.tar\\.gz$", "", file_gz)), ".rds")))
  }

  res <- regen_from_parsed(d, out_path)

  ## popgen summaries go to their OWN files -- diagnostics about the simulation,
  ## not inputs to the benchmark bundle
  if (!is.null(d$popgen)) {
    pg <- d$popgen
    ibd <- ibd_from_grm(res$GRM, d$env)
    pg$summary <- cbind(pg$summary, as.data.table(ibd))
    stem <- basename(sub("\\.tgz$|\\.tar\\.gz$", "", file_gz))
    saveRDS(pg, file.path(POPGEN_DIR, paste0(stem, ".rds")))
    fwrite(pg$summary, file.path(POPGEN_DIR, paste0(stem, "_summary.csv")))
    message("  wrote popgen: ", stem, ".rds  (Fst=", signif(pg$summary$fst_ntrl_wc, 3),
            ", delta_SA=", signif(pg$summary$delta_SA, 3),
            ", la_cor2=", signif(pg$summary$la_cor2, 3), ")")
  }

  invisible(res$out_path)
}

## ---- driver (same CLI shape as regen_sim_data.R) ---------------------
a   <- commandArgs(trailingOnly = TRUE)
V   <- if (length(a) >= 1) a[1] else "2"
cc  <- if (length(a) >= 2) a[2] else "1"
env <- if (length(a) >= 3) a[3] else "1"
chr <- if (length(a) >= 4) a[4] else "1"
tag <- if (length(a) >= 5) a[5] else "bgs"
envs <- if (env == "all") as.character(1:10) else env
chrs <- if (chr == "all") as.character(1:10) else chr

for (e in envs) for (ch in chrs) {
  stem <- sprintf("adapt_%s_chr%s_V%s_c%s_env%s", tag, ch, V, cc, e)
  out  <- file.path(OUT_DIR, paste0(stem, ".rds"))
  raw  <- file.path(RAW_DIR, paste0(stem, ".tgz"))
  if (file.exists(out))  { message("skip (exists): ", basename(out)); next }
  if (!file.exists(raw)) { message("!! missing raw: ", basename(raw));  next }
  message(sprintf("V%s_c%s_env%s_chr%s (%s)", V, cc, e, ch, tag))
  tryCatch(process_file(raw, out),
           error = function(err) message("  !! FAILED ", basename(out), ": ", conditionMessage(err)))
}
