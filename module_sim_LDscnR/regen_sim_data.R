## =====================================================================
## module_sim_LDscnR / regen_sim_data.R
##
## Regenerate the per-file simulated-data bundles FROM the parsed genotypes,
## self-contained and using only exported LDscnR functions -- the sim counterpart
## of module_sticklebacks_LDscnR/regen_3sp_data.R. One file = one chromosome (each
## is an independent Nemo sim with its own structure), matching how the benchmark
## pools 10 files per (V, c, env) and runs EMMAX per-file.
##
## Recomputes per file: compute_LD_decay (corr) -> compute_ld_w -> GRM -> EMMAX
## (+GC) -> LFMM (K=5, GC), and SAVES the GRM + its marker set so the structured
## null's surrogate EMMAX runs on the identical kinship as the observed emx_p
## (engine coherence). nobgs only (bgs is broken; co-author to reimplement).
##
## GRM basis (GRM_METHOD): "complexity_chain" is the validated sim default
## (ld_complexity_reduction rho=0.5 -> ld_prune_and_eMLG, keeps ~45% on these
## weak-LD sims). "ld_w_threshold" is the 3sp-style direct local-LD filter and is
## UNDER EVALUATION for the sims -- on 3sp the chain over-corrected, so we are
## checking whether an ld_w threshold is the consistent choice for both datasets;
## see module_sticklebacks_LDscnR + the sim GRM PR-AUC test.
##
## Run from the LDscnR-paper root (heavy; per-file, review before running):
##   Rscript module_sim_LDscnR/regen_sim_data.R  V c env chr   [tag]     # one file
##   Rscript module_sim_LDscnR/regen_sim_data.R  V c all all   [tag]     # a whole (V,c) cell
## defaults: tag = nobgs.
## =====================================================================

suppressMessages({ library(data.table); library(SNPRelate); library(LEA)
                   devtools::load_all("/Users/petrikem/gitlab/LDscnR", quiet = TRUE) })

IN_DIR  <- Sys.getenv("SIM_IN",  "/Volumes/Nemo/Nemo_sim/parsed_sim_data2")
OUT_DIR <- Sys.getenv("SIM_OUT", "/Volumes/Nemo/Nemo_sim/regen_sim_data")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

## ---- canonical settings ----------------------------------------------
RHO_GRID   <- c(seq(0.05, 0.95, by = 0.05), 0.99)          # ld_w columns
## ld_w computed IN PLACE (ld_w_rho) from the edge lists built for the decay fit,
## which are then dropped (keep_el = FALSE) -- no edge lists saved or recomputed.
DECAY_ARGS <- list(min_maf_decay = 0.1, q = 0.95, n_sub_bg = 5000, n_win_decay = 5,
                   overlap = 0.5, max_SNPs_decay = Inf, prob_robust = 0.95,
                   max_pairs = 5000, ld_method = "corr", n_strata = 20, keep_el = FALSE,
                   slide = 1000, rho_targets = c(0.99), cores = 1, ld_w_rho = RHO_GRID)
GRM_METHOD <- Sys.getenv("SIM_GRM", "complexity_chain")    # or "ld_w_threshold"
CR_RHO     <- 0.5                                           # ld_complexity_reduction rho
PRUNE_ARGS <- list(ld_w_col = "ld_w_095", ld_w_threshold = 0.025, score_threshold = 0.80,
                   min_r2 = 0.2, distance_threshold = 5e5, compute_unflagged_eMLG = FALSE)
GRM_LDW_THRESHOLD <- 0.02                                  # used only if GRM_METHOD == "ld_w_threshold"
LFMM_K     <- 5L

regen_file <- function(V, cc, env, chr, tag = "nobgs") {
  in_path  <- sprintf("%s/adapt_%s_chr%s_V%s_c%s_env%s.rds", IN_DIR,  tag, chr, V, cc, env)
  out_path <- sprintf("%s/adapt_%s_chr%s_V%s_c%s_env%s.rds", OUT_DIR, tag, chr, V, cc, env)
  d <- readRDS(in_path); GTs <- d$GTs; map <- data.table::as.data.table(d$map); env_ind <- d$env
  n <- nrow(GTs)

  gds_path <- tempfile(fileext = ".gds"); on.exit(unlink(gds_path), add = TRUE)
  gds <- create_gds_from_geno(geno = GTs, map = map, gds_path)

  ## LD-decay (corr; keep edge list) + per-SNP ld_w
  LD_decay <- do.call(compute_LD_decay, c(list(gds = gds, el_data_folder = NULL), DECAY_ARGS))
  ld_ws    <- LD_decay$ld_ws[map$marker, , drop = FALSE]                # computed in place
  ld95_col <- if ("rho_0.95" %in% colnames(ld_ws)) "rho_0.95" else "0.95"
  map[, ld_w_095 := ld_ws[, ld95_col]]

  ## GRM basis (shared kinship for observed + structured-null EMMAX)
  if (GRM_METHOD == "ld_w_threshold") {
    grm_markers <- map$marker[which(map$ld_w_095 < GRM_LDW_THRESHOLD)]
    stage1 <- NULL
  } else {
    stage1 <- ld_complexity_reduction(map = map, LD_decay = LD_decay, rho = CR_RHO)
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

  saveRDS(list(GTs = GTs, map = map, env = env_ind, LD_decay = LD_decay, ld_ws = ld_ws,
               GRM = GRM, grm_markers = grm_markers, grm_method = GRM_METHOD,
               complexity_reduction = if (is.null(stage1)) NULL else list(stage1 = stage1),
               emx_gif = gif),
          out_path)
  message("  wrote ", basename(out_path))
  invisible(out_path)
}

## ---- driver ----------------------------------------------------------
a   <- commandArgs(trailingOnly = TRUE)
V   <- if (length(a) >= 1) a[1] else "2"
cc  <- if (length(a) >= 2) a[2] else "1"
env <- if (length(a) >= 3) a[3] else "1"
chr <- if (length(a) >= 4) a[4] else "1"
tag <- if (length(a) >= 5) a[5] else "nobgs"
envs <- if (env == "all") as.character(1:10) else env
chrs <- if (chr == "all") as.character(1:10) else chr
for (e in envs) for (ch in chrs) {
  out <- sprintf("%s/adapt_%s_chr%s_V%s_c%s_env%s.rds", OUT_DIR, tag, ch, V, cc, e)
  if (file.exists(out)) { message("skip (exists): ", basename(out)); next }
  message(sprintf("V%s_c%s_env%s_chr%s (%s)", V, cc, e, ch, tag))
  tryCatch(regen_file(V, cc, e, ch, tag),
           error = function(err) message("  !! FAILED ", basename(out), ": ", conditionMessage(err)))
}
