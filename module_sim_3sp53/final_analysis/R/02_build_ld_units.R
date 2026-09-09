## =============================================================================
## final_analysis/R/02_build_ld_units.R
##
## Phenotype-blind Stage 1: GDS -> LD decay -> ld_complexity_reduction() ->
## GCTA kinship on the stage-1-pruned representatives. No association test
## here (that depends on the kinship this stage produces -- conflating the
## two makes "rebuild the kinship" and "rerun the association" impossible to
## decide independently, exactly module_sim_3sp53/R/02_bundle.R's own
## reasoning, reused here since it was never dataset-specific).
##
## Reads final_analysis/R/01_parse_nemo.R's corrected output (PATHS$parsed,
## the NEW "_final_v1" root -- never the old module_sim_3sp53 parsed bundles).
##
## Marker order invariant (Phase 3 gate 2's own requirement): GTs columns and
## map rows stay 1:1-aligned by `marker` from this stage's input through to
## its output -- asserted, not assumed.
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(SNPRelate); library(digest)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "R", "00_config.R"))
STAGE <- "02_build_ld_units"

build_ld_units <- function(tag, cell, rep, env, force = FALSE) {
  combo_id <- sprintf("%s_%s_rep%d_env%d", tag, cell, rep, env)
  parsed_file <- file.path(PATHS$parsed, sprintf("nemo_%s_rep%d_%s_env%d.rds", tag, rep, cell, env))
  if (!file.exists(parsed_file)) stop("R/01_parse_nemo.R has not produced: ", basename(parsed_file),
                                       " -- parse this combination first.")

  INPUTS <- unname(parsed_file)
  PARAMS <- list(decay_args = DECAY_ARGS, rho_grid = RHO_GRID, cr_rho = CR_RHO,
                 grm_basis = GRM_BASIS, grm_method = GRM_METHOD,
                 seed_bundle = SEEDS[["bundle"]], seed_clusters = SEEDS[["clusters"]])
  if (!force && !stage_stale(STAGE, INPUTS, PARAMS, target = combo_id)) {
    return(readRDS(file.path(stage_dir(STAGE, combo_id), "ld_units.rds")))
  }

  say("=== %s: %s/%s/rep%d/env%d ===\n\n", STAGE, tag, cell, rep, env)

  say("[1] reading %s\n", basename(parsed_file))
  b <- readRDS(parsed_file)
  GTs <- b$GTs; map <- b$map; env_dt <- b$env
  stopifnot("GTs columns and map rows must correspond 1:1" = identical(colnames(GTs), map$marker),
            "GTs rows and env rows must correspond 1:1" = nrow(GTs) == nrow(env_dt))
  say("    %d individuals x %s markers ; %d QTN\n", nrow(GTs), format(ncol(GTs), big.mark = ","), sum(map$true_QTN))

  say("\n[2] GDS\n")
  gds_path <- file.path(PATHS$out, "gds_tmp", paste0("sim_", combo_id, ".gds"))
  dir.create(dirname(gds_path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(gds_path)) unlink(gds_path)
  gds <- create_gds_from_geno(geno = GTs, map = map, gds_path)
  on.exit({ try(SNPRelate::snpgdsClose(gds), silent = TRUE); unlink(gds_path) }, add = TRUE)
  say("    %.1f MB\n", file.size(gds_path) / 1e6)

  say("\n[3] LD decay (n_win_decay=%d, seed %d)\n", DECAY_ARGS$n_win_decay, SEEDS[["bundle"]])
  set.seed(SEEDS[["bundle"]])
  LD_decay <- do.call(compute_LD_decay, c(list(gds = gds, ld_w_rho = RHO_GRID, seed = SEEDS[["bundle"]]), DECAY_ARGS))
  ld_ws <- LD_decay$ld_ws[map$marker, , drop = FALSE]
  ld95 <- if ("rho_0.95" %in% colnames(ld_ws)) "rho_0.95" else "0.95"
  map[, ld_w_095 := ld_ws[, ld95]]
  say("    %d chromosome(s) ; ld_w matrix %s x %d\n", nrow(LD_decay$decay_sum), format(nrow(ld_ws), big.mark = ","), ncol(ld_ws))

  say("\n[4] Stage-1 clustering (ld_complexity_reduction, rho=%.2f)\n", CR_RHO)
  set.seed(SEEDS[["clusters"]])
  stage1 <- ld_complexity_reduction(map = map, LD_decay = LD_decay, rho = CR_RHO, gds = gds)
  cl <- as.data.table(stage1$clusters)
  nl <- if ("n_loci" %in% names(cl)) cl$n_loci else lengths(cl$members)
  say("    %s Stage-1 units ; median size %.2f ; SIZE_FLOOR=%d retains %d units\n",
      format(nrow(cl), big.mark = ","), median(nl), SIZE_FLOOR, sum(nl >= SIZE_FLOOR))

  say("\n[5] kinship: basis=%s, method=%s\n", GRM_BASIS, GRM_METHOD)
  grm_markers <- switch(GRM_BASIS,
    stage1_pruned = unique(na.omit(stage1$pruned)),
    stop("unknown GRM_BASIS: ", GRM_BASIS))
  say("    %s of %s markers (%.1f%%)\n", format(length(grm_markers), big.mark = ","),
      format(nrow(map), big.mark = ","), 100 * length(grm_markers) / nrow(map))
  GRM <- SNPRelate::snpgdsGRM(gds, snp.id = grm_markers, method = GRM_METHOD, verbose = FALSE, autosome.only = FALSE)$grm
  say("    %d x %d GRM ; mean diagonal %.4f\n", nrow(GRM), ncol(GRM), mean(diag(GRM)))

  stopifnot("map/GTs correspondence broken after Stage 1" = identical(colnames(GTs), map$marker))

  OUT_DIR <- stage_dir(STAGE, combo_id)
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  out <- list(GTs = GTs, map = map, env = env_dt, LD_decay = LD_decay, stage1 = stage1,
             GRM = GRM, grm_markers = grm_markers,
             settings = list(tag = tag, cell = cell, rep = rep, env = env, params = PARAMS))
  OUT <- file.path(OUT_DIR, "ld_units.rds")
  saveRDS(out, OUT)
  write_receipt(STAGE, inputs = INPUTS, params = PARAMS, outputs = OUT, target = combo_id)
  say("\n[6] wrote %s (%.1f MB)\n", OUT, file.size(OUT) / 1e6)
  out
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 4) stop("Usage: Rscript R/02_build_ld_units.R <tag> <cell> <rep> <env>")
  build_ld_units(args[1], args[2], as.integer(args[3]), as.integer(args[4]))
}
