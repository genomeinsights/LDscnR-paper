## module_sim/R/regen_stats.R
## Regenerate the per-file summary statistics for the simulated data FROM SCRATCH,
## coherently and self-contained, now that we know exactly what the OR pipeline
## needs. Recomputes -- from the parsed genotype matrix (the exact final GTs) --
##   LD-decay -> ld_w -> LD-complexity-reduced GRM -> EMMAX (+ GC) -> LFMM (K=5, GC)
## and, crucially, SAVES THE GRM (and its pruned marker set) so the structured
## null's surrogate EMMAX runs on the identical kinship as the observed emx_p.
##
## Uses LDscnR functions only (compute_LD_decay / compute_ld_w /
## ld_complexity_reduction / ld_prune_and_eMLG / emmax) + LEA for LFMM. One file =
## one chromosome (each is an independent sim with its own GRM), matching how the
## OR pipeline pools 10 files and runs EMMAX per-file.
##
## Run (per your convention: reviewed before running):
##   Rscript module_sim/R/regen_stats.R  V c env chr   [tag]     # one file
##   Rscript module_sim/R/regen_stats.R  V c all all   [tag]     # a whole (V,c) cell
## defaults: tag = nobgs.  Reads parsed files, writes augmented ones to OUT_DIR.

suppressMessages({ library(data.table); library(SNPRelate); library(LEA)
                   devtools::load_all("/Users/petrikem/gitlab/LDscnR", quiet = TRUE) })

IN_DIR  <- Sys.getenv("SIM_IN",  "/Volumes/Nemo/Nemo_sim/parsed_sim_data2")
OUT_DIR <- Sys.getenv("SIM_OUT", "/Volumes/Nemo/Nemo_sim/regen_sim_data")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

## ---- canonical settings (from the original Parse_sim_data.R) ----------------
RHO_GRID   <- c(seq(0.05, 0.95, by = 0.05), 0.99)      # ld_w columns
DECAY_ARGS <- list(min_maf_decay = 0.1, q = 0.95, n_sub_bg = 5000, n_win_decay = 5,
                   overlap = 0.5, max_SNPs_decay = Inf, prob_robust = 0.95,
                   max_pairs = 5000, ld_method = "corr", n_strata = 20, keep_el = TRUE,
                   slide = 1000, rho_targets = c(0.99), cores = 1)
CR_RHO     <- 0.5                                       # ld_complexity_reduction rho
PRUNE_ARGS <- list(ld_w_col = "ld_w_095", ld_w_threshold = 0.025,
                   score_threshold = 0.80, min_r2 = 0.2, distance_threshold = 5e5)
LFMM_K     <- 5L

regen_file <- function(V, cc, env, chr, tag = "nobgs") {
  in_path  <- sprintf("%s/adapt_%s_chr%s_V%s_c%s_env%s.rds", IN_DIR,  tag, chr, V, cc, env)
  out_path <- sprintf("%s/adapt_%s_chr%s_V%s_c%s_env%s.rds", OUT_DIR, tag, chr, V, cc, env)
  d <- readRDS(in_path); GTs <- d$GTs; map <- data.table::as.data.table(d$map); env_ind <- d$env
  n <- nrow(GTs)

  gds_path <- tempfile(fileext = ".gds"); on.exit(unlink(gds_path), add = TRUE)
  gds <- create_gds_from_geno(geno = GTs, map = map, gds_path)

  ## LD-decay (keep edge list for the complexity reduction) + per-SNP ld_w
  ## (MAF filtering is via min_maf_decay in DECAY_ARGS; the branch compute_LD_decay
  ## reads allele frequencies from the gds, so no `maf` argument is passed)
  LD_decay <- do.call(compute_LD_decay, c(list(gds = gds, el_data_folder = NULL), DECAY_ARGS))
  ld_ws    <- compute_ld_w(rho = RHO_GRID, LD_decay)[map$marker, , drop = FALSE]  # rows named by marker
  ld95_col <- if ("rho_0.95" %in% colnames(ld_ws)) "rho_0.95" else "0.95"
  map[, ld_w_095 := ld_ws[, ld95_col]]                  # feeds ld_prune_and_eMLG

  ## LD-complexity-reduced GRM (the EMMAX kinship) -- the exact chain used
  stage1 <- ld_complexity_reduction(map = map, LD_decay = LD_decay, rho = CR_RHO)
  result <- do.call(ld_prune_and_eMLG, c(list(GTs = GTs, stage1 = stage1), PRUNE_ARGS))
  pruned_markers <- result$pruned
  message(sprintf("  GRM pruning: %d / %d markers kept", length(pruned_markers), nrow(map)))
  GRM <- snpgdsGRM(gds, snp.id = pruned_markers, method = "GCTA",
                   verbose = FALSE, autosome.only = FALSE)$grm

  ## EMMAX on the pruned GRM + genomic control if gif > 1.1
  emx <- emmax(env_ind$env, GTs, K = GRM)
  map[, emx_p := emx$pval][, emx_F := emx$F]
  gif <- stats::median(map$emx_F) / stats::qf(0.5, 1, n - 2, lower.tail = FALSE)
  if (gif > 1.1) { map[, emx_F := emx_F / gif]
                   map[, emx_p := stats::pf(emx_F, 1, n - 2, lower.tail = FALSE)] }

  ## LFMM (K=5) with genomic control
  tmp <- tempfile(); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  write.lfmm(GTs, file.path(tmp, "geno.lfmm")); write.env(env_ind$env, file.path(tmp, "grad.env"))
  proj <- lfmm2(file.path(tmp, "geno.lfmm"), file.path(tmp, "grad.env"), K = LFMM_K)
  pv   <- suppressWarnings(lfmm2.test(proj, file.path(tmp, "geno.lfmm"), file.path(tmp, "grad.env"),
                                      genomic.control = TRUE, full = TRUE))
  map[, lfmm_p := pv$pvalues][, lfmm_F := pv$fscores / pv$gif]

  saveRDS(list(GTs = GTs, map = map, env = env_ind, LD_decay = LD_decay, ld_ws = ld_ws,
               GRM = GRM, pruned_markers = pruned_markers,           # <- NEW: kinship kept
               emx_gif = gif),
          out_path)
  message("  wrote ", basename(out_path))
  invisible(out_path)
}

## ---- driver ----------------------------------------------------------------
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
  tryCatch(regen_file(V, cc, e, ch, tag),                # one bad file must not kill the batch
           error = function(err) message("  !! FAILED ", basename(out), ": ", conditionMessage(err)))
}
