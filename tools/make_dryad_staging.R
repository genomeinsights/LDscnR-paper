## =============================================================================
## tools/make_dryad_staging.R
##
## Builds the Dryad data folder (copies only; nothing in this repository or in
## LDscnR-NEMO / LD-pruning-paper is moved or modified). Dryad publishes under
## CC0 while the stickleback genotypes (Fang et al. 2021) are CC-BY, so no
## individual-level genotype matrix is copied: the Stage-1 bundles are replaced
## by slim copies without them, and every staged .rds is checked afterwards.
##
## Usage (from the repository root):
##   Rscript tools/make_dryad_staging.R [OUT_DIR]
##   default OUT_DIR: ../LDscnR_manuscript/deposit/dryad_data
## Environment: LDSCNR_NEMO (LDscnR-NEMO checkout, default ../LDscnR-NEMO),
##              LDSCNR_LDPP (LD-pruning-paper checkout, default ../LD-pruning-paper)
## =============================================================================
suppressMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
ROOT <- normalizePath(Sys.getenv("LDSCNR_PAPER_ROOT", "."), mustWork = TRUE)
OUT  <- if (length(args)) args[1] else file.path(ROOT, "..", "LDscnR_manuscript", "deposit", "dryad_data")
NEMO <- Sys.getenv("LDSCNR_NEMO", file.path(ROOT, "..", "LDscnR-NEMO"))
LDPP <- Sys.getenv("LDSCNR_LDPP", file.path(ROOT, "..", "LD-pruning-paper"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
OUT <- normalizePath(OUT)

copy_tree <- function(from, to, exclude = character()) {
  fs <- list.files(from, recursive = TRUE, all.files = FALSE)
  fs <- fs[!basename(fs) %in% c(".DS_Store", ".gitkeep")]
  for (x in exclude) fs <- fs[!grepl(x, fs)]
  for (f in fs) {
    dest <- file.path(to, f)
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    stopifnot(file.copy(file.path(from, f), dest, overwrite = TRUE, copy.date = TRUE))
  }
  length(fs)
}

## ---- simulation ---------------------------------------------------------------
cat("[simulation]\n")
n <- copy_tree(file.path(ROOT, "module_sim_3sp53", "final_analysis", "results"),
               file.path(OUT, "simulation", "results"))
cat(sprintf("  results: %d files\n", n))
## Nemo inputs: the canonical parameter set plus what make_prod.sh needs to
## regenerate the 1,540 .ini files and maps (R/, ini/ templates, the macOS patch).
ni <- file.path(OUT, "simulation", "nemo_inputs")
n <- copy_tree(file.path(NEMO, "params_3spC_53cM"), file.path(ni, "params_3spC_53cM"))
for (f in c("make_prod.sh", "PATCHES.md", "nemo-macos-bigpop.patch", "README.md"))
  if (file.exists(file.path(NEMO, f))) file.copy(file.path(NEMO, f), file.path(ni, f), overwrite = TRUE)
n <- n + copy_tree(file.path(NEMO, "R"), file.path(ni, "R")) + copy_tree(file.path(NEMO, "ini"), file.path(ni, "ini"))
cat(sprintf("  nemo_inputs: %d files\n", n))

## ---- empirical: module outputs except the full bundle -------------------------
slim_bundle <- function(b) {
  ## Figure 1B-C: r^2 among the markers of the chromosome-9 window, computed
  ## here from the genotypes so the genotypes need not be shipped.
  W_CHR <- "Chr9"; W0 <- 2732880; W1 <- 2747880
  fig1 <- NULL
  if (!is.null(b$stage1$map_snp) && any(b$stage1$map_snp$Chr == W_CHR)) {
    ms <- as.data.table(b$stage1$map_snp)[Chr == W_CHR & Pos >= W0 & Pos <= W1][order(Pos)]
    if (nrow(ms)) {
      G <- as.matrix(b$GTs[, ms$marker]); storage.mode(G) <- "double"
      R2 <- suppressWarnings(cor(G, use = "pairwise.complete.obs")^2)
      fig1 <- list(chr = W_CHR, from = W0, to = W1, markers = ms[, .(marker, Chr, Pos, CL_id)], r2 = R2)
    }
  }
  ## Sample metadata reduced to what the analyses use: a neutral row index (rows
  ## are aligned with eco and GRM), ecotype, population and lineage labels, and
  ## the locality grouping the three-spined permutation null holds fixed, as an
  ## anonymous code (L1, L2, ... in order of first appearance). Individual IDs,
  ## locality names and coordinates are left in the cited Zenodo records.
  ph <- as.data.table(b$pheno)
  pheno <- data.table(index = seq_len(nrow(ph)), ecotype = ph$ecotype, population = ph$pop_ID,
                      lineage = ph$lineage,
                      locality_group = paste0("L", match(ph$pop_locality, unique(ph$pop_locality))))
  list(n_markers = ncol(b$GTs), n_individuals = nrow(b$GTs),
       map = b$map, eco = b$eco, pheno = pheno,
       LD_decay = b$LD_decay, stage1 = b$stage1,
       GRM = b$GRM, grm_markers = b$grm_markers, settings = b$settings,
       fig1_window = fig1,
       note = paste("Slim copy of out/02_bundle/bundle.rds: the genotype matrix (GTs) and the",
                    "duplicate top-level ld_ws are removed; n_markers = ncol(GTs); pheno reduced to",
                    "index, ecotype, population, lineage and an anonymous locality code.",
                    "Genotypes: Zenodo 10.5281/zenodo.18466057 (CC-BY)."))
}
for (panel in c("3sp", "9sp")) {
  cat(sprintf("[empirical_%s]\n", panel))
  mod <- file.path(ROOT, paste0("module_", panel), "out")
  dest <- file.path(OUT, paste0("empirical_", panel))
  n <- copy_tree(mod, dest, exclude = c("^02_bundle/", "^run.*\\.log$", "^rerun_v091_logs/"))
  cat(sprintf("  module outputs: %d files\n", n))
  b <- readRDS(file.path(mod, "02_bundle", "bundle.rds"))
  s <- slim_bundle(b); rm(b); invisible(gc())
  saveRDS(s, file.path(dest, "bundle_slim.rds"))
  cat(sprintf("  bundle_slim.rds: %.0f MB (n_markers = %s)\n",
              file.size(file.path(dest, "bundle_slim.rds")) / 1e6, format(s$n_markers, big.mark = ",")))
  rm(s); invisible(gc())
}

## ---- Table 6 inputs (LD-pruning-paper cache), genotype-level parts removed ----
cat("[empirical_3sp/goldilocks_cache]\n")
gc_dir <- file.path(OUT, "empirical_3sp", "goldilocks_cache")
dir.create(gc_dir, recursive = TRUE, showWarnings = FALSE)
src <- file.path(LDPP, "cache_3sp")
file.copy(file.path(src, c("map_only.rds", "K_eco.rds")), gc_dir, overwrite = TRUE)
file.copy(file.path(LDPP, "results_3sp", "3sp_cmin_loci.csv"), gc_dir, overwrite = TRUE)
## stage2_Chr*_d1e+05.rds: partitions kept, per-individual eMLG consensus matrices removed
for (f in sort(list.files(src, "^stage2_Chr[0-9]+_d1e\\+05\\.rds$", full.names = TRUE))) {
  x <- readRDS(f)
  if (!is.null(x$res$eMLG)) x$res$eMLG <- NULL
  saveRDS(x, file.path(gc_dir, basename(f)))
}
cat(sprintf("  %d files (GTs_s1units.rds and the eMLG matrices are NOT included)\n",
            length(list.files(gc_dir))))

## ---- manifest -----------------------------------------------------------------
fs <- list.files(OUT, recursive = TRUE)
man <- data.table(file = fs, bytes = file.size(file.path(OUT, fs)),
                  sha256 = vapply(file.path(OUT, fs), digest::digest, "", algo = "sha256", file = TRUE))
fwrite(man, file.path(OUT, "MANIFEST.tsv"), sep = "\t")
cat(sprintf("\n%d files, %.2f GB -> %s\n", nrow(man), sum(man$bytes) / 1e9, OUT))
