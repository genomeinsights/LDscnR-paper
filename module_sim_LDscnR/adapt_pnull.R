## =====================================================================
## module_sim_LDscnR / adapt_pnull.R
##
## Adapter: turn the nulls pipeline's output into the two inputs that
## analyse_one_dataset.R consumes.
##
##   run_sim_nulls.R  ->  cell_<id>.rds        context: markers, map, ld_ws, decay_sum
##                        pnull_<eng>_<ty>_<id>.rds   p_obs + P_surr for one engine x null
##
##   this script      ->  panel_<id>.rds       GTs, map, ld_ws, decay_sum
##                        pvals_<id>_<eng>_<ty>.rds   p_obs, p_perm (named/aligned)
##
## Two things it has to fix, and both are the reason this is a script rather than
## a one-liner:
##
## 1. THE CONTEXT HAS NO GENOTYPES. run_sim_nulls.R deliberately leaves GTs out --
##    it stores ld_ws once per cell instead of repeating it per engine x null, and
##    GTs would dwarf that. But ld_edges() needs genotypes to measure pairwise r^2
##    between candidate markers: decay_sum sets the THRESHOLD, GTs supply the
##    values. So the adapter re-pools GTs from the bundles, using the identical
##    rule run_sim_nulls.R::pool_cell() uses (files ordered by chromosome number,
##    markers prefixed R<i>_), and then CHECKS the result against the context's
##    `markers` rather than trusting that the rule still matches.
##
## 2. NAMES ARE STORED ONCE, NOT PER VECTOR. p_obs is unnamed and P_surr is a
##    bare markers-x-B matrix; the names live in context$markers. The pipeline
##    downstream requires names, and mis-ordered names are the one error it
##    cannot detect, so they are reattached here and verified.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/adapt_pnull.R <pnull_file.rds> [outdir]
##   Rscript module_sim_LDscnR/adapt_pnull.R results/nulls/pnull_emmax_env_orth_V2_c1_env2.rds
## Env vars:
##   SIM_DATA  bundle dir the nulls were built from
##             (default /Volumes/Nemo/Nemo_sim/regen_sim_data_nobgs)
##   TAG       bundle tag (default nobgs)
##
## The panel is shared by every engine x null of a cell, so it is written once and
## reused; only the pvals file differs per null type. That is the same split
## analyse_one_dataset.R is built around: panel = the study, pvals = one method.
## =====================================================================
suppressMessages({ library(data.table); library(LDscnR) })
`%||%` <- function(a, b) if (is.null(a)) b else a   # base R only gained this in 4.4

a <- commandArgs(trailingOnly = TRUE)
if (!length(a)) stop("usage: adapt_pnull.R <pnull_file.rds> [outdir]")
PNULL <- a[1]
OUT   <- if (length(a) >= 2) a[2] else dirname(PNULL)
SIM_DATA <- Sys.getenv("SIM_DATA", "/Volumes/Nemo/Nemo_sim/regen_sim_data_nobgs")
TAG      <- Sys.getenv("TAG", "nobgs")
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)

pn <- readRDS(PNULL)
for (f in c("p_obs", "P_surr", "cell_id", "engine", "null_type"))
  if (is.null(pn[[f]])) stop("`", basename(PNULL), "` is missing `", f, "` -- is this a run_sim_nulls.R output?")

ctx_file <- file.path(dirname(PNULL), pn$context %||% paste0("cell_", pn$cell_id, ".rds"))
if (!file.exists(ctx_file))
  stop("context file not found: ", ctx_file, "\n  (run_sim_nulls.R writes it once per cell, alongside the pnull files)")
ctx <- readRDS(ctx_file)
markers <- ctx$markers
cat(sprintf("[1] %s\n    cell %s | engine %s | null %s | B = %d | %d markers\n",
            basename(PNULL), pn$cell_id, pn$engine, pn$null_type, ncol(pn$P_surr), length(markers)))

## ---- shape checks before anything is reconstructed --------------------------
if (length(pn$p_obs) != length(markers))
  stop(sprintf("p_obs has %d values but the context lists %d markers", length(pn$p_obs), length(markers)))
if (nrow(pn$P_surr) != length(markers))
  stop(sprintf("P_surr has %d rows but the context lists %d markers", nrow(pn$P_surr), length(markers)))
if (nrow(ctx$map) != length(markers))
  stop(sprintf("context map has %d rows but %d markers", nrow(ctx$map), length(markers)))
if (!identical(as.character(ctx$map$marker), as.character(markers)))
  stop("context map$marker is not in the same order as context markers")

## ---- 1. re-pool GTs, mirroring run_sim_nulls.R::pool_cell() -----------------
cell <- ctx$cell
V <- as.character(cell$V); CC <- as.character(cell$c); ENV <- as.character(cell$env)
tg <- as.character(cell$tag %||% TAG)
files <- list.files(SIM_DATA, full.names = TRUE,
  pattern = sprintf("^adapt_%s_chr[0-9]+_V%s_c%s_env%s[.]rds$", tg, V, CC, ENV))
files <- files[order(as.integer(sub(".*_chr([0-9]+)_.*", "\\1", basename(files))))]
if (!length(files)) stop("no bundles for ", pn$cell_id, " in ", SIM_DATA, " -- set SIM_DATA")
if (length(files) != (cell$n_files %||% length(files)))
  stop(sprintf("cell was built from %d files but %d are present in %s",
               cell$n_files, length(files), SIM_DATA))

panel_f <- file.path(OUT, sprintf("panel_%s.rds", pn$cell_id))
if (file.exists(panel_f)) {
  cat(sprintf("[2] panel: reusing %s\n", basename(panel_f)))
} else {
  cat(sprintf("[2] panel: re-pooling GTs from %d bundles\n", length(files))); flush.console()
  gts <- vector("list", length(files))
  for (i in seq_along(files)) {
    d <- readRDS(files[i]); m <- as.data.table(d$map)
    G <- d$GTs; colnames(G) <- paste0("R", i, "_", m$marker)
    gts[[i]] <- G
  }
  GTs <- do.call(cbind, gts); rm(gts); invisible(gc())
  ## The check that makes the reconstruction trustworthy: every marker the nulls
  ## were computed on must be present, and GTs is reordered to the context's
  ## order rather than assumed to already match it.
  miss <- setdiff(markers, colnames(GTs))
  if (length(miss))
    stop(sprintf("%d marker(s) in the null are absent from the bundles (e.g. %s) -- SIM_DATA does not match the data the nulls were built from",
                 length(miss), paste(utils::head(miss, 3), collapse = ", ")))
  GTs <- GTs[, markers, drop = FALSE]

  ld_ws <- ctx$ld_ws
  if (is.null(rownames(ld_ws))) rownames(ld_ws) <- markers
  if (!identical(rownames(ld_ws), as.character(markers))) {
    if (!all(markers %in% rownames(ld_ws))) stop("context ld_ws does not carry every marker")
    ld_ws <- ld_ws[markers, , drop = FALSE]
  }
  saveRDS(list(GTs = GTs, map = as.data.table(ctx$map), ld_ws = ld_ws,
               decay_sum = ctx$decay_sum, cell = pn$cell_id,
               n_ind = ctx$n_ind %||% nrow(GTs), source = "run_sim_nulls.R context"),
          panel_f)
  cat(sprintf("    %d markers x %d individuals -> %s (%.0f MB)\n",
              ncol(GTs), nrow(GTs), basename(panel_f), file.size(panel_f) / 1e6))
}

## ---- 2. reattach names to the p-values --------------------------------------
p_obs <- stats::setNames(as.numeric(pn$p_obs), markers)
P_surr <- pn$P_surr
rownames(P_surr) <- markers
pv_f <- file.path(OUT, sprintf("pvals_%s_%s_%s.rds", pn$cell_id, pn$engine, pn$null_type))
saveRDS(list(p_obs = p_obs, p_perm = P_surr,
             basis = pn$null_type, engine = pn$engine, B = ncol(P_surr),
             cell = pn$cell_id, source = basename(PNULL)), pv_f)
cat(sprintf("[3] p-values: %d observed + %d surrogates -> %s\n",
            length(p_obs), ncol(P_surr), basename(pv_f)))

cat(sprintf("\nNext:\n  Rscript module_sim_LDscnR/analyse_one_dataset.R %s %s\n", panel_f, pv_f))
