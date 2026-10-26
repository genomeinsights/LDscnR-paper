## =============================================================================
## R/00_fetch_inputs.R
##
## FETCH AND VERIFY THE EMPIRICAL INPUTS (three- and nine-spined stickleback).
##
## The manuscript's empirical analyses read four files:
##   module_3sp/data/3sp_data.RData, module_3sp/data/lfmm_F.rds,
##   module_9sp/data/9sp_data.rds,   module_9sp/data/lfmm_F.rds
## They are derived from the data of Fang et al. (2021, Mol Biol Evol 38:4205,
## doi:10.1093/molbev/msab144; raw data Zenodo doi:10.5281/zenodo.4722879) and
## are archived, byte for byte, in the LD-scaling Zenodo deposit
## (doi:10.5281/zenodo.18466057, CC-BY-4.0) as 3sp.zip and 9sp.zip.
##
## This script downloads those two archives (or uses local copies), checks each
## archive against the MD5 that Zenodo publishes, extracts the four files and
## checks each against the SHA-256 recorded in module_*/data/PROVENANCE.csv.
## It writes to an output folder (default rebuild_check/) and never touches
## module_*/data/ unless --install is given, in which case verified files are
## copied there only if no file of that name exists yet.
##
## Usage, from the repository root:
##   Rscript R/00_fetch_inputs.R                       # download (~6.1 GB) and verify
##   Rscript R/00_fetch_inputs.R --zip-dir=PATH        # use local 3sp.zip / 9sp.zip
##   Rscript R/00_fetch_inputs.R --out=DIR --install
## Environment: LDSCNR_PAPER_ROOT (repository root; default: working directory).
## Needs R packages digest and jsonlite.
## =============================================================================
args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "(=|$)"), args, value = TRUE)
  if (!length(hit)) return(default)
  v <- sub(paste0("^--", name, "=?"), "", hit[1])
  if (nzchar(v)) v else TRUE
}
ROOT    <- normalizePath(Sys.getenv("LDSCNR_PAPER_ROOT", getwd()), mustWork = TRUE)
OUT     <- file.path(ROOT, opt("out", "rebuild_check"))
ZIP_DIR <- opt("zip-dir", NULL)
INSTALL <- isTRUE(opt("install", FALSE))
stopifnot("run from the LDscnR-paper root or set LDSCNR_PAPER_ROOT" =
            dir.exists(file.path(ROOT, "module_3sp")) && dir.exists(file.path(ROOT, "module_9sp")))

RECORD <- "18466057"   # LD-scaling deposit holding the exact analysis inputs
FILES <- data.frame(
  archive = c("3sp.zip", "3sp.zip", "9sp.zip", "9sp.zip"),
  member  = c("3sp/3sp_data.RData", "3sp/lfmm_F.rds", "9sp/9sp_data.rds", "9sp/lfmm_F.rds"),
  module  = c("module_3sp", "module_3sp", "module_9sp", "module_9sp"),
  stringsAsFactors = FALSE)
FILES$file <- basename(FILES$member)

sha256 <- function(f) digest::digest(f, algo = "sha256", file = TRUE)
md5    <- function(f) unname(tools::md5sum(f))

## ---- 1. expected hashes -------------------------------------------------------
meta <- jsonlite::fromJSON(sprintf("https://zenodo.org/api/records/%s", RECORD))
zen <- meta$files
zen$md5 <- sub("^md5:", "", zen$checksum)
cat(sprintf("Zenodo %s: %s (%s)\n", RECORD, meta$metadata$title, meta$metadata$license$id))
prov <- do.call(rbind, lapply(unique(FILES$module), function(m) {
  p <- utils::read.csv(file.path(ROOT, m, "data", "PROVENANCE.csv"), stringsAsFactors = FALSE)
  p$module <- m; p }))

## ---- 2. archives: local copy or download, then MD5 ----------------------------
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
zip_path <- character()
for (a in unique(FILES$archive)) {
  z <- zen[zen$key == a, ]
  stopifnot("archive not listed in the Zenodo record" = nrow(z) == 1)
  f <- if (!is.null(ZIP_DIR)) file.path(ZIP_DIR, a) else file.path(OUT, a)
  if (is.null(ZIP_DIR) && !file.exists(f)) {
    cat(sprintf("downloading %s (%.2f GB) ...\n", a, z$size / 1e9))
    options(timeout = max(3600, getOption("timeout")))
    utils::download.file(z$links$self, f, mode = "wb", quiet = FALSE)
  }
  stopifnot("archive missing" = file.exists(f))
  got <- md5(f)
  cat(sprintf("  %-8s md5 %s %s Zenodo\n", a, got, if (identical(got, z$md5)) "==" else "!="))
  if (!identical(got, z$md5)) stop(a, ": MD5 does not match the Zenodo record -- corrupt or wrong file.")
  zip_path[a] <- f
}

## ---- 3. extract and check SHA-256 against PROVENANCE.csv ----------------------
res <- FILES
res$sha256 <- res$expected <- NA_character_
for (i in seq_len(nrow(FILES))) {
  dest <- file.path(OUT, FILES$module[i])
  utils::unzip(zip_path[[FILES$archive[i]]], files = FILES$member[i], exdir = dest, junkpaths = TRUE)
  f <- file.path(dest, FILES$file[i])
  res$sha256[i]   <- sha256(f)
  res$expected[i] <- prov$sha256[prov$module == FILES$module[i] & prov$file == FILES$file[i]]
}
res$identical <- res$sha256 == res$expected
for (i in seq_len(nrow(res)))
  cat(sprintf("  %-10s %-15s %s  %s\n", res$module[i], res$file[i], substr(res$sha256[i], 1, 16),
              if (res$identical[i]) "identical to PROVENANCE.csv" else "** DIFFERS from PROVENANCE.csv **"))
utils::write.csv(res[, c("module", "file", "archive", "member", "sha256", "expected", "identical")],
                 file.path(OUT, "verification.csv"), row.names = FALSE)
if (!all(res$identical)) stop("one or more files differ from the inputs the manuscript used.")

## ---- 4. optional install into module_*/data/ ----------------------------------
if (INSTALL) for (i in seq_len(nrow(res))) {
  target <- file.path(ROOT, res$module[i], "data", res$file[i])
  if (file.exists(target)) { cat("  exists, left untouched:", target, "\n"); next }
  file.copy(file.path(OUT, res$module[i], res$file[i]), target)
  cat("  installed:", target, "\n")
}
cat(sprintf("all %d inputs verified; report: %s\n", nrow(res), file.path(OUT, "verification.csv")))
