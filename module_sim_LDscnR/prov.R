## =====================================================================
## module_sim_LDscnR / prov.R
##
## write_prov(): stamp an output file with the parameter values it was actually
## produced with.
##
## WHY. Several scripts here take their output path from Sys.getenv("OUT"), so
## one script writes many differently-parameterised variants: engine_x_statistic.R
## produced both engine_x_statistic_bgs.csv and exs_full_bgs.csv off the same
## line. The grid (cell, tag, env, ...) survives as COLUMNS in the data, so those
## two are distinguishable. The SCALARS do not: SIM_DATA, ALPHA, LFMM_K and CORES
## appear nowhere in the output, and no committed driver records them. A file can
## therefore have a tracked producer and still not be regenerable, because the
## invocation is unrecorded -- which is worse than a missing producer, since it
## passes a "does a tracked script write this?" audit while failing in practice.
##
## SIDECAR, NOT A HEADER. The natural fix is a comment line at the top of the CSV,
## which travels with the file when it is copied out of the repo. These CSVs are
## read by fread()/read.csv() in dozens of places, though, and a leading comment
## would break those readers. A sidecar buys safety at the cost of being separable;
## write_prov() therefore also records enough to re-identify an orphaned copy
## (git commit, output path, grid) rather than relying on the pairing.
##
## Usage, immediately after the write:
##   source("module_sim_LDscnR/prov.R")   # from the repo root
##   fwrite(R, OUTF); write_prov(OUTF, list(SIM_DATA = SIM, CELLS = CELLS, ...))
## =====================================================================

write_prov <- function(path, params = list()) {
  if (!nzchar(path)) return(invisible(path))
  gitv <- function(args) tryCatch(suppressWarnings(system2("git", args, stdout = TRUE, stderr = FALSE)),
                                  error = function(e) character(0))
  head  <- gitv(c("rev-parse", "--short", "HEAD"))
  dirty <- length(gitv(c("status", "--porcelain"))) > 0
  a  <- commandArgs()
  sc <- sub("^--file=", "", grep("^--file=", a, value = TRUE))
  fmt <- function(v) paste(as.character(v), collapse = ",")
  lines <- c(
    "# provenance stamp -- see module_sim_LDscnR/prov.R",
    sprintf("output    %s", path),
    sprintf("script    %s", if (length(sc)) sc[1] else "(interactive)"),
    sprintf("git       %s%s", if (length(head)) head[1] else "(unknown)",
            if (isTRUE(dirty)) "  WORKING TREE DIRTY -- stamp may not match any commit" else ""),
    sprintf("when      %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    sprintf("R         %s", as.character(getRversion())),
    sprintf("LDscnR    %s", tryCatch(as.character(utils::packageVersion("LDscnR")),
                                     error = function(e) "(not installed)")),
    "params",
    if (length(params)) vapply(names(params), function(k) sprintf("  %-9s %s", k, fmt(params[[k]])),
                               character(1)) else "  (none recorded)")
  writeLines(lines, paste0(path, ".prov"))
  invisible(path)
}
