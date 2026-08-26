## =====================================================================
## regen_popgen_only.R -- rebuild the popgen/BGS summaries for bundles that
## already exist, WITHOUT touching the bundles themselves.
##
## Why not just re-run the pipeline: the driver skips a file whose bundle exists,
## so it would never reach the popgen stage; and forcing it would rebuild the
## bundle too. The decay fit subsamples pairs, so a rebuilt bundle differs from
## the one on disk (different ld_ws, GRM and p-values), which would silently
## invalidate every panel/pvals file already delivered from it.
##
## The popgen stage is deterministic given the raw tarball, so re-deriving it is
## safe and reproduces the previous numbers exactly, plus whatever columns have
## been added since (here: pi_bp). The GRM for the IBD statistic is READ from the
## existing bundle rather than recomputed, so the summaries stay tied to the
## bundle they describe.
##
##   Rscript regen_popgen_only.R <V> <c> <env|all> <chr|all> [tag]
## =====================================================================
Sys.setenv(SIM_RAW = Sys.getenv("SIM_RAW", "/Volumes/Nemo/Nemo_sim/Nemo_out_nobgs"))
## sourced from this script's own directory, so it works from anywhere
.here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
suppressMessages(source(file.path(.here, "parse_and_regen_sim_data.R")))

a   <- commandArgs(trailingOnly = TRUE)
V   <- if (length(a) >= 1) a[1] else "2"
CC  <- if (length(a) >= 2) a[2] else "1"
ENV <- if (length(a) >= 3) a[3] else "all"
CHR <- if (length(a) >= 4) a[4] else "all"
TG  <- if (length(a) >= 5) a[5] else "nobgs"
BUNDLES <- Sys.getenv("SIM_BUNDLES", file.path(SIM_ROOT, "regen_sim_data_nobgs"))
envs <- if (ENV == "all") as.character(1:10) else ENV
chrs <- if (CHR == "all") as.character(1:10) else CHR

for (e in envs) for (ch in chrs) {
  stem <- sprintf("adapt_%s_chr%s_V%s_c%s_env%s", TG, ch, V, CC, e)
  raw  <- file.path(RAW_DIR, paste0(stem, ".tgz"))
  bun  <- file.path(BUNDLES, paste0(stem, ".rds"))
  if (!file.exists(raw) || !file.exists(bun)) { message("skip (missing): ", stem); next }
  tmp <- file.path(TMP_ROOT, paste0("pg_", stem))
  ok <- tryCatch({
    d <- parse_raw(raw, tmp)
    if (is.null(d$popgen)) stop("popgen stage returned nothing")
    pg  <- d$popgen
    ibd <- ibd_from_grm(readRDS(bun)$GRM, d$env)      # GRM from the EXISTING bundle
    pg$summary <- cbind(pg$summary, as.data.table(ibd))
    saveRDS(pg, file.path(POPGEN_DIR, paste0(stem, ".rds")))
    fwrite(pg$summary, file.path(POPGEN_DIR, paste0(stem, "_summary.csv")))
    message(sprintf("  %s: pi_bp %.4g, SNP density %.0f/Mb", stem,
                    pg$summary$bgs_pi_bp_mean, pg$summary$bgs_snp_density))
    TRUE
  }, error = function(err) { message("  !! ", stem, ": ", conditionMessage(err)); FALSE })
  unlink(tmp, recursive = TRUE); gc()
}
