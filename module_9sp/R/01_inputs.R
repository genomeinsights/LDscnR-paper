## =============================================================================
## module_9sp/R/01_inputs.R
##
## RESOLVE AND HASH EVERY RAW INPUT. Computes nothing, writes no analysis output.
## Same purpose and structure as module_3sp/R/01_inputs.R -- copied, not
## re-invented, because it is infrastructure, not a modelling choice.
##
## No ecopeaks/recmap/MANIFEST.md cross-check here: 9sp has neither an external
## EcoPeaks-equivalent reference nor a cross-repository manifest to check
## against (see 00_config.R's header). This stage is correspondingly shorter,
## not because it does less checking of what IS declared, but because there is
## less declared to check.
## =============================================================================
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_9sp"), "R", "00_config.R"))
STAGE <- "01_inputs"

INPUTS <- c(raw_9sp = PATHS$raw_9sp)
if (identical(LFMM_SOURCE, "inherit")) INPUTS <- c(INPUTS, lfmm_F = PATHS$lfmm_F)

PARAMS <- list(maf_keep = MAF_KEEP, lfmm_source = LFMM_SOURCE)

say("=== %s ===\n\n", STAGE)

## ---- 1. do they exist, and what are they really? ----------------------------
say("[1] resolving %d declared inputs\n", length(INPUTS))
res <- rbindlist(lapply(names(INPUTS), function(nm) {
  p <- INPUTS[[nm]]
  ex <- file.exists(p)
  real <- if (ex) normalizePath(p, mustWork = FALSE) else NA_character_
  data.table(name = nm, declared = p, resolved = real,
             is_link = ex && (!identical(real, p) || nzchar(Sys.readlink(p))),
             bytes = if (ex) file.size(p) else NA_real_, exists = ex)
}))
for (i in seq_len(nrow(res))) say("    %-22s %-8s %10s %s\n", res$name[i],
    if (res$exists[i]) "OK" else "MISSING",
    if (is.na(res$bytes[i])) "--" else format(res$bytes[i], big.mark = ","),
    if (isTRUE(res$is_link[i])) paste("-> ", res$resolved[i]) else "")

missing <- res[exists == FALSE]
if (nrow(missing)) {
  say("\n[!] %d input(s) missing. Nothing downstream can be trusted; stopping.\n", nrow(missing))
  stop("missing inputs: ", paste(missing$name, collapse = ", "))
}

## ---- 2. hash the resolved bytes ---------------------------------------------
say("\n[2] hashing (SHA-256 over resolved paths)\n")
res[, sha256 := vapply(resolved, function(f) { t0 <- Sys.time(); h <- sha(f)
  say("    %-22s %s  (%.1fs)\n", basename(f), substr(h, 1, 16),
      as.numeric(difftime(Sys.time(), t0, units = "secs"))); h }, "")]

## ---- 2b. the copied inputs must still match what was recorded at copy time ----
say("\n[2b] copied inputs against data/PROVENANCE.csv\n")
if (file.exists(PATHS$provenance)) {
  prov <- fread(PATHS$provenance)
  for (i in seq_len(nrow(prov))) {
    r <- res[name == prov$name[i]]
    if (!nrow(r)) next
    say("    %-14s %s\n", prov$name[i],
        if (identical(r$sha256[1], prov$sha256[i])) "matches the recorded hash" else
        "** DIFFERS FROM PROVENANCE.csv -- the local copy has changed since it was recorded **")
  }
} else say("    PROVENANCE.csv absent -- copies cannot be verified\n")

## ---- 3. record ---------------------------------------------------------------
dir.create(stage_dir(STAGE), recursive = TRUE, showWarnings = FALSE)
fwrite(res[, .(name, declared, resolved, bytes, sha256)],
       file.path(stage_dir(STAGE), "inputs.csv"))
write_receipt(STAGE, inputs = unname(res$resolved), params = PARAMS,
              outputs = file.path(stage_dir(STAGE), "inputs.csv"))
say("\n[3] wrote %s\n", file.path(stage_dir(STAGE), "inputs.csv"))
say("    receipt: %s\n", receipt_path(STAGE))
say("\n    Nothing was computed. The next stage is 02_bundle.R, which is the expensive one.\n")
