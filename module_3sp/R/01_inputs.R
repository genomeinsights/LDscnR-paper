## =============================================================================
## module_3sp/R/01_inputs.R
##
## RESOLVE AND HASH EVERY RAW INPUT. Computes nothing, writes no analysis output.
##
## WHY IT IS A STAGE RATHER THAN A PREAMBLE. The pipeline it replaces reached across
## five locations, two of them through symlinks that no path grep could see, and its
## own producer had been deleted for six hours before anyone noticed. Three separate
## errors this week came from running a sound check against the wrong version of
## something. A stage that does nothing but pin what "the inputs" currently ARE turns
## that from a thing to remember into a thing on disk.
##
## It is also the answer to a question a referee can ask: which bytes produced this?
## After this stage runs, the answer is a table.
##
## FAILS LOUDLY AND EARLY. A missing input here costs seconds; the same input missing
## in stage 02 costs however long the decay fit had been running.
## =============================================================================
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_3sp"), "R", "00_config.R"))
STAGE <- "01_inputs"

INPUTS <- c(
  raw_3sp = PATHS$raw_3sp,
  recmap  = PATHS$recmap,
  setNames(file.path(PATHS$ecopeaks, ECOPEAK_BEDS), paste0("ecopeak_", sub("\\.bed$", "", ECOPEAK_BEDS)))
)
if (identical(LFMM_SOURCE, "inherit")) INPUTS <- c(INPUTS, lfmm_F = PATHS$lfmm_F)

PARAMS <- list(maf_keep = MAF_KEEP, lfmm_source = LFMM_SOURCE, ecopeak_beds = ECOPEAK_BEDS)

say("=== %s ===\n\n", STAGE)

## ---- 1. do they exist, and what are they really? ----------------------------
## resolve symlinks: 3sp_LDscnR_data.rds and cache_3sp are both links to other trees,
## and a hash of a link is not a hash of the bytes.
say("[1] resolving %d declared inputs\n", length(INPUTS))
res <- rbindlist(lapply(names(INPUTS), function(nm) {
  p <- INPUTS[[nm]]
  ex <- file.exists(p)
  real <- if (ex) normalizePath(p, mustWork = FALSE) else NA_character_
  ## was: !identical(real, normalizePath(p, mustWork=FALSE)) -- `real` IS
  ## normalizePath(p, mustWork=FALSE), computed two lines up, so that comparison was always
  ## FALSE (a value compared to itself) and the whole first half of the OR was vestigial;
  ## is_link worked only via the second half, Sys.readlink(). Fixed to compare against the
  ## UNRESOLVED declared path `p`, which is the actual symlink test (independent audit,
  ## ldscnr-26/fe, 2026-09-04 -- caught in the one stage whose entire job is precision about
  ## what's real).
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
## The 887 MB raw file dominates; this is the slow part of an otherwise instant stage,
## and it is the whole point of the stage.
say("\n[2] hashing (SHA-256 over resolved paths)\n")
res[, sha256 := vapply(resolved, function(f) { t0 <- Sys.time(); h <- sha(f)
  say("    %-22s %s  (%.1fs)\n", basename(f), substr(h, 1, 16), 
      as.numeric(difftime(Sys.time(), t0, units = "secs"))); h }, "")]

## ---- 2b. the copied inputs must still match what was recorded at copy time ----
## Checked against data/PROVENANCE.csv, NOT against the originals in another repository.
## The point of copying them in was to stop depending on that repository; verifying
## against it would have put the dependency straight back.
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

## ---- 3. cross-check against the manifests other sessions maintain ------------
## 3sp_data/MANIFEST.md is ldscnr-fe's record of the same untracked inputs. Where the
## two overlap they must agree; a disagreement means one of us is holding a different
## file than we think, which is the failure this stage exists to catch.
man <- path.expand("~/gitlab/LDscnR-paper/3sp_data/MANIFEST.md")
if (file.exists(man)) {
  txt <- readLines(man, warn = FALSE)
  hits <- 0L; disagree <- character()
  for (i in seq_len(nrow(res))) {
    b <- basename(res$resolved[i])
    ln <- grep(b, txt, fixed = TRUE, value = TRUE)
    if (!length(ln)) next
    hits <- hits + 1L
    if (!any(grepl(res$sha256[i], ln, fixed = TRUE))) disagree <- c(disagree, b)
  }
  say("\n[3] cross-checked %d input(s) against 3sp_data/MANIFEST.md\n", hits)
  if (length(disagree)) {
    say("    [!] HASH DISAGREEMENT: %s\n", paste(disagree, collapse = ", "))
    say("    The manifest and this stage are describing different bytes. Resolve before\n")
    say("    running stage 02 -- do not assume one is stale.\n")
  } else if (hits) say("    all agree\n")
} else say("\n[3] 3sp_data/MANIFEST.md absent -- no cross-check available\n")

## ---- 4. record ---------------------------------------------------------------
dir.create(stage_dir(STAGE), recursive = TRUE, showWarnings = FALSE)
fwrite(res[, .(name, declared, resolved, bytes, sha256)],
       file.path(stage_dir(STAGE), "inputs.csv"))
write_receipt(STAGE, inputs = unname(res$resolved), params = PARAMS,
              outputs = file.path(stage_dir(STAGE), "inputs.csv"))
say("\n[4] wrote %s\n", file.path(stage_dir(STAGE), "inputs.csv"))
say("    receipt: %s\n", receipt_path(STAGE))
say("\n    Nothing was computed. The next stage is 02_bundle.R, which is the expensive one.\n")
