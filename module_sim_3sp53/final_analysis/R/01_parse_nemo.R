## =============================================================================
## final_analysis/R/01_parse_nemo.R
##
## The corrected NEMO parser. CLAUDE_REANALYSIS_INSTRUCTIONS.md, Phase 1:
## "Nemo locus identifiers are one-based... Remove the current unconditional
## +1 transformation before joining loci to the type-specific reference-map
## rows."
##
## THE BUG THIS FIXES, and how it was confirmed (four independent lines of
## evidence, not just a docs read):
##   1. NEMO's own C++ source (nemo-release/src/servicenotifiers.cc:487-488)
##      writes the .map file's locus label AND its standalone "locus" column
##      as `loc_list[k]+1`, where loc_list[k] is the internal 0-based index.
##      i.e. NEMO's OWN OUTPUT FILE is already 1-based; the 0-based value
##      never reaches disk.
##   2. Empirical: on bgs/V0.5_c2/rep1/env1, offset=0 gives local_adapt_r2 =
##      0.374 (6/6 QTN with a valid effect) vs offset=+1's 0.251 (5/6, one
##      QTN loses its effect to an NA neighbour) -- the corrected join
##      explains 49% more variance in the very phenotype these QTN are
##      supposed to drive.
##   3. This exact bug was already found and fixed once before, in a
##      DIFFERENT module (module_sim_LDscnR/parse_and_regen_sim_data.R,
##      commit ee01f81, 2026-08-25) -- verified there against Nemo's own
##      chromosome/position columns: corrected indexing agreed 21/21 quant
##      loci on chromosome (Spearman 0.9997 on position); the old +1 offset
##      agreed 19/21 and left quant.101 unmapped entirely. That fix was later
##      silently reverted by an unrelated "sync" commit (007d431) and never
##      ported to module_sim/module_sim_3sp53's own parser, which is why it
##      is still present in both today.
##   4. ~/gitlab/LDscnR-NEMO/README.md documents it directly: "Nemo numbers
##      loci from 1, not 0... quant.101 exists when quanti_loci=101, and no
##      0-based scheme can produce that."
##
## Genotypes are NEVER affected by this bug (GTs columns are selected by
## nemo_marker identity, not by idx) -- only which reference-map row (Chr,
## Pos, type, allelic_values) each marker gets, i.e. every downstream "truth"
## quantity: QTN identity, distance-to-QTN, Va, precision/recall.
##
## Design: parse_nemo_run(tag, cell, rep, env, offset=0) does the full parse
## AND the assertion suite in one call, returning list(GTs, map, env, qa) --
## `qa` is a data.table, one row per assertion, machine-readable. `offset` is
## exposed ONLY so 01_parse_nemo_qa.R can call this same code path with
## offset=1 to reproduce the OLD (wrong) join for the three-run comparison;
## no caller in the real pipeline should ever pass offset=1.
## =============================================================================
suppressMessages({library(data.table)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "R", "00_config.R"))

## ---- assertion bookkeeping ---------------------------------------------------
## Every check is a hard stop on failure (per "Add hard assertions... stop
## immediately on a failed invariant") -- `.qa_rows` accumulates the PASSED
## checks (and is returned alongside the parse) so a human sees exactly what
## was verified, not just that nothing crashed.
.assert_qa <- function(qa_rows, label, cond, detail = "") {
  if (!isTRUE(cond)) {
    stop(sprintf("[FAILED ASSERTION] %s%s", label, if (nzchar(detail)) paste0(" -- ", detail) else ""), call. = FALSE)
  }
  qa_rows[[length(qa_rows) + 1]] <- data.table(check = label, status = "PASS", detail = detail)
  qa_rows
}

## ---- the join, parameterised by offset so the QA script can reproduce the ---
## OLD (wrong) behaviour for comparison. offset=0 is the only value the real
## pipeline ever uses.
.join_nemo_to_refmap <- function(nemo_map, refmap, offset) {
  nemo_map <- copy(nemo_map)
  nemo_map[, idx := idx_raw + offset]
  refmap <- copy(refmap)
  refmap[, indx := .I]

  idx_ntrl  <- refmap[type == "ntrl"][nemo_map[type == "ntrl",  idx], indx]
  idx_quant <- refmap[type == "QTN" ][nemo_map[type == "quant", idx], indx]
  idx_delet <- refmap[type == "delet"][nemo_map[type == "delet", idx], indx]
  refmap[idx_ntrl,  nemo_marker := nemo_map[type == "ntrl",  marker]]
  refmap[idx_quant, nemo_marker := nemo_map[type == "quant", marker]]
  refmap[idx_delet, nemo_marker := nemo_map[type == "delet", marker]]
  list(refmap = refmap, nemo_map = nemo_map)
}

parse_nemo_run <- function(tag, cell, rep, env, offset = 0L, keep_prefilter = FALSE) {
  qa <- list()
  stopifnot("offset must be 0 (corrected) or 1 (old, QA-only)" = offset %in% c(0L, 1L))

  run_dir  <- file.path(PATHS$raw_nemo, sprintf("adapt_%s_chr%d_%s_env%d", tag, rep, cell, env))
  geno_dir <- file.path(run_dir, "GENO")
  recmap_rds <- file.path(PATHS$raw_recmap_dir, sprintf("rec_map%d.rds", rep))
  env_txt    <- file.path(PATHS$raw_env_dir, sprintf("env_%d.txt", env))
  qa <- .assert_qa(qa, "raw inputs exist",
                   dir.exists(geno_dir) && file.exists(recmap_rds) && file.exists(env_txt),
                   sprintf("geno_dir=%s recmap=%s env=%s", geno_dir, recmap_rds, env_txt))

  files <- list.files(geno_dir, full.names = TRUE)
  map_file  <- files[grepl("\\.map$", files)]
  geno_file <- files[grepl("snp_geno", files)]
  qa <- .assert_qa(qa, "exactly one .map and one snp_geno file",
                   length(map_file) == 1L && length(geno_file) == 1L,
                   sprintf("%d map file(s), %d geno file(s)", length(map_file), length(geno_file)))

  map_nemo <- fread(map_file)
  GTs_raw  <- fread(geno_file)
  n_raw_individuals <- nrow(GTs_raw)
  qa <- .assert_qa(qa, "raw sample count is 320 before subsampling",
                   n_raw_individuals == 320L, sprintf("got %d", n_raw_individuals))

  nemo_map <- data.table(marker = map_nemo$trait.locus,
                         do.call(rbind, strsplit(map_nemo$trait.locus, ".", fixed = TRUE)))
  setnames(nemo_map, c("V1", "V2"), c("type", "idx_raw"))
  nemo_map[, idx_raw := as.numeric(idx_raw)]
  sample_info <- GTs_raw[, .(pop, ID)]
  GTs_raw <- as.matrix(GTs_raw[, 6:ncol(GTs_raw), with = FALSE])

  qa <- .assert_qa(qa, "each genotype column has exactly one NEMO map record",
                   ncol(GTs_raw) == nrow(nemo_map) && !anyDuplicated(nemo_map$marker),
                   sprintf("%d genotype columns, %d map records, %d duplicate marker names",
                           ncol(GTs_raw), nrow(nemo_map), anyDuplicated(nemo_map$marker)))

  refmap0 <- readRDS(recmap_rds)
  n_type <- refmap0[, .N, by = type]
  qa <- .assert_qa(qa, "type-specific idx is an integer in [1, n_loci_of_type]",
                   nemo_map[, all(idx_raw + offset == round(idx_raw + offset)) &&
                            all(sapply(unique(type), function(tt) {
                              rng <- range(idx_raw[type == tt] + offset)
                              lim <- n_type[type == c(ntrl="ntrl", quant="QTN", delet="delet")[tt], N]
                              rng[1] >= 1 && rng[2] <= lim
                            }))],
                   sprintf("offset=%d", offset))

  joined <- .join_nemo_to_refmap(nemo_map, refmap0, offset)
  refmap <- joined$refmap

  qa <- .assert_qa(qa, "every used NEMO map record has exactly one reference-map match within its type",
                   sum(!is.na(refmap$nemo_marker)) == nrow(nemo_map) &&
                     !anyDuplicated(refmap[!is.na(nemo_marker), nemo_marker]),
                   sprintf("%d of %d nemo records matched, %d duplicate refmap assignments",
                           sum(!is.na(refmap$nemo_marker)), nrow(nemo_map),
                           anyDuplicated(refmap[!is.na(nemo_marker), nemo_marker])))

  ## Reference-map structural check (offset-independent -- this is a property
  ## of rec_map<rep>.rds itself, not of the join): the unfiltered reference
  ## map must carry exactly 100 Chr1 QTN + 1 Chr2 QTN, i.e. quant.101 IS the
  ## single Chr2 QTN, in every reference map.
  refmap_qtn_by_chr <- refmap0[type == "QTN", .N, by = Chr][order(Chr)]
  if (offset == 0L) {
    ok <- nrow(refmap_qtn_by_chr) == 2L && refmap_qtn_by_chr[Chr == 1, N] == 100L && refmap_qtn_by_chr[Chr == 2, N] == 1L
    if (!ok) stop(sprintf(
      "[STOP POINT] reference map rec_map%d.rds does not carry 100 Chr1 QTN + 1 Chr2 QTN as expected -- got %s. Do not proceed without Petri's input.",
      rep, paste(sprintf("Chr%s=%d", refmap_qtn_by_chr$Chr, refmap_qtn_by_chr$N), collapse = ", ")), call. = FALSE)
    qa <- .assert_qa(qa, "reference map: 100 Chr1 QTN + 1 Chr2 QTN (quant.101 is the single Chr2 QTN)", TRUE,
                     sprintf("Chr1=%d, Chr2=%d", refmap_qtn_by_chr[Chr==1,N], refmap_qtn_by_chr[Chr==2,N]))
  }

  ## Per-run join check, ALWAYS applicable regardless of which QTN a given raw
  ## run happens to have saved (each run's own genotyper output only carries a
  ## small, variable subset of the 101 reference QTN -- confirmed directly:
  ## 2 to 14 of 101 across a 6-run/6-cell survey, with locus 101 itself present
  ## in only 1 of those 6 -- so a check hardcoded to "quant.101 must be in
  ## every run" is mis-specified, not a real invariant). What must ALWAYS be
  ## true: whichever QTN a run DOES contain, the join puts each one on the
  ## same chromosome the standalone reference map assigns to that index.
  if (offset == 0L) {
    quant_idx_here <- unique(nemo_map[type == "quant", idx_raw + offset])
    ref_qtn_only <- refmap0[type == "QTN"]
    expected_chr <- ref_qtn_only[quant_idx_here, Chr]
    got_chr <- as.integer(sub("Chr", "", refmap[match(sprintf("quant.%d", quant_idx_here), nemo_marker), Chr]))
    qa <- .assert_qa(qa, "every QTN locus present in this run lands on refmap0's chromosome for that index",
                     length(quant_idx_here) > 0 && all(expected_chr == got_chr),
                     sprintf("%d QTN checked (idx %s), includes idx 101 (the Chr2 QTN): %s",
                             length(quant_idx_here), paste(range(quant_idx_here), collapse="-"), 101 %in% quant_idx_here))
  }

  map <- refmap[nemo_marker %in% colnames(GTs_raw)]
  GTs <- GTs_raw[, map$nemo_marker]
  map[, Pos := bp]
  map[, Chr := paste0("Chr", Chr)]
  map[, marker := paste0(Chr, ":", Pos)]
  colnames(GTs) <- map$marker

  ord <- order(map$Chr, map$bp)
  map <- map[ord]; GTs <- GTs[, ord]
  dup <- duplicated(map$marker)
  if (any(dup)) { map <- map[!dup]; GTs <- GTs[, !dup] }

  qa <- .assert_qa(qa, "joined genotype columns and map rows remain in identical order",
                   identical(colnames(GTs), map$marker), "")
  qa <- .assert_qa(qa, "Chr/Pos/type/allelic_values come from the matched reference row",
                   all(map$Chr %in% c("Chr1", "Chr2")) && all(map$type %in% c("ntrl", "QTN", "delet")) &&
                     all(is.finite(map$Pos)),
                   "")
  qa <- .assert_qa(qa, "genotype values are valid diploid dosages",
                   all(GTs %in% c(0L, 1L, 2L, NA)), sprintf("range [%s, %s]", min(GTs, na.rm=TRUE), max(GTs, na.rm=TRUE)))
  qa <- .assert_qa(qa, "sample metadata remain aligned (nrow(GTs) == nrow(sample_info))",
                   nrow(GTs) == nrow(sample_info), "")

  n_by_chr_type_prefilter <- map[, .N, by = .(Chr, type)][order(Chr, type)]

  ## ---- environment ------------------------------------------------------------
  env_raw  <- scan(env_txt, what = character(), quiet = TRUE)
  env_raw  <- strsplit(env_raw, "}{", fixed = TRUE)[[1]]
  env_vals <- as.numeric(gsub("}}", "", gsub("{{", "", env_raw, fixed = TRUE), fixed = TRUE))
  env_grid <- data.table(expand.grid(x = 1:GRID_SIDE, y = GRID_SIDE:1))
  env_grid[, pop := .I]
  qa <- .assert_qa(qa, "environmental surface size matches the grid",
                   length(env_vals) == nrow(env_grid), sprintf("%d values vs %dx%d grid", length(env_vals), GRID_SIDE, GRID_SIDE))
  env_grid[, env := env_vals]
  env <- env_grid[match(sample_info$pop, env_grid$pop)]
  qa <- .assert_qa(qa, "every individual's env value matches its sampled population",
                   all(!is.na(env$env)) && all(env$pop == sample_info$pop), "")

  ## ---- subsample + MAF (counts reported before/after, per assertion list) -----
  keep_inds <- seq(1, nrow(GTs), by = SUBSAMPLE_STEP)
  GTs_sub <- GTs[keep_inds, ]; env_sub <- env[keep_inds]
  qa <- .assert_qa(qa, "subsampled individual count is 160",
                   nrow(GTs_sub) == 160L, sprintf("got %d (raw %d, step %d)", nrow(GTs_sub), n_raw_individuals, SUBSAMPLE_STEP))

  maf <- colSums(GTs_sub) / nrow(GTs_sub) / 2
  map[, MAF := pmin(maf, 1 - maf)]
  n_monomorphic <- sum(map$MAF == 0, na.rm = TRUE)
  keep_snps <- map$MAF > MAF_KEEP
  map_post <- map[keep_snps]; GTs_post <- GTs_sub[, keep_snps]
  n_by_chr_type_postfilter <- map_post[, .N, by = .(Chr, type)][order(Chr, type)]
  n_qtn_lost_to_maf <- sum(map$type == "QTN") - sum(map_post$type == "QTN")

  counts <- data.table(
    stage = c("post-join (pre-dup-removal already applied above)", "post-subsample (n_monomorphic among these)", sprintf("post-MAF>%.2f", MAF_KEEP)),
    n_markers = c(nrow(map), nrow(map), nrow(map_post)),
    n_QTN = c(sum(map$type == "QTN"), n_monomorphic, sum(map_post$type == "QTN")))

  map_out <- map_post[, .(Chr, Pos, marker, type, allelic_values, MAF)]
  map_out[, true_QTN := type == "QTN"]

  out <- list(GTs = GTs_post, map = map_out, env = env_sub,
             qa = rbindlist(qa),
             counts_by_chr_type = list(prefilter = n_by_chr_type_prefilter, postfilter = n_by_chr_type_postfilter),
             n_qtn_lost_to_maf = n_qtn_lost_to_maf,
             source = list(geno_dir = geno_dir, recmap = recmap_rds, env_file = env_txt, offset_used = offset))
  ## ADDITIVE only -- existing callers (run_combo.R, always keep_prefilter=FALSE)
  ## get an IDENTICAL return shape, so this never invalidates the parse
  ## stage's cached bundle/receipt. For callers that need the pre-MAF-filter
  ## genotypes/map (e.g. R/08_bgs_validation.R checking whether MAF>0.10
  ## filtering itself is masking BGS's effect on rare-variant diversity):
  ## GTs_sub/map here are post-subsample (160 ind), post-dup-removal, WITH
  ## the MAF column already attached, but before the MAF>0.10 cut.
  if (keep_prefilter) { out$GTs_prefilter <- GTs_sub; out$map_prefilter <- map }
  out
}
