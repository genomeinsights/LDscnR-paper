## ---------------------------------------------------------------------------
## extend_delet_map.R -- build params_V4 from params_V3 by adding deleterious loci
##
## Scenario B: 1000 deleterious loci per chromosome instead of 100. Everything
## else is carried over untouched:
##
##   map_ntrl<x>.txt, map_QTN<x>.txt, allelic_values_<x>.txt  copied byte-for-byte
##   the 100 existing delet positions                          kept, 900 added
##
## Leaving the ntrl/QTN maps identical is what keeps the bgs and nobgs arms
## paired: the nobgs .ini declares no deleterious trait at all, so its genetic
## map must be the same file the bgs arm uses for its neutral markers.
##
## New positions follow the rule in recombination_map_with_LRR.R: uniform in
## PHYSICAL space, avoiding every position already on the map. Uniform-in-bp is
## what makes low-recombination regions absorb more mutations per Morgan, which
## is the whole source of the BGS/recombination correlation.
##
## Chr1 and Chr2 share one draw, as they do now, so the pair stays matched on
## everything except QTN count.
##
## Usage:  Rscript extend_delet_map.R          (writes params_V4/ and rds/)
## ---------------------------------------------------------------------------

suppressMessages(library(data.table))

## Resolve paths relative to this script, so the folder can be moved or handed on.
script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) normalizePath(dirname(f)) else getwd()
}
NEMO <- Sys.getenv("BGS_NEMO_DIR", dirname(script_dir()))   # R/ -> NEMO/

V3_DIR  <- Sys.getenv("BGS_V3_DIR",
  "/Users/petrikem/gitlab/LD-scaling-genome-scans/Nemo_v3/chromosome_maps_500kb")
V3_RDS  <- Sys.getenv("BGS_V3_RDS",
  "/Users/petrikem/gitlab/LD-scaling-genome-scans/Nemo_v3/chromosome_maps_500kb_rds")
OUT_DIR <- Sys.getenv("BGS_OUT_DIR", file.path(NEMO, "params_V4"))
OUT_RDS <- file.path(OUT_DIR, "rds")

N_DELET_TOTAL <- 1000L    # per chromosome (was 100)
SEED          <- 20260825L

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_RDS, recursive = TRUE, showWarnings = FALSE)

## Nemo map format: "{{p1, p2, ...} {p1, p2, ...}}", one brace group per chromosome,
## positions ascending.
write_nemo_map <- function(dt, file) {
  blocks <- vapply(split(dt, dt$Chr), function(ch)
    paste0("{", paste(sort(ch$pos_nemo), collapse = ", "), "}"), character(1))
  write.table(paste0("{", paste(blocks, collapse = ""), "}"), file = file,
              col.names = FALSE, row.names = FALSE, quote = FALSE)
}

for (x in 1:10) {

  rd <- readRDS(file.path(V3_RDS, sprintf("rec_map%d.rds", x)))
  setorder(rd, Chr, bp)

  ## one draw shared by both chromosomes, as in the existing maps
  ref      <- rd[Chr == 1]
  taken    <- unique(rd$bp)
  n_new    <- N_DELET_TOTAL - ref[type == "delet", .N]
  bp_range <- range(ref$bp)

  set.seed(SEED + x)
  new_bp <- integer(0)
  while (length(new_bp) < n_new) {
    cand   <- as.integer(runif(2L * n_new, bp_range[1], bp_range[2]))
    cand   <- setdiff(unique(cand), c(taken, new_bp))
    new_bp <- c(new_bp, cand)
  }
  new_bp <- sort(new_bp[seq_len(n_new)])

  ## pos_nemo is piecewise-linear in bp (pos_nemo = trunc(cM * max(bp)/2000) and
  ## cM is linear between map anchors), so interpolating it directly is exact up
  ## to the truncation. Inside a zero-recombination block pos_nemo is constant,
  ## and new loci there correctly inherit that same value.
  new_rows <- rbindlist(lapply(1:2, function(ch) {
    src <- rd[Chr == ch]
    data.table(Chr = ch, type = "delet", bp = new_bp,
               cM       = approx(src$bp, src$cM,       xout = new_bp, rule = 2)$y,
               rec_rate = approx(src$bp, src$rec_rate, xout = new_bp, rule = 2)$y,
               pos_nemo = trunc(approx(src$bp, src$pos_nemo, xout = new_bp, rule = 2)$y),
               allelic_values = NA_real_)
  }))

  rd_new <- rbind(rd, new_rows, fill = TRUE)
  setorder(rd_new, Chr, bp)

  ## The params_V3 ntrl map already carries 88 duplicated bp on Chr1 (an artifact
  ## of unique()-before-trunc() in recombination_map_with_LRR.R). It is harmless --
  ## two markers at one map position -- and left alone, since removing them would
  ## change ntrl_loci and break comparability with the existing nobgs runs. What
  ## must hold is that the NEW deleterious positions are unique and collide with
  ## nothing already on the map.
  stopifnot(rd_new[type == "delet", .N, by = Chr]$N == N_DELET_TOTAL,
            !anyDuplicated(new_bp),
            !any(new_bp %in% taken),
            rd_new[Chr == 1, all(diff(pos_nemo) >= 0)])

  write_nemo_map(rd_new[type == "delet"], file.path(OUT_DIR, sprintf("map_delet%d.txt", x)))
  saveRDS(rd_new, file.path(OUT_RDS, sprintf("rec_map%d.rds", x)))

  for (f in c(sprintf("map_ntrl%d.txt", x), sprintf("map_QTN%d.txt", x),
              sprintf("allelic_values_%d.txt", x)))
    file.copy(file.path(V3_DIR, f), file.path(OUT_DIR, f), overwrite = TRUE)

  ## the dispersal matrices and environmental optima live one level up from the
  ## chromosome maps and are not chromosome-specific; copy once
  if (x == 1) {
    aux <- dirname(V3_DIR)
    for (f in c(sprintf("disp_mat_%s.txt", c("1", "1.5", "2")),
                sprintf("env_%d.txt", 1:10)))
      if (!file.copy(file.path(aux, f), file.path(OUT_DIR, f), overwrite = TRUE))
        stop("missing auxiliary input: ", file.path(aux, f))
  }

  ## deleterious loci per Morgan, the quantity that sets BGS strength
  d <- rd_new[Chr == 1 & type == "delet"]
  M <- diff(range(rd_new[Chr == 1, pos_nemo])) * 3.3e-5 / 100
  cat(sprintf("chr%-2d  delet %d/chr  map %.1f cM  %.0f loci/Morgan  %d in zero-recomb blocks\n",
              x, nrow(d), M * 100, nrow(d) / M,
              d[, sum(duplicated(pos_nemo) | duplicated(pos_nemo, fromLast = TRUE))]))
}

cat("\nwritten to", OUT_DIR, "\n")
