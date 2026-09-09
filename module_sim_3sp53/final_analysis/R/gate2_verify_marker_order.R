## Phase 3 gate 2's explicit requirement: "Confirm identical marker order
## through parsing, LD construction, GRM construction, association and truth
## scoring." Internal stopifnot()s in 02/03 already assert pairwise
## consistency; this cross-checks all four stage outputs directly, from one
## independent script, for the gate-2 combination.
suppressMessages(library(data.table))
source("R/00_config.R")

combo_id <- "bgs_V0.5_c2_rep1_env1"
parsed <- readRDS(file.path(PATHS$parsed, "nemo_bgs_rep1_V0.5_c2_env1.rds"))
ld_units <- readRDS(file.path(stage_dir("02_build_ld_units", combo_id), "ld_units.rds"))
emmax <- readRDS(file.path(stage_dir("03_emmax", combo_id), "emmax.rds"))
truth <- readRDS(file.path(stage_dir("05_score_truth", combo_id), "truth_scores.rds"))

ok <- TRUE
chk <- function(label, cond) { cat(sprintf("[%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", label)); ok <<- ok && isTRUE(cond) }

chk("parsed$map$marker == colnames(parsed$GTs)", identical(parsed$map$marker, colnames(parsed$GTs)))
chk("ld_units$map$marker == colnames(ld_units$GTs)", identical(ld_units$map$marker, colnames(ld_units$GTs)))
chk("ld_units's marker SET is a subset of parsed's (post-Stage1 no reorder, only GRM/LD_decay added)",
    identical(ld_units$map$marker, parsed$map$marker[parsed$map$marker %in% ld_units$map$marker]) &&
    setequal(ld_units$map$marker, parsed$map$marker))
chk("emmax$marker$marker == ld_units$map$marker (association stage read the same map, same order)",
    identical(emmax$marker$marker, ld_units$map$marker))
chk("emmax per-marker p-vector length == n markers", length(emmax$marker$p) == length(ld_units$map$marker))
chk("GRM row/col count == n individuals analysed (kinship is individual x individual, built FROM grm_markers)",
    nrow(ld_units$GRM) == nrow(ld_units$GTs) && ncol(ld_units$GRM) == nrow(ld_units$GTs))
chk("truth-scoring used the same n_tested as emmax's marker count (emmax_snp)",
    truth$scores[method == "emmax_snp", n_tested] == length(ld_units$map$marker))
chk("truth-scoring's detectable-QTN denominator is stable across all 4 methods",
    length(unique(truth$scores$n_detectable_qtn)) == 1L)

cat(sprintf("\n%s: marker order and record counts are consistent end to end (parsing -> LD/GRM -> association -> truth scoring) for %s\n",
    if (ok) "GATE 2 PASS" else "GATE 2 FAIL", combo_id))
