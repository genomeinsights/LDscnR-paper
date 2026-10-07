## =============================================================================
## final_analysis/R/24_single_floor_survivors.R
##
## PK (2026-10-07): which floor do "single-floor" regions come from? Every
## region at every floor of the density_half_double scheme (1/2/4) is tiered
## against the OTHER two floors with R/23's rule (shared significant core_snp
## + physical overlap). If single-floor survivors were concentrated at, and
## worse at, the lowest floor, that would support small units being noisier.
## Reads results/floor_stability_density_half_double_region_details.rds (R/22).
## =============================================================================
suppressMessages(library(data.table))
MODULE_ROOT <- path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis")
Sys.setenv(FLOOR_SCHEME = "density_half_double")
source(file.path(MODULE_ROOT, "R", "helpers_floor_stability.R"))
d <- readRDS(file.path(MODULE_ROOT, "results", "floor_stability_density_half_double_region_details.rds"))
d <- as.data.table(d)
out <- d[, {
  fl <- sort(unique(floor)); res <- list()
  for (F in c(1L, 2L, 4L)) {
    base <- .SD[floor == F]; if (!nrow(base)) next
    others <- setdiff(c(1L, 2L, 4L), F)
    rec <- sapply(others, function(o) match_region_sets(base, .SD[floor == o], use_core_snp = TRUE)$recovered)
    if (is.null(dim(rec))) rec <- matrix(rec, nrow = nrow(base))
    res[[length(res) + 1]] <- data.table(floor = F, TP = base$TP, n_markers = base$n_markers,
                                         n_floors = 1L + rowSums(rec))
  }
  rbindlist(res)
}, by = .(combo_id, method)]
cat("Regions found at only ONE of the floors 1/2/4, by the floor they occur at:\n")
print(out[n_floors == 1, .(n = .N, precision = round(mean(TP), 3), median_markers = as.numeric(median(n_markers))), by = .(floor)][order(floor)])
print(out[n_floors == 1, .(n = .N, precision = round(mean(TP), 3)), by = .(method, floor)][order(method, floor)])
cat("\nPrecision by floor x number of floors found at:\n")
print(dcast(out[, .(p = sprintf("%.2f (%d)", mean(TP), .N)), by = .(floor, n_floors)], floor ~ n_floors, value.var = "p"))
## rep-clustered bootstrap (10 replicate maps) of single-floor precision by floor, and floor-4 minus floor-1
out[, rep := as.integer(sub(".*_rep(\\d+)_.*", "\\1", combo_id))]
s1 <- out[n_floors == 1]
set.seed(1); B <- 2000; reps <- sort(unique(s1$rep))
bt <- t(replicate(B, { r <- sample(reps, replace = TRUE)
  x <- rbindlist(lapply(r, function(k) s1[rep == k]))
  p <- x[, mean(TP), by = floor][order(floor)]$V1; c(p, p[3] - p[1], p[3] - p[2]) }))
ci <- apply(bt, 2, quantile, c(0.025, 0.975), na.rm = TRUE)
colnames(ci) <- c("floor1", "floor2", "floor4", "f4_minus_f1", "f4_minus_f2"); print(round(ci, 3))
cat("\nOverall precision of ALL regions at each floor:\n"); print(out[, .(n = .N, precision = round(mean(TP), 3)), by = floor][order(floor)])

single <- out[n_floors == 1, .(n = .N, precision = mean(TP), median_markers = as.numeric(median(n_markers))), by = floor][order(floor)]
single[, `:=`(ci_lo = ci[1, 1:3], ci_hi = ci[2, 1:3])]
fwrite(single, file.path(MODULE_ROOT, "results", "floor_stability_density_half_double_single_floor.tsv"), sep = "\t")
fwrite(data.table(contrast = c("floor4_minus_floor1", "floor4_minus_floor2"),
                  diff = c(single$precision[3] - single$precision[1], single$precision[3] - single$precision[2]),
                  ci_lo = ci[1, 4:5], ci_hi = ci[2, 4:5]),
       file.path(MODULE_ROOT, "results", "floor_stability_density_half_double_single_floor_contrasts.tsv"), sep = "\t")
fwrite(dcast(out[, .(n = .N, precision = mean(TP)), by = .(floor, n_floors)], floor ~ n_floors, value.var = c("n", "precision")),
       file.path(MODULE_ROOT, "results", "floor_stability_density_half_double_by_floor_and_support.tsv"), sep = "\t")
