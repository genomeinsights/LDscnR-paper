## =============================================================================
## module_3sp/R/14_floor_tiers.R
##
## Floor-stability tiers for the three-spined panel (PK, 2026-10-07): regions
## at the reported marker-density floor (8) tiered by recovery at 0.5x and 2x
## that floor (4 and 16), observed only -- see helpers_floor_tiers.R. EcoPeak
## enrichment is then tested per tier with the same within-chromosome rotation
## null as 03_EMMAX.R. The simulation counterpart is
## module_sim_3sp53/final_analysis (FLOOR_SCHEME = density_half_double).
## =============================================================================
suppressMessages({library(data.table); library(LDscnR)})
devtools::load_all("~/gitlab/LDscnR")
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_3sp"), "R", "00_config.R"))
source(file.path(PATHS$module, "R", "helpers_floor_tiers.R"))
STAGE <- "14_floor_tiers"
say("=== %s ===\n\n", STAGE)

BUNDLE_PATH <- file.path(PATHS$out, "02_bundle", "bundle.rds")
SCAN_PATH <- file.path(PATHS$out, "03_EMMAX", "scan.rds")
b <- readRDS(BUNDLE_PATH)
sc <- readRDS(SCAN_PATH)
df <- density_floor(b$map)
say("[0] %.0f markers per Mb -> density floor %d (config SIZE_FLOOR = %d)\n", df$markers_per_mb, df$floor, SIZE_FLOOR)
stopifnot("SIZE_FLOOR must equal the marker-density rule" = df$floor == SIZE_FLOOR)

t0 <- Sys.time()
ft <- floor_tiers(b$GTs, b$map, b$stage1, b$GRM, b$LD_decay, b$eco, SIZE_FLOOR,
                  UNIT_REPR, REGION_ASSEMBLY, ALPHA)
say("[1] floors %s fitted in %.1f min\n", paste(ft$floors, collapse = "/"),
    as.numeric(difftime(Sys.time(), t0, units = "mins")))

## the canonical floor must reproduce 03_EMMAX.R exactly
for (st in c("unit", "simes")) {
  a <- ft[[st]]$tests$canon; r <- sc[[if (st == "unit") "consensus" else "simes"]]$test
  stopifnot(identical(a$units$significant, r$units$significant), nrow(a$regions) == nrow(r$regions))
}
say("[2] canonical floor reproduces 03_EMMAX.R (significant units and region count, both statistics)\n")

KG <- copy(sc$ecopeaks); LEN <- copy(sc$chrom_lengths)
rot_row <- function(reg, label) {
  if (!nrow(reg)) return(data.table(set = label, n = 0L))
  rr <- ld_region_rotation(reg[, .(Chr, from, to)], KG, LEN, scheme = ROTATION_SCHEME,
                           n_rotations = N_ROTATIONS, seed = 1L)
  data.table(set = label, n = nrow(reg), ecopeak = rr$observed, null_mean = rr$null_mean, fold = rr$fold, p = rr$p)
}
summ <- rbindlist(lapply(c("unit", "simes"), function(st) {
  reg <- ft[[st]]$regions
  rbindlist(list(
    rot_row(reg, "all"), rot_row(reg[tier >= 2], ">=2/3"), rot_row(reg[tier == 3], "3/3"),
    rot_row(reg[tier == 2], "2/3 only"), rot_row(reg[tier == 1], "1/3 only")
  ))[, statistic := st]
}), fill = TRUE)
print(summ)

## single-floor survivors: regions at each floor (4/8/16) found at none of the other two
single <- rbindlist(lapply(c("unit", "simes"), function(st) {
  af <- ft[[st]]$all_floors
  rbindlist(lapply(c("lo", "canon", "hi"), function(bs)
    rot_row(af[base == bs & n_floors == 1], sprintf("only at floor %d", ft$floors[[bs]]))))[, statistic := st]
}), fill = TRUE)
say("\n[3] single-floor survivors by floor:\n"); print(single)

dir.create(stage_dir(STAGE), recursive = TRUE, showWarnings = FALSE)
regions <- rbindlist(lapply(c("unit", "simes"), function(st) ft[[st]]$regions), fill = TRUE)
fwrite(regions, file.path(stage_dir(STAGE), "regions_by_tier.tsv"), sep = "\t")
fwrite(summ, file.path(stage_dir(STAGE), "ecopeak_by_tier.tsv"), sep = "\t")
fwrite(single, file.path(stage_dir(STAGE), "ecopeak_single_floor.tsv"), sep = "\t")
fwrite(rbindlist(lapply(c("unit", "simes"), function(st) ft[[st]]$all_floors), fill = TRUE),
       file.path(stage_dir(STAGE), "regions_all_floors.tsv"), sep = "\t")
write_receipt(STAGE, inputs = c(BUNDLE_PATH, SCAN_PATH),
              params = list(floors = ft$floors, alpha = ALPHA, unit_repr = UNIT_REPR,
                            n_rotations = N_ROTATIONS, rotation_scheme = ROTATION_SCHEME),
              outputs = file.path(stage_dir(STAGE), c("regions_by_tier.tsv", "ecopeak_by_tier.tsv", "ecopeak_single_floor.tsv", "regions_all_floors.tsv")))
say("\nwrote %s\n", stage_dir(STAGE))
