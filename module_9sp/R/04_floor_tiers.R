## =============================================================================
## module_9sp/R/04_floor_tiers.R
##
## Floor-stability tiers for the nine-spined panel (PK, 2026-10-07): regions
## at the reported marker-density floor (12) tiered by recovery at 0.5x and 2x
## that floor (6 and 24), observed only, with the lineage-residualised
## phenotype used by 03_EMMAX.R. Shared logic: module_3sp/R/helpers_floor_tiers.R.
## No EcoPeak-equivalent reference exists for this panel, so only tier counts
## are reported.
## =============================================================================
suppressMessages({library(data.table); library(LDscnR)})
devtools::load_all("~/gitlab/LDscnR")
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_9sp"), "R", "00_config.R"))
source(path.expand("~/gitlab/LDscnR-paper/module_3sp/R/helpers_floor_tiers.R"))
STAGE <- "04_floor_tiers"
say("=== %s ===\n\n", STAGE)

BUNDLE_PATH <- file.path(PATHS$out, "02_bundle", "bundle.rds")
SCAN_PATH <- file.path(PATHS$out, "03_EMMAX", "scan.rds")
b <- readRDS(BUNDLE_PATH)
sc <- readRDS(SCAN_PATH)
df <- density_floor(b$map)
say("[0] %.0f markers per Mb -> density floor %d (config SIZE_FLOOR = %d)\n", df$markers_per_mb, df$floor, SIZE_FLOOR)
stopifnot("SIZE_FLOOR must equal the marker-density rule" = df$floor == SIZE_FLOOR)

covar_obs <- b$pheno[[EMMAX_COVAR]]
eco_resid <- as.numeric(stats::resid(stats::lm(b$eco ~ covar_obs)))   # as in 03_EMMAX.R

t0 <- Sys.time()
ft <- floor_tiers(b$GTs, b$map, b$stage1, b$GRM, b$LD_decay, eco_resid, SIZE_FLOOR,
                  UNIT_REPR, REGION_ASSEMBLY, ALPHA)
say("[1] floors %s fitted in %.1f min\n", paste(ft$floors, collapse = "/"),
    as.numeric(difftime(Sys.time(), t0, units = "mins")))

for (st in c("unit", "simes")) {
  a <- ft[[st]]$tests$canon; r <- sc[[if (st == "unit") "consensus" else "simes"]]$test
  stopifnot(identical(a$units$significant, r$units$significant), nrow(a$regions) == nrow(r$regions))
}
say("[2] canonical floor reproduces 03_EMMAX.R (significant units and region count, both statistics)\n")

summ <- rbindlist(lapply(c("unit", "simes"), function(st) {
  ft[[st]]$regions[, .(n = .N), by = .(statistic, tier)][order(tier)]
}))
print(summ)
single <- rbindlist(lapply(c("unit", "simes"), function(st)
  ft[[st]]$all_floors[n_floors == 1, .(n_single_floor = .N), by = .(statistic, floor)]))
say("\n[3] single-floor survivors by floor:\n"); print(single)

dir.create(stage_dir(STAGE), recursive = TRUE, showWarnings = FALSE)
fwrite(rbindlist(lapply(c("unit", "simes"), function(st) ft[[st]]$regions), fill = TRUE),
       file.path(stage_dir(STAGE), "regions_by_tier.tsv"), sep = "\t")
fwrite(rbindlist(lapply(c("unit", "simes"), function(st) ft[[st]]$all_floors), fill = TRUE),
       file.path(stage_dir(STAGE), "regions_all_floors.tsv"), sep = "\t")
fwrite(summ, file.path(stage_dir(STAGE), "tier_counts.tsv"), sep = "\t")
write_receipt(STAGE, inputs = c(BUNDLE_PATH, SCAN_PATH),
              params = list(floors = ft$floors, alpha = ALPHA, unit_repr = UNIT_REPR, emmax_covar = EMMAX_COVAR),
              outputs = file.path(stage_dir(STAGE), c("regions_by_tier.tsv", "regions_all_floors.tsv", "tier_counts.tsv")))
say("\nwrote %s\n", stage_dir(STAGE))
