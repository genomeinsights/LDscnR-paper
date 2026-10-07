## =============================================================================
## module_3sp/R/15_marker_postfilter.R
##
## The size floor as a POST-FILTER on single-SNP scans (PK, 2026-10-07): the
## floor is not specific to LD aggregation. Marker-wise LFMM (and, for
## comparison, EMMAX) discoveries are taken at their own genome-wide BH
## threshold (q < 0.05, never recomputed), and at each floor only discoveries
## belonging to a Stage-1 cluster of at least that many markers are kept.
## Floor 1 keeps every discovery (singletons included) and is the
## unrestricted marker-wise comparator.
##
## Region assembly is the package's own (ld_outlier_test(assembly =
## "stage2_discovered")), so regions are built exactly as for the LD-aggregated
## analyses: every eligible cluster containing a discovery is passed in as a
## "significant unit" (p = 0; all others p = 1, so BH marks exactly these).
## EcoPeak enrichment uses the same within-chromosome rotation null as
## 03_EMMAX.R. Floors match the one-factor sensitivity table (stage 10) plus 1
## and 2.
## =============================================================================
suppressMessages({library(data.table); library(LDscnR)})
devtools::load_all("~/gitlab/LDscnR")
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_3sp"), "R", "00_config.R"))
STAGE <- "15_marker_postfilter"
say("=== %s ===\n\n", STAGE)

BUNDLE_PATH <- file.path(PATHS$out, "02_bundle", "bundle.rds")
SCAN_PATH <- file.path(PATHS$out, "03_EMMAX", "scan.rds")
LFMM_PATH <- file.path(PATHS$out, "04_lfmm", "lfmm.rds")
b <- readRDS(BUNDLE_PATH); sc <- readRDS(SCAN_PATH); lf <- readRDS(LFMM_PATH)
FLOORS <- c(1L, 2L, 4L, 8L, 12L, 16L, 24L, 48L)

stopifnot(length(lf$lfmm_p) == nrow(b$map))
pm_emmax <- emmax_fast(emmax_setup(b$GTs, b$GRM), b$eco)
sig <- list(
  LFMM  = b$map$marker[p.adjust(lf$lfmm_p, "BH") <= ALPHA],
  EMMAX = b$map$marker[p.adjust(pm_emmax, "BH") <= ALPHA])
say("[0] marker-wise BH q <= %.2f discoveries: LFMM %d, EMMAX %d\n", ALPHA, length(sig$LFMM), length(sig$EMMAX))

KG <- copy(sc$ecopeaks); LEN <- copy(sc$chrom_lengths)
one <- function(engine, fl) {
  u <- LDscnR:::.ld_outlier_units(b$stage1, b$map, fl)
  hit <- vapply(u$members, function(mm) any(mm %chin% sig[[engine]]), logical(1))
  n_snps <- sum(sig[[engine]] %chin% unlist(u$members[hit], use.names = FALSE))
  if (!any(hit)) return(data.table(engine = engine, floor = fl, snps_kept = 0L, clusters = 0L, regions = 0L))
  tt <- ld_outlier_test(b$stage1, b$map, as.numeric(!hit), statistic = "unit", size_floor = fl,
                        alpha = ALPHA, assembly = "stage2_discovered", GTs = b$GTs, LD_decay = b$LD_decay,
                        score_threshold = REGION_ASSEMBLY$score_threshold,
                        distance_threshold = REGION_ASSEMBLY$distance_threshold)
  stopifnot(identical(tt$units$significant, hit))
  reg <- as.data.table(tt$regions)
  rr <- ld_region_rotation(reg[, .(Chr, from, to)], KG, LEN, scheme = ROTATION_SCHEME,
                           n_rotations = N_ROTATIONS, seed = 1L)
  data.table(engine = engine, floor = fl, snps_kept = n_snps, clusters = sum(hit), regions = nrow(reg),
             ecopeak = rr$observed, frac = rr$observed / nrow(reg), null_mean = rr$null_mean,
             fold = rr$fold, p = rr$p)
}
t0 <- Sys.time()
res <- rbindlist(lapply(names(sig), function(e) rbindlist(lapply(FLOORS, function(f) one(e, f)), fill = TRUE)), fill = TRUE)
say("[1] %d engine x floor rows in %.1f min\n", nrow(res), as.numeric(difftime(Sys.time(), t0, units = "mins")))
print(res, digits = 3)

dir.create(stage_dir(STAGE), recursive = TRUE, showWarnings = FALSE)
fwrite(res, file.path(stage_dir(STAGE), "marker_postfilter_ecopeak.tsv"), sep = "\t")
write_receipt(STAGE, inputs = c(BUNDLE_PATH, SCAN_PATH, LFMM_PATH),
              params = list(floors = FLOORS, alpha = ALPHA, n_rotations = N_ROTATIONS,
                            rotation_scheme = ROTATION_SCHEME, region_assembly = REGION_ASSEMBLY),
              outputs = file.path(stage_dir(STAGE), "marker_postfilter_ecopeak.tsv"))
say("\nwrote %s\n", stage_dir(STAGE))
