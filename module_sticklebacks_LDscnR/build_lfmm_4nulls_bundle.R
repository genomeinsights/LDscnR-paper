## =====================================================================
## module_sticklebacks_LDscnR / build_lfmm_4nulls_bundle.R
##
## Build the SELF-CONTAINED input bundle for the cluster B=100 LFMM four-null run
## (run_lfmm_4nulls_3sp.R). Extends the earlier single-null portable bundle with the
## three fields the other nulls need: the LD-independent GRM (genetic MVN null), the
## per-locality label (regional permutation null), and the GPS coordinates (spatial
## MVN-kernel null). Everything the run needs travels in one .rds -- no raw VCF, no
## external volumes on the cluster.
##
## Run from the LDscnR-paper root:
##   Rscript module_sticklebacks_LDscnR/build_lfmm_4nulls_bundle.R
## Writes /Users/petrikem/3sp_lfmm_perm/data/3sp_lfmm_input_4nulls.rds
## =====================================================================
suppressMessages({ library(data.table) })

BND  <- "module_sticklebacks_LDscnR/data/3sp_LDscnR_data.rds"
RAW  <- "~/gitlab/LD-scaling-genome-scans/empirical_data/3sp/3sp_data.RData"
LAT  <- "module_sticklebacks_LDscnR/data/3sp_latent_basis.rds"   # shared top-K PC basis
OUT  <- "/Users/petrikem/3sp_lfmm_perm/data/3sp_lfmm_input_4nulls.rds"

d  <- readRDS(BND); map <- as.data.table(d$map)
if (!file.exists(LAT)) stop("run build_latent_basis_3sp.R first (shared latent basis missing)")
pc_latent <- readRDS(LAT)
e  <- new.env(); load(path.expand(RAW), envir = e); ph <- as.data.table(e$pheno_3sp)
stopifnot(nrow(ph) == nrow(d$GTs),
          all(d$eco == as.integer(ph$ecotype == "Marine")))

## per-individual keys and the pop/locality lookup tables the four generators use
pop         <- ph$pop_ID
pop_ecotype <- unique(ph[, .(pop = pop_ID, ecotype)])
pop_eco_loc <- unique(ph[, .(pop = pop_ID, ecotype, loc = pop_locality)])
coords      <- as.matrix(ph[, .(GPS_N_updated, GPS_E_updated)])
stopifnot(!anyDuplicated(pop_ecotype$pop), !anyDuplicated(pop_eco_loc$pop))

out <- list(
  GTs         = d$GTs,            # 117 x 790,578 genotype matrix
  ld_ws       = d$ld_ws,          # per-marker LD weight matrix (for ld_cscore)
  marker      = map$marker,       # marker ids, ordered as GTs columns
  eco         = d$eco,            # observed marine(1)/freshwater(0) ecotype
  pop         = pop,             # per-individual population id
  loc         = ph$pop_locality,  # per-individual locality (4 strata) -- regional perm
  coords      = coords,           # per-individual GPS (N,E) -- spatial kernel
  GRM         = d$GRM,            # LD-independent kinship (117x117) -- K-MVN (genetic) null
  pc_latent   = pc_latent,        # top-K genotype PCs (U, d, K) -- latent null (LFMM home field)
  pop_ecotype = pop_ecotype,      # one ecotype per pop -- global perm
  pop_eco_loc = pop_eco_loc)      # pop x ecotype x locality -- regional perm

if (!dir.exists(dirname(OUT))) dir.create(dirname(OUT), recursive = TRUE)
saveRDS(out, OUT)
cat(sprintf("wrote %s\n  GTs %dx%d ; GRM %dx%d ; latent U %dx%d ; %d pops ; %d localities (%s)\n  size %.0f MB\n",
            OUT, nrow(out$GTs), ncol(out$GTs), nrow(out$GRM), ncol(out$GRM),
            nrow(pc_latent$U), ncol(pc_latent$U),
            nrow(pop_ecotype), length(unique(out$loc)), paste(sort(unique(out$loc)), collapse=","),
            file.info(OUT)$size / 1e6))
