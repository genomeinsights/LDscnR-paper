## =============================================================================
## module_9sp/R/legacy_source.R
##
## THE CODE THAT ACTUALLY PRODUCED data/lfmm_F.rds, VERBATIM. Extracted,
## unedited beyond this header and the section markers below, from the legacy
## repository:
##   ~/gitlab/LD-scaling-genome-scans/empirical_data/R/9sp_sticlebacks.R
##   (parsing: lines 34-51; LFMM: lines 125-146; read 2026-09-09)
##
## WHY THIS FILE EXISTS. Matches module_3sp/R/legacy_source.R: 00_config.R's
## LFMM_SOURCE comment says lfmm_F.rds is "INHERITED, not computed by this
## pipeline". This file is the code that actually computed it, kept for
## provenance -- module_9sp has no current script that regenerates LFMM
## p-values (03_EMMAX.R and the figures use EMMAX consensus/Simes only; see
## PERMUTATION_VS_STRUCTURE.md), so unlike module_3sp there is no numbered
## stage here to point back at this file from, but the same "which bytes
## produced this" question applies to data/lfmm_F.rds either way.
##
## THIS FILE IS NOT EXECUTABLE HERE AND IS NOT MEANT TO BE. It assumes the
## legacy repository's working directory, several of its own R/ helpers
## sourced at that script's top (only the two actually needed for parsing --
## create_gds_from_geno, read_gds_ids -- are reproduced below), and external
## LEA/SNPRelate calls that write to ./tmp/ and run for a long time.
## module_9sp's own current pipeline (01_inputs.R onward) is how this
## dataset's genotype/map/GRM handling is reproduced today; this file answers
## "how was lfmm_F.rds first produced", not "how do I reproduce it now".
##
## VERBATIM MEANS VERBATIM, INCLUDING ITS QUIRKS. The LFMM branch below reads
## `data_9sp$ecotype_bin`, but the parsing block above it already ran
## `rm(data_9sp)`; the LFMM branch also ends in `q("no")`, which quits R
## entirely. In the original script these blocks were guarded by
## `if (!file.exists(...))` and evidently run as separate interactive steps
## over more than one R session, not straight through top to bottom. Nothing
## below has been "fixed" to make it runnable end-to-end -- that would no
## longer be the code that actually ran.
##
## SCOPE, DELIBERATELY NARROW. Only the two things module_9sp actually
## inherits from the legacy script: data parsing (raw bundle -> MAF-filtered
## GDS) and the LFMM analysis that produced lfmm_F.rds. NOT reproduced here,
## because module_9sp does not inherit any of it: the legacy script's own
## LD-decay estimation, single-linkage pruning, its own EMMAX call (a
## different GRM and a different emmax() from what 03_EMMAX.R uses), the
## C-score/OR machinery, and its Manhattan plotting.
## =============================================================================

## -----------------------------------------------------------------------------
## Supporting utility, copied from LD-scaling-genome-scans/R/ld_decay.R. That
## file is mostly the LD-decay machinery module_9sp does NOT inherit; these two
## generic GDS helpers just happened to live in it, and the parsing block below
## calls both.
## -----------------------------------------------------------------------------
create_gds_from_geno <- function(geno, map, gds_path) {
  stopifnot(ncol(geno) == nrow(map))
  snpgdsCreateGeno(
    gds_path,
    genmat         = t(round(geno)),     # SNPRelate expects SNP x sample if snpfirstdim=TRUE
    sample.id      = paste0("ind_", seq_len(nrow(geno))), ## sample name not important
    snp.id         = map$marker,
    snp.chromosome = map$Chr,
    snp.position   = map$Pos,
    snpfirstdim    = TRUE
  )
  snpgdsOpen(gds_path)
}

read_gds_ids <- function(gds) {
  list(
    snp_id  = read.gdsn(index.gdsn(gds, "snp.id")),
    snp_chr = read.gdsn(index.gdsn(gds, "snp.chromosome")),
    snp_pos = read.gdsn(index.gdsn(gds, "snp.position"))
  )
}

## -----------------------------------------------------------------------------
## PARSING -- verbatim, 9sp_sticlebacks.R lines 34-51.
## -----------------------------------------------------------------------------
if(!file.exists("./empirical_data")) message("Please download data form Zenondo")

data_9sp <-readRDS("./empirical_data/9sp/9sp_data.rds") ## contains SNP_res_9sp, GTs_9sp and pheno_9sp
GTs_9sp <- data_9sp$GT
map_9sp <- data_9sp$map
pheno_9sp <- data_9sp$pheno
eco_9sp <- data_9sp$ecotype_bin
##filter by maf
keep <- map_9sp$maf>0.1
GTs_9sp <- GTs_9sp[,map_9sp$maf>0.1]
map_9sp <- map_9sp[maf>0.1]
rm(data_9sp)
gc()
if(any(ls() %in% "gds_9sp")) snpgdsClose(gds_9sp); unlink(gds_9sp)
gds_9sp <- create_gds_from_geno(geno = GTs_9sp, map=map_9sp,"gds_9sp.gds")

# get id's for later use
ids <- read_gds_ids(gds_9sp)

## -----------------------------------------------------------------------------
## LFMM -- verbatim, 9sp_sticlebacks.R lines 125-146. This is the code that
## produced data/lfmm_F.rds (SHA-256 recorded in data/PROVENANCE.csv). Depends
## on `data_9sp` from the parsing block above -- see the quirks note in the
## file header: as written, this branch assumes `data_9sp` is still in scope,
## which the parsing block's own `rm(data_9sp)` above contradicts.
## -----------------------------------------------------------------------------
#### ==== Latent factor mixed model analyses (LFMM) ==== ####
#lfmm takes very long time
if(!file.exists("./empirical_data/9sp/lfmm_F.rds")){
  library(LEA)

  eco_bin  <- as.numeric(as.factor(data_9sp$ecotype_bin))
  write.lfmm(GTs_9sp, "./tmp/genotypes.lfmm")
  write.env(eco_bin, "./tmp/gradients.env")
  project = NULL

  project = lfmm2("./tmp/genotypes.lfmm", "./tmp/gradients.env", K=data_9sp$pheno[,length(unique(pop_locality))]) # K is number of populations
  pv = lfmm2.test(project, "./tmp/genotypes.lfmm", "./tmp/gradients.env",genomic.control = TRUE,full = TRUE)

  saveRDS(pv$f,"./empirical_data/9sp/lfmm_F.rds") # only F-value is needed
  q("no")
}else{
  map_9sp[,lfmm_F:=readRDS("./empirical_data/9sp/lfmm_F.rds")[keep]]
  lfmm_gif = map_9sp[,median(lfmm_F)/qf(0.5,1,nrow(GTs_9sp),lower.tail = FALSE)] ## inflation factor
  map_9sp[,lfmm_F_GC:=lfmm_F/lfmm_gif]  ## genomic control
  map_9sp[,lfmm_P:=pf(lfmm_F_GC,1,nrow(GTs_9sp),lower.tail = FALSE)]
  map_9sp[,lfmm_q:=p.adjust(lfmm_P,"fdr")]
}
