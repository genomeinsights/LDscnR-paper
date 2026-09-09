## =============================================================================
## module_3sp/R/legacy_source.R
##
## THE CODE THAT ACTUALLY PRODUCED data/lfmm_F.rds, VERBATIM. Extracted,
## unedited beyond this header and the section markers below, from the legacy
## repository:
##   ~/gitlab/LD-scaling-genome-scans/empirical_data/R/3sp_sticklebacks.R
##   (parsing: lines 34-48; LFMM: lines 121-140; read 2026-09-09)
##
## WHY THIS FILE EXISTS. 00_config.R's LFMM_SOURCE comment and 04_lfmm.R's own
## header used to say "no script in any repository produces [lfmm_F.rds]".
## That was true, and it is exactly the question a reader -- or a referee --
## would ask: which bytes produced this, and what code ran on them. This file
## is that code, kept for provenance. It does not make lfmm_F.rds
## reproducible (see the caveat below); it makes the history visible.
##
## THIS FILE IS NOT EXECUTABLE HERE AND IS NOT MEANT TO BE. It assumes the
## legacy repository's working directory, several of its own R/ helpers
## sourced at that script's top (only the two actually needed for parsing --
## create_gds_from_geno, read_gds_ids -- are reproduced below), and external
## LEA/SNPRelate calls that write to ./tmp/ and run for a long time.
## module_3sp's own current pipeline (01_inputs.R onward) is how this
## dataset's genotype/map/GRM handling is reproduced today; this file answers
## "how was lfmm_F.rds first produced", not "how do I reproduce it now".
##
## SCOPE, DELIBERATELY NARROW. Only the two things module_3sp actually
## inherits from the legacy script: data parsing (raw bundle -> MAF-filtered
## GDS) and the LFMM analysis that produced lfmm_F.rds. NOT reproduced here,
## because module_3sp does not inherit any of it: the legacy script's own
## LD-decay estimation, single-linkage pruning, its own EMMAX call (a
## different GRM and a different emmax() from what 03_EMMAX.R uses), the
## C-score/OR machinery, and its Manhattan plotting.
## =============================================================================

## -----------------------------------------------------------------------------
## Supporting utility, copied from LD-scaling-genome-scans/R/ld_decay.R. That
## file is mostly the LD-decay machinery module_3sp does NOT inherit; these two
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
## PARSING -- verbatim, 3sp_sticklebacks.R lines 34-48.
## -----------------------------------------------------------------------------
if(length(list.files("./empirical_data/"))==0) message("Please download data form Zenond")

load("./empirical_data/3sp/3sp_data.RData") ## contains SNP_res_3sp, GTs_3sp and pheno_3sp

##filter by maf
keep <- map_3sp$maf>0.1
GTs_3sp <- GTs_3sp[,map_3sp$maf>0.1]
map_3sp <- map_3sp[maf>0.1]
#if(any(ls() %in% "gds_3sp")) snpgdsClose(gds_3sp); unlink(gds_3sp)

gds_3sp <- create_gds_from_geno(geno = GTs_3sp, map=map_3sp,"gds_3sp.gds")


# get id's for later use
ids <- read_gds_ids(gds_3sp)

## -----------------------------------------------------------------------------
## LFMM -- verbatim, 3sp_sticklebacks.R lines 121-140. This is the code that
## produced data/lfmm_F.rds (SHA-256 recorded in data/PROVENANCE.csv).
## -----------------------------------------------------------------------------
#### ==== Latent factor mixed model analyses (LFMM) ==== ####
#lfmm takes very long time
if(!file.exists("./empirical_data/3sp/lfmm_F.rds")){
  library(LEA)
  phe <- as.numeric(as.factor(pheno_3sp$ecotype))

  write.lfmm(GTs_3sp, "./tmp/genotypes.lfmm")
  write.env(phe, "./tmp/gradients.env")
  project = NULL

  project = lfmm2("./tmp/genotypes.lfmm", "./tmp/gradients.env", K=7) # K is user defined
  pv = lfmm2.test(project, "./tmp/genotypes.lfmm", "./tmp/gradients.env",genomic.control = TRUE,full = TRUE)

  saveRDS(pv$f,"./empirical_data/3sp/lfmm_F.rds") # only F-value is needed

}else{
  map_3sp[,lfmm_F:=readRDS("./empirical_data/3sp/lfmm_F.rds")[keep]]
  map_3sp[,lfmm_P:=pf(lfmm_F,1,115,lower.tail = FALSE)]
  map_3sp[,lfmm_q:=p.adjust(lfmm_P,"fdr")]
}
