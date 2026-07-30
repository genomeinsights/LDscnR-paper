######################################################
##  Parsing raw data from simulations
##    Only shown for example data
######################################################

#----------------------------------------------------------
#  Load packages
#----------------------------------------------------------

library(data.table)
library(parallel)
library(SNPRelate)
library(LEA)
library(factoextra)
library(LDscnR)
library(ggplot2)
library(wesanderson)

source("./R_LDscnR/compute_ld_structure.R")
source("./R_LDscnR/compute_ld_w.R")
#----------------------------------------------------------
# Helper: files and fixed parameters
#----------------------------------------------------------
gz_files <- list.files("./Nemo_sim/Nemo_out_bgs/",full.names = "TRUE")
done <- list.files("./parsed_sim_data/",full.names = "TRUE")


# remove files already done
gz_files <- gz_files[!gsub(".tgz",".rds",basename(gz_files)) %in% basename(done)]

## BATCH SELECTOR -- the one line to change between terminal sessions.
## Run this script as-is for the "bgs" batch in one terminal, and with
## "adapt_bgs" swapped for "adapt_nobgs" in a second terminal -- get_tmp_dir()
## below then automatically routes each batch to its own tmp folder
## (./tmp_bgs vs ./tmp_nobgs), so the two sessions never collide.
gz_files <- gz_files[grepl("V2_c2_",gz_files) & grepl("env3.tgz",gz_files) & grepl("adapt_nobgs",gz_files)]
# file_gz <- gz_files[2]
# file_gz <- "./Nemo_sim/Nemo_out_bgs//adapt_bgs_chr1_V0.5_c1_env10.tgz"

outfolder <- "./parsed_sim_data"
keep_inds <- seq(1,320,by=2) # keep one from each population

min_maf = 0.05

#----------------------------------------------------------
# Helper: unpack / clean up a tarball
#----------------------------------------------------------
focal_QTN <- function(map,GTs, qtn_col = "true_QTN") {


  #ch = "1"
  rbindlist(lapply(unique(map$Chr), function(ch) {

    # map rows on this chromosome, restricted to SNPs present in GDS
    map_chr <- map[Chr == ch]

    markers_chr <- map_chr$marker

    qtn_chr <- map_chr[get(qtn_col) == TRUE, marker]

    if (length(qtn_chr) == 0L) {
      return(data.table(
        Chr = ch,
        marker = markers_chr,
        max_LD_with_QTN = 0,
        focal_QTN = NA_character_,
        bp_to_focal_QTN = NA_real_
      ))
    }


    r2 <- cor(GTs[, qtn_chr], GTs[,markers_chr],use = "pairwise.complete.obs")^2

    max_ld <- apply(r2, 2, max, na.rm = TRUE)
    focal_idx <- max.col(t(r2), ties.method = "first")
    focal_qtn <- qtn_chr[focal_idx]

    pos_vec <- map_chr$Pos
    names(pos_vec) <- map_chr$marker

    data.table(
      Chr = ch,
      marker = markers_chr,
      max_LD_with_QTN = max_ld,
      focal_QTN = focal_qtn,
      bp_to_focal_QTN = abs(pos_vec[markers_chr] - pos_vec[focal_qtn])
    )

  }), fill = TRUE)
}

get_va = function(map,GTs,qtn_rows){
  qtn_rows <- map[,which(type=="QTN")]

  if (length(qtn_rows) == 0L) return(numeric(0))  ## sapply(integer(0), ...) would otherwise
  ## silently return list() instead of numeric(0)

  va <- sapply(qtn_rows, function(i){
    a  <- map$allelic_values[i]               # effect size α
    gt <- GTs[, i, drop=TRUE]                 # genotype 0/1/2
    p  <- mean(gt)/2                          # allele frequency
    2 * p * (1 - p) * a^2                     # Va formula
  })

}

#----------------------------------------------------------
# LD-prune the marker set used to estimate the GRM, at a threshold
# tied to each chromosome's own fitted LD-decay curve (rho = 0.5,
# i.e. the r^2 halfway between short-range max and background LD)
# rather than an arbitrary fixed window -- addresses the same
# "proximal contamination" concern standard LD-pruning-before-kinship
# guidance is based on, but decay-curve-aware.
#
# Reuses LD_decay$by_chr[[ch]]$el -- the SAME precomputed edge list
# used for compute_ld_w() -- instead of computing pairwise
# correlations from scratch. This covers every marker on the
# chromosome (compute_LD_decay() was run with max_SNPs_decay = Inf),
# so nothing is missed; a pair not present in `el` fell outside its
# own slide-SNP window and is treated as below threshold (r2 ~ 0),
# rather than being recomputed here.
#
# Still targets COMPLETE linkage, not single linkage, for the same
# reason as before (avoiding chained, only-loosely-related markers
# merging via a long path of gradual decay) -- but done in two cheap
# stages instead of one expensive dense hclust:
#   1) connected components on the SPARSE thresholded graph
#      (single-linkage groups; cheap, same machinery as
#      LD_igraph_components())
#   2) WITHIN each component (far smaller than the whole chromosome),
#      refine to true complete-linkage sub-clusters via hclust on a
#      small dense sub-matrix built from `el`'s already-available r2
#      values (missing pairs = 0, i.e. not clustered together)
#
# Representative-SNP selection within each final sub-cluster is
# deterministic (not random, for reproducibility): the marker most
# correlated with the sub-cluster's mean genotype -- needs raw GTs,
# but only for these small groups, never for a whole chromosome at
# once. This assumes consistent allele-coding direction in GTs
# (dosages coded relative to the same reference allele); if that
# doesn't hold, a genuinely central marker could show a spuriously
# low or negative correlation with the mean.
#----------------------------------------------------------
prune_for_grm <- function(GTs, map, LD_decay, rho = 0.5, cores = 1) {

  chr_levels <- unique(map$Chr)

  kept_by_chr <- mclapply(chr_levels, function(ch) {

    chr_map <- map[Chr == ch]
    chr_markers <- chr_map$marker
    if (length(chr_markers) < 2) return(chr_markers)

    ds <- LD_decay$decay_sum[Chr == ch]
    if (nrow(ds) != 1) {
      stop("Expected exactly one decay_sum row for chr '", ch, "', found ", nrow(ds),
           ". Check LD_decay$decay_sum$Chr labels.")
    }
    r2_th <- ld_from_rho(b = ds$b, c = ds$c, rho = rho)

    el <- LD_decay$by_chr[[ch]]$el
    if (is.null(el)) {
      stop("LD_decay$by_chr[['", ch, "']]$el is NULL -- compute_LD_decay() must be run with ",
           "keep_el = TRUE (or a valid el_data_folder) for prune_for_grm() to reuse its edges.")
    }
    if (is.character(el)) el <- fread(el, showProgress = FALSE)  ## el_data_folder mode stores a path
    #plot(LD_decay,type="chr")
    edges_th <- el[r2 >= r2_th]
    if (nrow(edges_th) == 0) return(chr_markers)  ## nothing exceeds threshold -- keep everyone

    ## Stage 1 (cheap): connected components from the sparse thresholded
    ## graph -- single-linkage groups, same machinery as LD_igraph_components()
    g <- graph_from_data_frame(edges_th[, .(SNP1, SNP2)], directed = FALSE, vertices = data.table(name = chr_markers))
    comp <- components(g)
    #table(chr_markers %in% names(unlist(comp$membership)))
    #length(unique(comp$membership))
    #table(table(comp$membership)==1)
    #length(comp$membership)

    #cl_id <- 2212
    per_component <- lapply(unique(comp$membership), function(cl_id) {
      # print(cl_id)
      cluster_markers <- names(comp$membership)[comp$membership == cl_id]
      if (length(cluster_markers) == 1) return(cluster_markers)

      ## Stage 2: within this (single-linkage) component, refine to TRUE
      ## complete-linkage sub-clusters -- a component can still be a
      ## "chain" where not every pair inside it clears r2_th
      sub_edges <- edges_th[SNP1 %in% cluster_markers & SNP2 %in% cluster_markers]
      m <- length(cluster_markers)
      R2_sub <- matrix(0, m, m, dimnames = list(cluster_markers, cluster_markers))
      diag(R2_sub) <- 1
      i1 <- match(sub_edges$SNP1, cluster_markers)
      i2 <- match(sub_edges$SNP2, cluster_markers)
      R2_sub[cbind(i1, i2)] <- sub_edges$r2
      R2_sub[cbind(i2, i1)] <- sub_edges$r2

      d      <- as.dist(1 - R2_sub)
      hc     <- hclust(d, method = "complete")
      sub_cl <- cutree(hc, h = 1 - r2_th)

      ## representative per true complete-linkage sub-cluster: marker
      ## most correlated with the sub-cluster's mean genotype
      # sc_id <- 14
      vapply(unique(sub_cl), function(sc_id) {
        #print(sc_id)
        sub_markers <- cluster_markers[sub_cl == sc_id]
        if (length(sub_markers) == 1) return(sub_markers)
        if (length(sub_markers) == 2) return(sub_markers[1])


        sub_gt  <- as.matrix(GTs[, sub_markers, drop = FALSE])
        mean_gt <- rowMeans(sub_gt, na.rm = TRUE)
        cors    <- cor(sub_gt, mean_gt, use = "pairwise.complete.obs")[, 1]
        sub_markers[which.max(cors)]
      }, character(1))
    })

    out <- unlist(per_component, use.names = FALSE)

  }, mc.cores = cores)

  unlist(kept_by_chr, use.names = FALSE)
}



unpack_sim <- function(file_gz,n_comp=0,tmp_dir) {

  if(!file.exists(tmp_dir)) dir.create(tmp_dir, recursive = TRUE)
  folder <- sub("\\.tar\\.gz$", "", file_gz)
  folder <- basename(folder)
  folder <- paste0(tmp_dir,"/",folder)

  if(!file.exists(folder)){
    system(paste0("tar -xvzf ", file_gz, " --strip-components=",n_comp, " -C ", tmp_dir, "/")) ## on a mac you need the strip components in order not to recreate the entire path again
  }

  return(folder)

}

extraxt_params <- function(base_name){
  params <- data.table(t(strsplit(base_name, "_", fixed = TRUE)[[1]]))[, 1:6]
  setnames(params, c("sim", "bgs","Chr", "V","c", "env"))
  params[,env:=gsub(".tgz","",env)]
}

#----------------------------------------------------------
# Derive a per-batch tmp folder from the file's OWN bgs/nobgs tag,
# so two terminal sessions filtering to different subsets (e.g.
# "_bgs_" vs "_nobgs_") never write to the same tmp folder at once.
# Reuses extraxt_params()'s existing filename parsing (it already
# isolates the bgs/nobgs token cleanly as `params$bgs`) rather than
# a separate ad hoc string search -- this is pure string manipulation
# on file_gz's name, no unpacking needed, so it's safe to call before
# unpack_sim(). If you want to split on a DIFFERENT distinguishing
# string instead (not bgs/nobgs), swap params$bgs for params$V,
# params$env, or a grepl()-based check on base_name directly.
#----------------------------------------------------------
get_tmp_dir <- function(file_gz, base_dir = "./tmp") {
  base_name <- basename(sub("\\.tgz$|\\.tar\\.gz$", "", file_gz))
  bgs_tag <- extraxt_params(base_name)$bgs
  paste0(base_dir, "_", bgs_tag)
}

cleanup_sim <- function(base_name,tmp_dir) {
  system(paste("rm -r", tmp_dir))
}

#----------------------------------------------------------
# Core parser
#----------------------------------------------------------

parse_raw_data <- function(file_gz,
                           outfolder,
                           min_maf  = 0.05,
                           side,
                           keep) {


  # ---------------------------- #
  # 1) Unpack and and prepare data files
  # ---------------------------- #

  tmp_dir <- get_tmp_dir(file_gz)  ## e.g. "./tmp_bgs" or "./tmp_nobgs" -- per-batch, no collision
  folder <- unpack_sim(file_gz,n_comp = 4,tmp_dir = tmp_dir) ## unpack data into the batch's own tmp folder
  base_name <- basename(folder)

  files <- list.files(tmp_dir, recursive = TRUE, full.names = TRUE)

  if(length(files)==0){
    folder <- unpack_sim(file_gz,n_comp = 0,tmp_dir = tmp_dir)
    files <- list.files(tmp_dir, recursive = TRUE, full.names = TRUE)
  }


  message("Working on ", base_name)

  params <- extraxt_params(base_name)
  ## Parse parameter info from folder name

  ## get files

  ## Check that simulation produced output
  if (!any(grepl("snp_geno", files))) {
    message("Simulation did not work\n\n")
    cleanup_sim(folder, tmp_dir)
    return(invisible(NULL))
  }


  # ---------------------------- #
  # 2) Sim summary & map
  # ---------------------------- #


  map_nemo <- fread(files[grep(".map", files, fixed = TRUE)])
  GTs <- fread(files[grep("snp_geno", files, fixed = TRUE)])


  nemo_map <- data.table(marker=map_nemo$trait.locus,do.call(rbind,strsplit(map_nemo$trait.locus,".",fixed=TRUE)))
  setnames(nemo_map, c("V1","V2"),c("type","idx"))
  #nemo_map[,table(type)]

  sample_info <- GTs[,.(pop,ID)]
  GTs <- as.matrix(GTs[1:.N,6:ncol(GTs),with=FALSE])

  #"./Nemo_sim/maps_500kb_with_allelic_values/chromosome_maps_500kb_rds/rec_map4.rds "
  # original recombinatin map
  map <- readRDS(paste0("./Nemo_sim/maps_500kb_with_allelic_values/chromosome_maps_500kb_rds/rec_map",gsub("chr","",params$Chr),".rds"))
  map[,indx := .I]
  #map[type=="QTN"]
  #map[,hist(allelic_values)]

  # parse nemo map


  nemo_map[,idx:=as.numeric(idx)+1] ## because nemo starts from 0
  #nemo_map[,table(type)]
  # indexes in Nemo map
  ntrl_idx   <- nemo_map[type=="ntrl",idx]
  quanti_idx <- nemo_map[type=="quant",idx]
  delet_idx <- nemo_map[type=="delet" ,idx]

  # indexes in original map
  indx_ntrl   <- map[type=="ntrl"][ntrl_idx,indx]
  indx_quanti <- map[type=="QTN"][quanti_idx,indx]
  indx_delet <- map[type=="delet"][delet_idx,indx]

  map[indx_ntrl,nemo_marker     := nemo_map[type=="ntrl",marker]]
  map[indx_quanti,nemo_marker   := nemo_map[type=="quant",marker]]
  map[indx_delet,nemo_marker   := nemo_map[type=="delet" ,marker]]
  #map[type=="QTN" & !is.na(nemo_marker)]
  #map[type=="QTN"]
  map <- map[nemo_marker %in% colnames(GTs)]
  GTs <- GTs[,map$nemo_marker]
  map[,Pos:=bp]
  map[,marker := paste(paste0("Chr",Chr),Pos,sep=":")]
  colnames(GTs) <- map$marker
  setorder(map,Chr,bp)
  map[,chr_type := ifelse(Chr==1,"QTN","ntrl")]
  map[,Chr:=paste0("Chr",Chr)]

  GTs <- GTs[keep_inds,map$marker]
  #map <- map[,.(Chr,Pos,marker,type)]

  ## ------------------------------------------------
  ## Add environmental values
  ## ------------------------------------------------

  x <- gsub("env","",params$env)
  env_raw <- scan(paste0("./Nemo_sim/env/env_",x,".txt"), what = character(), quiet = TRUE)
  env_raw <- strsplit(env_raw,c("}{"),fixed=TRUE)[[1]]
  env_vals <- as.numeric(gsub("}}","",gsub("{{","",env_raw,fixed = TRUE),fixed = TRUE))
  side=48
  env <- data.table(expand.grid(x = 1:side, y = side:1))
  env[, pop := 1:(side * side)]
  env[, env := env_vals]
  new_order <- match(sample_info$pop,env$pop)
  env_ind <- env[new_order]
  env_ind <- env_ind[keep_inds]
  env_ind[,indx:=.I]


  # ggplot(env_ind[keep_inds], aes(x,y,col=env)) +
  #   geom_point() +
  #   scale_color_viridis_c(option="turbo")

  remove <- map[,duplicated(marker),by=Chr][,which(V1)]

  if(length(remove)>0){
    map <- map[-remove,]
    GTs <- GTs[,-remove]
  }

  maf <- colSums(GTs)/nrow(GTs)/2
  map[,MAF:=ifelse(maf<0.5,maf,1-maf)]
  keep_snps <- map$MAF>min_maf
  GTs <- GTs[,keep_snps]
  map <- map[MAF>min_maf,]
  qtn_rows <- map[,which(type=="QTN")]
  message("## ------------------------------------------------")
  message("Data contains ",nrow(map)," SNPs")
  message("## ------------------------------------------------")
  #cor(env_ind$env,GTs[,qtn_rows])^2

  ## compute Va
  #i <- qtn_rows[1]

  ## always create Va/sum_Va/p_Va so downstream code never hits a
  ## missing-column error, even when this simulation has zero QTNs
  ## (e.g. lost during the MAF filter above, or none segregating at
  ## all -- expect this less once simulations are conditioned on
  ## >=1 QTN, but handle it gracefully in the meantime)
  map[, `:=`(Va = NA_real_, sum_Va = NA_real_, p_Va = NA_real_)]

  if (length(qtn_rows) == 0L) {
    message("No QTNs present in this simulation -- skipping Va/p_Va computation.")
  } else {
    ## assign
    map[qtn_rows, Va := get_va(map,GTs)]

    ## total Va per chromosome
    map[qtn_rows, sum_Va := sum(Va)]

    ## proportion
    map[qtn_rows, p_Va := Va / sum_Va]
    #map[qtn_rows]
  }

  map[type == "QTN" ,true_QTN := TRUE]

  true_qtns <- map[which(true_QTN == TRUE), marker]

  fq <- focal_QTN(map, GTs,qtn_col = "true_QTN")

  map[fq, on = "marker",
      `:=`(
        focal_QTN = i.focal_QTN,
        bp_to_focal_QTN = i.bp_to_focal_QTN,
        max_LD_with_QTN = i.max_LD_with_QTN
      )]

  ## collect data
  map[type == "QTN", `:=`(
    focal_QTN = marker,
    max_LD_with_QTN = 1,
    bp_to_focal_QTN = 0
  )]

  map[chr_type=="ntrl",bp_to_focal_QTN:=max(Pos)]
  map[chr_type=="ntrl",max_LD_with_QTN:=0]

  message("Running LFMM ")


  #library(LEA)
  write.lfmm(GTs, file.path(tmp_dir, "genotypes.lfmm"))
  write.env(env_ind$env, file.path(tmp_dir, "gradients.env"))
  #write.env(env_null, file.path(tmp_dir, "gradients.env"))
  project = NULL

  project = lfmm2(file.path(tmp_dir, "genotypes.lfmm"), file.path(tmp_dir, "gradients.env"), K=5) # K is user defined
  pv = suppressWarnings(lfmm2.test(project, file.path(tmp_dir, "genotypes.lfmm"), file.path(tmp_dir, "gradients.env"),genomic.control = TRUE,full = TRUE))
  map[,lfmm_p:=pv$pvalues] ## add to map
  map[,lfmm_F:=pv$fscores/pv$gif] ## add to map

  message("Running EMX ")
  #plot(map$indx,-log10(p.adjust(pv$pvalues,"fdr")),pch=ifelse(map$type=="QTN",3,20),cex=ifelse(map$type=="QTN",3,1), main="all SNPs | LFMM")

  #plot(pv$fscores,pch=ifelse(map$type=="QTN",3,20),cex=ifelse(map$type=="QTN",3,1),col=ifelse(map$type=="QTN","firebrick4","grey"))
  on.exit(unlink(gds))
  on.exit(rm(gds))
  on.exit(unlink(gds_path))

  gds_path = tempfile(fileext = ".gds")
  gds <- create_gds_from_geno(geno = GTs,map=map,gds_path)

  ## LD_decay needs to come BEFORE the GRM step now, since prune_for_grm()
  ## converts rho=0.5 to an actual r^2 threshold using each chromosome's
  ## own fitted decay curve (a/b/c) -- moved up from its old position
  ## after EMMAX, nothing else here depends on GRM/EMMAX having run first.
  maf_vec <- setNames(map$MAF, map$marker)  ## needed for min_maf_decay to actually filter (was a no-op without it)
  LD_decay <- compute_LD_decay(gds,
                               maf = maf_vec,
                               min_maf_decay = 0.1,
                               el_data_folder = NULL,
                               q = 0.95,
                               n_sub_bg = 5000,
                               n_win_decay = 5,
                               overlap = 0.5,
                               max_SNPs_decay = Inf,
                               prob_robust = 0.95,
                               max_pairs = 5000,
                               ld_method = "r",
                               n_strata = 20,
                               keep_el = TRUE,
                               slide = 1000,
                               rho_targets = c(0.99),
                               cores = 1)
  # plot(LD_decay)
  #LD_decay$by_chr$Chr1
  # plot(LD_decay,type="chr",chr="Chr1")
  #LD_decay$decay_sum
  a <- LD_decay$by_chr$Chr1$decay_sum$a
  b <- LD_decay$by_chr$Chr1$decay_sum$b
  c <- LD_decay$by_chr$Chr1$decay_sum$c
  map[chr_type=="QTN", rho_d := a * bp_to_focal_QTN / (a * bp_to_focal_QTN + 1)]
  map[chr_type=="QTN", ld_rel := (max_LD_with_QTN - b) / (c - b)]
  map[chr_type=="QTN", ld_rel := pmin(pmax(ld_rel, 0), 1)]

  ## LD-prune the marker set used to estimate the GRM (rho=0.5 relative
  ## to each chromosome's own decay curve), keeping EMMAX itself testing
  ## the FULL marker set as before -- only the GRM's underlying markers
  ## are restricted, via snpgdsGRM(snp.id=...) on the SAME gds object.
  pruned_markers <- prune_for_grm(GTs, map, LD_decay, rho = 0.5, cores = 1)

  message("GRM pruning: ", length(pruned_markers), " / ", nrow(map), " markers kept")

  GRM <- snpgdsGRM(gds,snp.id = pruned_markers,method = "GCTA",verbose = FALSE,autosome.only = FALSE)$grm
  #GRM_full <- snpgdsGRM(gds,method = "GCTA",verbose = FALSE,autosome.only = FALSE)$grm


  ## EMMAX does not expect a file in 012 format so the maximum likelihood genotypes can be used (without rounding)
  emx <- emmax(env_ind$env,GTs,K = GRM)
  map[,emx_p:=emx$pval] ## add to map

  map[,emx_F:=emx$F] ## add to map

  ## genomic control if gif>1.1
  emx_gif = map[,median(emx_F)/qf(0.5,1,nrow(GTs)-2,lower.tail = FALSE)] ## inflation factor

  map[,emx_F:=if(emx_gif>1.1) emx_F/emx_gif else emx_F,]  ## genomic control
  if(emx_gif>1.1) map[,emx_p := pf(emx_F,df1=1,df2=nrow(GTs)-2,lower.tail=FALSE)]

  # ld_ws <- compute_ld_w(rho=c(seq(0.05,0.95,by=0.05),0.99),LD_decay)
  # dt <- data.table(q_ldw=ecdf(ld_ws[,"0.95"])(ld_ws[,"0.95"]),val=-log10(p.adjust(map$lfmm_p,"fdr")),max_LD_with_QTN=map$max_LD_with_QTN,Chr=map$Chr)
  # ggplot(dt[], aes(q_ldw,val,col=max_LD_with_QTN))+
  #   geom_point() +
  #   scale_color_viridis_c(option="turbo") +
  #   facet_grid(Chr~.)

  LD_decay$by_chr$`Chr1`$el <-NULL
  LD_decay$by_chr$`Chr2`$el <-NULL
  out_file <- paste0(file.path(outfolder,paste(params, collapse = "_")), ".rds")
  saveRDS(
    list(
      GTs      = GTs,
      map      = map,
      env      = env_ind,
      LD_decay = LD_decay,
      ld_ws = ld_ws
    ),
    out_file
  )

  message("removing temporary data\n\n")
  cleanup_sim(base_name, tmp_dir)

  invisible(NULL)

}

f <- file_gz <- gz_files[7]
lapply(gz_files, function(f) {

  tmp_dir_f <- get_tmp_dir(f)
  if(file.exists(tmp_dir_f)) system(paste("rm -r", tmp_dir_f))

  parse_raw_data(
    file_gz   = f,
    outfolder = "./parsed_sim_data",
    min_maf   = 0.05,
    side      = 48,
    keep      = keep_inds
  )
})

# -------------------------------
# for diagnostics only
# -------------------------------

if(FALSE){
  map[,indx:=.I]
  par(mfcol=c(2,1))
  rho = "0.9"
  map[, q_ldw := ecdf(ld_ws[,rho])(ld_ws[,rho])]
  map[, ld_w := ld_ws[,rho]]
  #map[,plot(ld_w)]
  map[, true_pos := rho_d <=0.99 & ld_rel > 0.25]
  map[chr_type=="ntrl",true_pos:=FALSE]
  map[, SNP_class := fifelse(true_pos, "SNPs near causal loci", "Neutral SNPs")]


  #map[,plot(MAF,ld_w)]
  plot(map$indx,-log10(p.adjust(map$lfmm_p,"fdr")),pch=ifelse(map$type=="QTN",3,20),cex=ifelse(map$type=="QTN",3,1),col=ifelse(map$true_pos,"firebrick4","grey"),
       main="all SNPs | LFMM")
  abline(h=1.3)
  qt = 0.8
  keep <-  map[,which(q_ldw>qt)]
  #map[keep & type=="QTN"]
  #hist(-log10(p.adjust(map$lfmm_p[keep],"fdr")))
  #hist(-log10(p.adjust(map$lfmm_p,"fdr"))[keep])

  #map[,plot(ld_ws,-log10(lfmm_p))]

  plot(map$indx[keep],-log10(p.adjust(map$lfmm_p[keep],"fdr")),
       pch=ifelse(map[keep]$type=="QTN",3,20),
       cex=ifelse(map[keep]$type=="QTN",3,1),col=ifelse(map[keep]$true_pos,"firebrick4","grey"),
       main="ld_w>ld_w[Q(95)] | LFMM",xlim=c(0,nrow(map)))
  abline(h=1.3)

  plot(map$indx,map$ld_w,pch=ifelse(map$type=="QTN",3,20),cex=ifelse(map$type=="QTN",3,1),col=ifelse(map$true_pos,"firebrick4","grey"))

  plot(map$indx,-log10(p.adjust(map$emx_p,"fdr")),pch=ifelse(map$type=="QTN",3,20),cex=ifelse(map$type=="QTN",3,1),col=ifelse(map$true_pos,"firebrick4","grey"),
       main="all SNPs | EMMAX")
  abline(h=1.3)

  plot(map$indx[keep],-log10(p.adjust(map$emx_p[keep],"fdr")),
       pch=ifelse(map[keep]$type=="QTN",3,20),
       cex=ifelse(map[keep]$type=="QTN",3,1),col=ifelse(map[keep]$true_pos,"firebrick4","grey"),
       main="ld_w>ld_w[Q(95)] | EMMAX",xlim=c(0,nrow(map)))
  abline(h=1.3)

  # map[, true_pos := rho_d <=0.99 & ld_rel > 0.25]
  # map[chr_type=="ntrl",true_pos:=FALSE]
  # map[, SNP_class := fifelse(true_pos, "SNPs near causal loci", "Neutral SNPs")]


  p1 <- ggplot(map,aes(q_ldw,lfmm_F,col=max_LD_with_QTN)) +
    geom_point(data=map[chr_type=="ntrl"],alpha=0.5,shape = 21, fill="grey40",size = 3.5, col = "black", stroke = 0.3) +
    #geom_point(data = map[!which(true_pos)],size = 1.5, col = "grey70", alpha = 0.35)+
    #geom_smooth(se=FALSE) +
    geom_point(data = map[which(true_pos)],
               aes(fill = max_LD_with_QTN),
               shape = 21, size = 3.5, col = "black", stroke = 0.3)+
    scale_shape_manual(values = c(20,3),name=NULL) +
    scale_fill_gradientn(
      colors = wes_palette("Zissou1", 100, type = "continuous"),name = expression("LD with QTN (" * italic(r)^2 * ")")
    ) +
    facet_grid(.~SNP_class)+
    scale_size_identity() +
    theme_bw(base_size = 22) +
    theme(strip.background = element_blank(),
          legend.background = element_blank(),
          legend.position = "inside",
          legend.position.inside = c(0.65,0.7))+
    ylab("LFMM association statistic") +
    xlab("Local LD rank quantile")
  p1

}

