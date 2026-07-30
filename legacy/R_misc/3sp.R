library(data.table)
library(parallel)
library(SNPRelate)
library(LDscnR)
source("./R/emmax.R")
load("../LD-scaling-genome-scans/empirical_data/3sp/3sp_data.RData") ## contains SNP_res_3sp, GTs_3sp and pheno_3sp

##filter by maf
keep <- map_3sp$maf>0.1
GTs_3sp <- GTs_3sp[,map_3sp$maf>0.1]
map_3sp <- map_3sp[maf>0.1]

map_3sp[,emx_F:=readRDS("../LD-scaling-genome-scans/empirical_data/3sp/emx_3sp.rds")$F] ## add to map
emx_gif = map_3sp[,median(emx_F)/qf(0.5,1,115,lower.tail = FALSE)] ## inflation factor
map_3sp[,emx_F_GC:=emx_F/emx_gif]  ## genomic control
map_3sp[,emx_p_GC:=pf(emx_F_GC,1,115,lower.tail = FALSE)] ## p-value
map_3sp[,emx_q:=p.adjust(emx_p_GC,"fdr")] ## fdr correction

map_3sp[,lfmm_F:=readRDS("../LD-scaling-genome-scans/empirical_data/3sp/lfmm_F.rds")[keep]]
map_3sp[,lfmm_P:=pf(lfmm_F,1,115,lower.tail = FALSE)]
map_3sp[,lfmm_q:=p.adjust(lfmm_P,"fdr")]



permute_ecotype_among_pops <- function(pop, pop_ecotype, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  pop_ecotype <- unique(pop_ecotype)

  stopifnot(all(c("pop", "ecotype") %in% names(pop_ecotype)))

  # One ecotype per population
  if (anyDuplicated(pop_ecotype$pop)) {
    stop("pop_ecotype must contain one ecotype value per population")
  }

  # Permute ecotype labels among populations
  pop_ecotype$ecotype_perm <- sample(pop_ecotype$ecotype)

  # Return permuted ecotype for each individual
  yperm <- pop_ecotype$ecotype_perm[match(pop, pop_ecotype$pop)]

  yperm
}

pop_ecotype <- unique(data.frame(
  pop = pheno_3sp$pop_ID,
  ecotype = as.numeric(as.factor(pheno_3sp$ecotype))
))

yperm <- permute_ecotype_among_pops(
  pop = pheno_3sp$pop_ID,
  pop_ecotype = pop_ecotype
)

dim(GTs_3sp)
gds_3sp <- create_gds_from_geno(geno = GTs_3sp, map=map_3sp,"gds_3sp.gds")

if(!file.exists("./3sp_data/ld_decay_3sp_MAF01.rds")){
  t1 <- Sys.time()
  ld_decay_3sp <- compute_LD_decay(
    gds=gds_3sp,
    el_data_folder = "./el_3sp/",
    n_win_decay = 20,
    max_SNPs_decay = Inf,
    keep_el = TRUE,
    slide=300,
    cores = 8,
    ld_method = "corr"
  )
  t2 <- Sys.time()
  saveRDS(ld_decay_3sp,"./3sp_data/ld_decay_3sp_MAF01.rds")
}else{
  ld_decay_3sp <- readRDS("./3sp_data/ld_decay_3sp_MAF01.rds")
}


library(MASS)
GRM <- snpgdsGRM(gds_3sp)$grm

ld_ws_3sp <- precalculate_ld_w(c(seq(0.75,0.95,by=0.05),0.99),ld_decay_3sp)

yperm <- permute_ecotype_among_pops(
  pop = pheno_3sp$pop_ID,
  pop_ecotype = pop_ecotype
)

emx_perm <- emmax(yperm,GTs_3sp,K = GRM) ## function lives in ./R/emmax.R
map_3sp[,emx_p_perm:=emx_perm$pval] ## add to map
map_3sp[,emx_p_F:=emx_perm$F] ## add to map

map_3sp[,ld_w := ld_ws_3sp[,"0.99"]]

map_3sp[,plot(-log10(p.adjust(emx_p_perm,"fdr")))]



map_3sp[,ld_w := ld_ws[,"0.99"]]
map_3sp[, indx := .I]

par(mfcol=c(2,1))
th = 0.2
map_3sp[,plot(indx,-log10(p.adjust(pf(lfmm_F,1,115,lower.tail = FALSE),"fdr")),cex=0.5,ylab="-log10(q)",col="grey40",main="LFMM | ld_w>0.2")]
map_3sp[ld_w>th,points(indx,-log10(p.adjust(pf(lfmm_F,1,115,lower.tail = FALSE),"fdr")),col="black",cex=1,pch=20)]
abline(h=1.3,col="steelblue")

map_3sp[,plot(indx,-log10(p.adjust(pf(emx_F,1,115,lower.tail = FALSE),"fdr")),cex=0.5,ylab="-log10(q)",col="grey40",main="EMMAX | ld_w>0.2",ylim=c(0,3))]
map_3sp[ld_w>th,points(indx,-log10(p.adjust(pf(emx_F,1,115,lower.tail = FALSE),"fdr")),col="black",cex=1,pch=20)]
abline(h=1.3,col="steelblue")

plot(ld_ws[,6],cex=0.5)

names(ld_decay_3sp$by_chr)

map_3sp[,ld_w := ld_ws[,"0.75"]]

par(mfcol=c(2,1))
qt = 0.9

map_joint_3sp[,plot(indx,-log10(p.adjust(pf(lfmm_F,1,115,lower.tail = FALSE),"fdr")),cex=0.5,ylab="-log10(q)",col="grey40",ylim=c(0,5),main="LFMM | ld_w>0.2")]
map_joint_3sp[ld_w>quantile(ld_w,qt,na.rm=TRUE),points(indx,-log10(p.adjust(pf(lfmm_F,1,115,lower.tail = FALSE),"fdr")),col="black",cex=1,pch=20)]
abline(h=1.3,col="steelblue")

plot(seq(0.8,0.999,by=0.001),sapply(seq(0.8,0.999,by=0.001),function(qt){
  map_joint_3sp[ld_w>quantile(ld_w,qt,na.rm=TRUE),cor(ld_w,lfmm_F)^2]
}),type="l")

which.max(sapply(seq(0.9,0.999,by=0.001),function(qt){
  map_joint_3sp[ld_w>quantile(ld_w,qt,na.rm=TRUE),cor(ld_w,lfmm_F)^2]
}))
seq(0.9,0.999,by=0.001)[89]

map_joint_3sp[,plot(ld_w,lfmm_F)]
qt = 0.9
abline(v=map_joint_3sp[,quantile(ld_w,qt,na.rm=TRUE)],col="red")
#map_3sp[,plot(ld_w)]
map_3sp[,quantile(ld_w,0.8,na.rm=TRUE)]
map_sub <- map_3sp[ld_w>quantile(ld_w,0.9,na.rm=TRUE)]
peaks <- map_sub[, .(Chr, marker, Pos)]
setorder(peaks, Chr, Pos)

max_gap = 1e6
peaks[, peak_cluster :=
        cumsum(is.na(Pos - shift(Pos)) |
                 (Pos - shift(Pos)) > max_gap),
      by = Chr]

peaks[, phy_uid := paste0(Chr, "_phy", peak_cluster)]

phy_cls <- split(peaks$marker, peaks$phy_uid)


LD_igraph_components <- function(GTs, map, markers, r2_th = 0.8, bp_th = Inf,full=FALSE) {
  gts <- as.matrix(GTs[,markers])
  storage.mode(gts) <- "double"

  stopifnot(all(colnames(gts) %in% map$marker))

  # order map to match genotype columns
  snp_map <- map[match(colnames(gts), marker)]

  # fast LD matrix

  R2 <- cor(gts, use = "pairwise.complete.obs")^2
  R2[is.na(R2)] <- 0

  # candidate SNP pairs above LD threshold
  idx <- which(R2 >= r2_th, arr.ind = TRUE)
  idx <- idx[idx[, 1] < idx[, 2], , drop = FALSE]

  edges <- data.table(
    from = colnames(gts)[idx[, 1]],
    to   = colnames(gts)[idx[, 2]],
    r2   = R2[idx]
  )

  # optional physical-distance filter
  if (is.finite(bp_th)) {
    pos <- snp_map$Pos
    chr <- snp_map$Chr

    edges[, `:=`(
      Chr_from = chr[match(from, colnames(gts))],
      Chr_to   = chr[match(to, colnames(gts))],
      dist_bp  = abs(pos[match(from, colnames(gts))] -
                       pos[match(to, colnames(gts))])
    )]

    edges <- edges[Chr_from == Chr_to & dist_bp <= bp_th]
  }

  g <- graph_from_data_frame(
    edges[, .(from, to)],
    directed = FALSE,
    vertices = data.table(name = colnames(gts))
  )

  comp <- components(g)

  clusters <- data.table(
    marker = names(comp$membership),
    CL_id = as.integer(comp$membership)
  )

  clusters[, n_loci := comp$csize[CL_id]]

  if(full){
    list(
      graph = g,
      edges = edges,
      clusters = clusters,
      components = comp
    )
  }else{
    clusters
  }


}

#colnames(GTs_3sp) <- map_3sp$marker
sl_cls <- rbindlist(parallel::mclapply(seq_along(phy_cls), function(cl) {
  print(cl)
  markers_cl <- phy_cls[[cl]]

  if (length(markers_cl) > 1) {

    out <- LD_igraph_components(
      GTs_3sp,
      map_3sp,
      markers = markers_cl,
      r2_th = 0.6,
      bp_th = 1e6,
      full = FALSE
    )

  } else {

    out <- data.table(
      marker = markers_cl,
      CL_id = 1L,
      n_loci = 1L
    )
  }

  out[, phy_cl := cl]

  # globally unique cluster name
  out[, CL_uid := paste0("phy", phy_cl, "_ld", CL_id)]

  out

}, mc.cores = 1))


map_joint_3sp <- map_3sp[sl_cls[n_loci>=1,], on = "marker",
                         `:=`(
                           CL_uid = i.CL_uid,
                           LD_component = i.CL_id,
                           n_loci = i.n_loci
                         )]


par(mfcol=c(2,1))
th = 0.2
map_joint_3sp[,plot(indx,-log10(p.adjust(pf(lfmm_F,1,115,lower.tail = FALSE),"fdr")),cex=0.5,ylab="-log10(q)",col="grey40",ylim=c(0,5),main="LFMM | ld_w>0.2")]
map_joint_3sp[n_loci>=1,points(indx,-log10(p.adjust(pf(lfmm_F,1,115,lower.tail = FALSE),"fdr")),col="black",cex=1,pch=20)]
abline(h=1.3,col="steelblue")

map_joint_3sp[,plot(indx,-log10(p.adjust(pf(emx_F,1,115,lower.tail = FALSE),"fdr")),cex=0.5,ylab="-log10(q)",col="grey40",main="EMMAX | ld_w>0.2",ylim=c(0,5))]
map_joint_3sp[n_loci>=1,points(indx,-log10(p.adjust(pf(emx_F,1,115,lower.tail = FALSE),"fdr")),col="black",cex=1,pch=20)]
abline(h=1.3,col="steelblue")

map_filt <- map_joint_3sp[]
map_filt[,log_p := -log10(p.adjust(pf(emx_F,1,115,lower.tail = FALSE),"fdr"))]
map_filt[,hist(log_p)]

map_filt[Chr=="Chr1" & log_p>1.3,table(CL_uid) ]

ggplot(map_filt[log_p>1.3],aes(Pos,log_p,col=CL_uid)) +
  geom_point() +
  facet_wrap(Chr~.,scale="free_x") +
  scale_color_manual(values=rep(col_vector,10)) +
  geom_hline(yintercept = 1.3)




GRM <- snpgdsGRM(gds_3sp)$grm

eco_bin  <- as.numeric(as.factor(pheno_3sp$ecotype))

# one row per population
pheno_3sp
pop_tab <- unique(pheno_3sp[, .(pop=pop_ID, region=pop_locality, ecotype)])

# permute ecotype among populations within region/drainage
pop_tab[, ecotype_perm := sample(ecotype), by = pop_locality]

# assign back to individuals
pheno_3sp[pop_tab, on = "pop_ID", ecotype_perm := i.ecotype_perm]
eco_perm  <- as.numeric(as.factor(pheno_3sp$ecotype_perm))
## EMMAX does not expect a file in 012 format so the maximum likelihood genotypes can be used (without rounding)
emx_perm <- emmax(eco_perm,GTs_3sp,K = GRM) ## function lives in ./R/emmax.R
map_3sp[,emx_F_perm:=emx_perm$F] ## add to map


th = 0.2
map_joint_3sp[,plot(indx,-log10(p.adjust(pf(emx_F_perm,1,115,lower.tail = FALSE),"fdr")),cex=0.5,ylab="-log10(q)",col="grey40",ylim=c(0,5),main="LFMM | ld_w>0.2")]
map_joint_3sp[n_loci>=1,points(indx,-log10(p.adjust(pf(emx_F_perm,1,115,lower.tail = FALSE),"fdr")),col="black",cex=1,pch=20)]
abline(h=1.3,col="steelblue")


ld_ws <- precalculate_ld_w(c(seq(0.75,0.95,by=0.05),0.99),ld_decay_3sp)


plot(ld_ws[,6],cex=0.5)

names(ld_decay_3sp$by_chr)

map_3sp[,ld_w := ld_ws[,"0.99"]]
map_3sp[,plot(ld_w)]
#map_3sp[,abline(h=quantile(ld_w,0.9,na.rm=TRUE),col="salmon")]
map_sub <- map_3sp[ld_w>quantile(ld_w,0.9,na.rm=TRUE)]

peaks <- map_sub[, .(Chr, marker, Pos)]
setorder(peaks, Chr, Pos)

max_gap = 1e6
peaks[, peak_cluster :=
        cumsum(is.na(Pos - shift(Pos)) |
                 (Pos - shift(Pos)) > max_gap),
      by = Chr]

peaks[, phy_uid := paste0(Chr, "_phy", peak_cluster)]

phy_cls <- split(peaks$marker, peaks$phy_uid)


LD_igraph_components <- function(GTs, map, markers, r2_th = 0.8, bp_th = Inf,full=FALSE) {
  gts <- as.matrix(GTs[,markers])
  storage.mode(gts) <- "double"

  stopifnot(all(colnames(gts) %in% map$marker))

  # order map to match genotype columns
  snp_map <- map[match(colnames(gts), marker)]

  # fast LD matrix

  R2 <- cor(gts, use = "pairwise.complete.obs")^2
  R2[is.na(R2)] <- 0

  # candidate SNP pairs above LD threshold
  idx <- which(R2 >= r2_th, arr.ind = TRUE)
  idx <- idx[idx[, 1] < idx[, 2], , drop = FALSE]

  edges <- data.table(
    from = colnames(gts)[idx[, 1]],
    to   = colnames(gts)[idx[, 2]],
    r2   = R2[idx]
  )

  # optional physical-distance filter
  if (is.finite(bp_th)) {
    pos <- snp_map$Pos
    chr <- snp_map$Chr

    edges[, `:=`(
      Chr_from = chr[match(from, colnames(gts))],
      Chr_to   = chr[match(to, colnames(gts))],
      dist_bp  = abs(pos[match(from, colnames(gts))] -
                       pos[match(to, colnames(gts))])
    )]

    edges <- edges[Chr_from == Chr_to & dist_bp <= bp_th]
  }

  g <- graph_from_data_frame(
    edges[, .(from, to)],
    directed = FALSE,
    vertices = data.table(name = colnames(gts))
  )

  comp <- components(g)

  clusters <- data.table(
    marker = names(comp$membership),
    CL_id = as.integer(comp$membership)
  )

  clusters[, n_loci := comp$csize[CL_id]]

  if(full){
    list(
      graph = g,
      edges = edges,
      clusters = clusters,
      components = comp
    )
  }else{
    clusters
  }


}

colnames(GTs_3sp) <- map_3sp$marker
sl_cls <- rbindlist(parallel::mclapply(seq_along(phy_cls), function(cl) {
  print(cl)
  markers_cl <- phy_cls[[cl]]

  if (length(markers_cl) > 1) {

    out <- LD_igraph_components(
      GTs_3sp,
      map_3sp,
      markers = markers_cl,
      r2_th = 0.8,
      bp_th = 1e6,
      full = FALSE
    )

  } else {

    out <- data.table(
      marker = markers_cl,
      CL_id = 1L,
      n_loci = 1L
    )
  }

  out[, phy_cl := cl]

  # globally unique cluster name
  out[, CL_uid := paste0("phy", phy_cl, "_ld", CL_id)]

  out

}, mc.cores = 1))


map_joint_3sp <- map_3sp[sl_cls[n_loci>=1,], on = "marker",
                         `:=`(
                           CL_uid = i.CL_uid,
                           LD_component = i.CL_id,
                           n_loci = i.n_loci
                         )]

peak_summary <- map_3sp[!is.na(CL_uid), .(
  start = min(Pos),
  end = max(Pos),
  width = max(Pos) - min(Pos),
  n_snps = .N,
  max_ld = max(ld_w,na.rm = TRUE),
  min_ld = min(ld_w,na.rm = TRUE),
  median_ld = median(ld_w,na.rm = TRUE)
), by = .(Chr, CL_uid)]

peak_summary[,table(n_snps>=10)]

map_joint_3sp[,CL_id := CL_uid]
eMLGs <- make_eMLGs(GTs_3sp,map_joint_3sp[!is.na(CL_id) & n_loci >=10,],cor_th = 0.8,l_min = 1,ncores=4)

eMLGs$map_eMLG[,hist(r2_eMLG)]


emx_eMLG <- emmax(eco_bin,eMLGs$eMLG,K = GRM) ## function lives in ./R/emmax.R
plot(-log10(emx_eMLG$pval))

map_joint_3sp[,plot(ld_w,col=ifelse(!is.na(CL_id) & n_loci>=2,"steelblue","salmon"))]

map_joint_3sp[Chr=="Chr1" & n_loci>=10,table(CL_id)]

ggplot(map_joint_3sp[Chr=="Chr1" & n_loci>=10,], aes(Pos,emx_F,col=CL_id)) +
  geom_point()

eMLGs$map_eMLG[CL_id=="phy2_ld1319"]

plot(round(eMLGs$eMLG[,"phy2_ld1319"]),eMLGs$eMLG[,"phy2_ld1319"])

cor.test(eco_bin,eMLGs$eMLG[,"phy2_ld1319"])

emx_eMLG$pval["phy2_ld1319"]

peak_summary[,hist(n_snps)]
map_joint_3sp[!duplicCL_id,hist(unique(CL_id))]
pruned_SNPs <- map_joint_3sp[,sample(marker,1),by=CL_id]$V1
GRM_pruned <- snpgdsGRM(gds_3sp,snp.id = pruned_SNPs)$grm

map_joint_3sp[!is.na(CL_id),length(unique(CL_id))]

eMLGs <- make_eMLGs(GTs_3sp,map_joint_3sp[!is.na(CL_id) & n_loci >=1],cor_th = 0.8,l_min = 1,ncores=4)

eMLGs$map_eMLG[,hist(r2_eMLG)]

yperm <- permute_ecotype_among_pops(
  pop = pheno_3sp$pop_ID,
  pop_ecotype = data.frame(
    pop = pheno_3sp$pop_ID,
    ecotype = as.numeric(as.factor(pheno_3sp$ecotype))
  )
)

emx_eMLG <- emmax(eco_bin,eMLGs$eMLG,K = GRM) ## function lives in ./R/emmax.R
plot(-log10(map_3sp$emx))

#emx_eMLG <- emmax(eco_bin,eMLGs$eMLG,K = GRM_pruned) ## function lives in ./R/emmax.R
plot(-log10(emx_eMLG$pval))

qqplot(-log10(ppoints(p_lm)),-log10(emx_eMLG$pval))
abline(0,1)


p_lm <- apply(eMLGs$eMLG,2,function(x) cor.test(x,eco_bin)$p.value)

qqplot(-log10(ppoints(p_lm)),-log10(p_lm))
abline(0,1)

map_joint_3sp[, indx := .I]

par(mfcol=c(2,1))
th = 0.3
map_joint_3sp[,plot(indx,-log10(p.adjust(pf(lfmm_F,1,115,lower.tail = FALSE),"fdr")),cex=0.5,ylab="-log10(q)",col="grey40",main="LFMM")]
map_joint_3sp[ld_w>th,points(indx,-log10(p.adjust(pf(lfmm_F,1,115,lower.tail = FALSE),"fdr")),col="black",cex=1,pch=20)]
abline(h=1.3,col="steelblue")

th = 0.3
map_joint_3sp[,plot(indx,-log10(p.adjust(pf(emx_F,1,115,lower.tail = FALSE),"fdr")),cex=0.5,ylab="-log10(q)",col="grey40",main="EMMAX")]
map_joint_3sp[ld_w>th,points(indx,-log10(p.adjust(pf(emx_F,1,115,lower.tail = FALSE),"fdr")),col="black",cex=1,pch=20)]
abline(h=1.3,col="steelblue")



cluster_tests <- map_joint_3sp[!is.na(CL_id) & ld_w>th, {
  z <- .SD[which.max(ld_w)]  # independent representative choice
  .(
    marker = z$marker,
    Pos = z$Pos,
    p_cluster = z$lfmm_P,
    ld_w = z$ld_w,
    n_loci = .N
  )
}, by = .(Chr, CL_id)]

cluster_tests[, q_cluster := p.adjust(p_cluster, "fdr")]




rep_snps <- map_joint_3sp[ld_w>th, .SD[which.max(ld_w)], by = CL_id]

rep_snps[, q_rep := p.adjust(lfmm_p, "fdr")]

map_joint_3sp[ld_w>0.1,plot(-log10(p.adjust(pf(lfmm_F,1,150,lower.tail = FALSE),"fdr")))]



tmp <- map_joint_3sp[,.(max_F = max(lfmm_F,na.rm=TRUE),mean_F = mean(lfmm_F,na.rm=TRUE)),by=CL_id]
#tmp[,plot(max_F,mean_F)]
#plot(-log10(p.adjust(pf(tmp$max_F,1,150,lower.tail = FALSE),"fdr")))
qqplot(-log10(ppoints(tmp$CL_id)),-log10(pf(tmp$max_F,1,150,lower.tail = FALSE)))
abline(0,1)


map_3sp[,emx_F:=emx$F] ## add to map



draws_p3_emx <- ld_rho_draws(gds_3sp,
                             ld_decay  = ld_decay_3sp,
                             F_vals     = map_3sp[,.(emx_F)], ## the F values are pre-calculated and must be provided, these will be used for getting F_prime.
                             q_vals     = NULL, ## q-values from original analyses can also be used but not necessary
                             n_draws    = 100,
                             stat_type  = "q",
                             rho        = NULL, ## grid of rho values used for ld_w
                             ld_ws      = ld_ws,
                             rho_d_lim  = list(min=0.75,max=0.95),
                             rho_ld_lim = list(min=0.75,max=0.95),
                             alpha_lim  = list(min=1.3,max=4), ## lowest alpha is 1/10^0.6=0.25; this range is defined in the -log10 scale!
                             lmin_lim   = list(min=1,max=20),
                             cores      = 8,
                             mode       = "joint"
)

map_3sp_manh <- add_consistency_to_map(map_3sp, consistency_obj = consistency_score(draws_p3_emx$draws))
map_3sp_manh <- add_ORs(gds_3sp, ld_decay_3sp, map_3sp_manh, stat="Joint_C", sign_th=0.2,sign_if="greater",mode="joint",rho_d=0.99, rho_ld=0.99,l_min=1)


layout <- prep_manhattan(
  map_3sp_manh[, .(
    Pos,
    Chr,
    marker,
    OR_id,
    emx_F,
    Joint_C
  )]
)

plot_manhattan_gg(
  layout,
  y_vars = c("Joint_C","emx_F"),
  y_labels = c("C","EMMAX"),
  titles = c("B Consistency score","C EMMAX"),
  thresholds = c(0.2,NA),
  col_var = "OR_id",
  ncol=1
)



keep <- map_3sp_manh[,which(Joint_C>0)]

map_3sp[,emx_F_perm:=emx_F]
x <- 1
region_perm_full <- list()

for(x in (length(region_perm_full)+1):100){
  print(x)
  yperm <- permute_ecotype_among_pops(
    pop = pheno_3sp$pop_ID,
    pop_ecotype = data.frame(
      pop = pheno_3sp$pop_ID,
      ecotype = as.numeric(as.factor(pheno_3sp$ecotype))
    )
  )

  emx_3sp_perm <- emmax(Y=yperm, X=GTs_3sp[,],K=GRM)
  map_3sp[,emx_F_perm := emx_3sp_perm$F]
  #map_3sp[,plot(emx_F_perm)]
  #ld_scan

  draws_p3_perm <- ld_rho_draws(gds_3sp,
                                ld_decay  = ld_decay_3sp,
                                F_vals     = map_3sp[,.(sample(emx_F_perm))], ## the F values are pre-calculated and must be provided, these will be used for getting F_prime.
                                q_vals     = NULL, ## q-values from original analyses can also be used but not necessary
                                n_draws    = 100,
                                stat_type  = "q",
                                rho        = NULL, ## grid of rho values used for ld_w
                                ld_ws      = ld_ws[,],
                                rho_d_lim  = list(min=0.75,max=0.95),
                                rho_ld_lim = list(min=0.75,max=0.95),
                                alpha_lim  = list(min=1.3,max=4), ## lowest alpha is 1/10^0.6=0.25; this range is defined in the -log10 scale!
                                lmin_lim   = list(min=1,max=20),
                                cores      = 8,
                                mode       = "joint"
  )

  draws_p3_perm$draws

  if(any(draws_p3_perm$draws$OR_size>0)){
    map_3sp_manh_perm <- add_consistency_to_map(map_3sp[,], consistency_obj = consistency_score(draws_p3_perm$draws))
    region_perm_full[[x]] <- map_3sp_manh_perm[,Joint_C]
  }
}

map_3sp_manh_perm <- add_consistency_to_map(map_3sp[,], consistency_obj = consistency_score(draws_p3_perm$draws))
map_3sp_manh_perm[is.na(Joint_C),Joint_C:=0]
map_3sp_manh_perm[,plot(ld_ws[,5],sample(emx_F_perm))]

region_perm
lengths(region_perm)
dt <- cbind(map_3sp_manh[,Joint_C],do.call(cbind,region_perm_full))

plot(dt[,4])
dt[is.na(dt)] <- 0
image(cor(dt[,1:10],use = "pair")^2)

tmp <- cor(dt[,1:10],use = "pair")^2
hist(tmp)

par(mfcol=c(2,1))
tmp <- dt[,3]-rowMeans(dt[,-(1:3)])
plot(tmp,main="Reference is permuted")

tmp <- dt[,1]-rowMeans(dt[,-1])
plot(tmp,main="Reference is observed")


qqplot(dt[,1],dt[,-1])

tmp[tmp<0] <- 0
plot(tmp)


perm_mean <- rowMeans(dt[,-(1:2)])

plot(ld_ws[,"0.95"],dt[,1])

plot(ld_ws[,"0.95"],map_3sp[,emx_F])
plot(ld_ws[,"0.95"],map_3sp[,emx_F_perm])

perm_sd <- apply(dt[,-(1:2)], 1, sd)

zC <- (dt[,2] - perm_mean) / perm_sd
plot(zC)


perm <- as.matrix(dt[,-1])
obs  <- dt[,1]

perm_mean <- rowMeans(perm)
obs_resid <- obs - perm_mean

apply(ld_ws,2,function(x) cor(x,perm_mean,use="pair")^2)

perm_max <- apply(perm,1,quantile,0.95)
perm_q_95 <- apply(perm,1,quantile,0.95)

png("tst.png",height = 12,width = 8,units = "in",res = 150)
par(mfcol=c(5,1))
plot(ld_ws[,"0.95"],cex=0.25)
plot(perm_q_95,cex=0.25)
plot(obs,cex=0.25)
plot(obs-perm_q_95,cex=0.25)
plot(dt[,6]-perm_q_95,cex=0.25)
abline(h=0.2)
plot(perm_q_95-ld_ws[,"0.95"],cex=0.25)
dev.off()
map_3sp[,residual_C := obs-perm_q_95]


png("tst.png",height = 8,width = 8,units = "in",res = 300)
plot_manhattan(map=map_3sp,
               gds_3sp,
               ld_decay_3sp,
               draws_p3_perm,
               ## for defining ORs
               sign_th=0.2,
               mode="joint",
               sign_if ="greater",
               rho_d = 0.9,
               rho_ld = 0.9,
               ## outlier regiosn are defined based on the first in y_vars
               y_vars = c("residual_C","emx_F"),
               y_labels = c("Residual C","emx_F"),
               titles = c(" "),
               thresholds = c(0.2,NA),
               col_var = "OR_id")
dev.off()


qqplot(perm_max,obs)
par(mfcol=c(1,3))
plot(obs,ld_ws[,"0.95"])
plot(obs,ld_ws[,"0.95"])

perm_resid <- sapply(seq_len(ncol(perm)), function(j) {
  perm[,j] - rowMeans(perm[,-j, drop=FALSE])
})

plot(rowMeans(perm_resid))

p_emp <- rowMeans(perm_resid >= obs_resid)
plot(-log10(p_emp))
#plot(dt[,2])

pairs(dt)

ggplot(dt)


map_3sp[,ld_w_99:=compute_ld_w(ld_decay_3sp,rho = 0.99)]
map_3sp[,plot(ld_w_99)]
map_3sp[is.na(ld_w_99),ld_w_99:=0]
rho_ld = 0.5
ch = "Chr1"
pruned_SNPs <- rbindlist(lapply(names(ld_decay_3sp$by_chr), function(ch){
  mp <- copy(map_3sp[Chr==ch])
  a_chr <- ld_decay_3sp$decay_sum[Chr==ch,a_pred]
  b_chr <- ld_decay_3sp$decay_sum[Chr==ch,b]
  c_chr <- ld_decay_3sp$decay_sum[Chr==ch,c_pred]


  ld_th <- ld_from_rho(b_chr, c_chr, rho = rho_ld)


  chr_obj <- copy(ld_decay_3sp$by_chr[[ch]])
  chr_obj$el <- fread(chr_obj$el,showProgress = FALSE)

  ed <- chr_obj$el[r2>ld_th,.(SNP1=paste(ch,pos1,sep=":"),SNP2=paste(ch,pos2,sep=":"))]
  g     <- igraph::graph_from_data_frame(ed, directed = FALSE)
  comps <- igraph::components(g)

  ors_chr <- split(names(comps$membership), comps$membership)
  ors_chr <- ors_chr[vapply(ors_chr, length, integer(1)) >= 1]

  blocks <- data.table(block_id=rep(paste0("OR_",ch,1:length(ors_chr)),sapply(ors_chr,length)),marker=unlist(ors_chr))
  mp[,block_id:=blocks$block_id[match(marker,blocks$marker)]]
  mp[is.na(block_id),block_id:=marker]

  core_snps <- mp[,.(core_snp=marker[which.max(ld_w_99)],n_snps=.N,mean_ld_w=mean(ld_w_99)),by=block_id]

  return(core_snps)
}))

ids <- .read_gds_ids(gds_3sp)
table(ids$snp_id %in% pruned_SNPs$core_snp)

pruned_SNPs$core_snp[!(pruned_SNPs$core_snp %in% ids$snp_id)]
pruned_SNPs

GRM <- snpgdsGRM(gds_3sp,method = "GCTA",snp.id = pruned_SNPs$core_snp,verbose = FALSE,autosome.only = FALSE)$grm

## the binary phenotype
eco_bin  <- as.numeric(as.factor(pheno_3sp$ecotype))

## EMMAX does not expect a file in 012 format so the maximum likelihood genotypes can be used (without rounding)
emx <- emmax(eco_bin,GTs_3sp,K = GRM) ## function lives in ./R/emmax.R

map_3sp[,emx_F:=emx$F] ## add to map
emx_gif = map_3sp[,median(emx_F)/qf(0.5,1,115,lower.tail = FALSE)] ## inflation factor
map_3sp[,emx_F_GC:=emx_F/emx_gif]  ## genomic control
map_3sp[,emx_p_GC:=pf(emx_F_GC,1,115,lower.tail = FALSE)] ## p-value
map_3sp[,emx_q:=p.adjust(emx_p_GC,"fdr")] ## fdr correction

library(LEA)
phe <- as.numeric(as.factor(pheno_3sp$ecotype))

write.lfmm(GTs_3sp, "./tmp/genotypes.lfmm")
write.env(phe, "./tmp/gradients.env")
project = NULL

project = lfmm2("./tmp/genotypes.lfmm", "./tmp/gradients.env", K=7) # K is user defined
pv = lfmm2.test(project, "./tmp/genotypes.lfmm", "./tmp/gradients.env",genomic.control = TRUE,full = TRUE)

saveRDS(pv$f,"./empirical_data/3sp/lfmm_F.rds") # only F-value is needed


ld_ws2 <- precalculate_ld_w(pmin(seq(0.75,1,by=0.05),0.99),ld_decay)
#plot(ld_ws2[,1],ld_ws[,1])


draws <- ld_rho_draws(gds,
                      ld_decay  = ld_decay,
                      F_vals     = map[,.(lfmm_F,emx_F)], ## the F values are pre-calculated and must be provided, these will be used for getting F_prime.
                      q_vals     = NULL, ## q-values from original analyses can also be used but not necessary
                      n_draws    = 100,
                      stat_type  = "q",
                      rho        = NULL, ## grid of rho values used for ld_w
                      ld_ws      = ld_ws,
                      rho_d_lim  = list(min=0.9,max=0.99),
                      rho_ld_lim = list(min=0.9,max=0.99),
                      alpha_lim  = list(min=0.6,max=4), ## lowest alpha is 1/10^0.6=0.25; this range is defined in the -log10 scale!
                      lmin_lim   = list(min=1,max=10),
                      cores      = cores,
                      mode       = "per_method"
)

LD_int_3sp  <- unlist(lapply(ld_struct_3sp$by_chr,function(x) x$ld_int))
map_3sp[,ld_w:=LD_int_3sp]
map_3sp[,emx_F:=readRDS("../LD-scaling-genome-scans/empirical_data/3sp/emx_3sp.rds")$F] ## add to map
emx_gif = map_3sp[,median(emx_F)/qf(0.5,1,115,lower.tail = FALSE)] ## inflation factor
map_3sp[,emx_F_GC:=emx_F/emx_gif]  ## genomic control
map_3sp[,emx_p_GC:=pf(emx_F_GC,1,115,lower.tail = FALSE)] ## p-value
map_3sp[,emx_q:=p.adjust(emx_p_GC,"fdr")] ## fdr correction

map_3sp[,lfmm_F:=readRDS("../LD-scaling-genome-scans/empirical_data/3sp/lfmm_F.rds")[keep]]
map_3sp[,lfmm_P:=pf(lfmm_F,1,115,lower.tail = FALSE)]
map_3sp[,lfmm_q:=p.adjust(lfmm_P,"fdr")]

lfmm_F_prime_full <- .compute_Fprime(map_3sp$lfmm_F,map_3sp$ld_w,n_rep=10,n_inds = nrow(GTs_3sp),full=TRUE)

plot(-log10(lfmm_F_prime_full[[1]]$q_prime),pch=20)
abline(h=1.31,lty=2,cex=2,col="salmon")

t1 <- Sys.time()
draws_3sp_int <- ld_rho_draws(gds=gds_3sp,
                              ld_struct  = ld_struct_3sp,
                              F_vals     = map_3sp[,..F_cols],
                              q_vals     = map_3sp[,..q_cols],
                              n_rho_w    = 1,
                              n_draws    = 25*25,
                              ld_w_int   = map_3sp$ld_w,
                              stat_type  = "q",
                              rho_w_lim  = NULL,
                              rho_d_lim  = list(min=0.5,max=0.99),
                              rho_ld_lim = list(min=0.5,max=0.99),
                              alpha_lim  = list(min=1.31,max=4),
                              lmin_lim   = list(min=1,max=10),
                              cores      = 8,
                              mode       = "joint"

)
t2 <- Sys.time()

map_3sp2 <- add_consistency_to_map(map_3sp, consistency_obj = consistency_score(draws_3sp_int$draws))

map_3sp2[,plot(Joint_C)]
map_3sp2[!is.na(Joint_C),]
ggplot(map_3sp2, aes(1:nrow(map_3sp2),Joint_C)) +
  geom_point()


ld_w_int <- compute_ld_summary(ld_structure=ld_struct,
                               method = c("ld_int"),
                               eps = 0.005,
                               d_window = derive_ld_radius(ld_struct$by_chr$Chr1$decay_sum$a, 1-rho_w),
                               shell_type = "median",
                               cores = cores)


draws_joint_3sp <- ld_rho_draws(gds = gds_3sp,
                                ld_struct  = ld_struct_3sp,
                                F_vals     = map_3sp[,.(emx_F_GC,lfmm_F)],
                                q_vals     = map_3sp[,.(emx_q,lfmm_q)],
                                n_inds     = nrow(GTs_3sp),
                                n_draws    = 100,
                                rho_d_lim  = list(min=0.5,max=max_rho),
                                rho_ld_lim = list(min=0.5,max=max_rho),
                                alpha_lim  = list(min=1.31,max=4),
                                lmin_lim   = list(min=1,max=10),
                                cores      = 8,
                                mode       = "joint"

)
joint <- rbindlist(draws_joint_3sp$draws)

draws_joint_3sp$draws$joint

plots_3sp <- lapply(names(draws_joint_3sp$draws),function(rho_name){
  map2 <- add_consistency_to_map(map_3sp, consistency_obj = consistency_score(draws_joint_3sp$draws[[rho_name]]))
  #map2 <- add_consistency_to_map(map2, consistency_obj = consistency_score(draws_per_method))
  ldw = unlist(ld_struct_3sp$ld_w_dt[,..rho_name])
  map2[,ld_w:=ldw]
  #map2[,plot(LD_AUC)]
  #plot(unlist(lapply(ld_struct$by_chr, function(x) x$LD_AUC$LD_AUC)))
  #map2[,plot(ld_w,max_LD_with_QTN )]
  #map2[,plot(Joint_C)]

  #map2[!is.na(Joint_C)]

  #map2[!is.na(Joint_C),plot(Joint_C)]
  sign_th = 0.2
  idx <- which(map2[,Joint_C>sign_th])
  el  <- get_el(gds_3sp, idx, slide_win_ld = -1,by_chr = TRUE)

  ORs_tbl <- detect_or(el,
                       q_vals=map2[,.(Joint_C)],
                       SNP_ids = map_3sp$marker,
                       SNP_chr = map_3sp$Chr,
                       ld_struct=ld_struct_3sp,
                       sign_th=sign_th,
                       sign_if = "greater",
                       rho_d=0.95,
                       rho_ld=0.95,
                       l_min = 1,
                       mode = "per_method",
                       ret_table = TRUE)


  #ORs_tbl[,table(do.call(rbind, strsplit(SNP,":"))[,1],OR_id)]
  map2 <- merge(map2,
                ORs_tbl,
                by.x = "marker",
                by.y = "SNP",
                all.x = TRUE)
  map2 <-  map2[order(as.numeric(gsub("Chr","",map2$Chr), map2$Pos))]
  map2[is.na(OR_id),OR_id:="ns"]
  #map2[,table(OR_id,Chr)]

  unique_ORs <- unique(na.omit(map2$OR_id))
  #unique_ORs <- unique_ORs[unique_ORs != "Chr12_OR10"]

  map2[, OR_factor := factor(OR_id, levels = sample(unique_ORs))]

  keep <- c(map2[,.(marker[which.min(Pos)]),by=Chr][,V1],map2[,.(marker[which.max(Pos)]),by=Chr][,V1])
  layout <- prep_manhattan(
    map2[Chr=="Chr7"][Joint_C>0 | marker %in% keep, .(
      bp = Pos,
      Chr,
      marker,
      Joint_C,
      lfmm_q = -log10(lfmm_q),
      ld_w,
      OR_factor
    )],spacer = 1e6
  )


  plot_manhattan_gg(
    layout,
    y_vars = c("ld_w","Joint_C"),
    y_labels = c("ld_w","Joint_C"),
    thresholds = c(NA,0.2),
    or_var = "OR_factor",
    ncol=2
  )+ggplot2::ggtitle(rho_name)

})

library(patchwork)

wrap_plots(plots_3sp, ncol = 2)
#(plots[[1]] | plots[[2]] | plots[[3]] | plots[[4]])
par(mfcol=c(2,1))
plot(ld_struct_3sp$ld_w_dt[grep("Chr7",SNP)]$ldw_rho_0.7,pch=20,cex=0.5,ylab = "ld_w")
plot(ld_struct_3sp$ld_w_dt[grep("Chr7",SNP)]$ldexc_rho_0.7,pch=20,cex=0.5,ylab = "ld_exc")

ld_struct_3sp$decay_summary[Chr=="Chr7",a/0.0009428298]
d_from_rho(a=2.449425e-05,rho=0.9,d0 = 0)
d_from_rho(a=2.449425e-05,rho=0.4,d0 = 0)
d_from_rho(a=0.0009428298,rho=0.95,d0 = 0)



plot(ld_struct_3sp$ld_w_dt[grep("Chr7",SNP)]$ldw_rho_0.7,ld_struct_3sp$ld_w_dt[grep("Chr7",SNP)]$ldw_rho_0.95)
plot(ld_struct_3sp$ld_w_dt[grep("Chr7",SNP)]$ldexc_rho_0.7,ld_struct_3sp$ld_w_dt[grep("Chr7",SNP)]$ldexc_rho_0.95)

plot(ld_struct_3sp$ld_w_dt$ldw_rho_0.7,ld_struct_3sp$ld_w_dt$ldexc_rho_0.7)
abline(0,1)
#SNP_ids == map_3sp$marker
#table(ld_struct_3sp_w1000$by_chr$Chr1$edges$SNP1 %in% map_3sp$marker)
#saveRDS(draws_joint,"draws_joint.rds")

# max_rho = 0.99
# t1 <- Sys.time()
# ld_struct_3sp_w1000 <-  compute_ld_structure(gds_3sp,
#                                          use         = "robust",
#                                          max_rho     = max_rho,
#                                          adapt_thin  = FALSE
# )
# t2 <- Sys.time()
# difftime(t2,t1)
# gc()
# t1 <- Sys.time()
# max_rho = 0.999
# ld_struct_3sp <-  compute_ld_structure(gds_3sp,
#                                        use         = "robust",
#                                        max_rho     = max_rho,
#                                        cores       = 8
# )
#ld_struct <- ld_struct_3sp
ld_struct$params
#gds <- gds_3sp
t1 <- Sys.time()
draws_3sp <- ld_rho_draws(gds = gds_3sp,
                          ld_struct  = ld_struct_3sp,
                          F_vals     = map_3sp[,.(emx_F_GC,lfmm_F)],
                          q_vals     = map_3sp[,.(emx_q,lfmm_q)],
                          n_inds     = nrow(GTs_3sp),
                          method     = "auc",
                          #n_rho      = 100,
                          n_or_draws = 100,
                          #rho_w_lim  = list(min=0,max=max_rho),
                          rho_d_lim  = list(min=0,max=max_rho),
                          rho_ld_lim = list(min=0,max=max_rho),
                          alpha_lim  = list(min=1.31,max=4),
                          lmin_lim   = list(min=1,max=30),
                          cores      = 8,
                          use        = "robust",
                          mode       = "joint"

)
t2 <- Sys.time()
difftime(t2,t1)
# draws_3sp_per_meth <- ld_rho_draws(
#   gds        = gds_3sp,
#   ld_struct  = ld_struct_3sp,
#   F_vals     = map_3sp[,.(emx_F_GC,lfmm_F)],
#   q_vals     = map_3sp[,.(emx_q,lfmm_q)],
#   n_inds     = n_inds,
#   SNP_ids    = SNP_ids,
#   n_rho      = 100,
#   n_or_draws = 25,
#   rho_w_lim  = list(min=0.75,max=max_rho),
#   rho_d_lim  = list(min=0.9,max=max_rho),
#   rho_ld_lim = list(min=0.5,max=max_rho),
#   alpha_lim  = list(min=1.31,max=4),
#   lmin_lim   = list(min=1,max=30),
#   cores      = 8,
#   mode       = "per_method"
# )

#draws_3sp
#saveRDS(draws_3sp,"draws_3sp.rds")
#draws_wide$draws[OR_size>0,hist(l)]
#tmp <- consistency_score(draws_3s)
#table(unique(unlist(draws_3sp_000$draws$OR)) %in% map_3sp$marker)

# table(tmp$consistency[,SNP] %in% map_3sp$marker)
# tmp$consistency$C[match(map_3sp$marker,tmp$consistency[,SNP])])


map2 <- add_consistency_to_map(map_3sp, consistency_obj = consistency_score(draws_3sp))
map2[,LD_AUC:= unlist(lapply(ld_struct_3sp$by_chr, function(x) x$LD_AUC$LD_AUC))]

#map2[,plot(lfmm_F,ld_w)]
#map2[,plot(Joint_C)]

#map2[!is.na(Joint_C),plot(Joint_C)]
sign_th = 0.2
#map2[,]
idx <- which(map2[,Joint_C>sign_th])
el  <- get_el(gds_3sp, idx, slide_win_ld = -1,by_chr = TRUE)


ORs_tbl <- detect_or(el,
                     q_vals=map2[,.(Joint_C)],
                     SNP_ids = map_3sp$marker,
                     SNP_chr = map_3sp$Chr,
                     ld_struct=ld_struct_3sp,
                     sign_th=sign_th,
                     sign_if = "greater",
                     rho_d=0.95,
                     rho_ld=0.95,
                     l_min = 1,
                     ret_table = TRUE)

#ORs_tbl[,table(do.call(rbind, strsplit(SNP,":"))[,1],OR_id)]
map2 <- merge(map2,
              ORs_tbl,
              by.x = "marker",
              by.y = "SNP",
              all.x = TRUE)
map2 <-  map2[order(as.numeric(gsub("Chr","",map2$Chr), map2$Pos))]
map2[is.na(OR_id),OR_id:="ns"]

map2[,table(OR_id,Chr)]

unique_ORs <- unique(na.omit(map2$OR_id))

map2[, OR_factor := factor(OR_id, levels = sample(unique_ORs))]

#map2[!is.na(OR_factor),.(max(max_LD_with_QTN),.N,focal_QTN[which.max(max_LD_with_QTN)]),by=OR_factor]
#map2[,plot(sqrt(LD_AUC),LD_AUC)]
tmp <- map2[, .(
  bp = Pos,
  Chr,
  marker,
  lfmm_q=-log10(lfmm_q),
  Joint_C,
  LD_AUC = sqrt(LD_AUC),
  OR_factor
)]
#tmp[Joint_C>0.05, Joint_C:=0.05]
layout <- prep_manhattan(
  tmp
)


plot_manhattan_gg(
  layout,
  y_vars = c("lfmm_q","Joint_C","LD_AUC"),
  y_labels = c("lfmm_q","Joint_C","LD_AUC"),
  thresholds = c(1.31,sign_th, NA),
  or_var = "OR_factor",
  ncol=1
)

##########
decay_sum <- ld_struct_3sp_w1000$by_chr$Chr1$summary


plot(seq(0.85,0.95,length.out=11),d_from_rho(a=decay_sum$robust$a,  d0 = decay_sum$robust$d0, rho = seq(0.85,0.95,length.out=11)))
par(mfcol=c(3,1))
plot(compute_ld_w(ld_struct_3sp_w1000,0.85),main="0.85")
plot(compute_ld_w(ld_struct_3sp_w1000,0.90),main="0.9")
plot(compute_ld_w(ld_struct_3sp_w1000,0.95),main="0.95")
gc()

