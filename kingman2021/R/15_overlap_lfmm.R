## Same Kingman-EcoPeak region overlap as R/12, but for the LFMM 3sp outlier regions
## at the manuscript's tight clustering (tau_C=0.05, l_min=3, rho_ld=0.60) rather than
## the module's earlier l_min=10 LFMM set. Region-level overlap in gasAcu1, within-
## chromosome rotation null (B=2000). Writes regions_tau0.05_lmin3_rho0.60_lfmm.csv,
## overlap_summary_lfmm.csv, overlap_detail_lfmm.csv.
suppressMessages({ library(data.table); library(LDscnR) })
P    <- path.expand("~/gitlab/LDscnR-paper/kingman2021")
LIFT <- file.path(P, "data", "liftover")
ROMAN <- c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV",
           "XVI","XVII","XVIII","XIX","XX","XXI")
rn <- function(x) match(sub("^chr", "", x), ROMAN)
B  <- 2000L; set.seed(1)

## ---- LFMM regions at tau=0.05, l_min=3, rho_ld=0.60 --------------------------------
d  <- readRDS(path.expand("~/gitlab/LDscnR-paper/module_sticklebacks_LDscnR/data/3sp_LDscnR_data.rds"))
m3 <- as.data.table(d$map)
C  <- ld_cscore(m3$lfmm_p, d$ld_ws, alpha = 0.05, qstar = seq(0, 0.95, by = 0.05)); names(C) <- m3$marker
uni <- names(C)[C > 0]
edges <- ld_edges(uni, d$GTs, m3[, .(marker, Chr, Pos)],
                  as.data.table(d$LD_decay$decay_sum), rho_ld = 0.60, dcap = 1e5)
regs <- ld_regions(names(C)[C >= 0.05], edges); regs <- regs[lengths(regs) >= 3L]
mpos <- stats::setNames(m3$Pos, m3$marker); mchr <- stats::setNames(as.character(m3$Chr), m3$marker)
L <- rbindlist(lapply(seq_along(regs), function(i) { r <- regs[[i]]
  data.table(method = "LFMM", region = i, Chr = unname(mchr[r[1]]),
             start = min(mpos[r]), end = max(mpos[r]), n_snp = length(r)) }))
L[, `:=`(span = end - start, chr_num = as.integer(gsub("Chr", "", Chr)))]
fwrite(L, file.path(P, "data", "regions_tau0.05_lmin3_rho0.60_lfmm.csv"))
cat(sprintf("[1] %d LFMM l_min=3 regions, span %.2f Mb (chroms: %s)\n",
            nrow(L), sum(L$span)/1e6, paste(sort(unique(L$Chr)), collapse = ",")))
rng <- m3[, .(lo = min(Pos), hi = max(Pos)), by = .(chr_num = as.integer(gsub("Chr","",Chr)))]; setkey(rng, chr_num)

## ---- overlap machinery (from R/09/R/12) --------------------------------------------
rd <- function(s) { Pk <- fread(file.path(LIFT, paste0("pv_",s,".bed")), header=FALSE,
                                col.names=c("chr","start","end","pv"))
  Pk[, chr_num := rn(chr)]; Pk[, p_snp := as.numeric(tstrsplit(pv,"\\|")[[1]])]
  Pk <- Pk[!is.na(chr_num), .(chr_num,start,end,p_snp)]; setkey(Pk,chr_num,start,end); Pk }
ovl <- function(A, Pk) { A2 <- A[, .(chr_num,start,end,id=.I)]
  o <- foverlaps(A2, Pk, by.x=c("chr_num","start","end"), type="any", nomatch=NULL)
  if (!nrow(o)) return(list(nhit=0L, bp=0)); o[, w := pmin(i.end,end)-pmax(i.start,start)]
  list(nhit=uniqueN(o$id), bp=sum(o$w)) }

SETS <- c("c155.specific","c155.sensitive","c150.specific","c150.sensitive")
summ <- rbindlist(lapply(SETS, function(s) {
  Pk <- rd(s); o <- ovl(L, Pk)
  r <- rng[J(L$chr_num)]; nh <- integer(B); nb <- numeric(B)
  for (b in seq_len(B)) { st <- r$lo + floor(runif(nrow(L))*pmax(1, r$hi-r$lo-L$span))
    on <- ovl(data.table(chr_num=L$chr_num, start=st, end=st+L$span), Pk)
    nh[b] <- on$nhit; nb[b] <- on$bp }
  data.table(set=s, n_peaks=nrow(Pk), n_regions=nrow(L), regions_hit=o$nhit,
             regions_hit_null=round(mean(nh),3), fold_region=round(o$nhit/mean(nh),3),
             p_region=(1+sum(nh>=o$nhit))/(B+1),
             bp_Mb=round(o$bp/1e6,3), fold_bp=round(o$bp/mean(nb),3), p_bp=(1+sum(nb>=o$bp))/(B+1)) }))
fwrite(summ, file.path(P,"data","overlap_summary_lfmm.csv"))
cat("\n=== LFMM l_min=3 region overlap vs Kingman EcoPeaks (rotation null B=2000) ===\n"); print(summ)

## ---- per-region detail (specific sets) ---------------------------------------------
det <- rbindlist(lapply(c(GlobalSpec="c155.specific", PacSpec="c150.specific"), function(s) {
  Pk <- rd(s)
  o <- foverlaps(L[,.(chr_num,start,end,Chr,region)], Pk, by.x=c("chr_num","start","end"),
                 type="any", nomatch=NULL)
  if (!nrow(o)) return(NULL)
  o[, `:=`(ov_start=pmax(i.start,start), ov_end=pmin(i.end,end))][, ov_kb := round((ov_end-ov_start)/1e3,2)]
  o[, .(Chr, region, reg_start=i.start, reg_end=i.end, p_snp, ov_kb)] }), idcol="set")
if (nrow(det)) { setorder(det, -ov_kb); fwrite(det, file.path(P,"data","overlap_detail_lfmm.csv")) }
cat(sprintf("\n[3] LFMM regions hitting a specific peak: Global=%d Pacific=%d (of %d)\n",
            uniqueN(det[set=="GlobalSpec", region]), uniqueN(det[set=="PacSpec", region]), nrow(L)))
