## module_sticklebacks/16_sim_machinery_full_3sp.R
## FULL sim-machinery run on 3sp with a RECONCILED engine: fast EMMAX with ONE GRM for
## BOTH observed and null. Engine = LD-pruned GCTA GRM (pruning.rds pruned_SNPs; the
## recipe direction -- pruning removes signal-carrying LD blocks so the GRM doesn't
## absorb the concentrated ecotype signal, unlike an all-markers GRM which over-corrects).
## Pipeline: per-SNP C (alpha=0.05) -> genetic structured null (y~MVN(0,K) _|_ ecotype)
## -> cluster (r2=0.1/0.5Mb) + l_min=10 region-FDR -> tau_C -> surviving loci.
## Run: Rscript module_sticklebacks/16_sim_machinery_full_3sp.R [B]
## Output: module_sticklebacks/sim_machinery_3sp.rds + fig_sim_machinery_3sp.png

suppressMessages({ library(LDscnR); library(data.table); library(ggplot2) })
source("module_sim/R/_config.R")
mod <- "/Users/petrikem/gitlab/LDscnR-paper/module_sticklebacks"
a <- commandArgs(trailingOnly = TRUE); B <- if (length(a) >= 1) as.integer(a[1]) else 100L
ALPHA <- 0.05; QSTAR <- seq(0, 0.95, by = 0.05); TAU <- seq(0.02, 1, by = 0.02)
R2LINK <- 0.1; DC <- 5e5; LMIN <- 10L
gcta_grm <- function(X){ p<-colMeans(X)/2; k<-p>0&p<1; X<-X[,k,drop=F]; p<-p[k]
  Z<-sweep(sweep(X,2,2*p,"-"),2,sqrt(2*p*(1-p)),"/"); tcrossprod(Z)/ncol(Z) }

sr <- readRDS(file.path(mod, "snp_stats_aligned.rds")); setDT(sr)
LDW <- readRDS("/Users/petrikem/gitlab/LDscnR-paper/3sp_data/ld_ws_3sp_MAF01.rds")[sr$marker, ]; RHO <- colnames(LDW)
decs <- as.data.table(readRDS(file.path(mod, "decay_sum_3sp.rds")))
e <- new.env(); load("/Users/petrikem/gitlab/LDscnR-paper/LFMM_3sp/data/3sp_data.RData", envir = e)
GTs <- e$GTs_3sp; colnames(GTs) <- e$map_3sp$marker; GTs <- GTs[, sr$marker]
eco <- as.integer(e$pheno_3sp$ecotype == "Marine"); n <- length(eco)

## reconciled engine: GCTA GRM from LOW-ld_w (LD-independent) markers, ld_w_095 < 0.05
## -- the recipe's neutral-background GRM; avoids absorbing the concentrated ecotype
## signal (recovers Chr4/Eda 0.53, Chr1 0.72 vs all-markers 0.09). Used for BOTH obs+null.
LDW_THRESH <- 0.05
keep_grm <- which(LDW[, "0.95"] < LDW_THRESH & is.finite(LDW[, "0.95"]))
K <- gcta_grm(GTs[, keep_grm]); pre <- fast_emmax_setup(GTs, K)
p_obs <- fast_emmax_p(pre, eco)
cat(sprintf("engine: low-ld_w GCTA GRM (ld_w<%.2f, %d markers); fast-EMMAX vs precomputed emx_p Pearson=%.3f\n",
            LDW_THRESH, length(keep_grm), cor(-log10(p_obs), -log10(sr$emx_p))))

eK <- eigen(K, symmetric = TRUE); Lv <- pmax(eK$values, 0); Vk <- eK$vectors
gen_surrogate <- function() { y <- as.numeric(Vk %*% (sqrt(Lv) * rnorm(n))); as.numeric(resid(lm(y ~ eco))) }
sparseC <- function(pv) { C <- cscore_count(pv, LDW, RHO, QSTAR, ALPHA); C[C > 0] }

cat(sprintf("3sp full run: %d SNPs, B=%d genetic surrogates, l_min=%d, r2=%.1f/%.1fMb\n", nrow(sr), B, LMIN, R2LINK, DC/1e6))
C_emx_full  <- cscore_count(p_obs,      LDW, RHO, QSTAR, ALPHA); C_obs  <- C_emx_full[C_emx_full > 0]
C_lfmm_full <- cscore_count(sr$lfmm_p,  LDW, RHO, QSTAR, ALPHA); C_lfmm <- C_lfmm_full[C_lfmm_full > 0]
set.seed(1); C_surr <- vector("list", B)
for (b in seq_len(B)) { C_surr[[b]] <- sparseC(fast_emmax_p(pre, gen_surrogate())); if (b %% 10 == 0) cat("surrogate", b, "/", B, "\n") }

uni <- unique(c(names(C_obs), names(C_lfmm), unlist(lapply(C_surr, names))))   # incl LFMM so it can cluster
cat(sprintf("union universe %d markers -> edge cache (r2=%.1f/%.1fMb)\n", length(uni), R2LINK, DC/1e6))
edges <- build_edge_cache(uni, sr[, .(marker, Chr, Pos)], GTs, decs, r2_link = R2LINK, dcap = DC); rm(GTs); gc()

n_reg <- function(Csp, tau, lmin) { mk <- names(Csp)[Csp >= tau]; if (!length(mk)) return(0L)
  sum(lengths(cluster_from_cache(mk, edges)) >= lmin) }
n_obs  <- vapply(TAU, function(t) n_reg(C_obs, t, LMIN), numeric(1))
n_null <- vapply(TAU, function(t) mean(vapply(C_surr, n_reg, numeric(1), tau = t, lmin = LMIN)), numeric(1))
FDR <- pmin(1, n_null / pmax(n_obs, 1))
tau_at <- function(q){ ok <- which(FDR <= q & n_obs > 0); if (length(ok)) TAU[min(ok)] else NA_real_ }
t05 <- tau_at(0.05); t10 <- tau_at(0.10)
res <- data.table(tau = TAU, n_obs = n_obs, n_null = round(n_null, 2), FDR = round(FDR, 3))
cat("\n=== region-FDR (l_min=10) vs tau_C ===\n"); print(res[tau %in% seq(0.05, 0.5, 0.05)])
cat(sprintf("\ntau_C: FDR<=0.10 -> %.2f | FDR<=0.05 -> %.2f\n", t10, t05))

## calibrated tau_C (EMMAX null); applied to EMMAX and to LFMM (LFMM has no fast null,
## so it inherits EMMAX's calibration -- directly and genomic-control-mapped)
tc <- if (!is.na(t05)) t05 else t10
CH <- setNames(sr$Chr, sr$marker)
get_loci <- function(Cvec, tau, tag) {
  mk <- names(Cvec)[Cvec >= tau]; reg <- if (length(mk)) cluster_from_cache(mk, edges) else list()
  reg <- reg[lengths(reg) >= LMIN]
  loci <- if (length(reg)) rbindlist(lapply(seq_along(reg), function(i) data.table(
    Chr = names(sort(table(CH[reg[[i]]]), decreasing = TRUE))[1], size = length(reg[[i]]), maxC = max(Cvec[reg[[i]]])))) else data.table(Chr=character(), size=integer(), maxC=numeric())
  cat(sprintf("\n=== %s: tau=%.3f, l_min=%d -> %d clusters across %d chr ===\n", tag, tau, LMIN, nrow(loci), if(nrow(loci)) uniqueN(loci$Chr) else 0L))
  if (nrow(loci)) print(loci[order(-size)][seq_len(min(.N, 12))])
  list(reg = reg, loci = loci) }

emx     <- get_loci(C_obs,  tc, "EMMAX @ EMMAX-calibrated tau_C")
lf_naive<- get_loci(C_lfmm, tc, "LFMM @ EMMAX tau_C (naive, same threshold)")
q_anchor<- mean(C_emx_full < tc); tc_lfmm <- as.numeric(quantile(C_lfmm_full, q_anchor))
lf_gc   <- get_loci(C_lfmm, tc_lfmm, sprintf("LFMM @ genomic-control-mapped tau_C=%.3f (EMMAX quantile %.4f)", tc_lfmm, q_anchor))

saveRDS(list(res = res, C_emx = C_emx_full, C_lfmm = C_lfmm_full, tau_c = tc, tau_lfmm = tc_lfmm,
             emx = emx$loci, lfmm_naive = lf_naive$loci, lfmm_gc = lf_gc$loci,
             engine = sprintf("lowldw_gcta_grm_<%.2f", LDW_THRESH), grm_n = length(keep_grm)),
        file.path(mod, "sim_machinery_3sp.rds"))

## two-panel Manhattan: EMMAX (its tau_C) and LFMM (GC-mapped tau_C)
chr_lev <- paste0("Chr", c(1:18,20,21)); sr[, Chr := factor(Chr, levels = chr_lev)]; setorder(sr, Chr, Pos)
clen <- sr[, .(mx = max(Pos)), by = Chr]; clen[, off := cumsum(shift(mx, fill = 0))]; xoff <- setNames(clen$off, as.character(clen$Chr))
sr[, gx := Pos + xoff[as.character(Chr)]]; ax <- clen[, .(center = off + mx/2, Chr)]
mkpanel <- function(Cfull, reg, tau, ttl) {
  d <- data.table(marker = names(Cfull), C = as.numeric(Cfull)); d <- d[sr[, .(marker, Chr, gx)], on = "marker"][C > 0]
  d[, surv := marker %in% unlist(reg)]
  ggplot() + geom_point(data = d[surv == FALSE], aes(gx, C), color = "grey70", size = 0.4) +
    geom_point(data = d[surv == TRUE], aes(gx, C, color = Chr), size = 1.2) +
    geom_hline(yintercept = tau, linetype = 2, color = "#D62828") +
    scale_x_continuous(breaks = ax$center, labels = sub("Chr", "", ax$Chr), expand = c(0.01, 0)) +
    scale_color_viridis_d(guide = "none") + labs(title = ttl, x = NULL, y = "C-score") +
    theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()) }
p <- mkpanel(C_emx_full, emx$reg, tc, sprintf("EMMAX (low-ld_w GRM), l_min=10, tau_C=%.3f -> %d loci", tc, nrow(emx$loci))) /
     mkpanel(C_lfmm_full, lf_gc$reg, tc_lfmm, sprintf("LFMM @ genomic-control-mapped tau_C=%.3f, l_min=10 -> %d loci (naive tau_C=%.3f gave %d)", tc_lfmm, nrow(lf_gc$loci), tc, nrow(lf_naive$loci)))
suppressMessages(library(patchwork))
ggsave(file.path(mod, "fig_sim_machinery_3sp.png"), p, width = 13, height = 7, dpi = 150)
cat("wrote fig_sim_machinery_3sp.png\n")
