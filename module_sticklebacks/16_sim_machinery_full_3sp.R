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

## reconciled engine: LD-pruned GCTA GRM, used for BOTH observed and null
pruned <- intersect(readRDS("/Users/petrikem/gitlab/LD-scaling-genome-scans/empirical_data/3sp/pruning.rds")$pruned_SNPs, colnames(GTs))
K <- gcta_grm(GTs[, pruned]); pre <- fast_emmax_setup(GTs, K)
p_obs <- fast_emmax_p(pre, eco)
cat(sprintf("engine: pruned GCTA GRM (%d markers); fast-EMMAX vs precomputed emx_p Pearson=%.3f\n",
            length(pruned), cor(-log10(p_obs), -log10(sr$emx_p))))

eK <- eigen(K, symmetric = TRUE); Lv <- pmax(eK$values, 0); Vk <- eK$vectors
gen_surrogate <- function() { y <- as.numeric(Vk %*% (sqrt(Lv) * rnorm(n))); as.numeric(resid(lm(y ~ eco))) }
sparseC <- function(pv) { C <- cscore_count(pv, LDW, RHO, QSTAR, ALPHA); C[C > 0] }

cat(sprintf("3sp full run: %d SNPs, B=%d genetic surrogates, l_min=%d, r2=%.1f/%.1fMb\n", nrow(sr), B, LMIN, R2LINK, DC/1e6))
C_obs <- sparseC(p_obs)
set.seed(1); C_surr <- vector("list", B)
for (b in seq_len(B)) { C_surr[[b]] <- sparseC(fast_emmax_p(pre, gen_surrogate())); if (b %% 10 == 0) cat("surrogate", b, "/", B, "\n") }

uni <- unique(c(names(C_obs), unlist(lapply(C_surr, names))))
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

## surviving loci at the calibrated tau_C
tc <- if (!is.na(t05)) t05 else t10
mk <- names(C_obs)[C_obs >= tc]; reg <- cluster_from_cache(mk, edges); reg <- reg[lengths(reg) >= LMIN]
CH <- setNames(sr$Chr, sr$marker)
loci <- rbindlist(lapply(seq_along(reg), function(i) data.table(
  Chr = names(sort(table(CH[reg[[i]]]), decreasing = TRUE))[1], size = length(reg[[i]]), maxC = max(C_obs[reg[[i]]]))))
cat(sprintf("\n=== EMMAX surviving loci at tau_C=%.2f, l_min=10 (%d clusters) ===\n", tc, nrow(loci)))
if (nrow(loci)) print(loci[order(-size)])
saveRDS(list(res = res, C_obs = C_obs, loci = loci, tau_c = tc, engine = "pruned_gcta_grm", pruned_n = length(pruned)),
        file.path(mod, "sim_machinery_3sp.rds"))

## Manhattan of C_obs with surviving-cluster members coloured
sr[, C := 0]; sr[names(C_obs), C := C_obs, on = "marker"]
survmk <- unlist(reg); sr[, surv := marker %in% survmk]
chr_lev <- paste0("Chr", c(1:18,20,21)); sr[, Chr := factor(Chr, levels = chr_lev)]; setorder(sr, Chr, Pos)
clen <- sr[, .(mx = max(Pos)), by = Chr]; clen[, off := cumsum(shift(mx, fill = 0))]; xoff <- setNames(clen$off, as.character(clen$Chr))
sr[, gx := Pos + xoff[as.character(Chr)]]; ax <- clen[, .(center = off + mx/2, Chr)]
p <- ggplot() +
  geom_point(data = sr[C > 0 & !surv], aes(gx, C), color = "grey70", size = 0.5) +
  geom_point(data = sr[surv == TRUE], aes(gx, C, color = Chr), size = 1.3) +
  geom_hline(yintercept = tc, linetype = 2, color = "#D62828") +
  scale_x_continuous(breaks = ax$center, labels = sub("Chr", "", ax$Chr), expand = c(0.01, 0)) +
  scale_color_viridis_d(guide = "none") +
  labs(title = sprintf("3sp full sim-machinery run (reconciled pruned-GRM engine): EMMAX C, l_min=10, tau_C=%.2f", tc),
       subtitle = sprintf("coloured = members of >=10-SNP clusters surviving FDR<=0.05; red dashed = tau_C; %d loci", nrow(loci)),
       x = "Chromosome", y = "C-score") + theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())
ggsave(file.path(mod, "fig_sim_machinery_3sp.png"), p, width = 12, height = 5, dpi = 150)
cat("wrote fig_sim_machinery_3sp.png\n")
