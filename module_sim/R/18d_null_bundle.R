## module_sim/R/18d_null_bundle.R
## Structured null that saves a REUSABLE BUNDLE so l_min never enters null GENERATION:
## the expensive part (surrogate phenotype -> fast EMMAX -> per-SNP C) is l_min-free;
## only the final cluster-and-count depends on l_min. From ONE run we save
##   - C_obs         : observed per-SNP C (sparse, C>0)
##   - C_surr        : list of B surrogate per-SNP C (sparse)
##   - edges         : edge cache over the UNION of every marker any surrogate lights up
##                     (+ observed) -> clustering any gate later needs no genotypes
##   - map           : marker/Chr/Pos/truth
## 18e_null_fdr.R then extracts FDR(tau_C, l_min) for ANY l_min as a pure lookup.
## Sim keeps the alpha sweep. Run: Rscript module_sim/R/18d_null_bundle.R [V c env] [B]

source("module_sim/R/_config.R")
a   <- commandArgs(trailingOnly = TRUE)
V   <- if (length(a) >= 1) a[1] else "2"
cc  <- if (length(a) >= 2) a[2] else "1"
env <- if (length(a) >= 3) a[3] else "1"
B   <- if (length(a) >= 4) as.integer(a[4]) else 100L
QSTAR <- seq(0, 0.95, by = 0.05); ALPHA_C <- c(0.001, 0.01, 0.05, 0.1)   # sim: alpha swept
RHO_LD <- 0.75; RHO_D <- 0.95; DCAP <- 5e5
gcta_grm <- function(X){ p<-colMeans(X)/2; k<-p>0&p<1; X<-X[,k,drop=F]; p<-p[k]
  Z<-sweep(sweep(X,2,2*p,"-"),2,sqrt(2*p*(1-p)),"/"); tcrossprod(Z)/ncol(Z) }

files <- list.files(SIM_DATA, full.names = TRUE,
                    pattern = sprintf("^adapt_bgs_chr[0-9]+_V%s_c%s_env%s\\.rds$", V, cc, env))
reps <- as.integer(sub(".*_chr([0-9]+)_.*", "\\1", basename(files))); files <- files[order(reps)]
Yobs <- COORDS <- NULL; PRE <- maps <- ldws <- mk_i <- gts <- decs <- vector("list", length(files))
for (i in seq_along(files)) { d <- readRDS(files[i]); m <- as.data.table(d$map)
  m[, marker := paste0("R", i, "_", marker)][, Chr := paste0("R", i, "_", Chr)]
  if (is.null(Yobs)) { Yobs <- d$env$env; COORDS <- cbind(d$env$x, d$env$y) }
  G <- d$GTs; colnames(G) <- m$marker
  PRE[[i]] <- fast_emmax_setup(G, gcta_grm(G)); gts[[i]] <- G
  lw <- d$ld_ws; rownames(lw) <- m$marker
  ds <- as.data.table(d$LD_decay$decay_sum); ds[, Chr := paste0("R", i, "_", Chr)]
  maps[[i]] <- m; ldws[[i]] <- lw; mk_i[[i]] <- m$marker; decs[[i]] <- ds }
map <- flag_true_positive_QTNs(rbindlist(maps, fill = TRUE))
LDW <- do.call(rbind, ldws)[map$marker, ]; GTs <- do.call(cbind, gts)[, map$marker]
decay_sum <- rbindlist(decs, fill = TRUE); RHO <- colnames(LDW); n <- length(Yobs)

Dm <- as.matrix(dist(COORDS)); l <- median(Dm[lower.tri(Dm)])
Ks <- exp(-0.5 * (Dm / l)^2); eK <- eigen(Ks, symmetric = TRUE); Lv <- pmax(eK$values, 0); Vk <- eK$vectors
gen_surrogate <- function() { y <- Vk %*% (sqrt(Lv) * rnorm(n)); as.numeric(resid(lm(y ~ Yobs))) }
emmax_pooled <- function(Yv) { pv <- numeric(0)
  for (i in seq_along(PRE)) pv <- c(pv, setNames(fast_emmax_p(PRE[[i]], Yv), mk_i[[i]]))
  pv[rownames(LDW)] }
sparseC <- function(pv) { C <- cscore_count(pv, LDW, RHO, QSTAR, ALPHA_C); C[C > 0] }   # keep only C>0

cat(sprintf("V%s_c%s_env%s: %d SNPs, B=%d surrogates -> null BUNDLE (spatial l=%.2f)\n", V, cc, env, nrow(map), B, l))
C_obs <- sparseC(emmax_pooled(Yobs))
set.seed(1); C_surr <- vector("list", B)
for (b in seq_len(B)) { C_surr[[b]] <- sparseC(emmax_pooled(gen_surrogate()))
  if (b %% 10 == 0) cat("surrogate", b, "/", B, "\n") }

uni <- unique(c(names(C_obs), unlist(lapply(C_surr, names))))
cat(sprintf("union universe (obs + all surrogates C>0): %d markers -> building edge cache\n", length(uni)))
edges <- build_edge_cache(uni, map, GTs, decay_sum, RHO_LD, RHO_D, DCAP)
saveRDS(list(C_obs = C_obs, C_surr = C_surr, edges = edges,
             map = map[, .(marker, Chr, Pos, true_pos_QTN)], B = B,
             params = list(V = V, c = cc, env = env, QSTAR = QSTAR, ALPHA_C = ALPHA_C,
                           RHO = RHO, RHO_LD = RHO_LD, RHO_D = RHO_D, DCAP = DCAP)),
        file.path(dir_data, sprintf("null_bundle_V%s_c%s_env%s.rds", V, cc, env)))
cat(sprintf("wrote null_bundle_V%s_c%s_env%s.rds (obs C>0=%d; extract any l_min with 18e)\n",
            V, cc, env, length(C_obs)))
