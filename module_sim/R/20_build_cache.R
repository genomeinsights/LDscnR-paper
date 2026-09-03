## module_sim/R/20_build_cache.R
## Build the AD-HOC CACHE (06a/06b pattern applied to the C-score chain): touch
## genotypes ONCE per (V,c,env), store everything the downstream C-score / tau_C /
## PR-AUC / null estimates need so they run WITHOUT re-pooling or re-clustering.
## Cache (git-ignored) per env: cache_V{V}_c{cc}_env{env}.rds =
##   map (marker,Chr,Pos,true_pos_QTN,emx_p,lfmm_p), LDW (pooled ld_ws, all rho),
##   decay_sum, th (r2min,dmax), qtab (precompute_QTN_LD over the candidate universe),
##   edges (build_edge_cache: r^2 graph among the universe for GTs-free clustering),
##   universe + the grid the universe was built with.
## Validates cluster_from_cache() == cluster_regions() on a test gate before saving.
## Run from LDscnR-paper/:  Rscript module_sim/R/20_build_cache.R [V c env...]
##   (env defaults to 1:5)

source("module_sim/R/_config.R")
a  <- commandArgs(trailingOnly = TRUE)
V  <- if (length(a) >= 1) a[1] else "2"
cc <- if (length(a) >= 2) a[2] else "1"
ENV <- if (length(a) >= 3) a[-(1:2)] else as.character(1:5)
QSTAR <- seq(0, 0.95, by = 0.05); ALPHA_C <- c(0.001, 0.01, 0.05, 0.1)   # universe = loosest grid
RHO_LD <- 0.75; RHO_D <- 0.95; DCAP <- 1e5; ALPHA_S <- c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1)

build_one <- function(env) {
  files <- list.files(SIM_DATA, full.names = TRUE,
                      pattern = sprintf("^adapt_bgs_chr[0-9]+_V%s_c%s_env%s\\.rds$", V, cc, env))
  reps <- as.integer(sub(".*_chr([0-9]+)_.*", "\\1", basename(files))); files <- files[order(reps)]
  maps <- gts <- ldws <- decs <- vector("list", length(files))
  for (i in seq_along(files)) { d <- readRDS(files[i]); m <- as.data.table(d$map)
    m[, `:=`(Chr = paste0("R", i, "_", Chr), marker = paste0("R", i, "_", marker))]
    G <- d$GTs; colnames(G) <- paste0("R", i, "_", colnames(G))
    lw <- d$ld_ws; rownames(lw) <- m$marker
    ds <- as.data.table(d$LD_decay$decay_sum); ds[, Chr := paste0("R", i, "_", Chr)]
    maps[[i]] <- m; gts[[i]] <- G; ldws[[i]] <- lw; decs[[i]] <- ds }
  map <- flag_true_positive_QTNs(rbindlist(maps, fill = TRUE))
  GTs <- do.call(cbind, gts); LDW <- do.call(rbind, ldws)[map$marker, ]
  decay_sum <- rbindlist(decs, fill = TRUE); RHO <- colnames(LDW)
  th <- score_thresholds(decay_sum, rho_r2 = 0.75, rho_d = 0.95, dmax_cap = DCAP)

  ## candidate universe = every marker any gate could pick (C>0 for either method, at
  ## the loosest grid) UNION the single-SNP alpha sets -> superset of all downstream gates
  Cf_e <- cscore_count(map$emx_p,  LDW, RHO, QSTAR, ALPHA_C)
  Cf_l <- cscore_count(map$lfmm_p, LDW, RHO, QSTAR, ALPHA_C)
  ss <- unlist(lapply(c("emx_p", "lfmm_p"), function(p)
    map$marker[p.adjust(map[[p]], "BH") < max(ALPHA_S)]))
  universe <- unique(c(names(Cf_e)[Cf_e > 0], names(Cf_l)[Cf_l > 0], ss))
  qtab  <- precompute_QTN_LD(GTs, map, universe, 2e6, cores = 4)
  edges <- build_edge_cache(universe, map, GTs, decay_sum, RHO_LD, RHO_D, DCAP)

  ## validate: cluster_from_cache reproduces cluster_regions on a mid gate
  gate <- names(Cf_e)[Cf_e >= 0.3]
  if (length(gate)) {
    a1 <- cluster_regions(gate, map, GTs, decay_sum = decay_sum, rho_ld = RHO_LD, rho_d = RHO_D, dcap_max = DCAP)
    a2 <- cluster_from_cache(gate, edges)
    key <- function(x) sort(vapply(x, function(v) paste(sort(v), collapse = "|"), ""))
    stopifnot(identical(key(a1), key(a2)))
  }
  cache <- list(map = map,          # FULL map -- diagnose_OR_classification needs type/Va/focal_QTN/...
                LDW = LDW, decay_sum = decay_sum, th = th, qtab = qtab, edges = edges,
                universe = universe, grid = list(RHO = RHO, QSTAR = QSTAR, ALPHA_C = ALPHA_C,
                RHO_LD = RHO_LD, RHO_D = RHO_D, DCAP = DCAP))
  f <- file.path(dir_data, sprintf("cache_V%s_c%s_env%s.rds", V, cc, env))
  saveRDS(cache, f)
  cat(sprintf("env%s: %d SNPs, universe=%d, qtab=%d, edges=%d chr | cluster check OK | %s\n",
              env, nrow(map), length(universe), nrow(qtab), length(edges), basename(f)))
}
for (e in ENV) build_one(e)
cat("cache build done\n")
