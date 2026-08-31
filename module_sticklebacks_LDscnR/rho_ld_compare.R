## Does 3sp work at rho_ld = 0.75? Compare against the published 0.60, on the
## CANONICAL (legacy) null bundles. The saved C-scores are rho_ld-independent --
## rho_ld enters only at ld_edges() -- so this needs no new EMMAX/LFMM scans.
##
## ld_from_rho(b, c, rho) = b + (c - b) * (1 - rho): LOWER rho_ld => HIGHER r^2
## link => stricter => more splitting. So 0.75 should MERGE relative to 0.60.
suppressMessages({ library(data.table); library(LDscnR) })
setwd("~/gitlab/LDscnR-paper")
RES <- "module_sticklebacks_LDscnR/results"
OUT <- "/private/tmp/claude-539526166/-Users-petrikem-gitlab-LDscnR/f5c2953d-f16b-4266-bda5-08c843e9b161/scratchpad"
TAU <- 0.05; LMIN <- 3L; DCAP <- 1e5; FDR <- 0.05
RHOS <- c(0.60, 0.75)

d <- readRDS("module_sticklebacks_LDscnR/data/3sp_LDscnR_data.rds")
map <- as.data.table(d$map); decay <- as.data.table(d$LD_decay$decay_sum)

legacy <- c(genetic = "null_uncapped_3sp.rds", global_perm = "null_popperm_3sp.rds",
            region_perm = "null_regionperm_3sp.rds", spatial = "null_spatial_3sp.rds",
            latent = "null_latent_3sp.rds")
nulls <- lapply(names(legacy), function(k) { x <- readRDS(file.path(RES, legacy[[k]]))
  class(x) <- "ld_null"; x$basis <- k; x$engine <- "EMMAX"; x })
names(nulls) <- names(legacy)

## what r^2 link does each rho imply, per chromosome?
cat("=== implied r^2 link per rho (from the decay fit) ===\n")
for (r in RHOS) {
  r2 <- ld_from_rho(decay$b, if ("c" %in% names(decay)) decay$c else rep(1, nrow(decay)), r)
  cat(sprintf("  rho_ld = %.2f -> r2 link: median %.3f  range %.3f - %.3f\n",
              r, median(r2), min(r2), max(r2)))
}

univ <- unique(unlist(lapply(nulls, `[[`, "universe"), use.names = FALSE))
INV <- c(21.40e6, 21.93e6); EDA <- c(12.79e6, 12.83e6)

res <- list(); regs <- list()
for (r in RHOS) {
  cat(sprintf("\n\n################ rho_ld = %.2f ################\n", r)); flush.console()
  t0 <- Sys.time()
  edges <- ld_edges(univ, d$GTs, map[, .(marker, Chr, Pos)], decay, rho_ld = r, dcap = DCAP)
  cat(sprintf("[edges] %d markers in %.0f s\n", length(univ), as.numeric(Sys.time()-t0, units="secs"))); flush.console()

  g <- ld_gate(nulls, edges, tau = TAU, l_min = LMIN)
  gg <- as.data.table(g)[, .(basis, B, obs_regions, med_regions, max_regions, ratio, pass)]
  gg[, rho_ld := r]; res[[length(res)+1L]] <- gg
  cat("\n-- gate --\n"); print(gg)

  s <- ld_region_scan(nulls$region_perm, edges, tau = TAU, l_min = LMIN, fdr = FDR)
  rr <- copy(s$regions); rr[, rho_ld := r]
  rr[, region := sprintf("%s:%.2f-%.2f", chr, lo/1e6, hi/1e6)]
  rr[chr == "Chr1" & lo < INV[2] & hi > INV[1], region := paste0(region, " (inversion)")]
  rr[chr == "Chr4" & lo < EDA[2] & hi > EDA[1], region := paste0(region, " (Eda)")]
  regs[[length(regs)+1L]] <- rr
  cat(sprintf("\n-- regional-permutation scan: %d regions, %d significant, max q_R = %.4f --\n",
              nrow(rr), sum(rr$sig), max(rr$q_R)))
  print(rr[order(-size), .(region, size, s_R = round(s_R,2), p_R = round(p_R,4), q_R = round(q_R,4), sig)], nrow = 40)
  cat(sprintf("[landmarks] Chr1 inversion regions: %d | Chr4 Eda regions: %d | Chr4 regions total: %d\n",
              nrow(rr[chr=="Chr1" & lo < INV[2] & hi > INV[1]]),
              nrow(rr[chr=="Chr4" & lo < EDA[2] & hi > EDA[1]]),
              nrow(rr[chr=="Chr4"])))
}
G <- rbindlist(res); R <- rbindlist(regs)
fwrite(G, file.path(OUT, "rho_gate_compare.csv")); fwrite(R, file.path(OUT, "rho_regions_compare.csv"))

cat("\n\n################ SIDE BY SIDE ################\n")
cat("\n-- gate: median surrogate regions / observed --\n")
print(dcast(G, basis ~ rho_ld, value.var = c("obs_regions","med_regions","pass")))
cat("\n-- regional-permutation scan --\n")
print(R[, .(n_regions = .N, n_sig = sum(sig), max_q = round(max(q_R),4),
            biggest = max(size), total_markers = sum(size)), by = rho_ld])
cat("\n[done]\n")
