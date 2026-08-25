## Acceptance test: the new package API must reproduce the paper's existing
## numbers exactly -- both the ad-hoc gate table (gate_background.csv) and the
## location-matched region p-values (region_emp_pvals_*.csv).
suppressMessages({ library(LDscnR); library(data.table) })
setwd("~/gitlab/LDscnR-paper")
RES <- "module_sticklebacks_LDscnR/results"
d <- readRDS("module_sticklebacks_LDscnR/data/3sp_LDscnR_data.rds")
map <- as.data.table(d$map); decay <- as.data.table(d$LD_decay$decay_sum)

files <- c(genetic      = "null_emmax_genetic_3sp.rds",
           global_perm  = "null_emmax_global_perm_3sp.rds",
           region_perm  = "null_emmax_region_perm_3sp.rds",
           spatial      = "null_emmax_spatial_3sp.rds",
           latent       = "null_emmax_latent_3sp.rds")
nulls <- lapply(files, function(f) { x <- readRDS(file.path(RES, f))
                                     class(x) <- "ld_null"; x })

## one shared edge graph over the union of universes, as ld_gate() documents
univ <- unique(unlist(lapply(nulls, `[[`, "universe"), use.names = FALSE))
cat(sprintf("[1] shared edge graph over %d markers\n", length(univ))); flush.console()
t0 <- Sys.time()
edges <- ld_edges(univ, d$GTs, map[, .(marker, Chr, Pos)], decay, rho_ld = 0.60, dcap = 5e5)
cat(sprintf("    built in %.1f s\n", as.numeric(Sys.time()-t0, units="secs"))); flush.console()

cat("\n########## A. ld_gate() vs the ad-hoc gate_background.csv ##########\n")
g <- ld_gate(nulls, edges, tau = 0.05, l_min = 3L)
print(g)
ref <- fread(file.path(RES, "gate_background.csv"))[set == "emmax_v2"]
key <- c("genetic (MVN on K)"="genetic", "global permutation"="global_perm",
         "regional permutation"="region_perm", "spatial (GP kernel)"="spatial",
         "latent (PC subspace)"="latent")
ref[, basis2 := key[basis]]
cmp <- merge(as.data.table(g)[, .(basis2 = basis, new_reg = med_regions, new_c = med_Cgt0,
                                  new_obs = obs_regions)],
             ref[, .(basis2, old_reg = med_regions, old_c = med_Cgt0, old_obs = obs_regions)],
             by = "basis2")
print(cmp)
cat("\nGATE MATCH (median regions):", identical(cmp$new_reg, as.numeric(cmp$old_reg)), "\n")
cat("GATE MATCH (median C>0)   :", identical(cmp$new_c,   as.numeric(cmp$old_c)),   "\n")
cat("GATE MATCH (observed)     :", identical(as.integer(cmp$new_obs), as.integer(cmp$old_obs)), "\n")

cat("\n########## B. ld_region_scan() vs region_emp_pvals_emmax_region_perm ##########\n")
s <- ld_region_scan(nulls$region_perm, edges, tau = 0.05, l_min = 3L, fdr = 0.05)
print(s)
old <- fread(file.path(RES, "region_emp_pvals_emmax_region_perm_tau0.05_lmin3_rho0.60.csv"))
new <- copy(s$regions)
setnames(new, c("chr","s_R"), c("Chr","score"))
m <- merge(old[, .(Chr, lo, hi, score, size, emp_p, emp_q)],
           new[, .(Chr, lo, hi, score, size, p_R, q_R)],
           by = c("Chr","lo","hi"), suffixes = c(".old",".new"))
cat(sprintf("\nmerged %d of %d observed regions (old) / %d (new)\n", nrow(m), nrow(old), nrow(new)))
m[, `:=`(dp = abs(emp_p - p_R), dq = abs(emp_q - q_R), ds = abs(score.old - score.new))]
print(m[order(-dp)][1:5, .(Chr, lo_Mb = round(lo/1e6,2), size.old, size.new,
                            emp_p, p_R, emp_q, q_R, dp, dq)])
cat("\nREGION COUNT MATCH:", nrow(m) == nrow(old) && nrow(m) == nrow(new), "\n")
cat("SCORE  MATCH (max abs diff):", max(m$ds), "\n")
cat("emp_p  MATCH (max abs diff):", max(m$dp), "\n")
cat("emp_q  MATCH (max abs diff):", max(m$dq), "\n")
cat("all q_R <= 0.0149 (framework's claim):", all(new$q_R <= 0.0149), "| max q_R =", max(new$q_R), "\n")
cat("n significant at q<0.05:", sum(new$sig), "of", nrow(new), "\n")
