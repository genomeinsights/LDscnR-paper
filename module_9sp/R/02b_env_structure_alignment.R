## =============================================================================
## module_9sp/R/02b_env_structure_alignment.R
##
## A pre-scan diagnostic shared by the three- and nine-spined stickleback
## panels. It asks two related but deliberately separate questions:
##
##   1. How much of the population-level environmental contrast is explained
##      by the first five axes of the genomic relationship matrix? This is the
##      same descriptive measure used for the simulations.
##   2. How much is explained by the sampling groups used to construct the
##      empirical null (geographic region in 3sp; lineage in 9sp)? This catches
##      discrete confounding that need not fall on the leading GRM axes.
##
## The same measurements are made for null phenotypes. No marker or Stage-1
## association test is run here; the script only reads the completed bundles.
## =============================================================================
suppressMessages(library(data.table))

ROOT <- path.expand("~/gitlab/LDscnR-paper")
OUT_DIR <- file.path(ROOT, "module_9sp", "out", "02b_env_structure_alignment")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

N_AXES <- 5L
N_DRAWS <- 1000L

.r2 <- function(y, X) {
  if (!length(y) || stats::var(y) == 0) return(NA_real_)
  summary(stats::lm(y ~ X))$r.squared
}

.adj_r2 <- function(y, X) {
  if (!length(y) || stats::var(y) == 0) return(NA_real_)
  summary(stats::lm(y ~ X))$adj.r.squared
}

.group_r2 <- function(y, group) {
  if (!length(y) || stats::var(y) == 0 || length(unique(group)) < 2L) return(NA_real_)
  summary(stats::lm(y ~ factor(group)))$r.squared
}

.pop_mean <- function(y, pop, pops) {
  as.numeric(tapply(y, pop, mean)[pops])
}

.prepare_panel <- function(bundle_path, species, structure_var) {
  b <- readRDS(bundle_path)
  ph <- b$pheno
  stopifnot(nrow(b$GRM) == nrow(ph), length(b$eco) == nrow(ph),
            structure_var %in% names(ph))

  pops <- unique(ph$pop_ID)
  pop_index <- match(ph$pop_ID, pops)
  Z <- matrix(0, nrow(ph), length(pops))
  Z[cbind(seq_len(nrow(ph)), pop_index)] <- 1
  Z <- sweep(Z, 2, colSums(Z), "/")
  K_pop <- crossprod(Z, b$GRM %*% Z)
  K_pop <- (K_pop + t(K_pop)) / 2
  eg_pop <- eigen(K_pop, symmetric = TRUE)
  n_axes <- min(N_AXES, length(pops) - 1L)
  axes <- eg_pop$vectors[, seq_len(n_axes), drop = FALSE]

  list(species = species, structure_var = structure_var, pheno = ph,
       eco = as.numeric(b$eco), GRM = b$GRM, pops = pops,
       pop_index = pop_index, axes = axes)
}

.observed_row <- function(x, phenotype, y_ind) {
  y_pop <- .pop_mean(y_ind, x$pheno$pop_ID, x$pops)
  group_ind <- x$pheno[[x$structure_var]]
  group_pop <- group_ind[match(x$pops, x$pheno$pop_ID)]
  data.table(
    species = x$species,
    phenotype = phenotype,
    structure_group = x$structure_var,
    n_individuals = length(y_ind),
    n_populations = length(y_pop),
    n_structure_groups = uniqueN(group_pop),
    grm_axes_r2_population = .r2(y_pop, x$axes),
    grm_axes_adjusted_r2_population = .adj_r2(y_pop, x$axes),
    structure_r2_population = .group_r2(y_pop, group_pop),
    structure_r2_individual_weighted = .group_r2(y_ind, group_ind)
  )
}

.null_row <- function(x, scheme, draw, phenotype, y_ind) {
  y_pop <- .pop_mean(y_ind, x$pheno$pop_ID, x$pops)
  group_ind <- x$pheno[[x$structure_var]]
  group_pop <- group_ind[match(x$pops, x$pheno$pop_ID)]
  data.table(
    species = x$species, scheme = scheme, draw = draw,
    phenotype = phenotype,
    grm_axes_r2_population = .r2(y_pop, x$axes),
    structure_r2_population = .group_r2(y_pop, group_pop),
    structure_r2_individual_weighted = .group_r2(y_ind, group_ind)
  )
}

cat("=== 02b_env_structure_alignment ===\n\n")
cat("[1] reading completed empirical bundles\n")
x3 <- .prepare_panel(
  file.path(ROOT, "module_3sp", "out", "02_bundle", "bundle.rds"),
  "Three-spined stickleback", "pop_locality"
)
x9 <- .prepare_panel(
  file.path(ROOT, "module_9sp", "out", "02_bundle", "bundle.rds"),
  "Nine-spined stickleback", "lineage"
)

## Raw habitat is the pre-scan diagnostic for both panels. The lineage-adjusted
## phenotype is also recorded for 9sp because it is the quantity actually tested.
y9_adjusted <- as.numeric(stats::resid(stats::lm(x9$eco ~ factor(x9$pheno$lineage))))
observed <- rbindlist(list(
  .observed_row(x3, "Raw habitat", x3$eco),
  .observed_row(x9, "Raw habitat", x9$eco),
  .observed_row(x9, "Lineage-adjusted habitat", y9_adjusted)
))

cat("[2] generating group-preserving null phenotypes\n")
null_rows <- vector("list", N_DRAWS * 4L)
ii <- 0L

## Three-spined stickleback: shuffle population labels within geographic region,
## exactly matching the primary null used by the association analysis.
popt3 <- unique(x3$pheno[, .(pop_ID, ecotype = x3$eco,
                              pop_locality = get(x3$structure_var))])
for (bb in seq_len(N_DRAWS)) {
  set.seed(bb)
  z <- copy(popt3)[, permuted := sample(ecotype), by = pop_locality]
  y <- z$permuted[match(x3$pheno$pop_ID, z$pop_ID)]
  ii <- ii + 1L
  null_rows[[ii]] <- .null_row(x3, "Within-region label permutation", bb,
                                "Raw habitat", y)
}

## Nine-spined stickleback: shuffle population labels within lineage. Record
## alignment for both the raw labels and the lineage-adjusted phenotype tested.
popt9 <- unique(x9$pheno[, .(pop_ID, ecotype = x9$eco,
                              lineage = get(x9$structure_var))])
for (bb in seq_len(N_DRAWS)) {
  set.seed(bb)
  z <- copy(popt9)[, permuted := sample(ecotype), by = lineage]
  y_raw <- z$permuted[match(x9$pheno$pop_ID, z$pop_ID)]
  y_adjusted <- as.numeric(stats::resid(
    stats::lm(y_raw ~ factor(x9$pheno$lineage))))
  ii <- ii + 1L
  null_rows[[ii]] <- .null_row(x9, "Within-lineage label permutation", bb,
                                "Raw habitat", y_raw)
  ii <- ii + 1L
  null_rows[[ii]] <- .null_row(x9, "Within-lineage label permutation", bb,
                                "Lineage-adjusted habitat", y_adjusted)
}

## Nine-spined sensitivity null: draw a continuous phenotype with covariance K
## and remove its linear association with the observed tested phenotype, exactly
## matching 03c_EMMAX_structurednull.R's surrogate generator.
eg9 <- eigen(x9$GRM, symmetric = TRUE)
values9 <- pmax(eg9$values, 0)
vectors9 <- eg9$vectors
for (bb in seq_len(N_DRAWS)) {
  set.seed(bb)
  s <- as.numeric(vectors9 %*% (sqrt(values9) * stats::rnorm(length(x9$eco))))
  y <- as.numeric(stats::resid(stats::lm(s ~ y9_adjusted)))
  ii <- ii + 1L
  null_rows[[ii]] <- .null_row(x9, "GRM-matched continuous null", bb,
                                "Lineage-adjusted habitat", y)
}
null_draws <- rbindlist(null_rows[seq_len(ii)])
stopifnot(
  nrow(null_draws) == 4L * N_DRAWS,
  uniqueN(null_draws, by = c("species", "scheme", "draw", "phenotype")) == nrow(null_draws),
  !anyNA(null_draws[, .(grm_axes_r2_population, structure_r2_population,
                        structure_r2_individual_weighted)])
)

## At the population level, permuting labels within the null-defining groups
## preserves the raw habitat composition of every group exactly. Assert this
## design property rather than relying on a visual check.
raw_group_null <- null_draws[
  phenotype == "Raw habitat" & grepl("label permutation$", scheme)
]
raw_group_obs <- observed[phenotype == "Raw habitat",
                          .(species, structure_r2_population)]
raw_group_null <- merge(raw_group_null, raw_group_obs, by = "species")
stopifnot(max(abs(raw_group_null$structure_r2_population.x -
                  raw_group_null$structure_r2_population.y)) < 1e-12)

.summarise_metric <- function(metric) {
  null_draws[, {
    obs <- observed[species == .BY$species & phenotype == .BY$phenotype,
                    get(metric)]
    stopifnot(length(obs) == 1L)
    x <- get(metric)
    q <- quantile(x, c(0.025, 0.975), na.rm = TRUE)
    list(observed = obs, null_mean = mean(x, na.rm = TRUE),
         null_median = median(x, na.rm = TRUE),
         null_sd = stats::sd(x, na.rm = TRUE),
         null_q025 = q[1], null_q975 = q[2],
         observed_within_null95 = obs >= q[1] - 1e-12 && obs <= q[2] + 1e-12,
         proportion_null_ge_observed = mean(x >= obs - 1e-12, na.rm = TRUE))
  }, by = .(species, phenotype, scheme)][, metric := metric]
}

null_summary <- rbindlist(lapply(
  c("grm_axes_r2_population", "structure_r2_population",
    "structure_r2_individual_weighted"),
  .summarise_metric
), use.names = TRUE)
setcolorder(null_summary, c("species", "phenotype", "scheme", "metric",
                            "observed", "null_mean", "null_median",
                            "null_sd", "null_q025", "null_q975",
                            "observed_within_null95",
                            "proportion_null_ge_observed"))

fwrite(observed, file.path(OUT_DIR, "alignment_observed.tsv"), sep = "\t")
fwrite(null_draws, file.path(OUT_DIR, "alignment_null_draws.tsv"), sep = "\t")
fwrite(null_summary, file.path(OUT_DIR, "alignment_null_summary.tsv"), sep = "\t")
saveRDS(list(observed = observed, null_draws = null_draws,
             null_summary = null_summary,
             settings = list(n_axes = N_AXES, n_draws = N_DRAWS)),
        file.path(OUT_DIR, "alignment.rds"))

cat("[3] observed alignment\n")
print(observed)
cat("\n[4] null comparison for the simulation-comparable GRM-axis measure\n")
print(null_summary[metric == "grm_axes_r2_population"])
cat("\nWrote outputs to", OUT_DIR, "\n")
