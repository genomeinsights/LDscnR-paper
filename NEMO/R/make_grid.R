## ---------------------------------------------------------------------------
## make_grid.R -- a one-replicate parameter grid, to check empirically which
## deleterious settings actually produce a detectable BGS signature.
##
## One cell (V2_c1 by default), one chromosome, one env, one replicate. Each
## deleterious variant gets its OWN burn-in, because the architecture has to be
## in place while diversity builds. All variants share a SINGLE nobgs run as the
## paired denominator -- nobgs has no deleterious trait, so it is the same run
## whatever the variant.
##
## Replication comes from windows, not runs: one chromosome gives ~59 windows of
## 500 kb, which is what carries the pi-vs-recombination contrast. With a single
## replicate the between-variant differences also contain realisation noise, so
## read the ranking and the magnitude, not small differences.
##
##   Rscript make_grid.R
##   BGS_GRID_CHR=5 BGS_GRID_GENS=4000 Rscript make_grid.R
## ---------------------------------------------------------------------------

suppressMessages(library(data.table))

script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) normalizePath(dirname(f)) else getwd()
}
NEMO <- Sys.getenv("BGS_NEMO_DIR", dirname(script_dir()))

CHR    <- as.integer(Sys.getenv("BGS_GRID_CHR",  "8"))
CELL   <- Sys.getenv("BGS_GRID_CELL", "V2_c1")
ENV    <- as.integer(Sys.getenv("BGS_GRID_ENV",  "1"))
GENS   <- as.integer(Sys.getenv("BGS_GRID_GENS", "10000"))
ROOT   <- Sys.getenv("BGS_GRID_ROOT", file.path(NEMO, "grid"))
TRIM   <- Sys.getenv("BGS_TRIM_MAF", "0")

## chr8 is the default because its Nemo map (26.7 cM) is closest to the mean over
## the ten sets, so its BGS strength is representative rather than extreme
## (chr4, 63 cM, is the weakest; chr5, 16 cM, the strongest).

## --- the grid --------------------------------------------------------------
## Spans U from 0.001 to 0.06 diploid per chromosome. The DFE is held fixed at
## lognormal(-4.6, 1.5) so that n_delet and u are the only things varying.
GRID <- data.table(
  n_delet = c(  100,   100,  1000,   1000,  1000),
  u       = c(3e-5,  1e-4,  5e-6, 1.5e-5,  3e-5))
GRID[, `:=`(U = 2 * n_delet * u,
            tag = sprintf("n%d_u%s", n_delet, format(u, scientific = TRUE, trim = TRUE)))]
GRID[, tag := gsub("[-+.]", "", tag)]

MEANLOG <- -4.6; SDLOG <- 1.5; H <- 0.5

## --- deleterious maps as nested subsets ------------------------------------
## params_V4 carries 1000 loci per chromosome, the first 100 of which are exactly
## the params_V3 set (extend_delet_map.R kept them and added 900). Taking the
## first n of the sorted positions therefore gives a nested series: every smaller
## variant is a subset of every larger one, so the variants differ in density and
## rate alone and not in which regions happen to carry deleterious loci.
rd <- readRDS(file.path(NEMO, "params_V4", "rds", sprintf("rec_map%d.rds", CHR)))
dir.create(file.path(ROOT, "params"), recursive = TRUE, showWarnings = FALSE)

for (n in unique(GRID$n_delet)) {
  blocks <- vapply(1:2, function(ch) {
    pos <- rd[Chr == ch & type == "delet"][order(bp)][seq_len(n), pos_nemo]
    paste0("{", paste(sort(pos), collapse = ", "), "}")
  }, character(1))
  writeLines(paste0("{", paste(blocks, collapse = ""), "}"),
             file.path(ROOT, "params", sprintf("map_delet%d_n%d.txt", CHR, n)))
}
for (f in c(sprintf("map_ntrl%d.txt", CHR), sprintf("map_QTN%d.txt", CHR),
            sprintf("allelic_values_%d.txt", CHR),
            sprintf("disp_mat_%s.txt", sub("^.*_c", "", CELL)),
            sprintf("env_%d.txt", ENV)))
  if (!file.copy(file.path(NEMO, "params_V4", f), file.path(ROOT, "params", f), overwrite = TRUE))
    stop("could not copy required input: ", f)

## --- ini generation --------------------------------------------------------
V <- sub("^V", "", strsplit(CELL, "_")[[1]][1])
C <- sub("^c", "", strsplit(CELL, "_")[[1]][2])
PARAMS <- "./params"   # run_pool.sh runs from grid/, and params/ sits inside it
SEED <- 110101L + 1000L * CHR + ENV

SAMPLE_PATCH <- paste0("{{", paste(c(
  197:200, 233:236, 245:248, 281:284, 293:296, 329:332, 341:344, 377:380,
  1079:1082, 1127:1130, 1175:1178, 1223:1226,
  1925:1928, 1961:1964, 1973:1976, 2009:2012, 2021:2024, 2057:2060, 2069:2072,
  2105:2108), collapse = ","), "}}")

fill <- function(tmpl, subs) {
  for (k in names(subs)) tmpl <- gsub(paste0("__", k, "__"), subs[[k]], tmpl, fixed = TRUE)
  tmpl
}
tmpl_burnin <- paste(readLines(file.path(NEMO, "ini", "burnin.ini.tmpl")), collapse = "\n")
tmpl_adapt  <- paste(readLines(file.path(NEMO, "ini", "adapt.ini.tmpl")),  collapse = "\n")

delet_block <- function(n, u) paste(
  sprintf("delet_dominance_mean          %s", H),
  sprintf("delet_effects_dist_param1     %s", MEANLOG),
  sprintf("delet_effects_dist_param2     %s", SDLOG),
  "delet_effects_distribution    lognormal",
  sprintf("delet_genetic_map             &%s/map_delet%d_n%d.txt", PARAMS, CHR, n),
  "delet_genetic_map_resolution  3.3e-5",
  sprintf("delet_loci                    %d", 2 * n),
  sprintf("delet_mutation_rate           %s", format(u, scientific = TRUE, trim = TRUE)),
  sep = "\n")

arm_subs <- function(bgs, n = NA, u = NA) list(
  DELET            = if (bgs) delet_block(n, u) else "",
  SEL_TRAIT        = if (bgs) "(quant, delet)"       else "quant",
  SEL_MODEL        = if (bgs) "(gaussian, direct)"   else "gaussian",
  STAT_DELET       = if (bgs) "adlt.delet "          else "",
  GENOTYPER_TRAITS = if (bgs) "(quant, delet, ntrl)" else "(quant, ntrl)",
  SOURCE_OVERRIDE  = if (bgs) "source_parameter_override delet_effects" else "")

dir.create(file.path(ROOT, "ini"), recursive = TRUE, showWarnings = FALSE)
burn <- adapt <- list()

emit <- function(tag, A) {
  base <- sprintf("%s_chr%d_%s", tag, CHR, CELL)
  writeLines(fill(tmpl_burnin, c(A, list(
    TAG = tag, CHR = CHR, V = V, C = C, PARAMS = PARAMS, GENERATIONS = GENS,
    SEED = SEED, ROOT = sprintf("./out/burnin_%s", base)))),
    file.path(ROOT, "ini", sprintf("burnin_%s.ini", base)))
  burn[[length(burn) + 1]] <<- data.frame(id = base,
    ini = sprintf("ini/burnin_%s.ini", base),
    expect = sprintf("out/burnin_%s/binary/%s_rep1_1.bin.bz2", base, base))

  aid <- sprintf("adapt_%s_chr%d_%s_env%d", tag, CHR, CELL, ENV)
  writeLines(fill(tmpl_adapt, c(A, list(
    TAG = tag, CHR = CHR, V = V, C = C, ENV = ENV, PARAMS = PARAMS, SEED = SEED,
    ROOT = sprintf("./out/%s", aid), SAMPLE_PATCH = SAMPLE_PATCH, TRIM_MAF = TRIM,
    SOURCE = sprintf("./out/burnin_%s/binary/%s_rep1", base, base)))),
    file.path(ROOT, "ini", sprintf("%s.ini", aid)))
  adapt[[length(adapt) + 1]] <<- data.frame(id = aid,
    ini = sprintf("ini/%s.ini", aid),
    expect = sprintf("out/%s/GENO/%s_1000_1.snp_geno", aid, aid), dep = base)
}

## the shared paired control, then one run per variant
emit("nobgs", arm_subs(FALSE))
for (i in seq_len(nrow(GRID))) emit(GRID$tag[i], arm_subs(TRUE, GRID$n_delet[i], GRID$u[i]))

write.table(do.call(rbind, burn),  file.path(ROOT, "manifest_burnin.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(do.call(rbind, adapt), file.path(ROOT, "manifest_adapt.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
fwrite(GRID, file.path(ROOT, "grid.tsv"), sep = "\t")

cat(sprintf("chr%d %s env%d | burn-in %d gens | trim_maf %s\n", CHR, CELL, ENV, GENS, TRIM))
print(GRID[, .(tag, n_delet, u, U_per_chr = U, load_pct = round(100 * (1 - exp(-2 * U)), 1))])
cat(sprintf("\n%d burn-ins + %d adapt runs -> %s\n", length(burn), length(adapt), ROOT))
