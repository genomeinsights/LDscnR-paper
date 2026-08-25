## Resolve paths relative to this script, so the folder can be moved or handed on.
script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) normalizePath(dirname(f)) else getwd()
}
NEMO <- Sys.getenv("BGS_NEMO_DIR", dirname(script_dir()))   # R/ -> NEMO/

## ---------------------------------------------------------------------------
## make_run.R -- generate every .ini for a background-selection re-run
##
## Two phases, as in the production bgs2 runs:
##
##   burn-in   10000 generations, flat optimum, one per (arm, chr, V, c)
##   adapt      1000 generations, env optimum, one per (arm, chr, V, c, env)
##
## The burn-in is shared by all 10 envs of its cell, which is a 10x saving over
## running 11000 generations per env. It also puts deleterious selection where
## it matters: the adapt phase is ~0.05 Ne, so the diversity landscape is
## inherited from the burn-in.
##
## bgs and nobgs differ ONLY by the deleterious block, and share the random seed
## per (chr, V, c, env), so B_obs = pi_bgs / pi_nobgs is a matched-pair contrast.
##
##   Rscript make_run.R
##   BGS_CELLS=V1_c1.5 BGS_CHRS=1,2 BGS_ENVS=1 Rscript make_run.R   # a subset
## ---------------------------------------------------------------------------

MOD    <- NEMO
ROOT   <- Sys.getenv("BGS_ROOT",   file.path(MOD, "run"))
PARAMS <- Sys.getenv("BGS_PARAMS", "../params_V4")   # relative to the working dir at run time
BURNIN <- as.integer(Sys.getenv("BGS_BURNIN_GENS", "10000"))
## 0 = keep every segregating site, drop monomorphic (see ini/adapt.ini.tmpl).
## Set to 0.01 to reproduce the bgs2 output exactly.
TRIM_MAF <- Sys.getenv("BGS_TRIM_MAF", "0")

splitenv <- function(v, d) strsplit(Sys.getenv(v, d), ",")[[1]]
CHRS <- as.integer(splitenv("BGS_CHRS", paste(1:10, collapse = ",")))
ENVS <- as.integer(splitenv("BGS_ENVS", paste(1:10, collapse = ",")))
ARMS <- splitenv("BGS_ARMS", "bgs,nobgs")
## default scope mirrors the reduced set the analysis pipeline now runs
## (SIM_CELLS=V0.5_c2,V1_c1.5,V2_c1 in Nemo_sim/pipeline/README.md)
CELLS <- splitenv("BGS_CELLS", "V0.5_c2,V1_c1.5,V2_c1")

## --- scenario B parameters -------------------------------------------------
## u = 1.5e-5 over 1000 loci/chromosome, h*s median 5.0e-3 (185x the drift
## barrier), reach h*s = 0.50 cM. Predicted mean B 0.65, Q1/Q5 pi contrast 73%,
## Haldane load 5.8%. See bgs_predict.R.
DELET <- paste(
  "delet_dominance_mean          0.5",
  "delet_effects_dist_param1     -4.6",
  "delet_effects_dist_param2     1.5",
  "delet_effects_distribution    lognormal",
  sprintf("delet_genetic_map             &%s/map_delet__CHR__.txt", PARAMS),
  "delet_genetic_map_resolution  3.3e-5",
  "delet_loci                    2000",
  "delet_mutation_rate           1.5e-5",
  sep = "\n")

## the 80 sampled patches, copied verbatim from the production bgs2 logs
SAMPLE_PATCH <- paste0("{{", paste(c(
  197:200, 233:236, 245:248, 281:284, 293:296, 329:332, 341:344, 377:380,
  1079:1082, 1127:1130, 1175:1178, 1223:1226,
  1925:1928, 1961:1964, 1973:1976, 2009:2012, 2021:2024, 2057:2060, 2069:2072,
  2105:2108), collapse = ","), "}}")

arm_subs <- function(arm) {
  bgs <- identical(arm, "bgs")
  list(TAG = arm,
       DELET             = if (bgs) DELET else "",
       SEL_TRAIT         = if (bgs) "(quant, delet)"        else "quant",
       SEL_MODEL         = if (bgs) "(gaussian, direct)"    else "gaussian",
       STAT_DELET        = if (bgs) "adlt.delet "           else "",
       GENOTYPER_TRAITS  = if (bgs) "(quant, delet, ntrl)"  else "(quant, ntrl)",
       SOURCE_OVERRIDE   = if (bgs) "source_parameter_override delet_effects" else "")
}

## Deterministic and identical across arms, so the paired runs see the same
## sequence of demographic and mutational events.
seed_of <- function(chr, cell, env)
  110101L + 1000L * chr + 100L * match(cell, CELLS) + env

fill <- function(tmpl, subs) {
  for (k in names(subs)) tmpl <- gsub(paste0("__", k, "__"), subs[[k]], tmpl, fixed = TRUE)
  tmpl
}

tmpl_burnin <- paste(readLines(file.path(MOD, "ini", "burnin.ini.tmpl")), collapse = "\n")
tmpl_adapt  <- paste(readLines(file.path(MOD, "ini", "adapt.ini.tmpl")),  collapse = "\n")

dir.create(file.path(ROOT, "ini"), recursive = TRUE, showWarnings = FALSE)
burn <- adapt <- list()

for (arm in ARMS) {
  A <- arm_subs(arm)
  for (chr in CHRS) for (cell in CELLS) {
    V <- sub("^V", "", strsplit(cell, "_")[[1]][1])
    C <- sub("^c", "", strsplit(cell, "_")[[1]][2])
    base <- sprintf("%s_chr%d_V%s_c%s", arm, chr, V, C)

    subs <- c(A, list(CHR = chr, V = V, C = C, PARAMS = PARAMS,
                      GENERATIONS = BURNIN, SEED = seed_of(chr, cell, 0L),
                      ROOT = sprintf("./out/burnin_%s", base)))
    f <- file.path(ROOT, "ini", sprintf("burnin_%s.ini", base))
    writeLines(fill(tmpl_burnin, subs), f)
    burn[[length(burn) + 1]] <- data.frame(
      id = base, ini = sprintf("ini/burnin_%s.ini", base),
      expect = sprintf("out/burnin_%s/binary/%s_rep1_1.bin.bz2", base, base))

    for (env in ENVS) {
      aid  <- sprintf("adapt_%s_chr%d_V%s_c%s_env%d", arm, chr, V, C, env)
      subs <- c(A, list(CHR = chr, V = V, C = C, ENV = env, PARAMS = PARAMS,
                        SEED = seed_of(chr, cell, env),
                        ROOT = sprintf("./out/%s", aid),
                        SAMPLE_PATCH = SAMPLE_PATCH, TRIM_MAF = TRIM_MAF,
                        SOURCE = sprintf("./out/burnin_%s/binary/%s_rep1", base, base)))
      writeLines(fill(tmpl_adapt, subs), file.path(ROOT, "ini", sprintf("%s.ini", aid)))
      adapt[[length(adapt) + 1]] <- data.frame(
        id = aid, ini = sprintf("ini/%s.ini", aid),
        expect = sprintf("out/%s/GENO/%s_1000_1.snp_geno", aid, aid),
        dep = base)
    }
  }
}

write.table(do.call(rbind, burn),  file.path(ROOT, "manifest_burnin.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(do.call(rbind, adapt), file.path(ROOT, "manifest_adapt.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat(sprintf("arms %s | chr %s | cells %s | envs %s\n",
            paste(ARMS, collapse = ","), paste(range(CHRS), collapse = "-"),
            paste(CELLS, collapse = ","), paste(range(ENVS), collapse = "-")))
cat(sprintf("%d burn-ins (%d gens) + %d adapt runs | trim_maf %s -> %s\n",
            length(burn), BURNIN, length(adapt), TRIM_MAF, ROOT))
