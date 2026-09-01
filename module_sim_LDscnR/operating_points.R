## =============================================================================
## module_sim_LDscnR / operating_points.R
##
## What do the sensitivity grid's LEVELS actually retain on this dataset? The
## grid's numbers (ld_w_threshold 0.0125/0.025/0.05, size floor 1/2/5) are only
## comparable across datasets if they land on comparable OPERATING POINTS, and
## the levels are on scales that are not portable.
##
## THE UNIVERSE MATTERS AND IS EASY TO GET WRONG. With compute_unflagged_eMLG =
## TRUE the stage-2 partition covers EVERY mapped marker, not just the
## ld_w-flagged ones -- verified here by sum(lengths(members)) == nrow(map)
## exactly. So "fraction of clusters with >= 2 markers" has two different right
## answers (47.7% over all clusters, 53.4% over flagged clusters only) and any
## cross-dataset comparison must say which. Both are reported.
##
## Cluster size here is in RAW MARKERS, not in stage-1 eMLGs. The identity check
## above is what establishes that, and it is the check to run on any dataset
## before comparing size distributions.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/operating_points.R
## Env: SIM_DATA, CELLS, TAGS, ENVS, FILES, CORES, OUT
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(parallel)})
SIM   <- Sys.getenv("SIM_DATA", "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5")
OUT   <- Sys.getenv("OUT", "module_sim_LDscnR/results/operating_points")
CELLS <- strsplit(Sys.getenv("CELLS", "V0.5_c1,V0.5_c2,V1_c1.5,V2_c1"), ",")[[1]]
TAGS  <- strsplit(Sys.getenv("TAGS", "nobgs,bgs"), ",")[[1]]
ENVS  <- as.integer(strsplit(Sys.getenv("ENVS", "1"), ",")[[1]])
FILES <- as.integer(strsplit(Sys.getenv("FILES", "1,2,3,4,5,6,7,8,9,10"), ",")[[1]])
CORES <- as.integer(Sys.getenv("CORES", "8"))
THR   <- c(0.0125, 0.025, 0.05)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

jobs <- CJ(CELL = CELLS, TAG = TAGS, ENV = ENVS, i = FILES, sorted = FALSE)
one <- function(k) {
  r <- jobs[k]
  f <- sprintf("%s/adapt_%s_chr%d_%s_env%d.rds", SIM, r$TAG, r$i, r$CELL, r$ENV)
  if (!file.exists(f)) return(NULL)
  x <- readRDS(f); m <- as.data.table(x$map)
  pr <- ld_prune_and_eMLG(GTs = x$GTs, stage1 = x$complexity_reduction$stage1,
          LD_decay = x$LD_decay, ld_w_col = "ld_w_095", ld_w_threshold = 0.025,
          score_threshold = 0.80, min_r2_rho = 0.5, distance_threshold = 1e5,
          compute_unflagged_eMLG = TRUE, cores = 1)
  g <- as.data.table(pr$groups)
  ## The identity that pins the size unit. If this fails, "members" is not raw
  ## markers and no size distribution built from it is comparable to one that is.
  stopifnot(sum(lengths(g$members)) == nrow(m))
  list(sz = data.table(CELL = r$CELL, TAG = r$TAG, ENV = r$ENV, i = r$i,
                       sz = lengths(g$members),
                       kind = fifelse(grepl("^U", g$group_id), "unflagged", "flagged")),
       hd = data.table(CELL = r$CELL, TAG = r$TAG, ENV = r$ENV, i = r$i,
                       n_map = nrow(m),
                       ret_0.0125 = mean(m$ld_w_095 >= THR[1], na.rm = TRUE),
                       ret_0.025  = mean(m$ld_w_095 >= THR[2], na.rm = TRUE),
                       ret_0.05   = mean(m$ld_w_095 >= THR[3], na.rm = TRUE)))
}
res <- Filter(Negate(is.null), mclapply(seq_len(nrow(jobs)), one, mc.cores = CORES))
sz <- rbindlist(lapply(res, `[[`, "sz")); hd <- rbindlist(lapply(res, `[[`, "hd"))

summ <- function(d, universe) data.table(universe = universe,
  clusters = nrow(d), markers = sum(as.numeric(d$sz)),
  pct_cl_ge2  = 100 * mean(d$sz >= 2),  pct_cl_ge5 = 100 * mean(d$sz >= 5),
  pct_cl_ge10 = 100 * mean(d$sz >= 10),
  pct_mk_ge2  = 100 * sum(as.numeric(d$sz[d$sz >= 2])) / sum(as.numeric(d$sz)),
  pct_mk_ge5  = 100 * sum(as.numeric(d$sz[d$sz >= 5])) / sum(as.numeric(d$sz)),
  median = as.double(median(d$sz)), q90 = as.double(quantile(d$sz, .9)),
  q95 = as.double(quantile(d$sz, .95)), q99 = as.double(quantile(d$sz, .99)))
S <- rbind(summ(sz, "all"), summ(sz[kind == "flagged"], "flagged_only"),
           summ(sz[kind == "unflagged"], "unflagged_only"))
## ld_w retention weighted by markers per file, not a mean of per-file fractions.
R <- data.table(ld_w_threshold = THR, pct_markers_retained = sapply(
  c("ret_0.0125","ret_0.025","ret_0.05"),
  function(nm) 100 * sum(hd[[nm]] * hd$n_map) / sum(hd$n_map)))
fwrite(S, file.path(OUT, "cluster_size_operating_points.csv"))
fwrite(R, file.path(OUT, "ldw_threshold_retention.csv"))
print(S); print(R)
